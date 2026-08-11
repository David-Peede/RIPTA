"""VCF processing and allele count extraction for RIPTA."""

import gzip
import multiprocessing as mp
import random
from collections import defaultdict
from itertools import chain, combinations

from ripta.stats import AlleleCountBlock, AlleleCountBlockResult

# Module-level globals populated by `_init_pool_worker` at pool startup and
# read by `_process_block_lines` when called with no explicit `shared`
# dictionary. These are used only by the multi-threaded path; the
# single-threaded path passes shared state via an explicit `shared` argument
# and never sets or reads these globals.
_pool_pop_indices = None
_pool_gp_inds = None
_pool_all_pops = None
_pool_all_pairs = None
_pool_sp_configs = None
_pool_n_seg_keys = None
_pool_fst_pairs = None
_pool_f2_pairs = None
_pool_f3_triples = None
_pool_f4_quads = None
_pool_block_size = None


def _build_pop_indices(config, all_pops):
    """Build a mapping from population labels to lists of VCF column indices.

    Args:
        config (dict): Validated configuration dictionary with 'pop_to_inds'
            (populated in config.py Phase 2) and 'idx_dicc'.
        all_pops (set): Set of population labels required by the requested
            statistics, as returned by `_collect_required_sets`. Only these
            populations are included in the returned mapping.

    Returns:
        pop_indices (dict): Mapping of population label strings to lists of
            0-based VCF column indices for each individual in that population.
    """
    # Retrieve the population-to-individuals mapping built in config.py Phase 2.
    pop_to_inds = config["pop_to_inds"]
    # Retrieve the individual-to-VCF-column-index mapping built in config.py Phase 4.
    idx_dicc = config["idx_dicc"]
    pop_indices = {}
    # Only populations actually needed by the requested statistics.
    for pop in all_pops:
        # Skip populations that are not present in the metadata mapping.
        if pop not in pop_to_inds:
            continue
        # Resolve each individual label to its VCF column index, skipping any
        # individual not present in the VCF header.
        indices = [idx_dicc[ind] for ind in pop_to_inds[pop] if ind in idx_dicc]
        if indices:
            pop_indices[pop] = indices
    return pop_indices


def _collect_required_sets(config):
    """Determine all populations, pairs, site pattern configs, n_seg keys, and f-stat keys.

    Args:
        config (dict): Validated configuration dictionary with populated configuration
            objects (sp_configs, pi_configs, dxy_configs, etc.).

    Returns:
        all_pops (set): Set of all unique population labels needed.
        all_pairs (set): Set of alphabetically sorted (pop_x, pop_y) tuples
            for all required between-population pairs.
        sp_configs (list): List of SitePatternConfig instances.
        n_seg_keys (set): Set of frozensets, one per unique population
            combination across all requested statistics.
        fst_pairs (set): Set of sorted (pop_x, pop_y) tuples for Fst
            per-site accumulation.
        f2_pairs (set): Set of sorted (pop_a, pop_b) tuples for f2
            per-site accumulation.
        f3_triples (set): Set of (pop_c, pop_a, pop_b) ordered tuples for
            f3 per-site accumulation.
        f4_quads (set): Set of (pop_a, pop_b, pop_c, pop_d) ordered tuples
            for f4 per-site accumulation.
    """
    all_pops = set()
    # Store pairs as sorted tuples for a consistent key ordering.
    all_pairs = set()
    n_seg_keys = set()

    def _register_pops(pops):
        """Register a set of population labels and their jointly-callable n_seg key."""
        all_pops.update(pops)
        n_seg_keys.add(frozenset(pops))

    # Collect from site pattern configs.
    sp_configs = list(config.get("sp_configs", []))
    for sp in sp_configs:
        _register_pops({sp.p1, sp.p2, sp.p3, sp.outgroup})

    # Collect from pairwise configs: pi, dxy, fst, D3.
    for cfg in config.get("pi_configs", []):
        _register_pops({cfg.pop})
    for cfg in config.get("dxy_configs", []):
        _register_pops({cfg.pop_x, cfg.pop_y})
    for cfg in config.get("fst_configs", []):
        _register_pops({cfg.pop_x, cfg.pop_y})
    for cfg in config.get("d3_configs", []):
        _register_pops({cfg.p1, cfg.p2, cfg.p3})

    # Collect from f-statistic configs: f2, f3, f3star, f4, f4ratio.
    for cfg_key in (
        "f2_configs",
        "f3_configs",
        "f3star_configs",
        "f4_configs",
        "f4ratio_configs",
    ):
        for cfg in config.get(cfg_key, []):
            pops = set()
            for field in ("pop_a", "pop_b", "pop_c", "pop_d", "pop_e"):
                val = getattr(cfg, field, "")
                if val:
                    pops.add(val)
            _register_pops(pops)

    # Build between-population pairs only for statistics that actually need
    # per-pair dxy accumulators: dxy itself and D3 (which reuses dxy sums).
    for cfg in config.get("dxy_configs", []):
        all_pairs.add(tuple(sorted([cfg.pop_x, cfg.pop_y])))
    for cfg in config.get("d3_configs", []):
        # D3 uses all three pairwise dxy values between p1, p2, and p3.
        for px, py in combinations(sorted({cfg.p1, cfg.p2, cfg.p3}), 2):
            all_pairs.add((px, py))

    # Collect Fst-specific population pairs for per-site accumulation.
    fst_pairs = set()
    for cfg in config.get("fst_configs", []):
        fst_pairs.add(tuple(sorted([cfg.pop_x, cfg.pop_y])))

    # Collect f2 population pairs for per-site accumulation.
    f2_pairs = set()
    for cfg in config.get("f2_configs", []):
        f2_pairs.add(tuple(sorted([cfg.pop_a, cfg.pop_b])))

    # Collect f3 population triples for per-site accumulation.
    # Key order is (pop_c, pop_a, pop_b) to match the f3 formula.
    f3_triples = set()
    for cfg in config.get("f3_configs", []):
        f3_triples.add((cfg.pop_c, cfg.pop_a, cfg.pop_b))
    # f3star uses the same f3 computation with the same population roles.
    for cfg in config.get("f3star_configs", []):
        f3_triples.add((cfg.pop_c, cfg.pop_a, cfg.pop_b))

    # Collect f4 population quads for per-site accumulation.
    # Key order is (pop_a, pop_b, pop_c, pop_d) matching f4(A, B; C, D).
    f4_quads = set()
    for cfg in config.get("f4_configs", []):
        f4_quads.add((cfg.pop_a, cfg.pop_b, cfg.pop_c, cfg.pop_d))
    # f4ratio requires two f4 computations per config entry.
    for cfg in config.get("f4ratio_configs", []):
        # Numerator f4: f4(pop_a, pop_d; pop_e, pop_c).
        f4_quads.add((cfg.pop_a, cfg.pop_d, cfg.pop_e, cfg.pop_c))
        # Denominator f4: f4(pop_a, pop_d; pop_b, pop_c).
        f4_quads.add((cfg.pop_a, cfg.pop_d, cfg.pop_b, cfg.pop_c))

    return (
        all_pops,
        all_pairs,
        sp_configs,
        n_seg_keys,
        fst_pairs,
        f2_pairs,
        f3_triples,
        f4_quads,
    )


def _init_block_accumulator(
    all_pops,
    all_pairs,
    sp_configs,
    n_seg_keys,
    fst_pairs,
    f2_pairs,
    f3_triples,
    f4_quads,
):
    """Create a fresh block accumulator with zeroed sums for all fields.

    Args:
        all_pops (set): Set of all population labels to track.
        all_pairs (set): Set of sorted (pop_x, pop_y) tuples for dxy sums.
        sp_configs (list): List of SitePatternConfig instances.
        n_seg_keys (set): Set of frozensets for n_seg tracking.
        fst_pairs (set): Set of sorted (pop_x, pop_y) tuples for Fst sums.
        f2_pairs (set): Set of sorted (pop_a, pop_b) tuples for f2 sums.
        f3_triples (set): Set of (pop_c, pop_a, pop_b) tuples for f3 sums.
        f4_quads (set): Set of (pop_a, pop_b, pop_c, pop_d) tuples for f4 sums.

    Returns:
        acc (dict): Mutable accumulator dict with zeroed sums.
    """
    return {
        "sum_aac": {pop: 0.0 for pop in all_pops},
        "sum_c": {pop: 0.0 for pop in all_pops},
        "sum_het_num": {pop: 0.0 for pop in all_pops},
        "sum_pi_diffs": {pop: 0.0 for pop in all_pops},
        "sum_pi_comps": {pop: 0.0 for pop in all_pops},
        "sum_dxy_diffs": {pair: 0.0 for pair in all_pairs},
        "sum_dxy_comps": {pair: 0.0 for pair in all_pairs},
        # Each SitePatternConfig maps to a defaultdict of pattern counts.
        "site_patterns": {sp: defaultdict(float) for sp in sp_configs},
        "n_seg": {key: 0 for key in n_seg_keys},
        # Number of sites that passed the non-biallelic and all-missing filters in this block.
        "site_count": 0,
        # Fst per-site numerator and denominator sums.
        "sum_fst_num": {pair: 0.0 for pair in fst_pairs},
        "sum_fst_num_unbiased": {pair: 0.0 for pair in fst_pairs},
        "sum_fst_denom_unbiased": {pair: 0.0 for pair in fst_pairs},
        # f2 per-site sums.
        "sum_f2": {pair: 0.0 for pair in f2_pairs},
        "sum_f2_unbiased": {pair: 0.0 for pair in f2_pairs},
        # f3 per-site sums.
        "sum_f3": {triple: 0.0 for triple in f3_triples},
        "sum_f3_unbiased": {triple: 0.0 for triple in f3_triples},
        # f4 per-site sums.
        "sum_f4": {quad: 0.0 for quad in f4_quads},
        # Large-sample pi per-site sum: sum_m 2 * p_m,x * (1 - p_m,x) (eq:pi-sample).
        "sum_pi_large": {pop: 0.0 for pop in all_pops},
        # Large-sample dxy per-site sum: sum_m (p_m,x*(1-p_m,y) + p_m,y*(1-p_m,x)) (eq:dxy-site-large).
        "sum_dxy_large": {pair: 0.0 for pair in all_pairs},
        # Large-sample Fst denominator per-site sum: sum_m D_m (eq:fst-D-large).
        "sum_fst_denom": {pair: 0.0 for pair in fst_pairs},
        # Large-sample half-heterozygosity for f3star denominator (eq:h-large).
        "sum_h_large": {pop: 0.0 for pop in all_pops},
        # Unbiased h_hat ratio-of-sums numerator for f3star denominator (eq:h-hat numerator).
        "sum_h_hat_num": {pop: 0.0 for pop in all_pops},
        # Unbiased h_hat ratio-of-sums denominator for f3star denominator (eq:h-hat denominator).
        "sum_h_hat_denom": {pop: 0.0 for pop in all_pops},
    }


def _compute_site_patterns(patterns, sp_config, site_aac, site_c, warnings, chrom, pos):
    """Update site pattern accumulators for one SitePatternConfig at one site.

    Args:
        patterns (defaultdict): Mutable defaultdict(float) of pattern counts
            for this SitePatternConfig, updated in place.
        sp_config (SitePatternConfig): The site pattern configuration.
        site_aac (dict): Per-population alternative allele counts at this site.
        site_c (dict): Per-population chromosome counts at this site.
        warnings (list): Mutable list of warning tuples, appended in place.
        chrom (str): Chromosome name from the VCF CHROM column.
        pos (str): Position string from the VCF POS column.

    Notes
    -----
    All site patterns are computed using derived allele frequencies for all
    four populations (P1, P2, P3, and the outgroup) per (eq:d-data) from
    RIPTA-equations. The derived allele is defined as the major allele in the
    outgroup population (the more common of the reference and alternative
    alleles); when the two allele frequencies are exactly tied at 0.5, the
    derived allele is chosen at random. Once the derived allele is fixed,
    each population's derived allele frequency is either the alternative
    allele frequency (when derived equals alternative) or its complement
    (when derived equals reference).

    Homogenized site patterns replace P2 with P3 in the frequency product;
    they are used in the denominator of fhom (eq:f-hom-data) from RIPTA-equations.
    """
    # Retrieve allele counts for the three focal populations.
    aac_p1 = site_aac.get(sp_config.p1, 0.0)
    c_p1 = site_c.get(sp_config.p1, 0)
    aac_p2 = site_aac.get(sp_config.p2, 0.0)
    c_p2 = site_c.get(sp_config.p2, 0)
    aac_p3 = site_aac.get(sp_config.p3, 0.0)
    c_p3 = site_c.get(sp_config.p3, 0)

    # Skip this site if any focal population has no called genotypes.
    if c_p1 == 0 or c_p2 == 0 or c_p3 == 0:
        return

    # Retrieve outgroup allele counts used to define the derived allele.
    aac_o = site_aac.get(sp_config.outgroup, 0.0)
    c_o = site_c.get(sp_config.outgroup, 0)

    # Skip this site if the outgroup has no called genotypes.
    if c_o == 0:
        return

    # Alternative allele frequency in the outgroup, used to choose the
    # major (derived) allele.
    freq_o_alt = aac_o / c_o

    # Decide whether the derived allele is the alternative allele (when the
    # alternative is the major allele in the outgroup) or the reference
    # allele (when the reference is the major allele). A tie at 0.5 is
    # broken by a uniform random draw.
    if freq_o_alt > 0.5:
        # Alternative allele is the minor allele in the outgroup; the
        # reference is the derived allele, so all derived allele frequencies
        # are complements of the alternative allele frequencies.
        derived_is_alt = False
    elif freq_o_alt < 0.5:
        # Alternative allele is the major allele in the outgroup and is
        # therefore the derived allele.
        derived_is_alt = True
    else:
        # Exact tie: choose the derived allele with a uniform random draw.
        derived_is_alt = random.random() < 0.5

    # Alternative allele frequencies for all four populations.
    freq_p1_alt = aac_p1 / c_p1
    freq_p2_alt = aac_p2 / c_p2
    freq_p3_alt = aac_p3 / c_p3

    # Derived allele frequencies: identical to the alternative allele
    # frequencies when the derived allele is the alternative, otherwise
    # their complement.
    if derived_is_alt:
        p1 = freq_p1_alt
        p2 = freq_p2_alt
        p3 = freq_p3_alt
        p_o = freq_o_alt
    else:
        p1 = 1.0 - freq_p1_alt
        p2 = 1.0 - freq_p2_alt
        p3 = 1.0 - freq_p3_alt
        p_o = 1.0 - freq_o_alt

    # Complement of the outgroup derived allele frequency, reused across all
    # six unpolarized site pattern products (eq:d-data) from RIPTA-equations.
    q_o = 1.0 - p_o

    # Standard site patterns (eq:d-data) from RIPTA-equations.
    # BBAA: P1=derived, P2=derived, P3=ancestral, O=ancestral.
    patterns["bbaa"] += p1 * p2 * (1 - p3) * q_o
    # BABA: P1=derived, P2=ancestral, P3=derived, O=ancestral.
    patterns["baba"] += p1 * (1 - p2) * p3 * q_o
    # ABBA: P1=ancestral, P2=derived, P3=derived, O=ancestral.
    patterns["abba"] += (1 - p1) * p2 * p3 * q_o
    # BAAA: P1=derived, P2=ancestral, P3=ancestral, O=ancestral.
    patterns["baaa"] += p1 * (1 - p2) * (1 - p3) * q_o
    # ABAA: P1=ancestral, P2=derived, P3=ancestral, O=ancestral.
    patterns["abaa"] += (1 - p1) * p2 * (1 - p3) * q_o
    # AABA: P1=ancestral, P2=ancestral, P3=derived, O=ancestral.
    patterns["aaba"] += (1 - p1) * (1 - p2) * p3 * q_o

    # Homogenized site patterns: P2 replaced by P3 in the frequency product
    # (eq:f-hom-data) from RIPTA-equations. Used only in the denominator of fhom.
    # ABBA_hom: P1=ancestral, P3=derived, P3=derived, O=ancestral.
    patterns["abba_hom"] += (1 - p1) * p3 * p3 * q_o
    # BABA_hom: P1=derived, P3=ancestral, P3=derived, O=ancestral.
    patterns["baba_hom"] += p1 * (1 - p3) * p3 * q_o


def _finalize_block(acc, chrom, block_idx, block_size):
    """Convert a mutable block accumulator into a frozen AlleleCountBlockResult.

    Args:
        acc (dict): Mutable accumulator dict with accumulated sums.
        chrom (str): Chrom name from the VCF CHROM column.
        block_idx (int): 0-based block index within this chrom.
        block_size (int): Block size in base pairs.

    Returns:
        result (AlleleCountBlockResult or None): The finalized block result,
            or None if the block contained no segregating sites.
    """
    # Skip empty blocks where no segregating sites were processed.
    if acc["site_count"] == 0:
        return None

    # Compute 1-based start and end positions for this block.
    block_start = block_idx * block_size + 1
    block_end = block_start + block_size - 1

    # Convert defaultdict site pattern accumulators to regular dicts.
    site_patterns = {
        sp: dict(patterns) for sp, patterns in acc["site_patterns"].items()
    }

    # Construct the frozen AlleleCountBlock from accumulated sums.
    block = AlleleCountBlock(
        sum_aac=dict(acc["sum_aac"]),
        sum_c=dict(acc["sum_c"]),
        sum_het_num=dict(acc["sum_het_num"]),
        sum_pi_diffs=dict(acc["sum_pi_diffs"]),
        sum_pi_comps=dict(acc["sum_pi_comps"]),
        sum_dxy_diffs=dict(acc["sum_dxy_diffs"]),
        sum_dxy_comps=dict(acc["sum_dxy_comps"]),
        site_patterns=site_patterns,
        n_seg=dict(acc["n_seg"]),
        sum_fst_num=dict(acc["sum_fst_num"]),
        sum_fst_num_unbiased=dict(acc["sum_fst_num_unbiased"]),
        sum_fst_denom_unbiased=dict(acc["sum_fst_denom_unbiased"]),
        sum_f2=dict(acc["sum_f2"]),
        sum_f2_unbiased=dict(acc["sum_f2_unbiased"]),
        sum_f3=dict(acc["sum_f3"]),
        sum_f3_unbiased=dict(acc["sum_f3_unbiased"]),
        sum_f4=dict(acc["sum_f4"]),
        sum_pi_large=dict(acc["sum_pi_large"]),
        sum_dxy_large=dict(acc["sum_dxy_large"]),
        sum_fst_denom=dict(acc["sum_fst_denom"]),
        sum_h_large=dict(acc["sum_h_large"]),
        sum_h_hat_num=dict(acc["sum_h_hat_num"]),
        sum_h_hat_denom=dict(acc["sum_h_hat_denom"]),
    )

    return AlleleCountBlockResult(
        chrom=chrom,
        block_start=block_start,
        block_end=block_end,
        block=block,
    )


def _produce_blocks(vcf_path, target_chroms, block_size):
    """Yield per-block groups of raw VCF data lines from one VCF file.

    Args:
        vcf_path (str): Path to the VCF file (plain text or gzip-compressed).
        target_chroms (set): Set of chrom names to include; chroms not in
            this set are silently skipped.
        block_size (int): Block size in base pairs used to compute the
            0-based block index `(pos - 1) // block_size`.

    Yields:
        block_tuple (tuple): A `(chrom, block_idx, block_lines)` tuple where
            `chrom` is the chrom name, `block_idx` is the 0-based block index,
            and `block_lines` is a list of raw VCF data line strings belonging
            to that block. Empty blocks are never yielded.

    Notes
    -----
    Groups consecutive VCF data lines that share the same `(chrom, block_idx)`
    into a single list, yielding one tuple per group. A chrom transition
    always starts a new block regardless of the block index on either side.
    Lines are passed through untouched: no strip, no full tab split, no
    per-site filtering — all biallelic, missingness, and FORMAT handling is
    performed downstream in `_process_block_lines`. The partial split
    `line.split('\\t', 2)` is used only to extract chrom and position for
    grouping.
    """
    # Select the appropriate file opener based on compression extension.
    opener = gzip.open if vcf_path.endswith(".gz") else open

    # Current grouping state: chrom, block index, and accumulated line list.
    current_chrom = None
    current_block_idx = None
    current_lines = []

    with opener(vcf_path, "rt") as vcf_handle:
        for line in vcf_handle:
            # Skip all VCF header and meta-information lines.
            if line.startswith("#"):
                continue

            # Partial tab split: only enough to extract chrom and position.
            parts = line.split("\t", 2)
            chrom = parts[0]

            # Skip chroms that are not in the requested target set.
            if chrom not in target_chroms:
                continue

            # Parse the position as an integer for block index computation.
            pos = int(parts[1])

            # Compute the 0-based block index from the genomic position.
            block_idx = (pos - 1) // block_size

            # A chrom transition always starts a new block; within the same
            # chrom, any change in block index also starts a new block.
            if chrom != current_chrom or block_idx != current_block_idx:
                # Yield the prior accumulated block if it held any lines.
                if current_lines:
                    yield (current_chrom, current_block_idx, current_lines)
                # Reset grouping state to start a new block.
                current_chrom = chrom
                current_block_idx = block_idx
                current_lines = []

            # Append the raw line untouched to the current block group.
            current_lines.append(line)

    # After end of file, yield the final accumulated block if non-empty.
    if current_lines:
        yield (current_chrom, current_block_idx, current_lines)


def _init_pool_worker(
    pop_indices,
    gp_inds,
    all_pops,
    all_pairs,
    sp_configs,
    n_seg_keys,
    fst_pairs,
    f2_pairs,
    f3_triples,
    f4_quads,
    block_size,
):
    """Populate module-level pool globals inside one worker process at pool startup.

    Args:
        pop_indices (dict): Mapping of population label to list of VCF column
            indices for individuals in that population.
        gp_inds (dict): Mapping of population label to set of individual IDs
            for populations using genotype probabilities.
        all_pops (set): Set of all unique population labels to extract.
        all_pairs (set): Set of alphabetically sorted (pop_x, pop_y) tuples
            for dxy sums.
        sp_configs (list): List of SitePatternConfig instances.
        n_seg_keys (set): Set of frozensets for n_seg tracking.
        fst_pairs (set): Set of alphabetically sorted (pop_x, pop_y) tuples
            for Fst per-site accumulation.
        f2_pairs (set): Set of alphabetically sorted (pop_a, pop_b) tuples
            for f2 per-site accumulation.
        f3_triples (set): Set of (pop_c, pop_a, pop_b) tuples for f3
            per-site accumulation.
        f4_quads (set): Set of (pop_a, pop_b, pop_c, pop_d) tuples for f4
            per-site accumulation.
        block_size (int): Block size in base pairs used for block index math
            and block-end computation.

    Returns:
        None: Side effect only; populates module-level `_pool_*` globals in
            the worker process.
    """
    # Declare the module-level globals so assignments populate them in the
    # worker process.
    global _pool_pop_indices
    global _pool_gp_inds
    global _pool_all_pops
    global _pool_all_pairs
    global _pool_sp_configs
    global _pool_n_seg_keys
    global _pool_fst_pairs
    global _pool_f2_pairs
    global _pool_f3_triples
    global _pool_f4_quads
    global _pool_block_size

    # Store the population-label-to-column-indices mapping for worker reads.
    _pool_pop_indices = pop_indices
    # Store the genotype-probability individuals mapping for worker reads.
    _pool_gp_inds = gp_inds
    # Store the set of all required population labels for worker reads.
    _pool_all_pops = all_pops
    # Store the set of required between-population pairs for dxy sums.
    _pool_all_pairs = all_pairs
    # Store the list of site pattern configuration objects.
    _pool_sp_configs = sp_configs
    # Store the set of frozenset keys used for the n_seg tracker.
    _pool_n_seg_keys = n_seg_keys
    # Store the set of population pairs requiring Fst per-site accumulation.
    _pool_fst_pairs = fst_pairs
    # Store the set of population pairs requiring f2 per-site accumulation.
    _pool_f2_pairs = f2_pairs
    # Store the set of population triples requiring f3 per-site accumulation.
    _pool_f3_triples = f3_triples
    # Store the set of population quads requiring f4 per-site accumulation.
    _pool_f4_quads = f4_quads
    # Store the block size in base pairs for downstream block index math.
    _pool_block_size = block_size


def _process_block_lines(block_tuple, shared=None):
    """Process one block of raw VCF lines into a finalized block result and warnings.

    Args:
        block_tuple (tuple): A `(chrom, block_idx, block_lines)` tuple
            produced by `_produce_blocks`. `block_lines` is a list of raw
            VCF data line strings, none of which have been stripped, fully
            split, or filtered.
        shared (dict or None): Optional dictionary containing all shared
            configuration required for per-line accumulation. When `None`
            (the multi-threaded default), the function reads shared
            configuration from module-level `_pool_*` globals populated by
            `_init_pool_worker`. When a dictionary is supplied (the
            single-threaded path), the function reads from that dictionary
            instead and does not touch module-level globals.

    Returns:
        result (AlleleCountBlockResult or None): The finalized block result,
            or None if the block contained no segregating sites after the
            biallelic and all-missing filters were applied.
        warnings (list): List of `(code, chrom, context, message)` warning
            tuples produced while processing this block.

    Notes
    -----
    For each VCF data line this function performs the full tab split,
    applies the biallelic SNP filter, parses the FORMAT field to
    locate the GP sub-field (per-line, no caching), extracts per-population
    allele counts through the hard-call or genotype-probability path, skips
    all-missing sites, and updates the block accumulator with per-population
    sums, per-pair dxy and Fst sums, n_seg counts, site pattern counts, and
    f2, f3, and f4 sums. After all lines in the block are processed, the
    accumulator is handed to `_finalize_block` to produce an
    `AlleleCountBlockResult`, which may be `None` if no segregating sites
    survived the filters.
    """
    # Resolve shared configuration from either the explicit dictionary
    # (single-threaded path) or the module-level pool globals (multi-threaded
    # path).
    if shared is not None:
        pop_indices = shared["pop_indices"]
        gp_inds = shared["gp_inds"]
        all_pops = shared["all_pops"]
        all_pairs = shared["all_pairs"]
        sp_configs = shared["sp_configs"]
        n_seg_keys = shared["n_seg_keys"]
        fst_pairs = shared["fst_pairs"]
        f2_pairs = shared["f2_pairs"]
        f3_triples = shared["f3_triples"]
        f4_quads = shared["f4_quads"]
        block_size = shared["block_size"]
    else:
        pop_indices = _pool_pop_indices
        gp_inds = _pool_gp_inds
        all_pops = _pool_all_pops
        all_pairs = _pool_all_pairs
        sp_configs = _pool_sp_configs
        n_seg_keys = _pool_n_seg_keys
        fst_pairs = _pool_fst_pairs
        f2_pairs = _pool_f2_pairs
        f3_triples = _pool_f3_triples
        f4_quads = _pool_f4_quads
        block_size = _pool_block_size

    # Unpack the block tuple into its three components.
    chrom, block_idx, block_lines = block_tuple

    # Per-block warnings list, returned alongside the finalized block result.
    warnings = []

    # Fresh accumulator with zeroed sums for every tracked key in this block.
    acc = _init_block_accumulator(
        all_pops,
        all_pairs,
        sp_configs,
        n_seg_keys,
        fst_pairs,
        f2_pairs,
        f3_triples,
        f4_quads,
    )

    for line in block_lines:
        # Parse the tab-separated VCF data fields.
        fields = line.strip().split("\t")

        # Position is needed only for warning context strings in this block;
        # chrom is already known from the block tuple.
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]

        # Check for biallelic SNP: single-character REF and ALT, no commas,
        # neither allele missing.
        if "," in alt or len(ref) != 1 or len(alt) != 1 or ref == "." or alt == ".":
            warnings.append(
                (
                    "V001",
                    chrom,
                    f"{chrom}:{pos}",
                    f"Non-biallelic site skipped at {chrom}:{pos}",
                )
            )
            continue

        # Parse the FORMAT field to locate the GP sub-field position.
        format_parts = fields[8].split(":")
        gp_idx = format_parts.index("GP") if "GP" in format_parts else None

        # ----------------------------------------------------------
        # Extract per-population allele counts at this site.
        # ----------------------------------------------------------
        site_aac = {}
        site_c = {}
        for pop in all_pops:
            # Populations with no individuals in the VCF get zero counts.
            if pop not in pop_indices:
                site_aac[pop] = 0.0
                site_c[pop] = 0
                continue

            indices = pop_indices[pop]
            n_aac = 0.0
            c = 0

            if pop in gp_inds:
                # GP path: expected allele counts from posterior genotype
                # probabilities (eq:h-hat) from RIPTA-equations.
                if gp_idx is not None:
                    for idx in indices:
                        sample = fields[idx]
                        # Skip individuals with missing genotypes.
                        if sample[0] == ".":
                            continue
                        # Extract GP sub-field for this individual.
                        gp_str = sample.split(":")[gp_idx]
                        # If the GP sub-field itself is missing, skip.
                        if gp_str == ".":
                            continue
                        # Parse the three GP values: P(hom_ref), P(het),
                        # P(hom_alt).
                        gp_vals = gp_str.split(",")
                        # Expected AAC = P(het) + 2 * P(hom_alt).
                        n_aac += float(gp_vals[1]) + 2 * float(gp_vals[2])
                        # Diploid: c = 2 per called individual.
                        c += 2
                # If GP field is absent from FORMAT at this site,
                # all individuals in this population are treated as
                # missing (c remains 0).
            else:
                # Hard-call path: count alleles from the GT field.
                for idx in indices:
                    gt = fields[idx][0:3]
                    # Skip individuals with any missing allele in GT.
                    if "." in gt:
                        continue
                    # Count alternative alleles ('1') in the genotype.
                    n_aac += gt.count("1")
                    # Chromosome count: total called alleles at this site.
                    c += gt.count("0") + gt.count("1")

            site_aac[pop] = n_aac
            site_c[pop] = c

        # V002: skip sites where every population has no called genotypes.
        if all(site_c[pop] == 0 for pop in all_pops):
            warnings.append(
                (
                    "V002",
                    chrom,
                    f"{chrom}:{pos}",
                    f"All individuals missing at {chrom}:{pos}",
                )
            )
            continue

        # ----------------------------------------------------------
        # Accumulate per-population sums for this segregating site.
        # ----------------------------------------------------------
        # Pre-compute per-site allele frequencies for Fst and
        # f-statistic accumulation loops below.
        site_freq = {}
        for pop in all_pops:
            aac = site_aac[pop]
            c = site_c[pop]
            # Per-site allele frequency, zero when no genotypes are called.
            site_freq[pop] = aac / c if c > 0 else 0.0
            # Sum of alternative allele counts.
            acc["sum_aac"][pop] += aac
            # Sum of chromosome counts.
            acc["sum_c"][pop] += c
            # Heterozygosity numerator: n_aac * (c - n_aac) (eq:h-hat).
            het_num = aac * (c - aac)
            acc["sum_het_num"][pop] += het_num
            # Pi differences: same as het_num (eq:pi-diffs).
            acc["sum_pi_diffs"][pop] += het_num
            # Pi comparisons: c choose 2 = c * (c - 1) / 2 (eq:pi-comps).
            acc["sum_pi_comps"][pop] += c * (c - 1) / 2
            # Large-sample pi per-site contribution: 2 * p_m,x * (1 - p_m,x) (eq:pi-sample).
            # Guard on c > 0: missing populations have site_freq 0.0 which would
            # incorrectly contribute a zero term to the per-site average.
            if c > 0:
                acc["sum_pi_large"][pop] += 2 * site_freq[pop] * (1 - site_freq[pop])
            # Large-sample half-heterozygosity per site for f3star denominator:
            # p_m,x * (1 - p_m,x) (eq:h-large). Distinct from sum_pi_large which
            # holds 2 * p * (1 - p); kept separate for semantic clarity.
            if c > 0:
                acc["sum_h_large"][pop] += site_freq[pop] * (1 - site_freq[pop])
            # Unbiased h_hat ratio-of-sums numerator: aac*(c - aac) (eq:h-hat numerator).
            # Semantically distinct from sum_het_num which is used in f2/f3 bias correction.
            acc["sum_h_hat_num"][pop] += het_num
            # Unbiased h_hat ratio-of-sums denominator: c*(c - 1) (eq:h-hat denominator).
            # When c < 2, the contribution is zero, which is the correct value.
            if c >= 2:
                acc["sum_h_hat_denom"][pop] += c * (c - 1)

        # ----------------------------------------------------------
        # Accumulate per-pair sums for between-population differences.
        # ----------------------------------------------------------
        for pair in all_pairs:
            px, py = pair
            aac_x = site_aac[px]
            c_x = site_c[px]
            aac_y = site_aac[py]
            c_y = site_c[py]
            # Dxy differences: cross-population allelic mismatches
            # (eq:dxy-diffs) from RIPTA-equations.
            acc["sum_dxy_diffs"][pair] += aac_x * (c_y - aac_y) + aac_y * (c_x - aac_x)
            # Dxy comparisons: cross-population allelic pairs
            # (eq:dxy-comps) from RIPTA-equations.
            acc["sum_dxy_comps"][pair] += c_x * c_y
            # Large-sample dxy: p_m,x*(1-p_m,y) + p_m,y*(1-p_m,x) (eq:dxy-site-large).
            # Guard on both counts > 0 so missing populations do not
            # contribute a zero term to the per-site average.
            if c_x > 0 and c_y > 0:
                freq_x = site_freq[px]
                freq_y = site_freq[py]
                acc["sum_dxy_large"][pair] += freq_x * (1 - freq_y) + freq_y * (
                    1 - freq_x
                )

        # ----------------------------------------------------------
        # Update n_seg: count jointly callable sites per frozenset.
        # ----------------------------------------------------------
        for key in n_seg_keys:
            if all(site_c[pop] > 0 for pop in key):
                acc["n_seg"][key] += 1

        # ----------------------------------------------------------
        # Update site pattern accumulators for each SitePatternConfig.
        # ----------------------------------------------------------
        for sp in sp_configs:
            _compute_site_patterns(
                acc["site_patterns"][sp],
                sp,
                site_aac,
                site_c,
                warnings,
                chrom,
                pos,
            )

        # ----------------------------------------------------------
        # Accumulate per-pair Fst numerator and denominator sums.
        # ----------------------------------------------------------
        for pair in fst_pairs:
            px, py = pair
            c_x = site_c[px]
            c_y = site_c[py]
            # Skip sites where either population has no called genotypes.
            if c_x == 0 or c_y == 0:
                continue
            # Large-sample Fst numerator: (p_x - p_y)^2 (eq:fst-N-large).
            freq_diff_sq = (site_freq[px] - site_freq[py]) ** 2
            acc["sum_fst_num"][pair] += freq_diff_sq
            # Per-site heterozygosity numerators for bias correction.
            het_num_x = site_aac[px] * (c_x - site_aac[px])
            het_num_y = site_aac[py] * (c_y - site_aac[py])
            # Unbiased h_hat/c terms (eq:fst-D-hat). When c < 2,
            # het_num is zero so the correction vanishes.
            h_hat_x_over_c = het_num_x / (c_x**2 * (c_x - 1)) if c_x >= 2 else 0.0
            h_hat_y_over_c = het_num_y / (c_y**2 * (c_y - 1)) if c_y >= 2 else 0.0
            # Unbiased h_hat terms for the denominator (eq:fst-hat).
            h_hat_x = het_num_x / (c_x * (c_x - 1)) if c_x >= 2 else 0.0
            h_hat_y = het_num_y / (c_y * (c_y - 1)) if c_y >= 2 else 0.0
            # Unbiased Fst numerator: (p_x - p_y)^2 - h_hat_x/c_x
            # - h_hat_y/c_y (eq:fst-D-hat).
            n_hat = freq_diff_sq - h_hat_x_over_c - h_hat_y_over_c
            acc["sum_fst_num_unbiased"][pair] += n_hat
            # Unbiased Fst denominator: N_hat + h_hat_x + h_hat_y
            # (eq:fst-hat).
            acc["sum_fst_denom_unbiased"][pair] += n_hat + h_hat_x + h_hat_y
            # Large-sample Fst denominator D_m:
            # p_x*(1-p_y) + p_y*(1-p_x) (eq:fst-D-large). The enclosing
            # loop already skipped sites where either c is zero.
            acc["sum_fst_denom"][pair] += site_freq[px] * (
                1 - site_freq[py]
            ) + site_freq[py] * (1 - site_freq[px])

        # ----------------------------------------------------------
        # Accumulate per-pair f2 sums.
        # ----------------------------------------------------------
        for pair in f2_pairs:
            pa, pb = pair
            c_a = site_c[pa]
            c_b = site_c[pb]
            # Skip sites where either population has no called genotypes.
            if c_a == 0 or c_b == 0:
                continue
            # Large-sample f2: (p_a - p_b)^2 (eq:f2-large).
            freq_diff_sq = (site_freq[pa] - site_freq[pb]) ** 2
            acc["sum_f2"][pair] += freq_diff_sq
            # Bias correction terms: h_hat_a/c_a and h_hat_b/c_b
            # (eq:f2). When c < 2, het_num is zero so the
            # correction vanishes.
            het_num_a = site_aac[pa] * (c_a - site_aac[pa])
            het_num_b = site_aac[pb] * (c_b - site_aac[pb])
            h_hat_a_over_c = het_num_a / (c_a**2 * (c_a - 1)) if c_a >= 2 else 0.0
            h_hat_b_over_c = het_num_b / (c_b**2 * (c_b - 1)) if c_b >= 2 else 0.0
            # Unbiased f2: (p_a - p_b)^2 - h_hat_a/c_a - h_hat_b/c_b
            # (eq:f2).
            acc["sum_f2_unbiased"][pair] += (
                freq_diff_sq - h_hat_a_over_c - h_hat_b_over_c
            )

        # ----------------------------------------------------------
        # Accumulate per-triple f3 sums.
        # ----------------------------------------------------------
        for triple in f3_triples:
            pc, pa, pb = triple
            c_c = site_c[pc]
            # Skip sites where any population has no called genotypes.
            if c_c == 0 or site_c[pa] == 0 or site_c[pb] == 0:
                continue
            # Large-sample f3: (p_c - p_a)(p_c - p_b) (eq:f3-large).
            p_c = site_freq[pc]
            f3_site = (p_c - site_freq[pa]) * (p_c - site_freq[pb])
            acc["sum_f3"][triple] += f3_site
            # Bias correction: h_hat_c/c_c since P_C appears in both
            # factors (eq:f3). When c_c < 2, correction vanishes.
            het_num_c = site_aac[pc] * (c_c - site_aac[pc])
            h_hat_c_over_c = het_num_c / (c_c**2 * (c_c - 1)) if c_c >= 2 else 0.0
            # Unbiased f3: (p_c - p_a)(p_c - p_b) - h_hat_c/c_c
            # (eq:f3).
            acc["sum_f3_unbiased"][triple] += f3_site - h_hat_c_over_c

        # ----------------------------------------------------------
        # Accumulate per-quad f4 sums.
        # ----------------------------------------------------------
        for quad in f4_quads:
            pa, pb, pc, pd = quad
            # Skip sites where any population has no called genotypes.
            if site_c[pa] == 0 or site_c[pb] == 0 or site_c[pc] == 0 or site_c[pd] == 0:
                continue
            # f4: (p_a - p_b)(p_c - p_d) — no bias correction needed
            # (eq:f4).
            acc["sum_f4"][quad] += (site_freq[pa] - site_freq[pb]) * (
                site_freq[pc] - site_freq[pd]
            )

        # Increment the segregating site counter for this block.
        acc["site_count"] += 1

    # Finalize the block after all lines have been accumulated.
    result = _finalize_block(acc, chrom, block_idx, block_size)

    return result, warnings


def _write_vcf_log(warnings, log_path_prefix, verbose):
    """Append VCF processing warnings to the gzipped TSV log file.

    Args:
        warnings (list): List of (code, chrom, context, message) tuples.
        log_path_prefix (str): Prefix for the log file path.
        verbose (bool): Whether to include human-readable messages in the
            MESSAGE column.

    Notes
    -----
    The log file and header are created by config.py during the config phase.
    This function appends warning rows only, using 'at' mode.
    """
    # Construct the log file path.
    log_path = f"{log_path_prefix}_ripta.log.gz"

    # Append warning rows to the existing log file (header written by config phase).
    with gzip.open(log_path, "at") as f:
        for code, chrom, context, message in warnings:
            if verbose:
                # Verbose rows include the human-readable MESSAGE column.
                f.write(f"VCF\t{code}\t{chrom}\t{context}\t{message}\n")
            else:
                # Non-verbose rows omit the MESSAGE column entirely.
                f.write(f"VCF\t{code}\t{chrom}\t{context}\n")


def process_vcfs(config):
    """Process VCF files and produce allele count blocks for all statistics.

    Args:
        config (dict): Validated configuration dictionary from read_config,
            with populated configuration objects (sp_configs, pi_configs, etc.),
            idx_dicc, and gp_inds.

    Returns:
        blocks (list): List of AlleleCountBlockResult objects, one per
            non-empty genomic block across all chroms, sorted by chrom
            order and block start position.

    Notes
    -----
    Builds a list of `(vcf_path, target_chroms)` pairs that uniformly
    represents both the single-VCF and multi-VCF input modes, then dispatches
    on `number_of_threads`. The single-threaded path builds an explicit
    `shared` dictionary once and calls `_process_block_lines` directly for
    every block yielded by `_produce_blocks`. The multi-threaded path builds
    a chained generator across all VCF items and dispatches blocks through
    `mp.Pool.imap_unordered` with `chunksize=1`; the pool initializer
    populates `_pool_*` globals in each worker process. Both paths share the
    same post-processing: warning log write, chrom-order sort, and
    `drop_last_block` filtering.
    """
    # Collect all required populations, pairs, configurations, and key sets.
    (
        all_pops,
        all_pairs,
        sp_configs,
        n_seg_keys,
        fst_pairs,
        f2_pairs,
        f3_triples,
        f4_quads,
    ) = _collect_required_sets(config)

    # Build mapping from population labels to VCF column indices, restricted
    # to populations actually required by the requested statistics.
    pop_indices = _build_pop_indices(config, all_pops)

    # Extract commonly used configuration values into local names.
    block_size = config["block_size"]
    drop_last_block = config["drop_last_block"]
    gp_inds = config["gp_inds"]

    # Resolve number_of_threads: convert the string 'max' to the CPU count.
    n_threads = config["number_of_threads"]
    if n_threads == "max":
        n_threads = mp.cpu_count()

    # Build the uniform list of (vcf_path, target_chroms) items for both
    # input modes. In single-VCF mode there is one item with all chroms;
    # in multi-VCF mode there is one item per chrom mapped to its own file.
    if isinstance(config["vcf_path"], str):
        vcf_items = [(config["vcf_path"], set(config["chroms"].keys()))]
    else:
        vcf_items = [(path, {chrom}) for chrom, path in config["vcf_path"].items()]

    all_blocks = []
    all_warnings = []

    if n_threads == 1:
        # Single-threaded path: build the explicit shared dictionary once
        # and pass it to every per-block call. Module-level `_pool_*` globals
        # are never set or read on this path.
        shared = {
            "pop_indices": pop_indices,
            "gp_inds": gp_inds,
            "all_pops": all_pops,
            "all_pairs": all_pairs,
            "sp_configs": sp_configs,
            "n_seg_keys": n_seg_keys,
            "fst_pairs": fst_pairs,
            "f2_pairs": f2_pairs,
            "f3_triples": f3_triples,
            "f4_quads": f4_quads,
            "block_size": block_size,
        }
        # Iterate over each VCF file and each block yielded from it, calling
        # the per-block worker directly with the explicit shared dictionary.
        for vcf_path, target_chroms in vcf_items:
            for block_tuple in _produce_blocks(
                vcf_path,
                target_chroms,
                block_size,
            ):
                result, warnings = _process_block_lines(block_tuple, shared)
                # Append the finalized block unless it contained no segregating sites.
                if result is not None:
                    all_blocks.append(result)
                # Extend the genome-wide warning list with this block's warnings.
                all_warnings.extend(warnings)
    else:
        # Multi-threaded path: chain the block generators across all VCF
        # items so one `imap_unordered` covers both input modes.
        block_generator = chain.from_iterable(
            _produce_blocks(vcf_path, target_chroms, block_size)
            for vcf_path, target_chroms in vcf_items
        )
        # Pool initializer arguments correspond one-to-one with the
        # parameters of `_init_pool_worker`, in the same order.
        initargs = (
            pop_indices,
            gp_inds,
            all_pops,
            all_pairs,
            sp_configs,
            n_seg_keys,
            fst_pairs,
            f2_pairs,
            f3_triples,
            f4_quads,
            block_size,
        )
        # Pool context manager handles shutdown on both success and
        # exception paths.
        with mp.Pool(
            n_threads,
            initializer=_init_pool_worker,
            initargs=initargs,
        ) as pool:
            try:
                # chunksize=1 is hardcoded because each block's per-line
                # accumulation is heavy relative to dispatch overhead and I
                # have found that chunksize=1 balances the load across workers.
                for result, warnings in pool.imap_unordered(
                    _process_block_lines,
                    block_generator,
                    chunksize=1,
                ):
                    # Append the finalized block unless it contained no segregating sites.
                    if result is not None:
                        all_blocks.append(result)
                    # Extend the genome-wide warning list with this block's warnings.
                    all_warnings.extend(warnings)
            except Exception:
                # Kill remaining workers before re-raising so the exception
                # trace reaches the caller.
                pool.terminate()
                raise

    # Append VCF processing warnings to the log file (always called; the
    # config phase created the file with the header row).
    _write_vcf_log(
        all_warnings,
        config["log_path_prefix"],
        config["verbose"],
    )

    # Sort blocks by the user-specified chrom order, then by block start,
    # so downstream consumers see a deterministic order.
    chrom_order = {chrom: idx for idx, chrom in enumerate(config["chroms"])}
    all_blocks.sort(
        key=lambda b: (chrom_order.get(b.chrom, len(chrom_order)), b.block_start)
    )

    # Post-processing: drop the last block per chrom when it extends past the
    # chrom length and drop_last_block is enabled. Only the highest-start block
    # per chrom is a candidate; all other blocks are retained.
    if drop_last_block:
        # Identify the highest block_start per chrom across all collected blocks.
        last_block_start_by_chrom = {}
        for block in all_blocks:
            prior = last_block_start_by_chrom.get(block.chrom)
            if prior is None or block.block_start > prior:
                last_block_start_by_chrom[block.chrom] = block.block_start
        # Retain a block unless it is the last for its chrom and extends past
        # the user-specified chrom length.
        retained_blocks = []
        for block in all_blocks:
            is_last = block.block_start == last_block_start_by_chrom[block.chrom]
            chrom_length = config["chroms"].get(block.chrom)
            if is_last and chrom_length is not None and block.block_end > chrom_length:
                continue
            retained_blocks.append(block)
        all_blocks = retained_blocks

    return all_blocks
