from dataclasses import dataclass


@dataclass(frozen=True)
class SitePatternConfig:
    """Population config for introgression statistics (D, Danc, Dplus, fhom, Dp, df, Dstar)."""

    stat_family: str = "site_patterns"
    p1: str = ""
    p2: str = ""
    p3: str = ""
    outgroup: str = ""
    stats: tuple = ()
    use_gp_p1: bool = False
    use_gp_p2: bool = False
    use_gp_p3: bool = False
    use_gp_outgroup: bool = False


@dataclass(frozen=True)
class PairwiseConfig:
    """Population config for pairwise and distance-based statistics (pi, dxy, fst, D3)."""

    stat_family: str = "pairwise"
    stat: str = ""
    pop: str = ""
    pop_x: str = ""
    pop_y: str = ""
    p1: str = ""
    p2: str = ""
    p3: str = ""
    unbiased: bool = True
    use_gp: bool = False
    use_gp_x: bool = False
    use_gp_y: bool = False
    use_gp_p1: bool = False
    use_gp_p2: bool = False
    use_gp_p3: bool = False


@dataclass(frozen=True)
class FStatConfig:
    """Population config for f-statistics (f2, f3, f3star, f4, f4ratio)."""

    stat_family: str = "fstat"
    stat: str = ""
    pop_a: str = ""
    pop_b: str = ""
    pop_c: str = ""
    pop_d: str = ""
    pop_e: str = ""
    unbiased: bool = True
    use_gp_a: bool = False
    use_gp_b: bool = False
    use_gp_c: bool = False
    use_gp_d: bool = False
    use_gp_e: bool = False


@dataclass(frozen=True)
class AlleleCountBlock:
    """Unified block container for all statistic families.

    Per-population fields use population label strings as keys. Per-population-pair
    fields use (pop_x, pop_y) sorted tuples as keys. Site pattern fields use
    SitePatternConfig instances as keys. Joint callable site counts use frozensets
    of population label strings as keys. Fst and f2 fields use (pop_x, pop_y)
    sorted tuples as keys. f3 fields use (pop_c, pop_a, pop_b) ordered tuples.
    f4 fields use (pop_a, pop_b, pop_c, pop_d) ordered tuples. The large-sample
    fields sum_pi_large, sum_dxy_large, and sum_fst_denom follow the same key
    conventions as their respective existing fields: sum_pi_large is keyed by
    population label strings (like sum_aac), and sum_dxy_large and sum_fst_denom
    are keyed by (pop_x, pop_y) sorted tuples (like sum_dxy_diffs and sum_fst_num).
    """

    # Sum of alternative allele counts across segregating sites in the block.
    sum_aac: dict
    # Sum of chromosome counts across segregating sites in the block.
    sum_c: dict
    # Sum of n_aac * (c - n_aac) — the numerator of the unbiased half-heterozygosity (eq:h-hat).
    sum_het_num: dict
    # Sum of within-population pairwise differences delta_m,x (eq:pi-diffs). Numerically identical to sum_het_num; kept separate for semantic clarity.
    sum_pi_diffs: dict
    # Sum of within-population pairwise comparisons kappa_m,x (eq:pi-comps).
    sum_pi_comps: dict
    # Sum of between-population pairwise differences delta_m,xy (eq:dxy-diffs).
    sum_dxy_diffs: dict
    # Sum of between-population pairwise comparisons kappa_m,xy (eq:dxy-comps).
    sum_dxy_comps: dict
    # Site pattern counts keyed by SitePatternConfig instances.
    site_patterns: dict
    # Counts of jointly callable sites keyed by frozenset of population labels.
    n_seg: dict
    # Large-sample Fst numerator per-site sum: sum_m (p_x - p_y)^2 (eq:fst-N-large).
    sum_fst_num: dict
    # Unbiased Fst numerator: sum_m N_hat_m (eq:fst-N-hat).
    sum_fst_num_unbiased: dict
    # Unbiased Fst denominator: sum_m D_hat_m (eq:fst-D-hat).
    sum_fst_denom_unbiased: dict
    # Large-sample f2: sum_m (p_a - p_b)^2 (eq:f2-large).
    sum_f2: dict
    # Unbiased f2: sum_m ((p_a - p_b)^2 - h_hat_a/c_a - h_hat_b/c_b) (eq:f2).
    sum_f2_unbiased: dict
    # Large-sample f3: sum_m (p_c - p_a)(p_c - p_b) (eq:f3-large).
    sum_f3: dict
    # Unbiased f3: sum_m ((p_c - p_a)(p_c - p_b) - h_hat_c/c_c) (eq:f3).
    sum_f3_unbiased: dict
    # f4: sum_m (p_a - p_b)(p_c - p_d) (eq:f4). No bias correction.
    sum_f4: dict
    # Large-sample pi per-site sum: sum_m 2 * p_m,x * (1 - p_m,x) (eq:pi-sample).
    sum_pi_large: dict
    # Large-sample dxy sum: sum_m (p_m,x * (1 - p_m,y) + p_m,y * (1 - p_m,x)) (eq:dxy-site-large).
    sum_dxy_large: dict
    # Large-sample Fst denominator per-site sum: sum_m p_m,x*(1-p_m,y) + p_m,y*(1-p_m,x) (eq:fst-D-large).
    sum_fst_denom: dict
    # Per-site large-sample half-heterozygosity sum for f3star denominator:
    # sum_m p_m,x * (1 - p_m,x) (eq:h-large). Follows the same naming convention
    # as sum_pi_large and sum_dxy_large. Distinct from sum_pi_large which stores
    # 2 * p * (1 - p); kept separate for semantic clarity.
    sum_h_large: dict
    # Unbiased half-heterozygosity ratio-of-sums numerator for f3star denominator:
    # sum_m aac_m,x * (c_m,x - aac_m,x) (numerator of eq:h-hat summed over sites).
    # Semantically distinct from sum_het_num which is used in f2/f3 bias correction.
    sum_h_hat_num: dict
    # Unbiased half-heterozygosity ratio-of-sums denominator for f3star denominator:
    # sum_m c_m,x * (c_m,x - 1) (denominator of eq:h-hat summed over sites).
    sum_h_hat_denom: dict


@dataclass(frozen=True)
class AlleleCountBlockResult:
    """Wraps an AlleleCountBlock with its genomic coordinates."""

    # Chromosome name from the VCF CHROM column.
    chrom: str
    # 1-based inclusive start position of this block.
    block_start: int
    # 1-based inclusive end position of this block.
    block_end: int
    # The computed allele count sums for this block.
    block: AlleleCountBlock
