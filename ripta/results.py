"""Genome-wide and windowed statistic compilation for RIPTA."""

import gzip

import numpy as np

from ripta.jackknife import weighted_block_jackknife

# Ordered list of the eight raw site pattern count columns reported in both
# the header and the data rows of the site pattern CSV.
_SITE_PATTERN_RAW_COLUMNS = (
    "BBAA",
    "BABA",
    "ABBA",
    "BAAA",
    "ABAA",
    "AABA",
    "ABBA_HOM",
    "BABA_HOM",
)
_SITE_PATTERN_RAW_KEYS = tuple(name.lower() for name in _SITE_PATTERN_RAW_COLUMNS)


def _open_results_csv(prefix, name, gzip_output):
    """Open a result CSV file for writing, choosing plain text or gzip.

    Args:
        prefix (str): The results path prefix from the configuration dictionary.
        name (str): Statistic-specific filename suffix (e.g. 'pi', 'site_patterns').
        gzip_output (bool): When true, append '.gz' and use gzip text mode;
            when false, write plain text.

    Returns:
        handle: A file handle opened in text-write mode that the caller is
            responsible for closing (use as a context manager).
    """
    # Build the base path including the .csv extension.
    base_path = f"{prefix}_{name}.csv"
    if gzip_output:
        # Gzip mode: append .gz and open in gzip text-write mode.
        return gzip.open(base_path + ".gz", "wt")
    # Plain text mode: open in standard text-write mode.
    return open(base_path, "w")


def compile_results(config, blocks):
    """Aggregate block-level allele count sums and write results CSVs.

    Args:
        config (dict): Validated configuration dictionary with populated configuration
            objects (sp_configs, pi_configs, dxy_configs, fst_configs,
            d3_configs, f2_configs, f3_configs, f3star_configs, f4_configs,
            f4ratio_configs) and global parameters.
        blocks (list): Flat list of AlleleCountBlockResult objects as returned
            by process_vcfs.

    Notes
    -----
    The output is based on config['windows']:
    - False: genome-wide jackknife mode (one CSV per statistic family with
      point estimate, SE, Z, and block count columns).
    - True: windowed point estimate mode (one CSV per statistic family with
      CHROM, WIND_START, WIND_END, and point estimate columns only).
    """
    if config["windows"]:
        _compile_windows(config, blocks)
    else:
        _compile_genome_wide(config, blocks)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _pair_key(pop_a, pop_b):
    """Returns the alphabetically sorted tuple key for a population pair.

    Args:
        pop_a (str): First population label.
        pop_b (str): Second population label.

    Returns:
        key (tuple): Alphabetically sorted (pop_x, pop_y) tuple.
    """
    return tuple(sorted([pop_a, pop_b]))


def _safe_ratio(numer, denom):
    """Return numer / denom, or NaN when denom is zero.

    Args:
        numer (float): Numerator value.
        denom (float): Denominator value.

    Returns:
        ratio (float): The ratio, or NaN if denom is zero.
    """
    if denom == 0:
        return float("nan")
    return numer / denom


def _clamp_bounded(theta):
    """Return theta if within [-1, 1], otherwise NaN.

    Args:
        theta (float): Raw point estimate.

    Returns:
        clamped (float): theta if -1 <= theta <= 1, otherwise NaN.

    Notes
    -----
    Applied to bounded ratios (D, Danc, Dplus, fhom, Dp, df, Dstar, D3, Fst,
    f4ratio) to suppress implausible windowed point estimates that fall
    outside the natural range of the statistic. Leave-one-out estimates are
    never clamped to avoid biasing jackknife variance estimates.
    """
    # NaN input passes through unchanged to propagate missing-data signals.
    if theta != theta:
        return theta
    # Clamp to NaN when outside the bounded interval.
    if theta < -1 or theta > 1:
        return float("nan")
    return theta


def _large_sample_dxy(aac_x, c_x, aac_y, c_y):
    """Compute large-sample dxy from allele frequencies.

    Args:
        aac_x (float): Alternative allele counts for population X.
        c_x (float): Chromosome counts for population X.
        aac_y (float): Alternative allele counts for population Y.
        c_y (float): Chromosome counts for population Y.

    Returns:
        dxy (float): Large-sample dxy; NaN when either chromosome count is zero.

    Notes
    -----
    Implements dxy = p_x * (1 - p_y) + p_y * (1 - p_x) (eq:dxy-site-large) using
    allele frequencies p_x = aac_x / c_x, p_y = aac_y / c_y.
    """
    # Return NaN when either population has no called genotypes
    if c_x == 0 or c_y == 0:
        return float("nan")
    # Compute allele frequencies.
    freq_x = aac_x / c_x
    freq_y = aac_y / c_y
    return freq_x * (1 - freq_y) + freq_y * (1 - freq_x)


def _sp_stats(stat_name, counts):
    """Compute a site pattern statistic from site pattern counts.

    Args:
        stat_name (str): Statistic name (D, Danc, Dplus, fhom, Dp, df, Dstar).
        counts (dict): Mapping of pattern name strings to float counts.

    Returns:
        theta (float): The point estimate; NaN when denominator is zero.

    Notes
    -----
    For Dp, this returns the SIGNED ratio (abba - baba) / (abba + baba + bbaa).

    Implements:
        D      (eq:d-data):       (abba - baba) / (abba + baba)
        Danc   (eq:d-anc-data):  (baaa - abaa) / (baaa + abaa)
        Dplus  (eq:d-plus-data): ((abba - baba) + (baaa - abaa)) /
                                  (abba + baba + baaa + abaa)
        fhom   (eq:f-hom-data):  (abba - baba) / (abba_hom - baba_hom)
        Dp     (eq:d-p-data):    (abba - baba) / (abba + baba + bbaa)  [signed]
        df     (eq:d-f-data):    (abba - baba) / (abba + baba + 2*bbaa)
        Dstar  (eq:d-star-data): (abba - baba) / (bbaa + abba - 2*baba)
    """
    # Extract pattern counts with zero defaults for missing patterns.
    abba = counts.get("abba", 0.0)
    baba = counts.get("baba", 0.0)

    if stat_name == "D":
        # D (eq:d-data).
        return _safe_ratio(abba - baba, abba + baba)

    if stat_name == "Danc":
        # Danc (eq:d-anc-data).
        baaa = counts.get("baaa", 0.0)
        abaa = counts.get("abaa", 0.0)
        return _safe_ratio(baaa - abaa, baaa + abaa)

    if stat_name == "Dplus":
        # Dplus (eq:d-plus-data).
        baaa = counts.get("baaa", 0.0)
        abaa = counts.get("abaa", 0.0)
        return _safe_ratio(
            (abba - baba) + (baaa - abaa),
            abba + baba + baaa + abaa,
        )

    if stat_name == "fhom":
        # fhom (eq:f-hom-data).
        abba_hom = counts.get("abba_hom", 0.0)
        baba_hom = counts.get("baba_hom", 0.0)
        return _safe_ratio(abba - baba, abba_hom - baba_hom)

    if stat_name == "Dp":
        # Dp (eq:d-p-data).
        bbaa = counts.get("bbaa", 0.0)
        return _safe_ratio(abba - baba, abba + baba + bbaa)

    if stat_name == "df":
        # df (eq:d-f-data).
        bbaa = counts.get("bbaa", 0.0)
        return _safe_ratio(abba - baba, abba + baba + 2 * bbaa)

    if stat_name == "Dstar":
        # Dstar (eq:d-star-data).
        bbaa = counts.get("bbaa", 0.0)
        return _safe_ratio(abba - baba, bbaa + abba - 2 * baba)

    return float("nan")


def _aggregate_blocks(blocks):
    """Sum all AlleleCountBlock fields across a list of block results.

    Args:
        blocks (list): List of AlleleCountBlockResult objects to aggregate.

    Returns:
        gw (dict): Mapping of field name to dict of (key -> aggregated value).
            Site patterns are nested: 'site_patterns' -> {sp_config -> {pattern -> total}}.
    """
    gw = {}

    for result in blocks:
        block = result.block

        # Aggregate all flat dictionary fields (per-population, per-pair, per-triple,
        # per-quad, and n_seg).
        for field in (
            "sum_aac",
            "sum_c",
            "sum_het_num",
            "sum_pi_diffs",
            "sum_pi_comps",
            "sum_dxy_diffs",
            "sum_dxy_comps",
            "sum_fst_num",
            "sum_fst_num_unbiased",
            "sum_fst_denom_unbiased",
            "sum_f2",
            "sum_f2_unbiased",
            "sum_f3",
            "sum_f3_unbiased",
            "sum_f4",
            "sum_pi_large",
            "sum_dxy_large",
            "sum_fst_denom",
            "sum_h_large",
            "sum_h_hat_num",
            "sum_h_hat_denom",
            "n_seg",
        ):
            block_dict = getattr(block, field)
            if field not in gw:
                gw[field] = {}
            for key, val in block_dict.items():
                gw[field][key] = gw[field].get(key, 0) + val

        # Aggregate site pattern counts (nested dictionary: sp_config -> pattern -> count).
        if "site_patterns" not in gw:
            gw["site_patterns"] = {}
        for sp_config, patterns in block.site_patterns.items():
            if sp_config not in gw["site_patterns"]:
                gw["site_patterns"][sp_config] = {}
            for pattern, count in patterns.items():
                gw["site_patterns"][sp_config][pattern] = (
                    gw["site_patterns"][sp_config].get(pattern, 0.0) + count
                )

    return gw


# ---------------------------------------------------------------------------
# Genome-wide jackknife mode
# ---------------------------------------------------------------------------


def _compile_genome_wide(config, blocks):
    """Compute genome-wide statistics with jackknife SE and write CSVs.

    Args:
        config (dict): Validated configuration dictionary.
        blocks (list): Flat list of AlleleCountBlockResult objects.
    """
    # Aggregate all block fields into genome-wide totals.
    gw = _aggregate_blocks(blocks)
    prefix = config["results_path_prefix"]
    # Output format flag: when true, every CSV is written as .csv.gz.
    gzip_output = config["gzip_output"]

    # ==================================================================
    # Site pattern statistics (D, Danc, Dplus, fhom, Dp, df, Dstar)
    # ==================================================================
    sp_configs = config.get("sp_configs", [])
    if sp_configs:
        # Determine the union of all requested statistics across all configurations,
        # preserving first-encountered order for column layout.
        all_sp_stats = []
        seen_stats = set()
        for sp in sp_configs:
            for stat in sp.stats:
                if stat not in seen_stats:
                    all_sp_stats.append(stat)
                    seen_stats.add(stat)

        # Build CSV header: population labels, raw pattern counts, N_SNPS,
        # then per-statistic point estimate, SE, Z, and block count columns.
        header = ["P1", "P2", "P3", "OUTGROUP"]
        header.extend(_SITE_PATTERN_RAW_COLUMNS)
        header.append("N_SNPS")
        for stat in all_sp_stats:
            header.extend([stat, f"{stat}_SE", f"{stat}_Z", f"{stat}_BLOCKS"])

        rows = []
        for sp in sp_configs:
            # Frozenset key for n_seg lookup.
            n_seg_key = frozenset({sp.p1, sp.p2, sp.p3, sp.outgroup})
            gw_sp = gw.get("site_patterns", {}).get(sp, {})

            # Build leave-one-out count dictionaries and weight array.
            loo_list = []
            m_js_list = []
            for result in blocks:
                block_n_seg = result.block.n_seg.get(n_seg_key, 0)
                # Skip blocks with no jointly callable sites.
                if block_n_seg == 0:
                    continue
                block_sp = result.block.site_patterns.get(sp, {})
                # Subtract this block from the genome-wide total for each site pattern.
                loo = {
                    pattern: gw_sp.get(pattern, 0.0) - block_sp.get(pattern, 0.0)
                    for pattern in gw_sp
                }
                loo_list.append(loo)
                m_js_list.append(block_n_seg)

            m_js = np.array(m_js_list, dtype=float)

            # Genome-wide N_SNPS for this site pattern configuration.
            gw_n_seg = gw.get("n_seg", {}).get(n_seg_key, 0)

            # Compute each requested stat.
            row = [sp.p1, sp.p2, sp.p3, sp.outgroup]
            # Append the eight raw site pattern counts in column order.
            for pattern_key in _SITE_PATTERN_RAW_KEYS:
                row.append(str(gw_sp.get(pattern_key, 0.0)))
            # Append the segregating site count column.
            row.append(str(gw_n_seg))
            for stat in all_sp_stats:
                # Genome-wide point estimate.
                theta = _sp_stats(stat, gw_sp)
                # Clamp bounded ratios to NaN when outside [-1, 1].
                theta = _clamp_bounded(theta)

                # Leave-one-out estimates (signed for Dp).
                theta_js = np.array([_sp_stats(stat, loo) for loo in loo_list])

                # Jackknife estimates.
                g, sigma, z = weighted_block_jackknife(theta, theta_js, m_js)
                row.extend([str(theta), str(sigma), str(z), str(g)])

            rows.append(",".join(row))

        # Write CSV.
        with _open_results_csv(prefix, "site_patterns", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # D3
    # ==================================================================
    d3_configs = config.get("d3_configs", [])
    if d3_configs:
        header = [
            "P1",
            "P2",
            "P3",
            "FORM",
            "N_SNPS",
            "D3",
            "D3_SE",
            "D3_Z",
            "D3_BLOCKS",
        ]
        rows = []
        for cfg in d3_configs:
            n_seg_key = frozenset({cfg.p1, cfg.p2, cfg.p3})
            gw_n_seg = gw.get("n_seg", {}).get(n_seg_key, 0)
            # Sorted pair keys for dxy lookup.
            pair_13 = _pair_key(cfg.p1, cfg.p3)
            pair_23 = _pair_key(cfg.p2, cfg.p3)

            if cfg.unbiased:
                # Unbiased D3 uses pixy-style dxy (eq:d-3-data + eq:dxy-diffs / eq:dxy-comps).
                gw_diffs_13 = gw.get("sum_dxy_diffs", {}).get(pair_13, 0.0)
                gw_comps_13 = gw.get("sum_dxy_comps", {}).get(pair_13, 0.0)
                gw_diffs_23 = gw.get("sum_dxy_diffs", {}).get(pair_23, 0.0)
                gw_comps_23 = gw.get("sum_dxy_comps", {}).get(pair_23, 0.0)
                # Genome-wide pixy dxy for each pair.
                d13 = _safe_ratio(gw_diffs_13, gw_comps_13)
                d23 = _safe_ratio(gw_diffs_23, gw_comps_23)
                theta = _clamp_bounded(_safe_ratio(d13 - d23, d13 + d23))

                # Leave-one-out estimates.
                theta_js_list = []
                m_js_list = []
                for result in blocks:
                    block_n_seg = result.block.n_seg.get(n_seg_key, 0)
                    if block_n_seg == 0:
                        continue
                    # Subtract this block's dxy sums from genome-wide.
                    loo_diffs_13 = gw_diffs_13 - result.block.sum_dxy_diffs.get(
                        pair_13, 0.0
                    )
                    loo_comps_13 = gw_comps_13 - result.block.sum_dxy_comps.get(
                        pair_13, 0.0
                    )
                    loo_diffs_23 = gw_diffs_23 - result.block.sum_dxy_diffs.get(
                        pair_23, 0.0
                    )
                    loo_comps_23 = gw_comps_23 - result.block.sum_dxy_comps.get(
                        pair_23, 0.0
                    )
                    # Leave-one-out pixy dxy.
                    d13_j = _safe_ratio(loo_diffs_13, loo_comps_13)
                    d23_j = _safe_ratio(loo_diffs_23, loo_comps_23)
                    theta_js_list.append(_safe_ratio(d13_j - d23_j, d13_j + d23_j))
                    m_js_list.append(block_n_seg)
            else:
                # Large-sample D3 uses per-site-averaged dxy (eq:dxy-genomewide-large).
                gw_dxy_large_13 = gw.get("sum_dxy_large", {}).get(pair_13, 0.0)
                gw_dxy_large_23 = gw.get("sum_dxy_large", {}).get(pair_23, 0.0)
                # Per-site averages: divide by the jointly-callable site count for each pair.
                n_seg_13 = gw.get("n_seg", {}).get(frozenset({cfg.p1, cfg.p3}), 0)
                n_seg_23 = gw.get("n_seg", {}).get(frozenset({cfg.p2, cfg.p3}), 0)
                d13 = _safe_ratio(gw_dxy_large_13, n_seg_13)
                d23 = _safe_ratio(gw_dxy_large_23, n_seg_23)
                theta = _clamp_bounded(_safe_ratio(d13 - d23, d13 + d23))

                # Leave-one-out estimates.
                theta_js_list = []
                m_js_list = []
                for result in blocks:
                    block_n_seg = result.block.n_seg.get(n_seg_key, 0)
                    if block_n_seg == 0:
                        continue
                    # Subtract this block's per-site dxy sum and n_seg for each pair.
                    loo_dxy_13 = gw_dxy_large_13 - result.block.sum_dxy_large.get(
                        pair_13, 0.0
                    )
                    loo_n_13 = n_seg_13 - result.block.n_seg.get(
                        frozenset({cfg.p1, cfg.p3}), 0
                    )
                    loo_dxy_23 = gw_dxy_large_23 - result.block.sum_dxy_large.get(
                        pair_23, 0.0
                    )
                    loo_n_23 = n_seg_23 - result.block.n_seg.get(
                        frozenset({cfg.p2, cfg.p3}), 0
                    )
                    d13_j = _safe_ratio(loo_dxy_13, loo_n_13)
                    d23_j = _safe_ratio(loo_dxy_23, loo_n_23)
                    theta_js_list.append(_safe_ratio(d13_j - d23_j, d13_j + d23_j))
                    m_js_list.append(block_n_seg)

            # Jackknife estimates.
            g, sigma, z_score = weighted_block_jackknife(
                theta,
                np.array(theta_js_list),
                np.array(m_js_list, dtype=float),
            )
            form_label = "unbiased" if cfg.unbiased else "large_sample"
            rows.append(
                f"{cfg.p1},{cfg.p2},{cfg.p3},{form_label},{gw_n_seg},"
                f"{theta},{sigma},{z_score},{g}"
            )

        # Write CSV.
        with _open_results_csv(prefix, "d3", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # pi
    # ==================================================================
    pi_configs = config.get("pi_configs", [])
    if pi_configs:
        header = ["POP", "FORM", "N_SNPS", "pi", "pi_SE", "pi_BLOCKS"]
        rows = []
        for cfg in pi_configs:
            pop = cfg.pop
            n_seg_key = frozenset({pop})
            gw_n_seg = gw.get("n_seg", {}).get(n_seg_key, 0)

            if cfg.unbiased:
                # Pixy pi (eq:pi-diffs / eq:pi-comps).
                gw_diffs = gw.get("sum_pi_diffs", {}).get(pop, 0.0)
                gw_comps = gw.get("sum_pi_comps", {}).get(pop, 0.0)
                theta = _safe_ratio(gw_diffs, gw_comps)

                # Leave-one-out estimates.
                theta_js_list = []
                m_js_list = []
                for result in blocks:
                    block_n_seg = result.block.n_seg.get(n_seg_key, 0)
                    if block_n_seg == 0:
                        continue
                    loo_diffs = gw_diffs - result.block.sum_pi_diffs.get(pop, 0.0)
                    loo_comps = gw_comps - result.block.sum_pi_comps.get(pop, 0.0)
                    theta_js_list.append(_safe_ratio(loo_diffs, loo_comps))
                    m_js_list.append(block_n_seg)
            else:
                # Large-sample pi: (1/M) * sum_m 2*p_m*(1-p_m) = sum_pi_large / n_seg (eq:pi-sample).
                gw_pi_large = gw.get("sum_pi_large", {}).get(pop, 0.0)
                theta = _safe_ratio(gw_pi_large, gw_n_seg)

                # Leave-one-out estimates.
                theta_js_list = []
                m_js_list = []
                for result in blocks:
                    block_n_seg = result.block.n_seg.get(n_seg_key, 0)
                    if block_n_seg == 0:
                        continue
                    # Remove this block's per-site sum and segregating site count.
                    loo_pi_large = gw_pi_large - result.block.sum_pi_large.get(pop, 0.0)
                    loo_n_seg = gw_n_seg - block_n_seg
                    theta_js_list.append(_safe_ratio(loo_pi_large, loo_n_seg))
                    m_js_list.append(block_n_seg)

            # Jackknife estimates.
            g, sigma, _z = weighted_block_jackknife(
                theta,
                np.array(theta_js_list),
                np.array(m_js_list, dtype=float),
            )
            form_label = "unbiased" if cfg.unbiased else "large_sample"
            rows.append(f"{pop},{form_label},{gw_n_seg},{theta},{sigma},{g}")

        # Write CSV.
        with _open_results_csv(prefix, "pi", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # dxy
    # ==================================================================
    dxy_configs = config.get("dxy_configs", [])
    if dxy_configs:
        header = ["POP_X", "POP_Y", "FORM", "N_SNPS", "dxy", "dxy_SE", "dxy_BLOCKS"]
        rows = []
        for cfg in dxy_configs:
            n_seg_key = frozenset({cfg.pop_x, cfg.pop_y})
            pair = _pair_key(cfg.pop_x, cfg.pop_y)
            gw_n_seg = gw.get("n_seg", {}).get(n_seg_key, 0)

            if cfg.unbiased:
                # Pixy dxy (eq:dxy-diffs / eq:dxy-comps).
                gw_diffs = gw.get("sum_dxy_diffs", {}).get(pair, 0.0)
                gw_comps = gw.get("sum_dxy_comps", {}).get(pair, 0.0)
                theta = _safe_ratio(gw_diffs, gw_comps)

                # Leave-one-out estimates.
                theta_js_list = []
                m_js_list = []
                for result in blocks:
                    block_n_seg = result.block.n_seg.get(n_seg_key, 0)
                    if block_n_seg == 0:
                        continue
                    loo_diffs = gw_diffs - result.block.sum_dxy_diffs.get(pair, 0.0)
                    loo_comps = gw_comps - result.block.sum_dxy_comps.get(pair, 0.0)
                    theta_js_list.append(_safe_ratio(loo_diffs, loo_comps))
                    m_js_list.append(block_n_seg)
            else:
                # Large-sample dxy: (1/M) * sum_m (p_x*(1-p_y)+p_y*(1-p_x)) = sum_dxy_large / n_seg (eq:dxy-genomewide-large).
                gw_dxy_large = gw.get("sum_dxy_large", {}).get(pair, 0.0)
                theta = _safe_ratio(gw_dxy_large, gw_n_seg)

                # Leave-one-out estimates.
                theta_js_list = []
                m_js_list = []
                for result in blocks:
                    block_n_seg = result.block.n_seg.get(n_seg_key, 0)
                    if block_n_seg == 0:
                        continue
                    # Remove this block's per-site sum and segregating site count.
                    loo_dxy_large = gw_dxy_large - result.block.sum_dxy_large.get(
                        pair, 0.0
                    )
                    loo_n_seg = gw_n_seg - block_n_seg
                    theta_js_list.append(_safe_ratio(loo_dxy_large, loo_n_seg))
                    m_js_list.append(block_n_seg)

            # Jackknife estimates.
            g, sigma, _z = weighted_block_jackknife(
                theta,
                np.array(theta_js_list),
                np.array(m_js_list, dtype=float),
            )
            form_label = "unbiased" if cfg.unbiased else "large_sample"
            rows.append(
                f"{cfg.pop_x},{cfg.pop_y},{form_label},{gw_n_seg},"
                f"{theta},{sigma},{g}"
            )

        # Write CSV.
        with _open_results_csv(prefix, "dxy", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # fst
    # ==================================================================
    fst_configs = config.get("fst_configs", [])
    if fst_configs:
        header = ["POP_X", "POP_Y", "FORM", "N_SNPS", "fst", "fst_SE", "fst_BLOCKS"]
        rows = []
        for cfg in fst_configs:
            n_seg_key = frozenset({cfg.pop_x, cfg.pop_y})
            pair = _pair_key(cfg.pop_x, cfg.pop_y)
            gw_n_seg = gw.get("n_seg", {}).get(n_seg_key, 0)

            if cfg.unbiased:
                # Unbiased Hudson Fst (eq:fst-N-hat / eq:fst-D-hat / eq:fst-hat).
                gw_num = gw.get("sum_fst_num_unbiased", {}).get(pair, 0.0)
                gw_den = gw.get("sum_fst_denom_unbiased", {}).get(pair, 0.0)
                theta = _clamp_bounded(_safe_ratio(gw_num, gw_den))

                # Leave-one-out estimates.
                theta_js_list = []
                m_js_list = []
                for result in blocks:
                    block_n_seg = result.block.n_seg.get(n_seg_key, 0)
                    if block_n_seg == 0:
                        continue
                    loo_num = gw_num - result.block.sum_fst_num_unbiased.get(pair, 0.0)
                    loo_den = gw_den - result.block.sum_fst_denom_unbiased.get(
                        pair, 0.0
                    )
                    theta_js_list.append(_safe_ratio(loo_num, loo_den))
                    m_js_list.append(block_n_seg)
            else:
                # Large-sample Fst: sum_fst_num / sum_fst_denom (eq:fst-N-large / eq:fst-D-large / eq:fst-large).
                gw_fst_num = gw.get("sum_fst_num", {}).get(pair, 0.0)
                gw_fst_denom = gw.get("sum_fst_denom", {}).get(pair, 0.0)
                theta = _clamp_bounded(_safe_ratio(gw_fst_num, gw_fst_denom))

                # Leave-one-out estimates.
                theta_js_list = []
                m_js_list = []
                for result in blocks:
                    block_n_seg = result.block.n_seg.get(n_seg_key, 0)
                    if block_n_seg == 0:
                        continue
                    # Remove this block's per-site numerator and denominator contributions.
                    loo_fst_num = gw_fst_num - result.block.sum_fst_num.get(pair, 0.0)
                    loo_fst_denom = gw_fst_denom - result.block.sum_fst_denom.get(
                        pair, 0.0
                    )
                    theta_js_list.append(_safe_ratio(loo_fst_num, loo_fst_denom))
                    m_js_list.append(block_n_seg)

            # Jackknife estimates.
            g, sigma, _z = weighted_block_jackknife(
                theta,
                np.array(theta_js_list),
                np.array(m_js_list, dtype=float),
            )
            form_label = "unbiased" if cfg.unbiased else "large_sample"
            rows.append(
                f"{cfg.pop_x},{cfg.pop_y},{form_label},{gw_n_seg},"
                f"{theta},{sigma},{g}"
            )

        # Write CSV.
        with _open_results_csv(prefix, "fst", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # f2
    # ==================================================================
    f2_configs = config.get("f2_configs", [])
    if f2_configs:
        header = ["POP_A", "POP_B", "FORM", "N_SNPS", "f2", "f2_SE", "f2_BLOCKS"]
        rows = []
        for cfg in f2_configs:
            n_seg_key = frozenset({cfg.pop_a, cfg.pop_b})
            pair = _pair_key(cfg.pop_a, cfg.pop_b)
            gw_n_seg = gw.get("n_seg", {}).get(n_seg_key, 0)

            # Select large-sample or unbiased sum field.
            field = "sum_f2_unbiased" if cfg.unbiased else "sum_f2"
            gw_sum = gw.get(field, {}).get(pair, 0.0)
            # f2 large-sample (eq:f2-large) or f2 unbiased (eq:f2).
            theta = _safe_ratio(gw_sum, gw_n_seg)

            # Leave-one-out estimates.
            theta_js_list = []
            m_js_list = []
            for result in blocks:
                block_n_seg = result.block.n_seg.get(n_seg_key, 0)
                if block_n_seg == 0:
                    continue
                loo_sum = gw_sum - getattr(result.block, field).get(pair, 0.0)
                loo_n = gw_n_seg - block_n_seg
                theta_js_list.append(_safe_ratio(loo_sum, loo_n))
                m_js_list.append(block_n_seg)

            # Jackknife estimates.
            g, sigma, _z = weighted_block_jackknife(
                theta,
                np.array(theta_js_list),
                np.array(m_js_list, dtype=float),
            )
            form_label = "unbiased" if cfg.unbiased else "large_sample"
            rows.append(
                f"{cfg.pop_a},{cfg.pop_b},{form_label},{gw_n_seg},"
                f"{theta},{sigma},{g}"
            )

        # Write CSV.
        with _open_results_csv(prefix, "f2", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # f3
    # ==================================================================
    f3_configs = config.get("f3_configs", [])
    if f3_configs:
        header = [
            "POP_C",
            "POP_A",
            "POP_B",
            "FORM",
            "N_SNPS",
            "f3",
            "f3_SE",
            "f3_Z",
            "f3_BLOCKS",
        ]
        rows = []
        for cfg in f3_configs:
            n_seg_key = frozenset({cfg.pop_a, cfg.pop_b, cfg.pop_c})
            triple = (cfg.pop_c, cfg.pop_a, cfg.pop_b)
            gw_n_seg = gw.get("n_seg", {}).get(n_seg_key, 0)

            # Select large-sample or unbiased sum field.
            field = "sum_f3_unbiased" if cfg.unbiased else "sum_f3"
            gw_sum = gw.get(field, {}).get(triple, 0.0)
            # f3 large-sample (eq:f3-large) or f3 unbiased (eq:f3).
            theta = _safe_ratio(gw_sum, gw_n_seg)

            # Leave-one-out estimates.
            theta_js_list = []
            m_js_list = []
            for result in blocks:
                block_n_seg = result.block.n_seg.get(n_seg_key, 0)
                if block_n_seg == 0:
                    continue
                loo_sum = gw_sum - getattr(result.block, field).get(triple, 0.0)
                loo_n = gw_n_seg - block_n_seg
                theta_js_list.append(_safe_ratio(loo_sum, loo_n))
                m_js_list.append(block_n_seg)

            # Jackknife estimates.
            g, sigma, z_score = weighted_block_jackknife(
                theta,
                np.array(theta_js_list),
                np.array(m_js_list, dtype=float),
            )
            form_label = "unbiased" if cfg.unbiased else "large_sample"
            rows.append(
                f"{cfg.pop_c},{cfg.pop_a},{cfg.pop_b},{form_label},{gw_n_seg},"
                f"{theta},{sigma},{z_score},{g}"
            )

        # Write CSV.
        with _open_results_csv(prefix, "f3", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # f3star
    # ==================================================================
    f3star_configs = config.get("f3star_configs", [])
    if f3star_configs:
        header = [
            "POP_C",
            "POP_A",
            "POP_B",
            "FORM",
            "N_SNPS",
            "f3star",
            "f3star_SE",
            "f3star_Z",
            "f3star_BLOCKS",
        ]
        rows = []
        for cfg in f3star_configs:
            n_seg_key = frozenset({cfg.pop_a, cfg.pop_b, cfg.pop_c})
            triple = (cfg.pop_c, cfg.pop_a, cfg.pop_b)
            gw_n_seg = gw.get("n_seg", {}).get(n_seg_key, 0)

            if cfg.unbiased:
                # Unbiased f3star (eq:f3-star): f3_hat / (2 * h_hat_pc).
                gw_f3 = gw.get("sum_f3_unbiased", {}).get(triple, 0.0)
                gw_h_hat_num = gw.get("sum_h_hat_num", {}).get(cfg.pop_c, 0.0)
                gw_h_hat_denom = gw.get("sum_h_hat_denom", {}).get(cfg.pop_c, 0.0)
                f3_val = _safe_ratio(gw_f3, gw_n_seg)
                # Ratio-of-sums h_hat_PC (eq:h-hat genome-wide).
                h_hat = _safe_ratio(gw_h_hat_num, gw_h_hat_denom)
                # f3-star (eq:f3-star).
                theta = _safe_ratio(f3_val, 2 * h_hat)

                # Leave-one-out estimates.
                theta_js_list = []
                m_js_list = []
                for result in blocks:
                    block_n_seg = result.block.n_seg.get(n_seg_key, 0)
                    if block_n_seg == 0:
                        continue
                    loo_f3 = gw_f3 - result.block.sum_f3_unbiased.get(triple, 0.0)
                    loo_n = gw_n_seg - block_n_seg
                    loo_h_num = gw_h_hat_num - result.block.sum_h_hat_num.get(
                        cfg.pop_c, 0.0
                    )
                    loo_h_denom = gw_h_hat_denom - result.block.sum_h_hat_denom.get(
                        cfg.pop_c, 0.0
                    )
                    f3_j = _safe_ratio(loo_f3, loo_n)
                    h_j = _safe_ratio(loo_h_num, loo_h_denom)
                    theta_js_list.append(_safe_ratio(f3_j, 2 * h_j))
                    m_js_list.append(block_n_seg)
            else:
                # Large-sample f3star: f3 / (2 * h_pc) where h_pc is per-site average.
                gw_f3 = gw.get("sum_f3", {}).get(triple, 0.0)
                gw_h_large = gw.get("sum_h_large", {}).get(cfg.pop_c, 0.0)
                f3_val = _safe_ratio(gw_f3, gw_n_seg)
                # Per-site average large-sample half-heterozygosity of P_C (eq:h-large).
                h_pc = _safe_ratio(gw_h_large, gw_n_seg)
                theta = _safe_ratio(f3_val, 2 * h_pc)

                # Leave-one-out estimates.
                theta_js_list = []
                m_js_list = []
                for result in blocks:
                    block_n_seg = result.block.n_seg.get(n_seg_key, 0)
                    if block_n_seg == 0:
                        continue
                    loo_f3 = gw_f3 - result.block.sum_f3.get(triple, 0.0)
                    loo_n = gw_n_seg - block_n_seg
                    loo_h = gw_h_large - result.block.sum_h_large.get(cfg.pop_c, 0.0)
                    f3_j = _safe_ratio(loo_f3, loo_n)
                    h_j = _safe_ratio(loo_h, loo_n)
                    theta_js_list.append(_safe_ratio(f3_j, 2 * h_j))
                    m_js_list.append(block_n_seg)

            # Jackknife estimates.
            g, sigma, z_score = weighted_block_jackknife(
                theta,
                np.array(theta_js_list),
                np.array(m_js_list, dtype=float),
            )
            form_label = "unbiased" if cfg.unbiased else "large_sample"
            rows.append(
                f"{cfg.pop_c},{cfg.pop_a},{cfg.pop_b},{form_label},{gw_n_seg},"
                f"{theta},{sigma},{z_score},{g}"
            )

        # Write CSV.
        with _open_results_csv(prefix, "f3star", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # f4
    # ==================================================================
    f4_configs = config.get("f4_configs", [])
    if f4_configs:
        header = [
            "POP_A",
            "POP_B",
            "POP_C",
            "POP_D",
            "N_SNPS",
            "f4",
            "f4_SE",
            "f4_Z",
            "f4_BLOCKS",
        ]
        rows = []
        for cfg in f4_configs:
            pops = {cfg.pop_a, cfg.pop_b, cfg.pop_c, cfg.pop_d}
            n_seg_key = frozenset(pops)
            quad = (cfg.pop_a, cfg.pop_b, cfg.pop_c, cfg.pop_d)
            gw_n_seg = gw.get("n_seg", {}).get(n_seg_key, 0)
            gw_sum = gw.get("sum_f4", {}).get(quad, 0.0)
            # f4 = (1/M) sum_m (p_a - p_b)(p_c - p_d) (eq:f4).
            theta = _safe_ratio(gw_sum, gw_n_seg)

            # Leave-one-out estimates.
            theta_js_list = []
            m_js_list = []
            for result in blocks:
                block_n_seg = result.block.n_seg.get(n_seg_key, 0)
                if block_n_seg == 0:
                    continue
                loo_sum = gw_sum - result.block.sum_f4.get(quad, 0.0)
                loo_n = gw_n_seg - block_n_seg
                theta_js_list.append(_safe_ratio(loo_sum, loo_n))
                m_js_list.append(block_n_seg)

            # Jackknife estimates.
            g, sigma, z_score = weighted_block_jackknife(
                theta,
                np.array(theta_js_list),
                np.array(m_js_list, dtype=float),
            )
            rows.append(
                f"{cfg.pop_a},{cfg.pop_b},{cfg.pop_c},{cfg.pop_d},"
                f"{gw_n_seg},{theta},{sigma},{z_score},{g}"
            )

        # Write CSV.
        with _open_results_csv(prefix, "f4", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # f4ratio
    # ==================================================================
    f4ratio_configs = config.get("f4ratio_configs", [])
    if f4ratio_configs:
        header = [
            "POP_A",
            "POP_B",
            "POP_C",
            "POP_D",
            "POP_E",
            "N_SNPS",
            "f4ratio",
            "f4ratio_SE",
            "f4ratio_Z",
            "f4ratio_BLOCKS",
        ]
        rows = []
        for cfg in f4ratio_configs:
            # Numerator f4: f4(pop_a, pop_d; pop_e, pop_c).
            num_quad = (cfg.pop_a, cfg.pop_d, cfg.pop_e, cfg.pop_c)
            # Denominator f4: f4(pop_a, pop_d; pop_b, pop_c).
            den_quad = (cfg.pop_a, cfg.pop_d, cfg.pop_b, cfg.pop_c)
            # Genome-wide raw sums; 1/M cancels in the ratio (eq:f4-ratio).
            gw_num = gw.get("sum_f4", {}).get(num_quad, 0.0)
            gw_den = gw.get("sum_f4", {}).get(den_quad, 0.0)
            # Point estimate of admixture proportion alpha.
            theta = _clamp_bounded(_safe_ratio(gw_num, gw_den))

            # Jackknife weight key: all five populations.
            n_seg_key = frozenset(
                {cfg.pop_a, cfg.pop_b, cfg.pop_c, cfg.pop_d, cfg.pop_e}
            )
            # Genome-wide N_SNPS for the five-population key.
            gw_n_seg = gw.get("n_seg", {}).get(n_seg_key, 0)
            # Leave-one-out estimates.
            theta_js_list = []
            m_js_list = []
            for result in blocks:
                block_n_seg = result.block.n_seg.get(n_seg_key, 0)
                if block_n_seg == 0:
                    continue
                # Leave-one-out: subtract this block from numerator and denominator independently.
                loo_num = gw_num - result.block.sum_f4.get(num_quad, 0.0)
                loo_den = gw_den - result.block.sum_f4.get(den_quad, 0.0)
                theta_js_list.append(_safe_ratio(loo_num, loo_den))
                m_js_list.append(block_n_seg)

            # Jackknife estimates.
            g, sigma, z_score = weighted_block_jackknife(
                theta,
                np.array(theta_js_list),
                np.array(m_js_list, dtype=float),
            )
            rows.append(
                f"{cfg.pop_a},{cfg.pop_b},{cfg.pop_c},{cfg.pop_d},"
                f"{cfg.pop_e},{gw_n_seg},"
                f"{theta},{sigma},{z_score},{g}"
            )

        # Write CSV.
        with _open_results_csv(prefix, "f4ratio", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")


# ---------------------------------------------------------------------------
# Windowed point-estimate mode
# ---------------------------------------------------------------------------


def _compile_windows(config, blocks):
    """Compute windowed point estimates and write CSVs.

    Args:
        config (dict): Validated configuration dictionary.
        blocks (list): Flat list of AlleleCountBlockResult objects.

    Notes
    -----
    Windows are non-overlapping base-pair ranges of size config['block_size']
    that reset per chrom. Each block is assigned to a window via:
        window_idx = (block_start - 1) // block_size
    Blocks in the same (chrom, window_idx) group are aggregated by summing
    their AlleleCountBlock fields. Point estimates are computed per window
    with no jackknife. Windows where the relevant n_seg sums to 0 are dropped
    and logged.
    """
    block_size = config["block_size"]
    prefix = config["results_path_prefix"]
    # Output format flag: when true, every CSV is written as .csv.gz.
    gzip_output = config["gzip_output"]
    log_path = f"{config['log_path_prefix']}_ripta.log.gz"

    # Sort order for chroms.
    chrom_order = {chrom: idx for idx, chrom in enumerate(config["chroms"])}

    # Group blocks by (chrom, window_idx) and aggregate.
    window_data = {}
    for result in blocks:
        window_idx = (result.block_start - 1) // block_size
        key = (result.chrom, window_idx)
        if key not in window_data:
            window_data[key] = []
        window_data[key].append(result)

    # Aggregate each window group.
    window_agg = {}
    for key, window_blocks in window_data.items():
        window_agg[key] = _aggregate_blocks(window_blocks)

    # Sorted window keys by chrom order then window index.
    sorted_keys = sorted(
        window_agg.keys(),
        key=lambda k: (chrom_order.get(k[0], len(chrom_order)), k[1]),
    )

    # Collect dropped window warnings for the log.
    dropped_warnings = []

    # ==================================================================
    # Site pattern statistics
    # ==================================================================
    sp_configs = config.get("sp_configs", [])
    if sp_configs:
        # Determine the union of all requested statistics across all configurations,
        # preserving first-encountered order for column layout.
        all_sp_stats = []
        seen_stats = set()
        for sp in sp_configs:
            for stat in sp.stats:
                if stat not in seen_stats:
                    all_sp_stats.append(stat)
                    seen_stats.add(stat)

        # Header order: window coordinates, N_SNPS, population labels, raw
        # site pattern counts, then per-statistic point estimate columns.
        header = [
            "CHROM",
            "WIND_START",
            "WIND_END",
            "N_SNPS",
            "P1",
            "P2",
            "P3",
            "OUTGROUP",
        ]
        header.extend(_SITE_PATTERN_RAW_COLUMNS)
        header.extend(all_sp_stats)

        rows = []
        # Iterate over every (chromosome, window) pair in the configured order.
        for chrom, window_idx in sorted_keys:
            # Aggregated allele count sums for this window.
            window_aggregate = window_agg[(chrom, window_idx)]
            # 1-based inclusive start and end positions of this window.
            wind_start = window_idx * block_size + 1
            wind_end = (window_idx + 1) * block_size
            for sp in sp_configs:
                # Frozenset key for jointly callable site count lookup.
                n_seg_key = frozenset({sp.p1, sp.p2, sp.p3, sp.outgroup})
                window_n_seg = window_aggregate.get("n_seg", {}).get(n_seg_key, 0)
                # Drop windows where no segregating site is callable for this configuration.
                if window_n_seg == 0:
                    dropped_warnings.append(
                        (
                            "R001",
                            chrom,
                            f"{chrom}:{wind_start}-{wind_end}",
                            f"Window dropped: n_seg=0 for site_patterns "
                            f"({sp.p1},{sp.p2},{sp.p3},{sp.outgroup})",
                        )
                    )
                    continue
                # Site pattern counts aggregated over this window.
                window_sp = window_aggregate.get("site_patterns", {}).get(sp, {})
                # Build the row in the column order.
                row = [
                    chrom,
                    str(wind_start),
                    str(wind_end),
                    str(window_n_seg),
                    sp.p1,
                    sp.p2,
                    sp.p3,
                    sp.outgroup,
                ]
                # Append the eight raw site pattern counts in column order.
                for pattern_key in _SITE_PATTERN_RAW_KEYS:
                    row.append(str(window_sp.get(pattern_key, 0.0)))
                # Compute and append each requested point estimate.
                for stat in all_sp_stats:
                    theta = _sp_stats(stat, window_sp)
                    # Clamp bounded ratios to NaN when outside [-1, 1].
                    theta = _clamp_bounded(theta)
                    row.append(str(theta))
                rows.append(",".join(row))

        # Write CSV.
        with _open_results_csv(prefix, "site_patterns", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # D3
    # ==================================================================
    d3_configs = config.get("d3_configs", [])
    if d3_configs:
        header = [
            "CHROM",
            "WIND_START",
            "WIND_END",
            "N_SNPS",
            "P1",
            "P2",
            "P3",
            "FORM",
            "D3",
        ]
        rows = []
        for chrom, window_idx in sorted_keys:
            # Aggregated allele count sums for this window.
            window_aggregate = window_agg[(chrom, window_idx)]
            # 1-based inclusive start and end positions of this window.
            wind_start = window_idx * block_size + 1
            wind_end = (window_idx + 1) * block_size
            for cfg in d3_configs:
                # Frozenset key for jointly callable site count lookup.
                n_seg_key = frozenset({cfg.p1, cfg.p2, cfg.p3})
                window_n_seg = window_aggregate.get("n_seg", {}).get(n_seg_key, 0)
                # Drop windows with no jointly callable sites for the three populations.
                if window_n_seg == 0:
                    dropped_warnings.append(
                        (
                            "R001",
                            chrom,
                            f"{chrom}:{wind_start}-{wind_end}",
                            f"Window dropped: n_seg=0 for D3 ({cfg.p1},{cfg.p2},{cfg.p3})",
                        )
                    )
                    continue
                # Sorted pair keys for dxy lookups (p1-p3 and p2-p3).
                pair_13 = _pair_key(cfg.p1, cfg.p3)
                pair_23 = _pair_key(cfg.p2, cfg.p3)
                if cfg.unbiased:
                    # Unbiased D3 uses pixy-style per-window dxy ratios.
                    d13 = _safe_ratio(
                        window_aggregate.get("sum_dxy_diffs", {}).get(pair_13, 0.0),
                        window_aggregate.get("sum_dxy_comps", {}).get(pair_13, 0.0),
                    )
                    d23 = _safe_ratio(
                        window_aggregate.get("sum_dxy_diffs", {}).get(pair_23, 0.0),
                        window_aggregate.get("sum_dxy_comps", {}).get(pair_23, 0.0),
                    )
                else:
                    # Large-sample windowed D3 uses per-site-averaged dxy within this window.
                    window_dxy_large_13 = window_aggregate.get("sum_dxy_large", {}).get(
                        pair_13, 0.0
                    )
                    window_dxy_large_23 = window_aggregate.get("sum_dxy_large", {}).get(
                        pair_23, 0.0
                    )
                    window_n_seg_13 = window_aggregate.get("n_seg", {}).get(
                        frozenset({cfg.p1, cfg.p3}), 0
                    )
                    window_n_seg_23 = window_aggregate.get("n_seg", {}).get(
                        frozenset({cfg.p2, cfg.p3}), 0
                    )
                    d13 = _safe_ratio(window_dxy_large_13, window_n_seg_13)
                    d23 = _safe_ratio(window_dxy_large_23, window_n_seg_23)
                # Windowed D3 = (d13 - d23) / (d13 + d23).
                theta = _clamp_bounded(_safe_ratio(d13 - d23, d13 + d23))
                form_label = "unbiased" if cfg.unbiased else "large_sample"
                rows.append(
                    f"{chrom},{wind_start},{wind_end},{window_n_seg},"
                    f"{cfg.p1},{cfg.p2},{cfg.p3},{form_label},{theta}"
                )

        # Write CSV.
        with _open_results_csv(prefix, "d3", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # pi
    # ==================================================================
    pi_configs = config.get("pi_configs", [])
    if pi_configs:
        header = ["CHROM", "WIND_START", "WIND_END", "N_SNPS", "POP", "FORM", "pi"]
        rows = []
        for chrom, window_idx in sorted_keys:
            # Aggregated allele count sums for this window.
            window_aggregate = window_agg[(chrom, window_idx)]
            # 1-based inclusive start and end positions of this window.
            wind_start = window_idx * block_size + 1
            wind_end = (window_idx + 1) * block_size
            for cfg in pi_configs:
                # Frozenset key for the single-population n_seg lookup.
                n_seg_key = frozenset({cfg.pop})
                window_n_seg = window_aggregate.get("n_seg", {}).get(n_seg_key, 0)
                # Drop windows with no callable sites for this population.
                if window_n_seg == 0:
                    dropped_warnings.append(
                        (
                            "R001",
                            chrom,
                            f"{chrom}:{wind_start}-{wind_end}",
                            f"Window dropped: n_seg=0 for pi ({cfg.pop})",
                        )
                    )
                    continue
                if cfg.unbiased:
                    # Pixy pi: sum of pairwise differences over comparisons.
                    theta = _safe_ratio(
                        window_aggregate.get("sum_pi_diffs", {}).get(cfg.pop, 0.0),
                        window_aggregate.get("sum_pi_comps", {}).get(cfg.pop, 0.0),
                    )
                else:
                    # Large-sample pi: per-site average of 2*p*(1-p) (eq:pi-sample).
                    window_pi_large = window_aggregate.get("sum_pi_large", {}).get(
                        cfg.pop, 0.0
                    )
                    theta = _safe_ratio(window_pi_large, window_n_seg)
                form_label = "unbiased" if cfg.unbiased else "large_sample"
                rows.append(
                    f"{chrom},{wind_start},{wind_end},{window_n_seg},"
                    f"{cfg.pop},{form_label},{theta}"
                )

        # Write CSV.
        with _open_results_csv(prefix, "pi", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # dxy
    # ==================================================================
    dxy_configs = config.get("dxy_configs", [])
    if dxy_configs:
        header = [
            "CHROM",
            "WIND_START",
            "WIND_END",
            "N_SNPS",
            "POP_X",
            "POP_Y",
            "FORM",
            "dxy",
        ]
        rows = []
        for chrom, window_idx in sorted_keys:
            # Aggregated allele count sums for this window.
            window_aggregate = window_agg[(chrom, window_idx)]
            # 1-based inclusive start and end positions of this window.
            wind_start = window_idx * block_size + 1
            wind_end = (window_idx + 1) * block_size
            for cfg in dxy_configs:
                # Frozenset key for jointly callable site count lookup.
                n_seg_key = frozenset({cfg.pop_x, cfg.pop_y})
                window_n_seg = window_aggregate.get("n_seg", {}).get(n_seg_key, 0)
                # Drop windows with no jointly callable sites for the pair.
                if window_n_seg == 0:
                    dropped_warnings.append(
                        (
                            "R001",
                            chrom,
                            f"{chrom}:{wind_start}-{wind_end}",
                            f"Window dropped: n_seg=0 for dxy ({cfg.pop_x},{cfg.pop_y})",
                        )
                    )
                    continue
                # Sorted pair key for dxy lookups.
                pair = _pair_key(cfg.pop_x, cfg.pop_y)
                if cfg.unbiased:
                    # Pixy dxy: cross-population mismatches over comparisons.
                    theta = _safe_ratio(
                        window_aggregate.get("sum_dxy_diffs", {}).get(pair, 0.0),
                        window_aggregate.get("sum_dxy_comps", {}).get(pair, 0.0),
                    )
                else:
                    # Large-sample dxy: per-site average of p_x*(1-p_y)+p_y*(1-p_x) (eq:dxy-site-large).
                    window_dxy_large = window_aggregate.get("sum_dxy_large", {}).get(
                        pair, 0.0
                    )
                    theta = _safe_ratio(window_dxy_large, window_n_seg)
                form_label = "unbiased" if cfg.unbiased else "large_sample"
                rows.append(
                    f"{chrom},{wind_start},{wind_end},{window_n_seg},"
                    f"{cfg.pop_x},{cfg.pop_y},{form_label},{theta}"
                )

        # Write CSV.
        with _open_results_csv(prefix, "dxy", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # fst
    # ==================================================================
    fst_configs = config.get("fst_configs", [])
    if fst_configs:
        header = [
            "CHROM",
            "WIND_START",
            "WIND_END",
            "N_SNPS",
            "POP_X",
            "POP_Y",
            "FORM",
            "fst",
        ]
        rows = []
        for chrom, window_idx in sorted_keys:
            # Aggregated allele count sums for this window.
            window_aggregate = window_agg[(chrom, window_idx)]
            # 1-based inclusive start and end positions of this window.
            wind_start = window_idx * block_size + 1
            wind_end = (window_idx + 1) * block_size
            for cfg in fst_configs:
                # Frozenset key for jointly callable site count lookup.
                n_seg_key = frozenset({cfg.pop_x, cfg.pop_y})
                window_n_seg = window_aggregate.get("n_seg", {}).get(n_seg_key, 0)
                # Drop windows with no jointly callable sites for the pair.
                if window_n_seg == 0:
                    dropped_warnings.append(
                        (
                            "R001",
                            chrom,
                            f"{chrom}:{wind_start}-{wind_end}",
                            f"Window dropped: n_seg=0 for fst ({cfg.pop_x},{cfg.pop_y})",
                        )
                    )
                    continue
                # Sorted pair key for fst sum lookups.
                pair = _pair_key(cfg.pop_x, cfg.pop_y)
                if cfg.unbiased:
                    # Unbiased Hudson Fst: numerator over denominator from per-site sums.
                    theta = _clamp_bounded(
                        _safe_ratio(
                            window_aggregate.get("sum_fst_num_unbiased", {}).get(
                                pair, 0.0
                            ),
                            window_aggregate.get("sum_fst_denom_unbiased", {}).get(
                                pair, 0.0
                            ),
                        )
                    )
                else:
                    # Large-sample Fst: sum_fst_num / sum_fst_denom (eq:fst-N-large / eq:fst-D-large).
                    window_fst_num = window_aggregate.get("sum_fst_num", {}).get(
                        pair, 0.0
                    )
                    window_fst_denom = window_aggregate.get("sum_fst_denom", {}).get(
                        pair, 0.0
                    )
                    theta = _clamp_bounded(
                        _safe_ratio(window_fst_num, window_fst_denom)
                    )
                form_label = "unbiased" if cfg.unbiased else "large_sample"
                rows.append(
                    f"{chrom},{wind_start},{wind_end},{window_n_seg},"
                    f"{cfg.pop_x},{cfg.pop_y},{form_label},{theta}"
                )

        # Write CSV.
        with _open_results_csv(prefix, "fst", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # f2
    # ==================================================================
    f2_configs = config.get("f2_configs", [])
    if f2_configs:
        header = [
            "CHROM",
            "WIND_START",
            "WIND_END",
            "N_SNPS",
            "POP_A",
            "POP_B",
            "FORM",
            "f2",
        ]
        rows = []
        for chrom, window_idx in sorted_keys:
            # Aggregated allele count sums for this window.
            window_aggregate = window_agg[(chrom, window_idx)]
            # 1-based inclusive start and end positions of this window.
            wind_start = window_idx * block_size + 1
            wind_end = (window_idx + 1) * block_size
            for cfg in f2_configs:
                # Frozenset key for jointly callable site count lookup.
                n_seg_key = frozenset({cfg.pop_a, cfg.pop_b})
                window_n_seg = window_aggregate.get("n_seg", {}).get(n_seg_key, 0)
                # Drop windows with no callable sites for the pair.
                if window_n_seg == 0:
                    dropped_warnings.append(
                        (
                            "R001",
                            chrom,
                            f"{chrom}:{wind_start}-{wind_end}",
                            f"Window dropped: n_seg=0 for f2 ({cfg.pop_a},{cfg.pop_b})",
                        )
                    )
                    continue
                # Sorted pair key for f2 sum lookups.
                pair = _pair_key(cfg.pop_a, cfg.pop_b)
                # Select large-sample or unbiased sum field based on the configuration.
                field = "sum_f2_unbiased" if cfg.unbiased else "sum_f2"
                # Windowed f2 = sum / n_seg.
                theta = _safe_ratio(
                    window_aggregate.get(field, {}).get(pair, 0.0),
                    window_n_seg,
                )
                form_label = "unbiased" if cfg.unbiased else "large_sample"
                rows.append(
                    f"{chrom},{wind_start},{wind_end},{window_n_seg},"
                    f"{cfg.pop_a},{cfg.pop_b},{form_label},{theta}"
                )

        # Write CSV.
        with _open_results_csv(prefix, "f2", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # f3
    # ==================================================================
    f3_configs = config.get("f3_configs", [])
    if f3_configs:
        header = [
            "CHROM",
            "WIND_START",
            "WIND_END",
            "N_SNPS",
            "POP_C",
            "POP_A",
            "POP_B",
            "FORM",
            "f3",
        ]
        rows = []
        for chrom, window_idx in sorted_keys:
            # Aggregated allele count sums for this window.
            window_aggregate = window_agg[(chrom, window_idx)]
            # 1-based inclusive start and end positions of this window.
            wind_start = window_idx * block_size + 1
            wind_end = (window_idx + 1) * block_size
            for cfg in f3_configs:
                # Frozenset key for jointly callable site count lookup.
                n_seg_key = frozenset({cfg.pop_a, cfg.pop_b, cfg.pop_c})
                window_n_seg = window_aggregate.get("n_seg", {}).get(n_seg_key, 0)
                # Drop windows with no callable sites for the triple.
                if window_n_seg == 0:
                    dropped_warnings.append(
                        (
                            "R001",
                            chrom,
                            f"{chrom}:{wind_start}-{wind_end}",
                            f"Window dropped: n_seg=0 for f3 "
                            f"({cfg.pop_c};{cfg.pop_a},{cfg.pop_b})",
                        )
                    )
                    continue
                # Ordered triple key with pop_c first.
                triple = (cfg.pop_c, cfg.pop_a, cfg.pop_b)
                # Select large-sample or unbiased sum field based on the configuration.
                field = "sum_f3_unbiased" if cfg.unbiased else "sum_f3"
                # Windowed f3 = sum / n_seg.
                theta = _safe_ratio(
                    window_aggregate.get(field, {}).get(triple, 0.0),
                    window_n_seg,
                )
                form_label = "unbiased" if cfg.unbiased else "large_sample"
                rows.append(
                    f"{chrom},{wind_start},{wind_end},{window_n_seg},"
                    f"{cfg.pop_c},{cfg.pop_a},{cfg.pop_b},{form_label},{theta}"
                )

        # Write CSV.
        with _open_results_csv(prefix, "f3", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # f3star
    # ==================================================================
    f3star_configs = config.get("f3star_configs", [])
    if f3star_configs:
        header = [
            "CHROM",
            "WIND_START",
            "WIND_END",
            "N_SNPS",
            "POP_C",
            "POP_A",
            "POP_B",
            "FORM",
            "f3star",
        ]
        rows = []
        for chrom, window_idx in sorted_keys:
            # Aggregated allele count sums for this window.
            window_aggregate = window_agg[(chrom, window_idx)]
            # 1-based inclusive start and end positions of this window.
            wind_start = window_idx * block_size + 1
            wind_end = (window_idx + 1) * block_size
            for cfg in f3star_configs:
                # Frozenset key for jointly callable site count lookup.
                n_seg_key = frozenset({cfg.pop_a, cfg.pop_b, cfg.pop_c})
                window_n_seg = window_aggregate.get("n_seg", {}).get(n_seg_key, 0)
                # Drop windows with no callable sites for the triple.
                if window_n_seg == 0:
                    dropped_warnings.append(
                        (
                            "R001",
                            chrom,
                            f"{chrom}:{wind_start}-{wind_end}",
                            f"Window dropped: n_seg=0 for f3star "
                            f"({cfg.pop_c};{cfg.pop_a},{cfg.pop_b})",
                        )
                    )
                    continue
                # Ordered triple key with pop_c first.
                triple = (cfg.pop_c, cfg.pop_a, cfg.pop_b)
                if cfg.unbiased:
                    # Unbiased windowed f3star (eq:f3-star).
                    f3_val = _safe_ratio(
                        window_aggregate.get("sum_f3_unbiased", {}).get(triple, 0.0),
                        window_n_seg,
                    )
                    # Ratio-of-sums h_hat_PC (eq:h-hat).
                    h_hat = _safe_ratio(
                        window_aggregate.get("sum_h_hat_num", {}).get(cfg.pop_c, 0.0),
                        window_aggregate.get("sum_h_hat_denom", {}).get(cfg.pop_c, 0.0),
                    )
                    # f3star = f3_hat / (2 * h_hat_pc).
                    theta = _safe_ratio(f3_val, 2 * h_hat)
                else:
                    # Large-sample windowed f3star (eq:f3-star).
                    f3_val = _safe_ratio(
                        window_aggregate.get("sum_f3", {}).get(triple, 0.0),
                        window_n_seg,
                    )
                    # Per-site average large-sample half-heterozygosity of P_C (eq:h-large).
                    h_pc = _safe_ratio(
                        window_aggregate.get("sum_h_large", {}).get(cfg.pop_c, 0.0),
                        window_n_seg,
                    )
                    # f3star = f3 / (2 * h_pc).
                    theta = _safe_ratio(f3_val, 2 * h_pc)
                form_label = "unbiased" if cfg.unbiased else "large_sample"
                rows.append(
                    f"{chrom},{wind_start},{wind_end},{window_n_seg},"
                    f"{cfg.pop_c},{cfg.pop_a},{cfg.pop_b},{form_label},{theta}"
                )

        # Write CSV.
        with _open_results_csv(prefix, "f3star", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # f4
    # ==================================================================
    f4_configs = config.get("f4_configs", [])
    if f4_configs:
        header = [
            "CHROM",
            "WIND_START",
            "WIND_END",
            "N_SNPS",
            "POP_A",
            "POP_B",
            "POP_C",
            "POP_D",
            "f4",
        ]
        rows = []
        for chrom, window_idx in sorted_keys:
            # Aggregated allele count sums for this window.
            window_aggregate = window_agg[(chrom, window_idx)]
            # 1-based inclusive start and end positions of this window.
            wind_start = window_idx * block_size + 1
            wind_end = (window_idx + 1) * block_size
            for cfg in f4_configs:
                # Frozenset key for jointly callable site count lookup.
                n_seg_key = frozenset({cfg.pop_a, cfg.pop_b, cfg.pop_c, cfg.pop_d})
                window_n_seg = window_aggregate.get("n_seg", {}).get(n_seg_key, 0)
                # Drop windows with no callable sites for the quad.
                if window_n_seg == 0:
                    dropped_warnings.append(
                        (
                            "R001",
                            chrom,
                            f"{chrom}:{wind_start}-{wind_end}",
                            f"Window dropped: n_seg=0 for f4 "
                            f"({cfg.pop_a},{cfg.pop_b};{cfg.pop_c},{cfg.pop_d})",
                        )
                    )
                    continue
                # Ordered quad key (pop_a, pop_b, pop_c, pop_d).
                quad = (cfg.pop_a, cfg.pop_b, cfg.pop_c, cfg.pop_d)
                # Windowed f4 = sum / n_seg.
                theta = _safe_ratio(
                    window_aggregate.get("sum_f4", {}).get(quad, 0.0),
                    window_n_seg,
                )
                rows.append(
                    f"{chrom},{wind_start},{wind_end},{window_n_seg},"
                    f"{cfg.pop_a},{cfg.pop_b},{cfg.pop_c},{cfg.pop_d},"
                    f"{theta}"
                )

        # Write CSV.
        with _open_results_csv(prefix, "f4", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # f4ratio
    # ==================================================================
    f4ratio_configs = config.get("f4ratio_configs", [])
    if f4ratio_configs:
        header = [
            "CHROM",
            "WIND_START",
            "WIND_END",
            "N_SNPS",
            "POP_A",
            "POP_B",
            "POP_C",
            "POP_D",
            "POP_E",
            "f4ratio",
        ]
        rows = []
        for chrom, window_idx in sorted_keys:
            # Aggregated allele count sums for this window.
            window_aggregate = window_agg[(chrom, window_idx)]
            # 1-based inclusive start and end positions of this window.
            wind_start = window_idx * block_size + 1
            wind_end = (window_idx + 1) * block_size
            for cfg in f4ratio_configs:
                # Frozenset key for the five-population n_seg lookup.
                n_seg_key = frozenset(
                    {cfg.pop_a, cfg.pop_b, cfg.pop_c, cfg.pop_d, cfg.pop_e}
                )
                window_n_seg = window_aggregate.get("n_seg", {}).get(n_seg_key, 0)
                # Drop windows with no callable sites for the five-population set.
                if window_n_seg == 0:
                    dropped_warnings.append(
                        (
                            "R001",
                            chrom,
                            f"{chrom}:{wind_start}-{wind_end}",
                            f"Window dropped: n_seg=0 for f4ratio "
                            f"({cfg.pop_a},{cfg.pop_b},{cfg.pop_c},"
                            f"{cfg.pop_d},{cfg.pop_e})",
                        )
                    )
                    continue
                # Numerator f4: f4(pop_a, pop_d; pop_e, pop_c).
                num_quad = (cfg.pop_a, cfg.pop_d, cfg.pop_e, cfg.pop_c)
                # Denominator f4: f4(pop_a, pop_d; pop_b, pop_c).
                den_quad = (cfg.pop_a, cfg.pop_d, cfg.pop_b, cfg.pop_c)
                # Windowed admixture proportion alpha; 1/M cancels in the ratio.
                alpha = _clamp_bounded(
                    _safe_ratio(
                        window_aggregate.get("sum_f4", {}).get(num_quad, 0.0),
                        window_aggregate.get("sum_f4", {}).get(den_quad, 0.0),
                    )
                )
                rows.append(
                    f"{chrom},{wind_start},{wind_end},{window_n_seg},"
                    f"{cfg.pop_a},{cfg.pop_b},{cfg.pop_c},{cfg.pop_d},"
                    f"{cfg.pop_e},{alpha}"
                )

        # Write CSV.
        with _open_results_csv(prefix, "f4ratio", gzip_output) as f:
            f.write(",".join(header) + "\n")
            for row_str in rows:
                f.write(row_str + "\n")

    # ==================================================================
    # Append dropped window warnings to the gzipped log file (header
    # written by config phase; always called even when list is empty).
    # ==================================================================
    with gzip.open(log_path, "at") as f:
        for code, chrom, context, message in dropped_warnings:
            if config["verbose"]:
                f.write(f"RESULTS\t{code}\t{chrom}\t{context}\t{message}\n")
            else:
                f.write(f"RESULTS\t{code}\t{chrom}\t{context}\n")
