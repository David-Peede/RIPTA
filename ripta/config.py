import gzip
import os

import yaml

from ripta.stats import FStatConfig, PairwiseConfig, SitePatternConfig

# Valid statistic names for site_pattern_statistics entries.
_VALID_SP_STATS = frozenset({"D", "Danc", "Dplus", "fhom", "Dp", "df", "Dstar"})

# Validation specifications for each statistics config section. Each spec defines
# the required population fields, GP boolean fields, whether an unbiased flag is
# expected, and whether a stats list is expected.
_STAT_SPECS = {
    "site_pattern_statistics": {
        "pop_fields": ("p1", "p2", "p3", "outgroup"),
        "gp_fields": ("use_gp_p1", "use_gp_p2", "use_gp_p3", "use_gp_outgroup"),
        "has_unbiased": False,
        "has_stats_list": True,
    },
    "D3": {
        "pop_fields": ("p1", "p2", "p3"),
        "gp_fields": ("use_gp_p1", "use_gp_p2", "use_gp_p3"),
        "has_unbiased": True,
        "has_stats_list": False,
    },
    "pi": {
        "pop_fields": ("pop",),
        "gp_fields": ("use_gp",),
        "has_unbiased": True,
        "has_stats_list": False,
    },
    "dxy": {
        "pop_fields": ("pop_x", "pop_y"),
        "gp_fields": ("use_gp_x", "use_gp_y"),
        "has_unbiased": True,
        "has_stats_list": False,
    },
    "fst": {
        "pop_fields": ("pop_x", "pop_y"),
        "gp_fields": ("use_gp_x", "use_gp_y"),
        "has_unbiased": True,
        "has_stats_list": False,
    },
    "f2": {
        "pop_fields": ("pop_a", "pop_b"),
        "gp_fields": ("use_gp_a", "use_gp_b"),
        "has_unbiased": True,
        "has_stats_list": False,
    },
    "f3": {
        "pop_fields": ("pop_a", "pop_b", "pop_c"),
        "gp_fields": ("use_gp_a", "use_gp_b", "use_gp_c"),
        "has_unbiased": True,
        "has_stats_list": False,
    },
    "f3star": {
        "pop_fields": ("pop_a", "pop_b", "pop_c"),
        "gp_fields": ("use_gp_a", "use_gp_b", "use_gp_c"),
        "has_unbiased": True,
        "has_stats_list": False,
    },
    "f4": {
        "pop_fields": ("pop_a", "pop_b", "pop_c", "pop_d"),
        "gp_fields": ("use_gp_a", "use_gp_b", "use_gp_c", "use_gp_d"),
        "has_unbiased": False,
        "has_stats_list": False,
    },
    "f4ratio": {
        "pop_fields": ("pop_a", "pop_b", "pop_c", "pop_d", "pop_e"),
        "gp_fields": ("use_gp_a", "use_gp_b", "use_gp_c", "use_gp_d", "use_gp_e"),
        "has_unbiased": False,
        "has_stats_list": False,
    },
}


def _read_metadata(meta_data_path):
    """Read the population metadata TSV and return a population-to-individuals mapping.

    Args:
        meta_data_path (str): Path to the tab-separated metadata file.

    Returns:
        pop_to_inds (dict): Mapping of population label strings to lists of individual
            ID strings.
        errors (list): List of (code, context, message) error tuples.

    Notes
    -----
    Emits error code C004 for malformed metadata rows, missing required IND or POP
    columns, and underlying file read errors.
    """
    errors = []
    pop_to_inds = {}
    # Track individual IDs already seen across all populations to detect duplicates.
    seen_inds = set()

    try:
        with open(meta_data_path, "r") as f:
            # Read the header line of the metadata file.
            header = f.readline().strip().split("\t")

            # Both IND and POP columns are required in the header.
            if "IND" not in header or "POP" not in header:
                errors.append(
                    (
                        "C004",
                        meta_data_path,
                        f"Metadata file must contain IND and POP columns, "
                        f"found: {', '.join(header)}",
                    )
                )
                return pop_to_inds, errors

            # Locate the column indices for IND and POP.
            ind_col = header.index("IND")
            pop_col = header.index("POP")

            # Parse each data row after the header.
            for line_num, line in enumerate(f, 2):
                # Strip leading and trailing whitespace from the data row.
                stripped = line.strip()
                # Skip blank lines.
                if not stripped:
                    continue
                # Split the stripped line into tab-separated fields.
                fields = stripped.split("\t")
                # Row must have enough columns to cover both IND and POP positions.
                min_cols = max(ind_col, pop_col) + 1
                if len(fields) < min_cols:
                    errors.append(
                        (
                            "C004",
                            f"{meta_data_path}:{line_num}",
                            f"Line {line_num} has {len(fields)} columns, "
                            f"need at least {min_cols}",
                        )
                    )
                    continue
                # Extract individual ID and population label.
                ind = fields[ind_col]
                pop = fields[pop_col]
                # Reject the row if this individual ID has already been registered.
                if ind in seen_inds:
                    errors.append(
                        (
                            "C004",
                            f"{meta_data_path}:{line_num}",
                            f"Duplicate individual ID '{ind}' on line {line_num}",
                        )
                    )
                    continue
                # Record the individual ID as seen across all populations.
                seen_inds.add(ind)
                # Accumulate individuals under their population label.
                if pop not in pop_to_inds:
                    pop_to_inds[pop] = []
                pop_to_inds[pop].append(ind)

    except Exception as error:
        errors.append(("C004", meta_data_path, f"Error reading metadata file: {error}"))

    return pop_to_inds, errors


def _validate_stat_entry(entry, entry_idx, section_name, spec, valid_pops):
    """Validate a single entry dict from a statistics config section.

    Args:
        entry (dict): One entry from the statistics section list.
        entry_idx (int): Zero-based index of this entry in the list.
        section_name (str): Name of the statistics section (e.g., 'pi', 'f4').
        spec (dict): Validation spec from _STAT_SPECS for this section.
        valid_pops (set): Set of population label strings found in metadata.

    Returns:
        errors (list): List of (code, context, message) error tuples.

    Notes
    -----
    Emits the following error codes: C005 for missing required fields or wrong field
    types; C006 when a population label is not present in the metadata POP column;
    C007 for an empty or invalid stats list in site_pattern_statistics.
    """
    errors = []
    # Label for error context strings.
    context = f"{section_name}[{entry_idx}]"

    # Validate each population field: must be present, a string, and in metadata.
    for pop_field in spec["pop_fields"]:
        if pop_field not in entry:
            errors.append(("C005", context, f"Missing required field '{pop_field}'"))
        elif not isinstance(entry[pop_field], str):
            errors.append(
                (
                    "C005",
                    f"{context}.{pop_field}",
                    f"Field '{pop_field}' must be a string, "
                    f"got {type(entry[pop_field]).__name__}",
                )
            )
        elif entry[pop_field] not in valid_pops:
            errors.append(
                (
                    "C006",
                    f"{context}.{pop_field}",
                    f"Population '{entry[pop_field]}' not found in metadata POP column",
                )
            )

    # Validate each GP boolean field: must be present and boolean.
    for gp_field in spec["gp_fields"]:
        if gp_field not in entry:
            errors.append(("C005", context, f"Missing required field '{gp_field}'"))
        elif not isinstance(entry[gp_field], bool):
            errors.append(
                (
                    "C005",
                    f"{context}.{gp_field}",
                    f"Field '{gp_field}' must be a boolean, "
                    f"got {type(entry[gp_field]).__name__}",
                )
            )

    # Validate the unbiased flag when required by this statistic family.
    if spec["has_unbiased"]:
        if "unbiased" not in entry:
            errors.append(("C005", context, "Missing required field 'unbiased'"))
        elif not isinstance(entry["unbiased"], bool):
            errors.append(
                (
                    "C005",
                    f"{context}.unbiased",
                    f"Field 'unbiased' must be a boolean, "
                    f"got {type(entry['unbiased']).__name__}",
                )
            )

    # Validate the stats list for site_pattern_statistics entries.
    if spec["has_stats_list"]:
        if "stats" not in entry:
            errors.append(("C005", context, "Missing required field 'stats'"))
        elif not isinstance(entry["stats"], list) or len(entry["stats"]) == 0:
            errors.append(
                ("C007", f"{context}.stats", "Field 'stats' must be a non-empty list")
            )
        else:
            # Check each stat name against the valid set.
            invalid = set(entry["stats"]) - _VALID_SP_STATS
            if invalid:
                errors.append(
                    (
                        "C007",
                        f"{context}.stats",
                        f"Invalid statistic names: {', '.join(sorted(invalid))}; "
                        f"valid names are: {', '.join(sorted(_VALID_SP_STATS))}",
                    )
                )

    return errors


def _read_vcf_header(vcf_path):
    """Read and return the #CHROM header line from a VCF file as a list of field strings.

    Args:
        vcf_path (str): Path to the VCF file (plain text or gzip-compressed).

    Returns:
        header (list): List of header field strings, or None if no header was found.
    """
    # Select the appropriate file opener based on file extension.
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            # The #CHROM line is the VCF column header.
            if line.startswith("#CHROM"):
                # Strip trailing whitespace and split into tab-separated field strings.
                return line.strip().split("\t")
            # Stop at the first non-meta-information line.
            if not line.startswith("##"):
                break
    return None


def _write_error_log(errors, log_path_prefix, verbose):
    """Write config validation errors to the gzipped TSV log file.

    Args:
        errors (list): List of (code, context, message) error tuples.
        log_path_prefix (str): Prefix for the log file path.
        verbose (bool): Whether to include human-readable messages in the log.
    """
    # Construct the log file path.
    log_path = f"{log_path_prefix}_ripta.log.gz"

    # Ensure the parent directory exists for the log file.
    parent = os.path.dirname(log_path)
    if parent:
        os.makedirs(parent, exist_ok=True)

    # Write the header and error rows as gzipped TSV.
    with gzip.open(log_path, "wt") as f:
        # Header includes MESSAGE only in verbose mode to keep column counts consistent.
        if verbose:
            f.write("CATEGORY\tCODE\tCHROM\tCONTEXT\tMESSAGE\n")
        else:
            f.write("CATEGORY\tCODE\tCHROM\tCONTEXT\n")
        for code, context, message in errors:
            # CHROM column is "." for non-VCF (configuration) errors per the log schema.
            if verbose:
                # Verbose rows include the human-readable MESSAGE column.
                f.write(f"CONFIG\t{code}\t.\t{context}\t{message}\n")
            else:
                # Non-verbose rows omit the MESSAGE column entirely.
                f.write(f"CONFIG\t{code}\t.\t{context}\n")


def read_config(config_path):
    """Read and validate the YAML config file, returning a populated configuration dictionary.

    Args:
        config_path (str): Path to the YAML configuration file.

    Returns:
        config (dict): Validated configuration dictionary with populated configuration
            objects, idx_dicc, and gp_inds.

    Raises:
        FileNotFoundError: If the configuration file does not exist.
        ValueError: If any configuration parameter fails validation.

    Notes
    -----
    Validation proceeds in phases: global parameters, metadata, statistics sections,
    VCF headers. All errors are collected, written to the log file, and then raised
    together so the user sees every problem in a single pass. The log file is written
    to {log_path_prefix}_ripta.log.gz.
    """
    # Collect all errors as (code, context, message) tuples.
    errors = []

    # Check that the config file itself exists before attempting to parse it.
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file {config_path} not found")

    # Parse the YAML file.
    try:
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
    except yaml.YAMLError as error:
        raise ValueError(f"Error parsing YAML file: {error}")

    # -----------------------------------------------------------------------
    # Phase 1: Validate global parameters.
    # -----------------------------------------------------------------------

    # Check that all required global keys are present.
    required_globals = {
        "chroms",
        "vcf_path",
        "meta_data_path",
        "results_path_prefix",
        "log_path_prefix",
        "block_size",
        "drop_last_block",
        "number_of_threads",
        "verbose",
        "windows",
    }
    # Identify which required keys are absent from the configuration dictionary.
    missing_globals = required_globals - set(config.keys())
    for key in sorted(missing_globals):
        errors.append(("C001", key, f"Missing required global parameter '{key}'"))

    # If any required global key is missing, log what we can and raise early
    # because downstream validation depends on these keys being present.
    if missing_globals:
        verbose = config.get("verbose", True)
        # Fall back to verbose=True when the user-supplied value is not a boolean.
        if not isinstance(verbose, bool):
            verbose = True
        log_prefix = config.get("log_path_prefix")
        if log_prefix:
            try:
                _write_error_log(errors, log_prefix, verbose)
            except Exception:
                pass
        raise ValueError(
            "Configuration validation failed:\n"
            + "\n".join(f"[{code}] {context}: {msg}" for code, context, msg in errors)
        )

    # Validate chroms: must be a non-empty dict of string keys to positive integers.
    chroms_ok = True
    if not isinstance(config["chroms"], dict) or len(config["chroms"]) == 0:
        errors.append(
            (
                "C002",
                "chroms",
                "chroms must be a non-empty dict mapping chrom_id strings to "
                "positive integers",
            )
        )
        chroms_ok = False
    else:
        for chrom_id, length in config["chroms"].items():
            if not isinstance(chrom_id, str):
                errors.append(
                    (
                        "C002",
                        f"chroms.{chrom_id}",
                        f"Chrom ID must be a string, got {type(chrom_id).__name__}",
                    )
                )
                chroms_ok = False
            if not isinstance(length, int) or length <= 0:
                errors.append(
                    (
                        "C002",
                        f"chroms.{chrom_id}",
                        f"Chrom length must be a positive integer, got {length}",
                    )
                )
                chroms_ok = False

    # Validate vcf_path: string (single-VCF mode) or dict matching chroms keys.
    vcf_path = config["vcf_path"]
    if isinstance(vcf_path, str):
        # Single-VCF mode — valid structure.
        pass
    elif isinstance(vcf_path, dict):
        # Multi-VCF mode — keys must exactly match chroms keys.
        if chroms_ok:
            chroms_keys = set(config["chroms"].keys())
            vcf_keys = set(vcf_path.keys())
            if chroms_keys != vcf_keys:
                missing_in_vcf = chroms_keys - vcf_keys
                extra_in_vcf = vcf_keys - chroms_keys
                parts = []
                if missing_in_vcf:
                    parts.append(
                        f"missing from vcf_path: "
                        f"{', '.join(sorted(missing_in_vcf))}"
                    )
                if extra_in_vcf:
                    parts.append(
                        f"extra in vcf_path: " f"{', '.join(sorted(extra_in_vcf))}"
                    )
                errors.append(
                    (
                        "C003",
                        "vcf_path",
                        f"vcf_path dict keys must exactly match chroms keys; "
                        f"{'; '.join(parts)}",
                    )
                )
    else:
        errors.append(
            (
                "C003",
                "vcf_path",
                f"vcf_path must be a string or dict, got {type(vcf_path).__name__}",
            )
        )

    # Validate meta_data_path: file must exist.
    meta_exists = os.path.exists(config["meta_data_path"])
    if not meta_exists:
        errors.append(
            (
                "C003",
                "meta_data_path",
                f"Metadata file {config['meta_data_path']} not found",
            )
        )

    # Validate results_path_prefix: parent directory must exist.
    results_parent = os.path.dirname(config["results_path_prefix"])
    if results_parent and not os.path.isdir(results_parent):
        errors.append(
            (
                "C003",
                "results_path_prefix",
                f"Parent directory {results_parent} does not exist",
            )
        )

    # Validate log_path_prefix: parent directory must exist.
    log_parent = os.path.dirname(config["log_path_prefix"])
    if log_parent and not os.path.isdir(log_parent):
        errors.append(
            ("C003", "log_path_prefix", f"Parent directory {log_parent} does not exist")
        )

    # Validate block_size: positive integer strictly less than every chrom length.
    if not isinstance(config["block_size"], int) or config["block_size"] <= 0:
        errors.append(
            (
                "C002",
                "block_size",
                f"block_size must be a positive integer, got {config['block_size']}",
            )
        )
    elif chroms_ok:
        # Check that block_size is strictly less than every chrom length.
        for chrom_id, length in config["chroms"].items():
            if isinstance(length, int) and config["block_size"] >= length:
                errors.append(
                    (
                        "C002",
                        "block_size",
                        f"block_size ({config['block_size']}) must be strictly less "
                        f"than chrom length for {chrom_id} ({length})",
                    )
                )

    # Validate drop_last_block: must be a boolean.
    if not isinstance(config["drop_last_block"], bool):
        errors.append(
            (
                "C002",
                "drop_last_block",
                f"drop_last_block must be a boolean, "
                f"got {type(config['drop_last_block']).__name__}",
            )
        )

    # Validate number_of_threads: positive integer or the string 'max'.
    n_threads = config["number_of_threads"]
    if not ((isinstance(n_threads, int) and n_threads > 0) or n_threads == "max"):
        errors.append(
            (
                "C002",
                "number_of_threads",
                f"number_of_threads must be a positive integer or 'max', "
                f"got {n_threads}",
            )
        )

    # Validate verbose: must be a boolean.
    if not isinstance(config["verbose"], bool):
        errors.append(
            (
                "C002",
                "verbose",
                f"verbose must be a boolean, got {type(config['verbose']).__name__}",
            )
        )

    # Validate windows: must be a boolean.
    if not isinstance(config["windows"], bool):
        errors.append(
            (
                "C002",
                "windows",
                f"windows must be a boolean, got {type(config['windows']).__name__}",
            )
        )

    # Validate gzip_output: optional global, defaults to False when absent.
    # When present, it must be a boolean. Controls whether result CSVs are
    # written as plain .csv or gzipped .csv.gz.
    if "gzip_output" not in config:
        config["gzip_output"] = False
    elif not isinstance(config["gzip_output"], bool):
        errors.append(
            (
                "C002",
                "gzip_output",
                f"gzip_output must be a boolean, got "
                f"{type(config['gzip_output']).__name__}",
            )
        )

    # -----------------------------------------------------------------------
    # Phase 2: Read and validate metadata.
    # -----------------------------------------------------------------------

    # Mapping of population label strings to lists of individual ID strings.
    pop_to_inds = {}
    # Set of population labels recognized by downstream statistic validation.
    valid_pops = set()
    # Flag indicating whether at least one population was successfully loaded.
    metadata_ok = False

    if meta_exists:
        # Read the metadata TSV and collect any C004 errors raised during parsing.
        pop_to_inds, meta_errors = _read_metadata(config["meta_data_path"])
        errors.extend(meta_errors)
        # Derive the valid population label set from the loaded metadata.
        valid_pops = set(pop_to_inds.keys())
        # Treat metadata as usable when at least one population was loaded.
        metadata_ok = len(pop_to_inds) > 0
        # Store the population-to-individuals mapping in the configuration
        # dictionary so io.py does not need to re-open the metadata file.
        if metadata_ok:
            config["pop_to_inds"] = pop_to_inds

    # -----------------------------------------------------------------------
    # Phase 3: Validate statistics sections.
    # -----------------------------------------------------------------------

    # Only validate statistics entries if metadata loaded successfully,
    # because population labels need to be checked against valid_pops.
    if metadata_ok:
        for section_name, spec in _STAT_SPECS.items():
            # Statistics sections are all optional.
            if section_name not in config:
                continue

            entries = config[section_name]
            # Each statistics section must be a list of config dicts.
            if not isinstance(entries, list):
                errors.append(
                    (
                        "C005",
                        section_name,
                        f"Statistics section '{section_name}' must be a list of "
                        f"config dicts, got {type(entries).__name__}",
                    )
                )
                continue

            # Validate each entry in the section.
            for entry_idx, entry in enumerate(entries):
                if not isinstance(entry, dict):
                    errors.append(
                        (
                            "C005",
                            f"{section_name}[{entry_idx}]",
                            f"Each entry must be a dict, "
                            f"got {type(entry).__name__}",
                        )
                    )
                    continue
                entry_errors = _validate_stat_entry(
                    entry,
                    entry_idx,
                    section_name,
                    spec,
                    valid_pops,
                )
                errors.extend(entry_errors)

    # -----------------------------------------------------------------------
    # Phase 4: Read VCF headers and build idx_dicc.
    # -----------------------------------------------------------------------

    # Initialize the individual-to-column-index mapping.
    idx_dicc = {}

    # Collect all populations referenced across all statistic configs.
    required_pops = set()
    for section_name, spec in _STAT_SPECS.items():
        if section_name not in config:
            continue
        entries = config[section_name]
        if not isinstance(entries, list):
            continue
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            for pop_field in spec["pop_fields"]:
                pop_label = entry.get(pop_field)
                if isinstance(pop_label, str):
                    required_pops.add(pop_label)

    # Collect all individual IDs that belong to required populations.
    required_inds = set()
    for pop in required_pops:
        if pop in pop_to_inds:
            required_inds.update(pop_to_inds[pop])

    # Determine which VCF file(s) to read headers from.
    vcf_paths_to_read = []
    if isinstance(config["vcf_path"], str):
        vcf_paths_to_read = [config["vcf_path"]]
    elif isinstance(config["vcf_path"], dict):
        vcf_paths_to_read = list(config["vcf_path"].values())

    # Read VCF headers and validate individual presence.
    headers = []
    for vcf in vcf_paths_to_read:
        if not os.path.exists(vcf):
            errors.append(("C010", vcf, f"VCF file {vcf} not found"))
            continue
        try:
            header = _read_vcf_header(vcf)
            if header is None:
                errors.append(
                    ("C010", vcf, f"No #CHROM header line found in VCF file {vcf}")
                )
            else:
                headers.append((vcf, header))
        except Exception as error:
            errors.append(("C010", vcf, f"Error reading VCF file {vcf}: {error}"))

    # Build idx_dicc from the first valid VCF header.
    if headers and required_inds:
        first_vcf, first_header = headers[0]
        # Pre-build a name-to-index mapping across all required individuals.
        header_index = {name: idx for idx, name in enumerate(first_header)}
        for ind in sorted(required_inds):
            col_idx = header_index.get(ind)
            if col_idx is not None:
                # Map each individual to their 0-based column index in the VCF.
                idx_dicc[ind] = col_idx
            else:
                errors.append(
                    (
                        "C011",
                        f"{first_vcf}:{ind}",
                        f"Individual '{ind}' not found in VCF header of {first_vcf}",
                    )
                )

        # In multi-VCF mode, verify all headers have consistent sample columns.
        if len(headers) > 1:
            # Sample columns start at index 9 in VCF format.
            first_samples = first_header[9:]
            for other_vcf, other_header in headers[1:]:
                other_samples = other_header[9:]
                if first_samples != other_samples:
                    errors.append(
                        (
                            "C012",
                            f"{first_vcf} vs {other_vcf}",
                            "VCF headers have inconsistent sample columns",
                        )
                    )

    # -----------------------------------------------------------------------
    # Phase 5: Build config objects and gp_inds.
    # -----------------------------------------------------------------------

    # Default containers used when metadata_ok is False; Phase 5 is skipped entirely
    # in that case because constructing statistic configuration objects requires a
    # validated population mapping.
    sp_configs = []

    if metadata_ok:
        # Build SitePatternConfig instances from site_pattern_statistics entries.
        if "site_pattern_statistics" in config and isinstance(
            config["site_pattern_statistics"], list
        ):
            for entry in config["site_pattern_statistics"]:
                if not isinstance(entry, dict):
                    continue
                try:
                    sp_configs.append(
                        SitePatternConfig(
                            p1=entry["p1"],
                            p2=entry["p2"],
                            p3=entry["p3"],
                            outgroup=entry["outgroup"],
                            stats=tuple(entry["stats"]),
                            use_gp_p1=entry["use_gp_p1"],
                            use_gp_p2=entry["use_gp_p2"],
                            use_gp_p3=entry["use_gp_p3"],
                            use_gp_outgroup=entry["use_gp_outgroup"],
                        )
                    )
                except (KeyError, TypeError):
                    # Validation errors for missing fields already recorded above.
                    continue

        # Build PairwiseConfig instances for D3, pi, dxy, and fst.
        pairwise_config_keys = {
            "D3": "d3_configs",
            "pi": "pi_configs",
            "dxy": "dxy_configs",
            "fst": "fst_configs",
        }
        for section_name, config_key in pairwise_config_keys.items():
            configs = []
            if section_name in config and isinstance(config[section_name], list):
                spec = _STAT_SPECS[section_name]
                for entry in config[section_name]:
                    if not isinstance(entry, dict):
                        continue
                    try:
                        # Build kwargs from the spec field lists.
                        kwargs = {"stat": section_name}
                        for pop_field in spec["pop_fields"]:
                            kwargs[pop_field] = entry[pop_field]
                        for gp_field in spec["gp_fields"]:
                            kwargs[gp_field] = entry[gp_field]
                        if spec["has_unbiased"]:
                            kwargs["unbiased"] = entry["unbiased"]
                        configs.append(PairwiseConfig(**kwargs))
                    except (KeyError, TypeError):
                        continue
            # Store the configuration list (empty if section was absent).
            config[config_key] = configs

        # Build FStatConfig instances for f2, f3, f3star, f4, and f4ratio.
        fstat_sections = ("f2", "f3", "f3star", "f4", "f4ratio")
        for section_name in fstat_sections:
            config_key = f"{section_name}_configs"
            configs = []
            if section_name in config and isinstance(config[section_name], list):
                spec = _STAT_SPECS[section_name]
                for entry in config[section_name]:
                    if not isinstance(entry, dict):
                        continue
                    try:
                        # Build kwargs from the spec field lists.
                        kwargs = {"stat": section_name}
                        for pop_field in spec["pop_fields"]:
                            kwargs[pop_field] = entry[pop_field]
                        for gp_field in spec["gp_fields"]:
                            kwargs[gp_field] = entry[gp_field]
                        if spec["has_unbiased"]:
                            kwargs["unbiased"] = entry["unbiased"]
                        configs.append(FStatConfig(**kwargs))
                    except (KeyError, TypeError):
                        continue
            # Store the configuration list (empty if section was absent).
            config[config_key] = configs
    else:
        # Metadata could not be loaded; default every configuration list to empty
        # so downstream phases and callers can rely on the keys being present.
        for empty_key in (
            "d3_configs",
            "pi_configs",
            "dxy_configs",
            "fst_configs",
            "f2_configs",
            "f3_configs",
            "f3star_configs",
            "f4_configs",
            "f4ratio_configs",
        ):
            config[empty_key] = []

    # Build gp_inds: populations flagged for GP use across any statistic config.
    gp_inds = {}
    for section_name, spec in _STAT_SPECS.items():
        if section_name not in config or not isinstance(config[section_name], list):
            continue
        for entry in config[section_name]:
            if not isinstance(entry, dict):
                continue
            # Check each (population field, GP flag) pair.
            for pop_field, gp_field in zip(spec["pop_fields"], spec["gp_fields"]):
                if entry.get(gp_field) is True:
                    pop_label = entry.get(pop_field)
                    if pop_label and pop_label in pop_to_inds:
                        # Add all individuals in this population to the GP set.
                        if pop_label not in gp_inds:
                            gp_inds[pop_label] = set()
                        gp_inds[pop_label].update(pop_to_inds[pop_label])

    # Add all derived objects to the configuration dictionary.
    config["sp_configs"] = sp_configs
    config["idx_dicc"] = idx_dicc
    config["gp_inds"] = gp_inds

    # -----------------------------------------------------------------------
    # Phase 6: Write log and raise if any errors were collected.
    # -----------------------------------------------------------------------

    if errors:
        # Determine verbose setting for log output; default to True if invalid.
        verbose = config.get("verbose", True)
        if not isinstance(verbose, bool):
            verbose = True

        # Attempt to write the error log before raising.
        log_prefix = config.get("log_path_prefix")
        if log_prefix:
            try:
                _write_error_log(errors, log_prefix, verbose)
            except Exception:
                pass

        # Raise with all collected errors.
        raise ValueError(
            "Configuration validation failed:\n"
            + "\n".join(f"[{code}] {context}: {msg}" for code, context, msg in errors)
        )

    # Always create the log file with a header row so that a clean run produces
    # a header-only file, distinguishing "no warnings" from "RIPTA failed to run".
    # Subsequent phases (io, results) append with 'at' mode.
    log_path = f"{config['log_path_prefix']}_ripta.log.gz"
    parent = os.path.dirname(log_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    if config["verbose"]:
        log_header = "CATEGORY\tCODE\tCHROM\tCONTEXT\tMESSAGE\n"
    else:
        log_header = "CATEGORY\tCODE\tCHROM\tCONTEXT\n"
    with gzip.open(log_path, "wt") as f:
        f.write(log_header)

    return config
