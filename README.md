![ripta_logo](./ripta_logo.png)

# Re-evaluating Introgression site Patterns Through Ancestral alleles

If you have any questions or feature requests, please leave a detailed issue. If you use `RIPTA` in your work, please cite doi: https://doi.org/10.1101/2022.12.02.518851. 



## `site_patterns.py`

```bash
usage: site_patterns.py [-h] -y YAML_PATH

arguments:
  -h, --help            Show this help message and exit.
  -y YAML_PATH, --yaml_path YAML_PATH
                        Path to the YAML configuration file.
```


## `YAML` Configuration File Parameters

Currently, all parameters are required and must be specified in the configuration file; however, this will likely change once I implement non-site pattern statistics and sliding-window capabilities.

- `vcf_path`: List of file paths to the VCF files to be analyzed.
- `contig_length`: List of contig (chromosome) lengths for each corresponding VCF file.
- `block_size`: Size of the blocks for the weighted-block Jackknife procedure (see the `appendix` directory for the exact equations used for the weighted-block Jackknife procedure written by Dr. Patterson).
- `meta_data_path`: Path to a tab-separated text file with two columns: Individuals and Population. The individual IDs must be identical to how they appear in all VCF headers.
- `results_path_prefix`: Path to write the results to i.e., `{results_path_prefix}_site_patterns.csv`.
- `log_path_prefix`: Path to write all the errors and warnings to i.e., `{log_path_prefix}_config_errors.log` and `{log_path_prefix}_vcf_processing_warnings.log`.
- `use_raf`: Boolean to denote filtering by the reference panel allele frequency (designed for analyzing imputed genomes).
- `raf_threshold`: Float to denote the reference panel allele frequency threshold for filtering (designed for analyzing imputed genomes).
- `number_of_threads`: Positive integer for the number of threads to use for parallelization; to use all available threads set this parameter to `'max'`.
- `compute_site_patterns`: Boolean to denote you want to compute site pattern statistics.
- `site_pattern_config_path`: Path to a tab-separated text file with eight columns: P1, P2, P3, P4, P1 Sampling, P2 Sampling, P3 Sampling, and P4 Sampling.
	- For the first four columns (P1, P2, P3, P4) specify the individual ID or population ID which appears in the meta data file.
	- For the last four columns (P1 Sampling, P2 Sampling, P3 Sampling, P4 Sampling) specify one of the three implemented sampling methods:
		- `pop_freq`: To denote computing allele frequencies for a population of individuals.
		- `ind_freq`: To denote computing allele frequencies for a single individual (i.e., 0 = homozygous ancestral, 0.5 = heterozygous, 1 = homozygous derived).
		- `ind_gp`: To denote computing allele frequencies for an individual from genotype probabilities.
- `site_pattern_stats`: List of site pattern statistics to compute currently the options are `'D'`, `'Danc'`, `'D+'`, `'fhom'`, `'fanc'`, and/or `'f+'`.


## Understanding Error and Warning Messages

Before computing any summary statistics, `RIPTA` will first validate the entire `YAML` configuration file and report all the detected errors before exiting. During summary statistic computation, `RIPTA` will log a warning for any position in the VCF file that is excluded from the analysis. Below is a list of all possible error and warning messages, their meanings, and suggested solutions.

### Error Messages

- `'CONFIG ERROR: VCF file {VCF} not found!'`
	- Indicates that no VCF file exists at the path specified by the `vcf_path` parameter.
	- Solution: Ensure that each VCF file path is correctly specified.
- `'CONFIG ERROR: Meta data file {META_DATA} not found!'`
	- Indicates that no meta data file exists at the path specified by the `meta_data_path` parameter.
	- Solution: Ensure the path for the meta data file is correctly specified.
- `'CONFIG ERROR: Site pattern configuration file {SITE_PATTERN_CONFIG} not found!'`
	- Indicates that no site pattern configuration file exists at the path specified by the `site_pattern_config_path` parameter.
	- Solution: Ensure the path for the site pattern configuration file is correctly specified.
- `'CONFIG ERROR: The vcf_path ({N_VCFS}) and contig_length ({N_CONTIG_LENGTHS}) paramters must have the same number of entries!'`
	- Indicates a mismatch between the number of VCF files specified by `vcf_path` and the number of contig lengths specified by `contig_length` parameters.
	- Solution: Ensure that the configuration file specifies equal numbers of VCF files and contig lengths.
- `CONFIG ERROR: The Block size must be smaller than the contig lengths!`
	- Indicates that the block size specified by the `block_size` parameter is larger than the specified contig length.
	- Solution: Ensure that the block size is smaller than each specified contig length.
- `CONFIG ERROR: The RAF threshold must be between 0 and 1!`
	- Indicates that the RAF threshold specified by the `raf_threshold` parameter is invalid.
	- Solution: Ensure that the RAF threshold is set to a value between 0 and 1.
- `CONFIG ERROR: The number of threads must be a positive integer or max!"`
	- Indicates an invalid number of threads as specified by the  `number_of_threads` parameter.
	- Solution: Ensure that the number of threads is a positive integer or `'max'` to use all available threads.
- `'CONFIG ERROR: The site_pattern_stats paramter must be specified when the compute_site_patterns parameters is set to True!'`
	- Indicates missing site pattern statistics specified by the `site_pattern_stats` parameter.
	- Solution: Ensure at least one valid site pattern statistic is specified.
- `'CONFIG ERROR: Invalid site pattern statistics ({INVALID_STATS}) specified!'`
	- Indicates that there are invalid site pattern statistics specified by the `site_pattern_stats` parameter.
	- Solution: Refer to the `YAML Configuration File Parameters` section for valid site pattern statistics.
- `'SITE PATTERN CONFIG ERROR: Line {LINE_NUMBER} has {N_COLUMNS} columns, expected 8!'`
	- Indicates that a line in the site pattern configuration file has an incorrect number of columns.
	- Solution: Refer to the `YAML Configuration File Parameters` section to correctly format the site pattern configuration file.
- `'SITE PATTERN CONFIG ERROR: Invalid sampling option "{SAMPLING_OPTION}" for P{X} "{POP_ID}" on line {LINE_NUMBER}, the sampling option must be either "pop_freq", "ind_freq", or "ind_gp"!'`
	- Indicates an invalid sampling option in the site pattern configuration file.
	- Solution: Refer to the `YAML Configuration File Parameters` section to correctly format the site pattern configuration file.
- `'SITE PATTERN CONFIG ERROR: Error reading site pattern configuration file - {ERROR}'`
	- Indicates that the site pattern configuration file could not be read.
	- Solution: Refer to the `YAML Configuration File Parameters` section to correctly format the site pattern configuration file.
- `'META DATA ERROR: Line {LINE_NUMBER} is invalid, expected at least 2 tab-separated columns!'`
	- Indicates that a line in the meta data file has an incorrect number of columns.
	- Solution: Refer to the `YAML Configuration File Parameters` section to correctly format the meta data file.
- `'META DATA ERROR: Error reading meta data file - {ERROR}'`
	- Indicates that the meta data file could not be read.
	- Solution: Refer to the `YAML Configuration File Parameters` section to correctly format the meta data file.
- `'META DATA ERROR: Required populations not found: {MISSING_POPULATIONS}!'`
	- Indicates missing individuals for populations required for site pattern statistics.
	- Solution: Ensure all individuals for all required populations are in the meta data file.
- `'META DATA ERROR: Required individuals not found: {MISSING_ INDIVIDUALS}!'`
	- Indicates missing individuals in the metadata file for required site pattern statistics.
	- Solution: Ensure all required individuals are in the meta data file.
- `'META DATA ERROR: No samples found in meta data file, please check the errors and ensure the file has at least 2 tab-separated columns where column 1 = individual ID and column 2 = population ID!'`
	- Indicates that none of the required individuals or populations to compute the site pattern statistics were found in the meta data file.
	- Solution: Refer to the `YAML Configuration File Parameters` section to correctly format the meta data file.
- `'VCF ERROR: No header line found in the VCF file {VCF} - (please ensure the VCF file is formatted properly, for help see sections 1.5-1.6 https://samtools.github.io/hts-specs/VCFv4.3.pdf)!'`
	- Indicates that no header line was found in the VCF file.
	- Solution: Ensure the VCF file includes a properly formatted header.
- `'VCF ERROR: Sample {MISSING_INDIVIDUAL} from population {POPULATION} not found in VCF header of {VCF}!'`
	- Indicates that an individual required for computing site pattern statistics is missing from the VCF file header.
	- Solution: Ensure that all required individuals are present in the VCF file header.
- `'VCF ERROR: No samples found for required population {POPULATION} in the VCF file {VCF}!'`
	- Indicates that all required individuals from a population are missing in the VCF file header.
	- Solution: Ensure that all required individuals are present in the VCF file header.
- `'VCF ERROR: Unable to open the VCF file {VCF} - {ERROR}'`
	- Indicates that the VCF file could not be open.
	- Solution: Ensure that the VCF file is not corrupted and if the VCF file ends with `.gz` that the VCF file is actually gzipped or bgzipped.
- `'VCF ERROR: Population mismatch between VCF files {VCF_X} and {VCF_Y} - please ensure that the header lines are consistent across all VCF files!'`
	- Indicates that the two VCF files have different headers.
	- Solution: Ensure that all VCF files have the same exact header.
- `'VCF ERROR: Sample mismatch in population {POPULATION} between VCF files {VCF_X} and {VCF_Y} - please ensure that the header lines are consistent across all VCF files!'`
	- Indicates that the two VCF files have different headers.
	- Solution: Ensure that all VCF files have the same exact header.
- `'CONFIG ERROR: Configuration file {YAML_PATH} not found!'`
	- Indicates that no configuration file exists at the file path specified by the `YAML_PATH` parameter.
	- Solution: Ensure the `YAML_PATH` argument is set correctly.
- `'YAML ERROR: Error parsing YAML file - {ERROR}'`
	- Indicates an issue in reading the `YAML` configuration file.
	- Solution: Ensure that the `YAML` configuration file file is correctly formatted.
- `'CONFIG ERROR: Missing the following required parameters {MISSING_REQUIRED_PARAMS}!'`
	- Indicates that not all required parameters were specified in the `YAML` configuration file.
	- Solution: Ensure that the `YAML` configuration file file is correctly formatted.
- `'LOGGING ERROR: Failed to write to log file {LOG_FILE} - {ERROR}'`
	- Indicates an inability to write to the path specified by the `log_path_prefix` parameter.
	- Solution: Ensure that `log_path_prefix` specifies a writable path.
- `CONFIG ERROR: Configuration file validation failed!`
	- Indicates one or more errors were raised during validation.
	- Solution: Correct all errors in the log file, then rerun `RIPTA`.


### Warning Messages

- `'WARNING: Position {POSITION} on chromosome {CHROMOSOME} is not bi-allelic, continuing to the next site...'`
	- Indicates the presence of a non-bi-allelic SNP.
- `'WARNING: Position {CURRENT_POSITION} on chromosome {CURRENT_CHROMOSOME} does not have the RAF flag in the INFO columns, continuing to the next site...`
	- Indicates the absence of the RAF flag in the INFO field for a position when the `use_raf` parameter is set to `True`.
- `'WARNING: Position {POSITION} on chromosome {CHROMOSOME} does meet the RAF threshold, continuing to the next site...'`
	- Indicates the position does not meet the specified RAF threshold when the `use_raf` parameter is set to `True`.
- `'WARNING: Position {POSITION} on chromosome {CHROMOSOME} has no GP field for {INDIVIDUAL}...'`
	- Indicates missing GP field information at a position when the `ind_gp` sampling option is specified.
- `'WARNING: Position {POSITION} on chromosome {CHROMOSOME} has no genotype information for {INDIVIDUAL}...`
	- Indicates missing genotype information for an individual at a specified position.



