# Import dependicies.
import argparse
from collections import namedtuple
import gzip
import multiprocessing as mp
import numpy as np
import os
import yaml


# Define the site pattern configuration named tuple.
SitePatternConfig = namedtuple('SitePatternConfig', 
                              ['p1', 'p2', 'p3', 'p4', 
                               'p1_info', 'p2_info', 'p3_info', 'p4_info'])

# Define the site pattern block named tuple.
SitePatternBlock = namedtuple('SitePatternBlock',
                              ['abba', 'baba', 'baaa', 'abaa',
                               'abba_hom', 'baba_hom', 'baaa_hom', 'abaa_hom', 'n_seg'])

# Define a function to validate the file paths in the configuration file.
def validate_paths(config):
    '''Validates all file paths in the configuration file.'''
    # Intialize a list to store the errors.
    errors = []
    
    # For every VCF.
    for vcf in config['vcf_path']:
        # If the VCF path does not exsist.
        if not os.path.exists(vcf):
            # Update the errors.
            errors.append(f'CONFIG ERROR: VCF file {vcf} not found!')
            
    # If the meta data file does not exsist.
    if not os.path.exists(config['meta_data_path']):
        # Update the errors.
        errors.append(f"CONFIG ERROR: Meta data file {config['meta_data_path']} not found!")
        
    # If the compute_site_patterns flag is true.
    if config.get('compute_site_patterns'):
        # If the site pattern configuration file does not exsist.
        if not os.path.exists(config['site_pattern_config_path']):
            # Update the errors.
            errors.append(f"CONFIG ERROR: Site pattern configuration file {config['site_pattern_config_path']} not found!")
    return errors

# Define a function to validate the numeric parameters.
def validate_numerics(config):
    '''Validates the number of VCF files, contig lengths, and block sizes.'''
    # Intialize a list to store the errors.
    errors = []
    
    # If the the number contig lengths does not match the number of VCF files.
    if len(config['contig_length']) != len(config['vcf_path']):
        # Update the errors.
        errors.append(f"CONFIG ERROR: The vcf_path ({len(config['vcf_path'])}) and contig_length ({len(config['contig_length'])}) paramters must have the same number of entries!")
        
    # If all the block sizes are not smaller than the respective contigs.
    elif not all(config['block_size'] < contig for contig in config['contig_length']):
        # Update the errors.
        errors.append("CONFIG ERROR: The Block size must be smaller than the contig lengths!")
        
    # If the alternative allele frequency option for the refernce panel is set to true, make sure the value is between 0 and 1.
    if config['use_raf']:
        if (config['raf_threshold'] < 0) or (config['raf_threshold'] > 1):
            # Update the errors.
            errors.append("CONFIG ERROR: The RAF threshold must be between 0 and 1!")
            
    # If the number of threads is not an integer greater than 0 or max.
    if not (isinstance(config['number_of_threads'], int) and config['number_of_threads'] > 0 or config['number_of_threads'] == 'max'):
        # Update the errors.
        errors.append("CONFIG ERROR: The number of threads must be a positive integer or max!")
    return errors

# Define a function to validate the site pattern statistics.
def validate_sp_stats(config):
    '''Validates the specified site pattern statistics.'''
    # Intialize a list to store the errors.
    errors = []
    
    # If the compute_site_patterns parameter is set to true.
    if config.get('compute_site_patterns'):
        # Intialize the implemented site pattern statistics.
        valid_stats = {'D', 'Danc', 'D+', 'fhom', 'fanc', 'f+'}
        # If no site pattern statistics are listed.
        if not config.get('site_pattern_stats'):
            # Update the errors.
            errors.append("CONFIG ERROR: The site_pattern_stats paramter must be specified when the compute_site_patterns parameters is set to True!")
        # Else, site pattern statistics were listed.
        else:
            # Determine if there are any invalid site pattern statistics.
            invalid_stats = set(config['site_pattern_stats']) - valid_stats
            # If there are invalide site pattern statistics.
            if invalid_stats:
                # Update the errors.
                errors.append(f"CONFIG ERROR: Invalid site pattern statistics ({', '.join(invalid_stats)}) specified!")
    return errors

# Define a function to validate the site pattern configuration file.
def validate_sp_config(config):
    '''
    Validates the site pattern configuration file and returns required samples/populations.
    First four columns should be the site pattern configurations and the last four columns
    are the sampling options.
    
    Args:
        config (dict): Configuration dictionary.
        
    Returns:
        Tuple of (list of named tuples of site pattern configurations, set of required populations, set of required individuals, set of individuals to use GP, list of errors).
    '''
    # Initialize lists and sets to store validation results.
    errors = []
    required_pops = set()
    required_inds = set()
    gp_inds = set()
    # Initialize a list to store the site pattern configurations.
    sp_configs =[]
    
    
    # Only proceed if we're computing site patterns.
    if not config.get('compute_site_patterns'):
        return sp_configs, required_pops, required_inds, errors
        
    # Try to open the site pattern configuration file.
    try:
        # Open the file.
        with open(config['site_pattern_config_path'], 'rt') as sp_data:
            # Process each line in the file.
            for line_num, line in enumerate(sp_data, 1):
                # Split the line by tabs.
                spline = line.strip().split('\t')
                
                # If there is not exactly eight columns.
                if len(spline) != 8:
                    # Update the errors.
                    errors.append(f'SITE PATTERN CONFIG ERROR: Line {line_num} has {len(spline)} columns, expected 8!')
                    # Continue to the next line.
                    continue
                    
                # Unpack the site pattern configuration info.
                p1, p2, p3, p4, p1_info, p2_info, p3_info, p4_info = spline
                # Intialize the valid sampling options.
                valid_sampling = {'pop_freq', 'ind_freq', 'ind_gp'}
                # Check if all sampling options are valid.
                is_sp_config_valid = all(sampling in valid_sampling for sampling in [p1_info, p2_info, p3_info, p4_info])
                
                # If all sampling methods are valid.
                if is_sp_config_valid:
                    # Create a named tuple for this site pattern configuration.
                    sp_config = SitePatternConfig(
                        p1=p1, p2=p2, p3=p3, p4=p4,
                        p1_info=p1_info, p2_info=p2_info, 
                        p3_info=p3_info, p4_info=p4_info
                    )
                    # Update the site pattern configuration list.
                    sp_configs.append(sp_config)
                    # For each group and its sampling method.
                    for group, sampling in zip([p1, p2, p3, p4], [p1_info, p2_info, p3_info, p4_info]):
                        # Add the group to required populations or individuals based on sampling method.
                        if sampling == 'pop_freq':
                            required_pops.add(group)
                        elif sampling == 'ind_gp':
                            gp_inds.add(group)
                            required_inds.add(group)
                        else:
                            required_inds.add(group)
                # Else, at least one sampling method is invalid.
                else:
                    # For each group and its sampling method.
                    for i, (group, sampling) in enumerate(zip([p1, p2, p3, p4], [p1_info, p2_info, p3_info, p4_info])):
                        # If the group doesn't have a valid sampling option.
                        if sampling not in valid_sampling:
                            # Update the errors.
                            errors.append(f'SITE PATTERN CONFIG ERROR: Invalid sampling option "{sampling}" for P{i + 1} "{group}" on line {line_num}, the sampling option must be either {", ".join(valid_sampling)}!')
                 
    # Handle file reading errors.
    except Exception as error:
        # Update the errors.
        errors.append(f'SITE PATTERN CONFIG ERROR: Error reading site pattern configuration file - {str(error)}')
    return sp_configs, required_pops, required_inds, gp_inds, errors

# Define a function to validate the meta data.
def validate_meta_data(config, required_pops=None, required_inds=None):
    '''
    Parse population meta data file and create the population dictionary.
    First column should be individual IDs, second column should be population IDs, separated by tabs.
    
    Args:
        config (dict): Configuration dictionary.
        required_pops (set): Set of required populations to validate.
        required_inds (set): Set of required individuals to validate.
        
    Returns:
        Tuple of (population dictionary, list of errors).
    '''
    # Initialize storage for errors and population dictionary.
    errors = []
    pop_dicc = {}
    # Initialize sets for tracking found populations and individuals.
    found_pops = set()
    found_inds = set()
    
    # Try to open and parse the meta data file.
    try:
        # Open the meta data.
        with open(config['meta_data_path'], 'rt') as pop_data:
            # Process each line in the file.
            for line_num, line in enumerate(pop_data, 1):
                try:
                    # Split line by tabs.
                    spline = line.strip().split('\t')
                    # If there are less than two columns.
                    if len(spline) < 2:
                        # Update the erros.
                        errors.append(f'META DATA ERROR: Line {line_num} is invalid, expected at least 2 tab-separated columns!')
                        # Continue 
                        continue
                        
                    # Extract individual and population IDs.
                    ind, pop = spline[0], spline[1]
                    
                    # Check if this individual/population is required.
                    if ((required_inds is None or ind in required_inds) or 
                        (required_pops is None or pop in required_pops)):
                        # Initialize the population sub-dictionary if it doesn't already exsists.
                        if pop not in pop_dicc:
                            pop_dicc[pop] = {'ind': [], 'idx': []}
                        # Add individual to population.
                        pop_dicc[pop]['ind'].append(ind)
                        # Update tracking sets.
                        found_pops.add(pop)
                        found_inds.add(ind)
                        
                # Handle line parsing errors.
                except Exception as error:
                    # Update the errors.
                    errors.append(f'META DATA ERROR: Error parsing line {line_num} - {str(error)}')
                    
    # Handle file reading errors.
    except Exception as error:
        # Update the errors.
        errors.append(f'META DATA ERROR: Error reading meta data file - {str(error)}')
        
    # Validate that all required populations were found.
    if required_pops:
        # Check for missing populations.
        missing_pops = required_pops - found_pops
        # If there are missing populations.
        if missing_pops:
            # Update the errors.
            errors.append(f'META DATA ERROR: Required populations not found: {", ".join(missing_pops)}!')
            
    # Validate that all required individuals were found.
    if required_inds:
        # Check for missing individuals.
        missing_inds = required_inds - found_inds
        # If there are missing individuals.
        if missing_inds:
            # Update the errors.
            errors.append(f'META DATA ERROR: Required individuals not found: {", ".join(missing_inds)}')
            
    # Validate that that the meta dat was found.
    if not pop_dicc:
        # Update the errors.
        errors.append('META DATA ERROR: No samples found in meta data file, please check the errors and ensure the file has at least 2 tab-separated columns where column 1 = individual ID and column 2 = population ID!')
    return pop_dicc, errors

# Define a function to validate the VCF files.
def validate_vcfs(config, pop_dicc):
    '''
    Validate samples in VCF files against population metadata and the VCF file integrity.
    
    
    Args:
        config (dict): Configuration dictionary.
        pop_dicc (dict): Population dictionary.
        
    Returns:
        Tuple of (updated population dictionary, list of errors).
    '''
    # Initialize a list to store the errors and population information.
    errors = []
    pop_info_list = []
    
    # For each VCF file.
    for vcf in config['vcf_path']:
        # Try to open the VCF file.
        try:
            # Flexibly open the VCF file.
            opener = gzip.open if vcf.endswith('.gz') else open
            with opener(vcf, 'rt') as vcf_data:
                # Initialize a variable for finding the header line.
                header = None
                # Iterate through the VCF file line by line.
                for line in vcf_data:
                    # If this is the header line.
                    if line.startswith('#CHROM'):
                        # Split line by tabs.
                        header = line.strip().split('\t')
                        # Stop iterating through the VCF file.
                        break
                    # Else-if this is not a meta information line.
                    elif not line.startswith('##'):
                        # Stop iterating through the VCF file.
                        break
                        
                # If we didn't find the header.
                if not header:
                    # Update the errors.
                    errors.append(f'VCF ERROR: No header line found in the VCF file {vcf} - (please ensure the VCF file is formatted properly, for help see sections 1.5-1.6 https://samtools.github.io/hts-specs/VCFv4.3.pdf)!')
                # Else, the VCF has a header.
                else:
                    # Create a copy of population dictionary for this VCF.
                    vcf_pop_dicc = pop_dicc.copy()
                    # For every required population.
                    for pop, pop_info in vcf_pop_dicc.items():
                        # Copy the samples and initialize a list to store their respective indices.
                        pop_samples = pop_info['ind'].copy()
                        indices = []
                        # For every individual.
                        for ind in pop_samples:
                            # If the individual is in the header.
                            if ind in header:
                                # Append their index.
                                indices.append(header.index(ind))
                            # Else.
                            else:
                                # Remove the individual from the population.
                                pop_info['ind'].remove(ind)
                                # Update the errors.
                                errors.append(f'VCF ERROR: Sample {ind} from population {pop} not found in VCF header of {vcf}!')

                        # Update the dictionary.
                        pop_info['idx'] = indices
                        # If no indices were found.
                        if not indices and (required_pops is None or pop in required_pops):
                            # Update the errors.
                            errors.append(f'VCF ERROR: No samples found for required population {pop} in the VCF file {vcf}!')

                    # Store population info for this VCF.
                    pop_info_list.append((vcf, vcf_pop_dicc))
                
        # Handle file reading errors.        
        except Exception as error:
            # Update the errors.
            errors.append(f'VCF ERROR: Unable to open the VCF file {vcf} - {str(error)}')
            
    # If the population info list isn't empty.
    if pop_info_list:
        # Unpack the first tuple.
        first_vcf, first_pop_dicc = pop_info_list[0]
        # For all other vcf files (if they exist).
        for c_vcf, c_pop_dicc in pop_info_list[1:]:
            # If the populations are not consistent across VCF files.
            if set(first_pop_dicc.keys()) != set(c_pop_dicc.keys()):
                # Update the errors.
                errors.append(f'VCF ERROR: Population mismatch between VCF files {first_vcf} and {c_vcf} - please ensure that the header lines are consistent across all VCF files!')
                # Continue to the next population.
                continue
                
            # For every population.
            for pop in first_pop_dicc:
                # If both the individual IDs and their respective indices aren't the same.
                if (set(first_pop_dicc[pop]['ind']) != set(c_pop_dicc[pop]['ind']) or 
                    first_pop_dicc[pop]['idx'] != c_pop_dicc[pop]['idx']):
                    # Update the errors.
                    errors.append(f'VCF ERROR: Sample mismatch in population {pop} between VCF files {first_vcf} and {c_vcf} - please ensure that the header lines are consistent across all VCF files!')
                    
        # If all VCFs are consistent.
        if not any('mismatch' in error for error in errors):
            # Use the first VCF's population dictionary.
            pop_dicc = first_pop_dicc
    return pop_dicc, errors

# Define a function to read and validate the configuration file.
def read_config(config_path):
    '''
    Read and validate the YAML configuration file.
    
    
    Args:
        config_path (str): Path to the YAML configuration file.
        
    Returns:
        Contents of YAML file in Python dictionary.
        
    Raises:
        ValueError: If a configuration paramter is invalid.
        FileNotFoundError: If the configuration file doesn't exist.
    '''
    # Initialize a list to store all errors.
    all_errors = []
    
    # If the configuration file does not exist.
    if not os.path.exists(config_path):
        # Raise FileNotFoundError.
        raise FileNotFoundError(f"CONFIG ERROR: Configuration file {config_path} not found!")
        
    # Try to read and parse the YAML.
    try:
        # Open the YAML.
        with open(config_path, 'r') as file:
            # Initialize the configuration file
            config = yaml.safe_load(file)
    # Handle the YAML error.
    except yaml.YAMLError as error:
        # Raise ValueError.
        raise ValueError(f"YAML ERROR: Error parsing YAML file - {error}")
        
    # Run all validation checks and collect errors.
    all_errors.extend(validate_paths(config))
    all_errors.extend(validate_numerics(config))
    all_errors.extend(validate_sp_stats(config))
    
    # Initialize the set of required parameters.
    required_params = {
        'vcf_path', 'contig_length', 'block_size', 'meta_data_path',
        'results_path_prefix', 'log_path_prefix', 'use_raf',
        'raf_threshold', 'number_of_threads', 'compute_site_patterns', 
    }
    # Determine if there are any missing parameters.
    missing_params = required_params - set(config.keys())
    # If there are missing required parameters.
    if missing_params:
        # Update the errors.
        all_errors.append(f"CONFIG ERROR: Missing the following required parameters {', '.join(missing_params)}!")
    
    # For every file path key
    for key in ['vcf_path', 'meta_data_path', 'results_path_prefix', 
                'log_path_prefix', 'site_pattern_config_path']:
        # If the key is in the configuration.
        if key in config:
            # If the key's value is a list.
            if isinstance(config[key], list):
                # Normalize all the file paths.
                config[key] = [os.path.normpath(path) for path in config[key]]
            # Else-if the value is a single path and not None.
            elif config[key] is not None:
                # Normalize the file path.
                config[key] = os.path.normpath(config[key])
    
    # Get required populations and individuals from site pattern configurations.
    sp_configs, required_pops, required_inds, gp_inds, sp_errors = validate_sp_config(config)
    # Add any site pattern configuration errors.
    all_errors.extend(sp_errors)
    # Validate the meta data and VCF files.
    raw_pop_dicc, meta_data_errors = validate_meta_data(config, required_pops, required_inds)
    final_pop_dicc, vcf_errors = validate_vcfs(config, raw_pop_dicc)
    # Add meta data and VCF validation errors.
    all_errors.extend(meta_data_errors)
    all_errors.extend(vcf_errors)
    
    # Create the config error log file.
    log_file = f"{config['log_path_prefix']}_config_errors.log"
    # Try to write the log file.
    try:
        # Ensure the directory exists.
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        # Write all errors to the log file.
        with open(log_file, 'w') as log_data:
            for i, error in enumerate(all_errors):
                log_data.write(f'[{i}] {error}\n')
    
    # Handle file writting errors.
    except Exception as error:
        # Add error logging failure to the errors list
        all_errors.append(f'LOGGING ERROR: Failed to write to log file {log_file} - {str(error)}')
    # If there are any validation errors.
    if all_errors:
        # Raise ValueError with all error messages
        raise ValueError('CONFIG ERROR: Configuration file validation failed!\n' + 
                        '\n'.join(f'[{i}] {error}' for i, error in enumerate(all_errors)))
    
    # Initialize a dictionary to store the focal indicies.
    idx_dicc = {}
    
    # Add the required population indices to the dictionary.
    for pop in required_pops:
        if pop in final_pop_dicc:
            idx_dicc[pop] = final_pop_dicc[pop]['idx']
    
    # Add the required individual indices to the dictionary.
    for ind in required_inds:
        for pop in final_pop_dicc.values():
            if ind in pop['ind']:
                ind_idx = pop['idx'][pop['ind'].index(ind)]
                idx_dicc[ind] = [ind_idx]
                break
    
    # Update the configuration dictionary with the individuals to compute using GP, site pattern configurations and indicies.
    config['gp_inds'] = gp_inds
    config['sp_configs'] = sp_configs
    config['idx_dicc'] = idx_dicc
    return config

# Define a function to extract the genotype information for a block.
def extract_block_info(config, block):
    '''
    Extracts the allele counts and missingness information for a given block.
    
    Args:
        config (dict): Configuration dictionary.
        block (list): List of VCF file lines.
        
    Returns:
        Tuple of (dictionary of allele counts, dictionary of missing genotype information, list of warnings).
    '''
    # Intialize a list to store warnings.
    warnings = []
    # Intialize a dictionary to store the allele counts.
    ac_dicc = {key: [] for key in config['idx_dicc']}
    # Intialize a dictionary to store the misisng information.
    is_missing_dicc = {key: [] for key in config['idx_dicc']}

    # For every line in the block.
    for line in block:
        # Split the line by tabs.
        spline = line.strip().split('\t')
        # Grab the refernce and alternative alleles.
        alleles = [spline[3], spline[4]]

        # If the site is not bi-allelic in the VCF file.
        if '.' in alleles or (len(alleles[0]) + len(alleles[1])) != 2:
            # Update the warnings.
            warnings.append(f'WARNING: Position {spline[1]} on chromosome {spline[0]} is not bi-allelic, continuing to the next site...\n')
            # Continue to the next site.
            continue
        # Else-if we are supposed to use an RAF threshold.
        elif config['use_raf']:
            # Extract the info fieled.
            info = spline[7]
            # If the RAF field is not in the INFO columns.
            if 'RAF' not in info:
                # Update the warnings.
                warnings.append(f'WARNING: Position {spline[1]} on chromosome {spline[0]} does not have the RAF flag in the INFO columns, continuing to the next site...\n')
                # Continue to the next site.
                continue
            # Else the RAF field is in the INFO column.
            else:
                # Split the INFO column.
                spinfo = info.split(';')
                # Find the RAF flag index.
                raf_idx = [spinfo.index(flag) for flag in spinfo if flag.startswith('RAF')]
                # Determine if the RAF flag index was found.
                raf_idx = raf_idx[0] if len(raf_idx) == 1 else -1
                # If the RAF flag was not found or it outside the threshold.
                if raf_idx > -1 and (float(spinfo[raf_idx][4:]) < config['raf_threshold'] or float(spinfo[raf_idx][4:]) > 1 - config['raf_threshold']):
                    # Update the warnings.
                    warnings.append(f'WARNING: Position {spline[1]} on chromosome {spline[0]} does meet the RAF threshold, continuing to the next site...\n')
                    # Continue to the next site.
                    continue

        # For every population/individual.
        for key, indicies in config['idx_dicc'].items():
            # Intialize allele counters.
            rac, aac = 0, 0

            # If this is a GP individual.
            if key in config['gp_inds']:
                # Determine the GP flag index if it exists.
                gp_idx = spline[8].split(':').index('GP') if 'GP' in spline[8] else -1
                # If the GP flag was not found.
                if gp_idx == -1:
                    # Update the warnings.
                    warnings.append(f'WARNING: Position {spline[1]} on chromosome {spline[0]} has no GP field for {key}...\n')
                # Else, the GP flag was found.
                else:
                    # For every index.
                    for idx in indicies:
                        # If the genotype is not missing.
                        if not spline[idx][0] == '.':
                            # Extract the genotype probabilities.
                            prob_hom_ref, prob_het, prob_hom_alt = spline[idx].split(':')[gp_idx].split(',')
                            # Update the allele counts.
                            rac += (float(prob_hom_ref) * 2) + float(prob_het)
                            aac += (float(prob_hom_alt) * 2) + float(prob_het)
            # Else, we are not using GP.
            else:
                # For every index.
                for idx in indicies:
                    # Update the allele counts.
                    rac += spline[idx][0:3].count('0')
                    aac += spline[idx][0:3].count('1')

            # Update the allele count dictionary.
            ac_dicc[key].append([rac, aac])
            # If no alleles were detected.
            if rac + aac == 0:
                # Update the missingness dictionary.
                is_missing_dicc[key].append(True)
                # Update the warnings.
                warnings.append(f'WARNING: Position {spline[1]} on chromosome {spline[0]} has no genotype information for {key}...\n')
            # Else, alleles were detected.
            else:
                # Update the missingness dictionary.
                is_missing_dicc[key].append(False)

    # For every population/individual.
    for key in config['idx_dicc']:
        # Convert the lists to numpy arrays.
        ac_dicc[key] = np.array(ac_dicc[key])
        is_missing_dicc[key] = np.array(is_missing_dicc[key])
    return ac_dicc, is_missing_dicc, warnings

# Define a function for the site pattern block worker.
def sp_block_worker(config, file, lock, block_queue, results_queue):
    '''
    Retrieves a block from the producer queue annd places the results for the block in the consumer queue.
    
    Args:
        config (dict): Configuration dictionary.
        file (TextIO): An open text object for writting warning messages to.
        lock (multiprocessing.Lock): A lock object to synchronize access to resources across processes.
        block_queue (multiprocessing.Queue): The producer queue to communicate between processes or threads.
        results_queue (multiprocessing.Queue): The consumer queue to communicate between processes or threads.
    '''
    # Ensure that all avaiable workers are ready for the next block in the block queue.
    while True:
        # Grab the block information from the queue.
        block = block_queue.get()
        # If we get the work done signal kill the worker.
        if not block:
            return
        # Extract the allele count information, missingness, information, and warnings.
        ac_dicc, is_missing_dicc, warnings = extract_block_info(config, block)
        # If there are warnings to write.
        if warnings:
            # Write to the warnings file with thread safe access.
            with lock:
                file.writelines(warnings)
        # Intialize a dictionary to store the results.
        block_dicc = {}
        # If there is no information in the current block.
        if all(value.size == 0 for value in ac_dicc.values()) and all(value.size == 0 for value in is_missing_dicc.values()):
            # For every site pattern configuration.
            for sp_config in config['sp_configs']:
                # Store the results as a named tuple.
                block_dicc[sp_config] = SitePatternBlock(
                    abba=0, baba=0, baaa=0, abaa=0,
                    abba_hom=0, baba_hom=0, baaa_hom=0, abaa_hom=0, n_seg=0,
                )
        # Else there is information in the block to do computations on.
        else:
            # For every site pattern configuration.
            for sp_config in config['sp_configs']:
                # Determine where all samples have valid genotypes.
                is_called = ~is_missing_dicc[sp_config.p1] & ~is_missing_dicc[sp_config.p2] & ~is_missing_dicc[sp_config.p3] & ~is_missing_dicc[sp_config.p4]
                # Extract the allele count information for the ingroup.
                ingroup_ac = ac_dicc[sp_config.p1][is_called] + ac_dicc[sp_config.p2][is_called] + ac_dicc[sp_config.p3][is_called]
                # Determine what sites are segregating.
                is_seg = (ingroup_ac[:, 0] != ingroup_ac.sum(axis=1)) & (ingroup_ac[:, 1] != ingroup_ac.sum(axis=1))
                # If there are segregating sites to do computations on.
                if is_seg.sum() != 0:
                    # Extract the allele counts.
                    p1_ac, p2_ac, p3_ac, p4_ac = ac_dicc[sp_config.p1][is_called][is_seg], ac_dicc[sp_config.p2][is_called][is_seg], ac_dicc[sp_config.p3][is_called][is_seg], ac_dicc[sp_config.p4][is_called][is_seg]
                    # Compute the alternative allele frequencies for the outgroup.
                    p4_aaf = np.argmax(p4_ac, axis=1)
                    # Polarize the ingroup samples.
                    p1_daf = np.where(p4_aaf == 1, 1 - (p1_ac[:, 1] / np.round(p1_ac.sum(axis=1))), p1_ac[:, 1] / np.round(p1_ac.sum(axis=1)))
                    p2_daf = np.where(p4_aaf == 1, 1 - (p2_ac[:, 1] / np.round(p2_ac.sum(axis=1))), p2_ac[:, 1] / np.round(p2_ac.sum(axis=1)))
                    p3_daf = np.where(p4_aaf == 1, 1 - (p3_ac[:, 1] / np.round(p3_ac.sum(axis=1))), p3_ac[:, 1] / np.round(p3_ac.sum(axis=1)))
                    # Store the results as a named tuple.
                    block_dicc[sp_config] = SitePatternBlock(
                        abba=np.nansum((1 - p1_daf) * p2_daf * p3_daf),
                        baba=np.nansum(p1_daf * (1 - p2_daf) * p3_daf),
                        baaa=np.nansum(p1_daf * (1 - p2_daf) * (1 - p3_daf)),
                        abaa=np.nansum((1 - p1_daf) * p2_daf * (1 - p3_daf)),
                        abba_hom=np.nansum((1 - p1_daf) * p3_daf * p3_daf),
                        baba_hom=np.nansum(p1_daf * (1 - p3_daf) * p3_daf),
                        baaa_hom=np.nansum(p1_daf * (1 - p3_daf) * (1 - p3_daf)),
                        abaa_hom=np.nansum((1 - p1_daf) * p3_daf * (1 - p3_daf)),
                        n_seg=is_seg.sum(),
                    )
                # Else, there are no segregating sites to do computations on.
                else:
                    # Store the results as a named tuple.
                    block_dicc[sp_config] = SitePatternBlock(
                        abba=0, baba=0, baaa=0, abaa=0,
                        abba_hom=0, baba_hom=0, baaa_hom=0, abaa_hom=0, n_seg=0,
                    )
        # Put the results in the results queue.
        results_queue.put(block_dicc)

# Define a function to process the VCF file.
def consumer_producer_vcf_processor(config):
    '''
    Process VCF files in blocks using a consumer producer architecture.
    
    Args:
        config (dict): Configuration dictionary.
        
    Returns:
        List of named tuples with block information.
    '''
    # Create the warnings log file.
    log_file = f"{config['log_path_prefix']}_vcf_processing_warnings.log"
    # Create a manager for creating a queue in a shared memory space and open the log file.
    with mp.Manager() as manager, open(log_file, 'w', buffering=1) as file:
        # Create a shared memory queue for the blocks and initialize results list.
        block_queue = manager.Queue()
        results_queue = manager.Queue()
        # Create the lock for thread safe file writting.
        lock = manager.Lock()
        # Determine the number of threads.
        n_threads = mp.cpu_count() if config['number_of_threads'] == 'max' else config['number_of_threads']
        
        # Spawn the block workers.
        block_workers = [
            mp.Process(target=sp_block_worker, args=(config, file, lock, block_queue, results_queue))
            for _ in range(n_threads)
        ]
        
        # Allow the block workers to begin.
        for worker in block_workers:
            # Start the worker process.
            worker.start()
        
        # For each VCF file in the configuration.
        for vcf, contig_length, in zip(config['vcf_path'], config['contig_length']):
            # Initialize block tracking variables.
            c_block_lines = []
            c_block_idx = None
            
            # Flexibly open the VCF file.
            opener = gzip.open if vcf.endswith('.gz') else open
            with opener(vcf, 'rt') as vcf_data:
                # Process each line in the VCF file.
                for line in vcf_data:
                    # Skip header lines.
                    if line.startswith('#'):
                        continue
                    
                    # Grab the position from the current line.
                    pos = int(line.strip().split('\t', 2)[1])
                    # Calculate which block this position belongs to.
                    block_idx = (pos - 1) // config['block_size']
                    
                    # If this is our first position or it belongs to a new block.
                    if c_block_idx is None:
                        # Initialize the first block.
                        c_block_idx = block_idx
                        c_block_lines = [line]
                    # Else if this position belongs to the same block.
                    elif block_idx == c_block_idx:
                        # Add the line to the current block.
                        c_block_lines.append(line)
                    # Else this position belongs to a later block.
                    else:
                        # Put the completed block in the queue.
                        block_queue.put(c_block_lines)
                        # Start the new block.
                        c_block_idx = block_idx
                        c_block_lines = [line]
                
                # If we have a final block to process.
                if c_block_lines:
                    # Put the final block in the queue.
                    block_queue.put(c_block_lines)
        
        # Send termination signal to workers.
        for _ in range(n_threads):
            # Put None in queue to signal completion.
            block_queue.put(None)
            
        # Wait for all workers to finish.
        for worker in block_workers:
            # Join the worker process.
            worker.join()
        
        # Intialize a list to store all of the blocked results.
        blocks = []
        # While thera are results to collect.
        while not results_queue.empty():
            # Collect the results.
            blocks.append(results_queue.get())
            
        return blocks

# Define a function to compute the standard error and corresponding Z-score from a weighted block jackknife procedure.
def weighted_block_jackknife_procedure(theta, n, theta_js,  m_js):
    '''
    Weighted block jackknife procedure to assess genome-wide significance.
    
    
    Args:
        theta (float): The observed genome-wide estimator.
        n (int): The total number of informative SNPs genome wide.
        theta_js (numpy.ndarray): Array of psuedovalues per block.
        m_js (numpy.ndarray): Array of weights associated with each psuedovalue per block.
    
    Returns:
        Tuple of (standard error, Z-score, and number of blocks with informative SNPs).
    '''
    # Determine what blocks have defined values.
    is_defined = ~np.isnan(theta_js)
    # Mask the blocks with undefined values.
    theta_js, m_js = theta_js[is_defined], m_js[is_defined]
    # Determine the number blocks.
    g = theta_js.size
    # If there are no valid blocks.
    if g == 0:
        return 'NULL', 'NULL', 0
    # Else there are valid blocks.
    else:
        # Computed the weighted squared differences.
        weighted_squared_diffs = np.fromiter(
            ((m_j / (n - m_j)) * ((theta - theta_j) ** 2) for theta_j, m_j in zip(theta_js, m_js)),
            dtype=np.float64, count=g,
        )
        # Compute the standard error of the blocked jackknife distribution.
        sigma = np.sqrt((1 / g) * np.sum(weighted_squared_diffs))
        return sigma, theta / sigma, g

# Define a function to compile the results.
def compile_results(config, blocks):
    '''
    Computes and writes the results to a CSV file
    
    
    Args:
        config (dict): Configuration dictionary.
        blocks (list): List of named tuples with block information, from consumer_producer_vcf_processor.
    '''
    # Intialize a dictionary to store the genomewide site pattern counts.
    gw_info = {
        sp_config: {'abba': 0, 'baba': 0, 'baaa': 0, 'abaa': 0, 'abba_hom': 0, 'baba_hom': 0, 'baaa_hom': 0, 'abaa_hom': 0, 'n_seg': 0} for sp_config in config['sp_configs']
    }
    # For every block.
    for block in blocks:
        # For all the configurations in this block.
        for key, values in block.items():
            # Update the genome wide values.
            gw_info[key]['abba'] += values.abba
            gw_info[key]['baba'] += values.baba
            gw_info[key]['baaa'] += values.baaa
            gw_info[key]['abaa'] += values.abaa
            gw_info[key]['abba_hom'] += values.abba_hom
            gw_info[key]['baba_hom'] += values.baba_hom
            gw_info[key]['baaa_hom'] += values.baaa_hom
            gw_info[key]['abaa_hom'] += values.abaa_hom
            gw_info[key]['n_seg'] += values.n_seg

    # Intialize a dictionary to store the jackknife information.
    jackknife_info = {sp_config: {'theta_js': [], 'm_js': []} for sp_config in config['sp_configs']}
    # For every block.
    for block in blocks:
        # For all the configurations in this block.
        for key, values in block.items():
            # If there are informative sites.
            if values.n_seg != 0:
                # Update the list with the jackknife information.
                jackknife_info[key]['theta_js'].append(SitePatternBlock(
                    abba=gw_info[key]['abba'] - values.abba,
                    baba=gw_info[key]['abba'] - values.baba,
                    baaa=gw_info[key]['abba'] - values.baaa,
                    abaa=gw_info[key]['abba'] - values.abaa,
                    abba_hom=gw_info[key]['abba_hom'] - values.abba_hom,
                    baba_hom=gw_info[key]['abba_hom'] - values.baba_hom,
                    baaa_hom=gw_info[key]['abba_hom'] - values.baaa_hom,
                    abaa_hom=gw_info[key]['abba_hom'] - values.abaa_hom,
                    n_seg=0,
                    ))
                jackknife_info[key]['m_js'].append(values.n_seg)

    # Intialize dummy variables to maintain the correct order when writting the output.
    write_D = False
    write_Danc = False
    write_Dplus = False
    write_fhom = False
    write_fanc = False
    write_fplus = False

    # If we are calculating Patterson's D.
    if 'D' in config['site_pattern_stats']:
        # Set the dummy variable to true.
        write_D = True
        # For every site pattern configuration.
        for sp_config in config['sp_configs']:
            # Compute the numerator and the denominator.
            numer = gw_info[sp_config]['abba'] - gw_info[sp_config]['baba']
            denom = gw_info[sp_config]['abba'] + gw_info[sp_config]['baba']
            # If the genome-wide value is undefined.
            if denom == 0:
                # Update the results.
                gw_info[sp_config]['D'] = 'NULL'
                gw_info[sp_config]['D_SE'] = 'NULL'
                gw_info[sp_config]['D_Z'] = 'NULL'
                gw_info[sp_config]['D_BLOCKS'] = 'NULL'
            # Else, there are values to do computations on.
            else:
                # Compute the observed value.
                D_theta = numer / denom
                # Compile the jackknife values.
                D_theta_js = np.array([
                    (theta_j.abba - theta_j.baba) / (theta_j.abba + theta_j.baba)
                    for theta_j in jackknife_info[sp_config]['theta_js']
                ])
                # Perform the weighted jackknife procedure.
                D_se, D_z, D_blocks = weighted_block_jackknife_procedure(
                    D_theta, gw_info[sp_config]['n_seg'], D_theta_js, np.array(jackknife_info[sp_config]['m_js']),
                )
                # Update the results.
                gw_info[sp_config]['D'] = D_theta
                gw_info[sp_config]['D_SE'] = D_se
                gw_info[sp_config]['D_Z'] = D_z
                gw_info[sp_config]['D_BLOCKS'] = D_blocks

    # If we are calculating D_anc.
    if 'Danc' in config['site_pattern_stats']:
        # Set the dummy variable to true.
        write_Danc = True
        # For every site pattern configuration.
        for sp_config in config['sp_configs']:
            # Compute the numerator and the denominator.
            numer = gw_info[sp_config]['baaa'] - gw_info[sp_config]['abaa']
            denom = gw_info[sp_config]['baaa'] + gw_info[sp_config]['abaa']
            # If the genome-wide value is undefined.
            if denom == 0:
                # Update the results.
                gw_info[sp_config]['Danc'] = 'NULL'
                gw_info[sp_config]['Danc_SE'] = 'NULL'
                gw_info[sp_config]['Danc_Z'] = 'NULL'
                gw_info[sp_config]['Danc_BLOCKS'] = 'NULL'
            # Else, there are values to do computations on.
            else:
                # Compute the observed value.
                Danc_theta = numer / denom
                # Compile the jackknife values.
                Danc_theta_js = np.array([
                    (theta_j.baaa - theta_j.abaa) / (theta_j.baaa + theta_j.abaa)
                    for theta_j in jackknife_info[sp_config]['theta_js']
                ])
                # Perform the weighted jackknife procedure.
                Danc_se, Danc_z, Danc_blocks = weighted_block_jackknife_procedure(
                    Danc_theta, gw_info[sp_config]['n_seg'], Danc_theta_js, np.array(jackknife_info[sp_config]['m_js']),
                )
                # Update the results.
                gw_info[sp_config]['Danc'] = Danc_theta
                gw_info[sp_config]['Danc_SE'] = Danc_se
                gw_info[sp_config]['Danc_Z'] = Danc_z
                gw_info[sp_config]['Danc_BLOCKS'] = Danc_blocks

    # If we are calculating D+.
    if 'D+' in config['site_pattern_stats']:
        # Set the dummy variable to true.
        write_Dplus = True
        # For every site pattern configuration.
        for sp_config in config['sp_configs']:
            # Compute the numerator and the denominator.
            numer = (gw_info[sp_config]['abba'] - gw_info[sp_config]['baba']) + (gw_info[sp_config]['baaa'] - gw_info[sp_config]['abaa'])
            denom = gw_info[sp_config]['abba'] + gw_info[sp_config]['baba'] + gw_info[sp_config]['baaa'] + gw_info[sp_config]['abaa']
            # If the genome-wide value is undefined.
            if denom == 0:
                # Update the results.
                gw_info[sp_config]['D+'] = 'NULL'
                gw_info[sp_config]['D+_SE'] = 'NULL'
                gw_info[sp_config]['D+_Z'] = 'NULL'
                gw_info[sp_config]['D+_BLOCKS'] = 'NULL'
            # Else, there are values to do computations on.
            else:
                # Compute the observed value.
                Dplus_theta = numer / denom
                # Compile the jackknife values.
                Dplus_theta_js = np.array([
                    ((theta_j.abba - theta_j.baba) + (theta_j.baaa - theta_j.abaa)) / (theta_j.abba + theta_j.baba + theta_j.baaa + theta_j.abaa)
                    for theta_j in jackknife_info[sp_config]['theta_js']
                ])
                # Perform the weighted jackknife procedure.
                Dplus_se, Dplus_z, Dplus_blocks = weighted_block_jackknife_procedure(
                    Dplus_theta, gw_info[sp_config]['n_seg'], Dplus_theta_js, np.array(jackknife_info[sp_config]['m_js']),
                )
                # Update the results.
                gw_info[sp_config]['D+'] = Dplus_theta
                gw_info[sp_config]['D+_SE'] = Dplus_se
                gw_info[sp_config]['D+_Z'] = Dplus_z
                gw_info[sp_config]['D+_BLOCKS'] = Dplus_blocks
    
    # If we are calculating f_hom.
    if 'fhom' in config['site_pattern_stats']:
        # Set the dummy variable to true.
        write_fhom = True
        # For every site pattern configuration.
        for sp_config in config['sp_configs']:
            # Compute the numerator and the denominator.
            numer = gw_info[sp_config]['abba'] - gw_info[sp_config]['baba']
            denom = gw_info[sp_config]['abba_hom'] - gw_info[sp_config]['baba_hom']
            # If the genome-wide value is undefined.
            if denom == 0:
                # Update the results.
                gw_info[sp_config]['fhom'] = 'NULL'
                gw_info[sp_config]['fhom_SE'] = 'NULL'
                gw_info[sp_config]['fhom_Z'] = 'NULL'
                gw_info[sp_config]['fhom_BLOCKS'] = 'NULL'
            # Else, there are values to do computations on.
            else:
                # Compute the observed value.
                fhom_theta = numer / denom
                # Compile the jackknife values.
                fhom_theta_js = np.array([
                    (theta_j.abba - theta_j.baba) / (theta_j.abba_hom - theta_j.baba_hom)
                    for theta_j in jackknife_info[sp_config]['theta_js']
                ])
                # Perform the weighted jackknife procedure.
                fhom_se, fhom_z, fhom_blocks = weighted_block_jackknife_procedure(
                    fhom_theta, gw_info[sp_config]['n_seg'], fhom_theta_js, np.array(jackknife_info[sp_config]['m_js']),
                )
                # Update the results.
                gw_info[sp_config]['fhom'] = fhom_theta
                gw_info[sp_config]['fhom_SE'] = fhom_se
                gw_info[sp_config]['fhom_Z'] = fhom_z
                gw_info[sp_config]['fhom_BLOCKS'] = fhom_blocks
    
    # If we are calculating f_anc.
    if 'fanc' in config['site_pattern_stats']:
        # Set the dummy variable to true.
        write_fanc = True
        # For every site pattern configuration.
        for sp_config in config['sp_configs']:
            # Compute the numerator and the denominator.
            numer = gw_info[sp_config]['baaa'] - gw_info[sp_config]['abaa']
            denom = gw_info[sp_config]['baaa_hom'] - gw_info[sp_config]['abaa_hom']
            # If the genome-wide value is undefined.
            if denom == 0:
                # Update the results.
                gw_info[sp_config]['fanc'] = 'NULL'
                gw_info[sp_config]['fanc_SE'] = 'NULL'
                gw_info[sp_config]['fanc_Z'] = 'NULL'
                gw_info[sp_config]['fanc_BLOCKS'] = 'NULL'
            # Else, there are values to do computations on.
            else:
                # Compute the observed value.
                fanc_theta = numer / denom
                # Compile the jackknife values.
                fanc_theta_js = np.array([
                    (theta_j.baaa - theta_j.abaa) / (theta_j.baaa_hom + theta_j.abaa_hom)
                    for theta_j in jackknife_info[sp_config]['theta_js']
                ])
                # Perform the weighted jackknife procedure.
                fanc_se, fanc_z, fanc_blocks = weighted_block_jackknife_procedure(
                    fanc_theta, gw_info[sp_config]['n_seg'], fanc_theta_js, np.array(jackknife_info[sp_config]['m_js']),
                )
                # Update the results.
                gw_info[sp_config]['fanc'] = fanc_theta
                gw_info[sp_config]['fanc_SE'] = fanc_se
                gw_info[sp_config]['fanc_Z'] = fanc_z
                gw_info[sp_config]['fanc_BLOCKS'] = fanc_blocks
                
    # If we are calculating f+.
    if 'f+' in config['site_pattern_stats']:
        # Set the dummy variable to true.
        write_fplus = True
        # For every site pattern configuration.
        for sp_config in config['sp_configs']:
            # Compute the numerator and the denominator.
            numer = (gw_info[sp_config]['abba'] - gw_info[sp_config]['baba']) + (gw_info[sp_config]['baaa'] - gw_info[sp_config]['abaa'])
            denom = (gw_info[sp_config]['abba_hom'] - gw_info[sp_config]['baba_hom']) + (gw_info[sp_config]['baaa_hom'] - gw_info[sp_config]['abaa_hom'])
            # If the genome-wide value is undefined.
            if denom == 0:
                # Update the results.
                gw_info[sp_config]['f+'] = 'NULL'
                gw_info[sp_config]['f+_SE'] = 'NULL'
                gw_info[sp_config]['f+_Z'] = 'NULL'
                gw_info[sp_config]['f+_BLOCKS'] = 'NULL'
            # Else, there are values to do computations on.
            else:
                # Compute the observed value.
                fplus_theta = numer / denom
                # Compile the jackknife values.
                fplus_theta_js = np.array([
                    ((theta_j.abba - theta_j.baba) + (theta_j.baaa - theta_j.abaa)) / ((theta_j.abba_hom - theta_j.baba_hom) + (theta_j.baaa_hom - theta_j.abaa_hom))
                    for theta_j in jackknife_info[sp_config]['theta_js']
                ])
                # Perform the weighted jackknife procedure.
                fplus_se, fplus_z, fplus_blocks = weighted_block_jackknife_procedure(
                    fplus_theta, gw_info[sp_config]['n_seg'], fplus_theta_js, np.array(jackknife_info[sp_config]['m_js']),
                )
                # Update the results.
                gw_info[sp_config]['f+'] = fplus_theta
                gw_info[sp_config]['f+_SE'] = fplus_se
                gw_info[sp_config]['f+_Z'] = fplus_z
                gw_info[sp_config]['f+_BLOCKS'] = fplus_blocks

    # Create the site pattern results file.
    results_file = f"{config['results_path_prefix']}_site_patterns.csv"
    # Open the results file.
    with open(results_file, 'w') as results_data:
        # Intialize a header list.
        header_list = [
            'P1', 'P2', 'P3', 'P4',
            'ABBA', 'BABA', 'BAAA', 'ABAA',
            'ABBA_HOM', 'BABA_HOM', 'BAAA_HOM', 'ABAA_HOM', 'N_SEG'
        ]
        # Add the specified site pattern statistics to the header.
        if write_D:
            header_list.extend(['D', 'D_SE', 'D_Z', 'D_BLOCKS'])
        if write_Danc:
            header_list.extend(['Danc', 'Danc_SE', 'Danc_Z', 'Danc_BLOCKS'])
        if write_Dplus:
            header_list.extend(['D+', 'D+_SE', 'D+_Z', 'D+_BLOCKS'])
        if write_fhom:
            header_list.extend(['fhom', 'fhom_SE', 'fhom_Z', 'fhom_BLOCKS'])
        if write_fanc:
            header_list.extend(['fanc', 'fanc_SE', 'fanc_Z', 'fanc_BLOCKS'])
        if write_fplus:
            header_list.extend(['f+', 'f+_SE', 'f+_Z', 'f+_BLOCKS'])
        # Write the header list to the results file.
        results_data.write(','.join(header_list)+'\n')
        # For every site pattern configuration.
        for sp_config in gw_info:
            # Intialize a list to store the results.
            sp_config_list = [
                sp_config.p1, sp_config.p2, sp_config.p3, sp_config.p4,
                str(gw_info[sp_config]['abba']), str(gw_info[sp_config]['baba']),
                str(gw_info[sp_config]['baaa']), str(gw_info[sp_config]['abaa']),
                str(gw_info[sp_config]['abba_hom']), str(gw_info[sp_config]['baba_hom']),
                str(gw_info[sp_config]['baaa_hom']), str(gw_info[sp_config]['abaa_hom']),
                str(gw_info[sp_config]['n_seg']),
            ]
            # Add the specified site pattern statistics.
            if write_D:
                sp_config_list.extend([
                    str(gw_info[sp_config]['D']), str(gw_info[sp_config]['D_SE']),
                    str(gw_info[sp_config]['D_Z']), str(gw_info[sp_config]['D_BLOCKS']),
                ])
            if write_Danc:
                sp_config_list.extend([
                    str(gw_info[sp_config]['Danc']), str(gw_info[sp_config]['Danc_SE']),
                    str(gw_info[sp_config]['Danc_Z']), str(gw_info[sp_config]['Danc_BLOCKS']),
                ])
            if write_Dplus:
                sp_config_list.extend([
                    str(gw_info[sp_config]['D+']), str(gw_info[sp_config]['D+_SE']),
                    str(gw_info[sp_config]['D+_Z']), str(gw_info[sp_config]['D+_BLOCKS']),
                ])
            if write_fhom:
                sp_config_list.extend([
                    str(gw_info[sp_config]['fhom']), str(gw_info[sp_config]['fhom_SE']),
                    str(gw_info[sp_config]['fhom_Z']), str(gw_info[sp_config]['fhom_BLOCKS']),
                ])
            if write_fanc:
                sp_config_list.extend([
                    str(gw_info[sp_config]['fanc']), str(gw_info[sp_config]['fanc_SE']),
                    str(gw_info[sp_config]['fanc_Z']), str(gw_info[sp_config]['fanc_BLOCKS']),
                ])
            if write_fplus:
                sp_config_list.extend([
                    str(gw_info[sp_config]['f+']), str(gw_info[sp_config]['f+_SE']),
                    str(gw_info[sp_config]['f+_Z']), str(gw_info[sp_config]['f+_BLOCKS']),
                ])
            # Write the results.
            results_data.write(','.join(sp_config_list)+'\n')
    return

# Intialize the command line options.
parser = argparse.ArgumentParser()
parser.add_argument(
    '-y', '--yaml_path', required=True, type=str, action='store',
    help='Path to the YAML configuration file.',
)
args = parser.parse_args()

# Load and validate the configuration file.
config = read_config(args.yaml_path)
# Extract the block information.
blocks = consumer_producer_vcf_processor(config)
# Compile the results.
compile_results(config, blocks)