# Import packages.
import argparse
import gzip
import random
import sys

# Intialize the command line options.
parser = argparse.ArgumentParser()
parser.add_argument(
    '-v', '--vcf_file', required=True,
    type=str, action='store',
    help='Path to the input vcf file.',
)
parser.add_argument(
    '-m', '--meta_data', required=True,
    type=str, action='store',
    help='Path to the tab delimited text file where the first column contains'\
    +' the sample name as it appears in the header of the input vcf file and the'\
    +' second column contains the population the sample belongs to.',
)
parser.add_argument(
    '-f', '--frequency', choices=('True','False'),
    required=True, action='store', type=str,
    help='If true calculate site patterns from derived allele frequencies,'\
    +' if false calculate site patterns by randomly sampling one chromosome.',
)
parser.add_argument(
    '-p1', '--p1_population', required=True,
    type=str, action='store',
    help='Population in the meta data file to use as P1 for comparisons'\
    +' (ie a potential recipient population).',
)
parser.add_argument(
    '-p2', '--p2_population', required=True,
    type=str, action='store',
    help='Population in the meta data file to use as P2 for comparisons'\
    +' (ie a potential recipient population).',
)
parser.add_argument(
    '-p3', '--p3_population', required=True,
    type=str, action='store',
    help='Population in the meta data file to use as P3 for comparisons'\
    +' (ie the source population).',
)
parser.add_argument(
    '-p4', '--p4_population', required=True,
    type=str, action='store',
    help='Population in the meta data file to use as P3 for comparisons'\
    +' (ie the outrgroup population used for polarization).',
)
parser.add_argument(
    '-cl', '--contig_length', required=True,
    type=int, action='store',
    help='Total length in base pairs of conting.'
)
parser.add_argument(
    '-bs', '--block_size', required=True,
    type=int, action='store',
    help='Block size to be sampled with replacemnt for builidng bootstrapped'\
    +' contigs.',
)
parser.add_argument(
    '-r', '--replicates', required=True,
    type=int, action='store',
    help='Number of bootstrapped replicates to perform.'\
)
parser.add_argument(
    '-p', '--path', required=False,
    action='store', default='./',
    help='Path for results and qc file (default = current working directory).',
)
args = parser.parse_args()

# [0] Intialize output files.

# Intailzie the output file prefix.
prefix = args.path+'{0}_{1}_{2}_{3}_'.format(args.p1_population, args.p2_population, args.p3_population, args.p4_population)
# Intialize the output files.
out_file = open(prefix+'bootstraps.txt', 'w')
log_file = open(prefix+'bootstraps_log.txt', 'w')


# [1] Extract the meta data.

# Intialize an empty dictionary to store the focal population meta data.
pop_dicc = {
    args.p1_population: {'IND': [], 'IDX': []}, args.p2_population: {'IND': [], 'IDX': []},
    args.p3_population: {'IND': [], 'IDX': []}, args.p4_population: {'IND': [], 'IDX': []},
}
# Open the meta data file.
with open(args.meta_data, 'r') as pop_data:
    # For ever line in the meta data file...
    for line in pop_data:
        # Split the line by tabs.
        spline = line.split()
        # Grab the current sample and its corresponding population.
        ind = spline[0]
        pop = spline[1]
        # If the current population is a focal population.
        if pop in list(pop_dicc.keys()):
            # Append the sample to the population dictionary.
            pop_dicc[pop]['IND'].append(ind)
        # Else...
        else:
            # Compose the error message.
            err = 'QC: {0} from {1} does not belong to a focal population...'.format(ind, pop)
            # Log the error.
            log_file.write(err+'\n')


# [2] Estimate site patterns.

# Initialize the freq_flag.
freq_flag = args.frequency == 'True'
# If the file is gzipped...
if args.vcf_file.endswith('.gz'):
    # Intialize the gzipped vcf file.
    input_file = gzip.open(args.vcf_file, 'rt')
# Else...
else:
    # Intialize the vcf file.
    input_file = open(args.vcf_file, 'rt')
# Open the input file.
with input_file as data:
    # For every line in the vcf file...
    for line in data:
        # If the current line is a part of the meta information...
        if line.startswith('##'):
            # Continue to the next line in the vcf file.
            continue
        # Else-if the current line is the header line...
        elif line.startswith('#'):
            # Intialize the fatal error tracker.
            fatal = False
            # Split the header line by tabs.
            header = line.split()
            # For every focal population...
            for key in pop_dicc.keys():
                # For every sample in the current population...
                for ind in pop_dicc[key]['IND']:
                    # If the sample appears in the header line...
                    if ind in header:
                        # Append the sample's column index w/in the vcf file.
                        pop_dicc[key]['IDX'].append(header.index(ind))
                    # Else...
                    else:
                        # Remove the sample from the sample list.
                        pop_dicc[key]['IND'].remove(ind)
                        # Compose the error message.
                        err = 'QC: {0} from {1} does not appear in the vcf header...'.format(ind, key)
                        # Log the error.
                        log_file.write(err+'\n')
                # If the current population has no samples in the vcf file...
                if len(pop_dicc[key]['IDX']) == 0:
                    # Compose the error message.
                    err = 'FATAL: no {0} samples appear in the vcf header...'.format(key)\
                    +' make sure samples in the meta data file are identical to how they'\
                    +' appear in the vcf header...'
                    # Log the error.
                    log_file.write(err+'\n')
                    # Record the fatal error.
                    fatal = True
            # If a fatal error occured...
            if fatal:
                # Stop all calculations.
                break
            # Else...
            else:
                # Intialize a dictinary to store site patterns by position.
                site_patterns = {}
        # Else...
        else:
            # Split the header line by tabs.
            spline = line.split()
            # Grab the refernce and alternative alleles.
            alleles = [spline[3], spline[4]]
            # If the site is monomorphic...
            if '.' in alleles:
                # Grab the chromosome and position.
                chrom = spline[0]
                pos = spline[1]
                # Compose the error message.
                err = 'QC: pos {0} on chrom {1} is invariant, continuing to the next site...'.format(pos, chrom)
                # Log the error.
                log_file.write(err+'\n')
                # Continue to the next line...
                continue
            # Else-if the site contains a structural or multiallelic variant...
            elif (len(alleles[0]) + len(alleles[1])) != 2:
                # Grab the chromosome and position.
                chrom = spline[0]
                pos = spline[1]
                # Compose the error message.
                err = 'QC: pos {0} on chrom {1} is not a bi-allelic snp, continuing to the next site...'.format(pos, chrom)
                # Log the error.
                log_file.write(err+'\n')
                # Continue to the next line...
                continue
            # Else-if site patterns are to be calculated from derived allele frequencies...
            elif freq_flag:
                # Grab the poistion.
                pos = int(spline[1])
                # Intialize a dictionary to store site patterns for this position.
                site_patterns[pos] = {
                    'ABBA': 0, 'ABBA_HOM': 0,
                    'BABA': 0, 'BABA_HOM': 0,
                    'BAAA': 0, 'BAAA_HOM': 0,
                    'ABAA': 0, 'ABAA_HOM': 0,
                }
                # Intialize a alternative allele frequency dictionary.
                freq_dicc = {}
                # For every focal population...
                for key in pop_dicc.keys():
                    # Intialize an alternative allele counter.
                    alt_allele_counter = 0
                    # For every individual in the population...
                    for idx in pop_dicc[key]['IDX']:
                        # Count the number of alternative alleles.
                        alt_allele_counter += spline[idx][0:3].count('1')
                    # Determine the alternative allele frequency.
                    freq_dicc[key] = float(alt_allele_counter) / (len(pop_dicc[key]['IDX']) * 2)
                # If the ancestral allele is the refernce allele...
                if freq_dicc[args.p4_population] == 0.0:
                    # Calculate site patterns.
                    site_patterns[pos]['ABBA'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p2_population] * freq_dicc[args.p3_population]
                    site_patterns[pos]['ABBA_HOM'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p3_population] * freq_dicc[args.p3_population]
                    site_patterns[pos]['BABA'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p2_population]) * freq_dicc[args.p3_population]
                    site_patterns[pos]['BABA_HOM'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p3_population]) * freq_dicc[args.p3_population]
                    site_patterns[pos]['BAAA'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p2_population]) * (1 - freq_dicc[args.p3_population])
                    site_patterns[pos]['BAAA_HOM'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p3_population]) * (1 - freq_dicc[args.p3_population])
                    site_patterns[pos]['ABAA'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p2_population] * (1 - freq_dicc[args.p3_population])
                    site_patterns[pos]['ABAA_HOM'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p3_population] * (1 - freq_dicc[args.p3_population])
                # Else-if the ancestral allele is the alternative allele...
                elif freq_dicc[args.p4_population] == 1.0:
                    # Calculate site patterns.
                    site_patterns[pos]['ABBA'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p2_population]) * (1 - freq_dicc[args.p3_population])
                    site_patterns[pos]['ABBA_HOM'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p3_population]) * (1 - freq_dicc[args.p3_population])
                    site_patterns[pos]['BABA'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p2_population] * (1 - freq_dicc[args.p3_population])
                    site_patterns[pos]['BABA_HOM'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p3_population] * (1 - freq_dicc[args.p3_population])
                    site_patterns[pos]['BAAA'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p2_population] * freq_dicc[args.p3_population]
                    site_patterns[pos]['BAAA_HOM'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p3_population] * freq_dicc[args.p3_population]
                    site_patterns[pos]['ABAA'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p2_population]) * freq_dicc[args.p3_population]
                    site_patterns[pos]['ABAA_HOM'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p3_population]) * freq_dicc[args.p3_population]
                # Else...
                else:
                    # Calculate site patterns.
                    site_patterns[pos]['ABBA'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p2_population] * freq_dicc[args.p3_population] * (1 - freq_dicc[args.p4_population])
                    site_patterns[pos]['ABBA_HOM'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p3_population] * freq_dicc[args.p3_population] * (1 - freq_dicc[args.p4_population])
                    site_patterns[pos]['BABA'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p2_population]) * freq_dicc[args.p3_population] * (1 - freq_dicc[args.p4_population])
                    site_patterns[pos]['BABA_HOM'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p3_population]) * freq_dicc[args.p3_population] * (1 - freq_dicc[args.p4_population])
                    site_patterns[pos]['BAAA'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p2_population]) * (1 - freq_dicc[args.p3_population]) * (1 - freq_dicc[args.p4_population])
                    site_patterns[pos]['BAAA_HOM'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p3_population]) * (1 - freq_dicc[args.p3_population]) * (1 - freq_dicc[args.p4_population])
                    site_patterns[pos]['ABAA'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p2_population] * (1 - freq_dicc[args.p3_population]) * (1 - freq_dicc[args.p4_population])
                    site_patterns[pos]['ABAA_HOM'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p3_population] * (1 - freq_dicc[args.p3_population]) * (1 - freq_dicc[args.p4_population])
            # Else...
            else:
                # Grab the poistion.
                pos = int(spline[1])
                # Build the dictionary.
                site_patterns[pos] = {}
                # For every P1 sample...
                for p1 in pop_dicc[args.p1_population]['IND']:
                    # For every P2 sample...
                    for p2 in pop_dicc[args.p2_population]['IND']:
                        # For every P3 sample...
                        for p3 in pop_dicc[args.p3_population]['IND']:
                            # For every P4 sample...
                            for p4 in pop_dicc[args.p4_population]['IND']:
                                # Intialize the quartet.
                                quartet = '{0}-{1}-{2}-{3}'.format(p1, p2, p3, p4)
                                # Fill the dictionary to store site patterns.
                                site_patterns[pos][quartet] = {
                                    'ABBA': 0, 'ABBA_HOM': 0,
                                    'BABA': 0, 'BABA_HOM': 0,
                                    'BAAA': 0, 'BAAA_HOM': 0,
                                    'ABAA': 0, 'ABAA_HOM': 0,
                                }
                # Intialize RNG values.
                rng_vals = [0, 2]
                # Randomly select a chromosome to sample.
                random_chr = random.choice(rng_vals)
                # For every P1 sample...
                for p1_idx in range(len(pop_dicc[args.p1_population]['IND'])):
                    # For every P2 sample...
                    for p2_idx in range(len(pop_dicc[args.p2_population]['IND'])):
                        # For every P3 sample...
                        for p3_idx in range(len(pop_dicc[args.p3_population]['IND'])):
                            # For every P4 sample...
                            for p4_idx in range(len(pop_dicc[args.p4_population]['IND'])):
                                # Extract the samples.
                                p1_ind = pop_dicc[args.p1_population]['IND'][p1_idx]
                                p2_ind = pop_dicc[args.p2_population]['IND'][p2_idx]
                                p3_ind = pop_dicc[args.p3_population]['IND'][p3_idx]
                                p4_ind = pop_dicc[args.p4_population]['IND'][p4_idx]
                                # Construct the quartet.
                                quartet = '{0}-{1}-{2}-{3}'.format(p1_ind, p2_ind, p3_ind, p4_ind)
                                # Randomly sample an allele for each focal sample.
                                p1 = spline[pop_dicc[args.p1_population]['IDX'][p1_idx]][random_chr]
                                p2 = spline[pop_dicc[args.p2_population]['IDX'][p2_idx]][random_chr]
                                p3 = spline[pop_dicc[args.p3_population]['IDX'][p3_idx]][random_chr]
                                p4 = spline[pop_dicc[args.p4_population]['IDX'][p4_idx]][random_chr]
                                # Determine the site pattern.
                                if ((p1 == '0') and (p2 == '1') and (p3 == '1') and (p4 == '0')):
                                    site_patterns[pos][quartet]['ABBA'] += 1
                                    site_patterns[pos][quartet]['ABBA_HOM'] += 1
                                elif ((p1 == '1') and (p2 == '0') and (p3 == '0') and (p4 == '1')):
                                    site_patterns[pos][quartet]['ABBA'] += 1
                                    site_patterns[pos][quartet]['ABBA_HOM'] += 1
                                elif ((p1 == '0') and (p2 == '0') and (p3 == '1') and (p4 == '0')):
                                    site_patterns[pos][quartet]['ABBA_HOM'] += 1
                                elif ((p1 == '1') and (p2 == '1') and (p3 == '0') and (p4 == '1')):
                                    site_patterns[pos][quartet]['ABBA_HOM'] += 1
                                elif ((p1 == '1') and (p2 == '0') and (p3 == '1') and (p4 == '0')):
                                    site_patterns[pos][quartet]['BABA'] += 1
                                elif ((p1 == '0') and (p2 == '1') and (p3 == '0') and (p4 == '1')):
                                    site_patterns[pos][quartet]['BABA'] += 1
                                elif ((p1 == '1') and (p2 == '0') and (p3 == '0') and (p4 == '0')):
                                    site_patterns[pos][quartet]['BAAA'] += 1
                                    site_patterns[pos][quartet]['BAAA_HOM'] += 1
                                elif ((p1 == '0') and (p2 == '1') and (p3 == '1') and (p4 == '1')):
                                    site_patterns[pos][quartet]['BAAA'] += 1
                                    site_patterns[pos][quartet]['BAAA_HOM'] += 1
                                elif ((p1 == '1') and (p2 == '1') and (p3 == '0') and (p4 == '0')):
                                    site_patterns[pos][quartet]['BAAA_HOM'] += 1
                                elif ((p1 == '0') and (p2 == '0') and (p3 == '1') and (p4 == '1')):
                                    site_patterns[pos][quartet]['BAAA_HOM'] += 1
                                elif ((p1 == '0') and (p2 == '1') and (p3 == '0') and (p4 == '0')):
                                    site_patterns[pos][quartet]['ABAA'] += 1
                                elif ((p1 == '1') and (p2 == '0') and (p3 == '1') and (p4 == '1')):
                                    site_patterns[pos][quartet]['ABAA'] += 1
                                else:
                                    continue


# [3] Perform bootstrapping.

# Intialize a header list.
header_list = [
    'P1', 'P2', 'P3', 'P4',
    'ABBA', 'BABA', 'BAAA', 'ABAA',
    'ABBA_HOM', 'BABA_HOM', 'BAAA_HOM', 'ABAA_HOM',
    'D', 'Danc', 'D+', 'fhom', 'fanc', 'f+', 'BS_REP',
]
# Write the header list to the results file.
out_file.write('\t'.join(header_list)+'\n')
# Determine the number of blocks needed to build a bootstrapped contig.
blocks = args.contig_length // args.block_size
# For every bootstrap replicate...
for rep in range(args.replicates):
    # Intialize a list to store start positions.
    starts = []
    # For every block.
    for _ in range(blocks):
        # Generate a starting position.
        starts.append(random.randint(1, args.contig_length-args.block_size))
    # Sort the starting positions.
    starts = sorted(starts)
    # Intialzie a dictioanry to store position weights.
    pos_weight_dicc = {}
    # For every start position...
    for start in starts:
        # For every position covered by the start and end position...
        for pos in range(start, start+args.block_size):
            # If the current position is already known...
            if pos in pos_weight_dicc:
                # Add an additional weight.
                pos_weight_dicc[pos] += 1
            # Else...
            else:
                # Append the position to the dictionary.
                pos_weight_dicc[pos] = 1
    # Determine which bootstrapped positions overlap with the observed positions.
    overlap_list = list(set(list(site_patterns.keys())) & set(list(pos_weight_dicc.keys())))
    # If site patterns are to be calculated from derived allele frequencies...
    if freq_flag:
        # Intialize a dictionary to store bootstrapped site patterns.
        bootstrap_rep = {
            'ABBA': 0, 'ABBA_HOM': 0,
            'BABA': 0, 'BABA_HOM': 0,
            'BAAA': 0, 'BAAA_HOM': 0,
            'ABAA': 0, 'ABAA_HOM': 0,
        }
        # If there are overlapping positions...
        if len(overlap_list) > 0:
            # For every overlapping position...
            for pos in overlap_list:
                # For every site pattern...
                for key in bootstrap_rep.keys():
                    # Update the bootstrapped dictionary.
                    bootstrap_rep[key] += (site_patterns[pos][key] * pos_weight_dicc[pos])
        # Calculate numerators and denonimators for detection metrics.
        d_num = (bootstrap_rep['ABBA'] - bootstrap_rep['BABA'])
        d_den = (bootstrap_rep['ABBA'] + bootstrap_rep['BABA'])
        danc_num = (bootstrap_rep['BAAA'] - bootstrap_rep['ABAA'])
        danc_den = (bootstrap_rep['BAAA'] + bootstrap_rep['ABAA'])
        dplus_num = d_num + danc_num
        dplus_den = d_den + danc_den
        # Calculate numerators and denonimators for quantification metrics.
        fhom_num = d_num
        fhom_den = (bootstrap_rep['ABBA_HOM'] - bootstrap_rep['BABA_HOM'])
        fanc_num = danc_num
        fanc_den = (bootstrap_rep['BAAA_HOM'] - bootstrap_rep['ABAA_HOM'])
        fplus_num = dplus_num
        fplus_den = fhom_den + fanc_den
        # If D is undefined...
        if d_den == 0:
            # Set the value to nan.
            bootstrap_rep['D'] = 'nan'
        # Else...
        else:
            # Calculate D.
            bootstrap_rep['D'] = d_num / float(d_den)
        # If Danc is undefined...
        if danc_den == 0:
            # Set the value to nan.
            bootstrap_rep['Danc'] = 'nan'
        # Else...
        else:
            # Calculate Danc.
            bootstrap_rep['Danc'] = danc_num / float(danc_den)
        # If D+ is undefined...
        if dplus_den == 0:
            # Set the value to nan.
            bootstrap_rep['D+'] = 'nan'
        # Else...
        else:
            # Calculate D+.
            bootstrap_rep['D+'] = dplus_num / float(dplus_den)
        # If fhom is undefined...
        if fhom_den == 0:
            # Set the value to nan.
            bootstrap_rep['fhom'] = 'nan'
        # Else...
        else:
            # Calculate fhom.
            bootstrap_rep['fhom'] = fhom_num / float(fhom_den)
        # If fanc is undefined...
        if fanc_den == 0:
            # Set the value to nan.
            bootstrap_rep['fanc'] = 'nan'
        # Else...
        else:
            # Calculate fanc.
            bootstrap_rep['fanc'] = fanc_num / float(fanc_den)
        # If f+ is undefined...
        if fplus_den == 0:
            # Set the value to nan.
            bootstrap_rep['f+'] = 'nan'
        # Else...
        else:
            # Calculate f+.
            bootstrap_rep['f+'] = fplus_num / float(fplus_den)
        # Intialize the results list.
        results_list = [
            args.p1_population, args.p2_population, args.p3_population, args.p4_population,
            str(bootstrap_rep['ABBA']), str(bootstrap_rep['BABA']),
            str(bootstrap_rep['BAAA']), str(bootstrap_rep['ABAA']),
            str(bootstrap_rep['ABBA_HOM']), str(bootstrap_rep['BABA_HOM']),
            str(bootstrap_rep['BAAA_HOM']), str(bootstrap_rep['ABAA_HOM']),
            str(bootstrap_rep['D']), str(bootstrap_rep['Danc']), str(bootstrap_rep['D+']),
            str(bootstrap_rep['fhom']), str(bootstrap_rep['fanc']), str(bootstrap_rep['f+']), str(rep),
        ]
        # Write the results list to the results file.
        out_file.write('\t'.join(results_list)+'\n')
    # Else...
    else:
        # Intialize a dictionary to store bootstrapped site patterns.
        bootstrap_rep = {}
        # For every P1 sample...
        for p1 in pop_dicc[args.p1_population]['IND']:
            # For every P2 sample...
            for p2 in pop_dicc[args.p2_population]['IND']:
                # For every P3 sample...
                for p3 in pop_dicc[args.p3_population]['IND']:
                    # For every P4 sample...
                    for p4 in pop_dicc[args.p4_population]['IND']:
                        # Intialize the quartet.
                        quartet = '{0}-{1}-{2}-{3}'.format(p1, p2, p3, p4)
                        # Fill the dictionary to store site patterns.
                        bootstrap_rep[quartet] = {
                            'ABBA': 0, 'ABBA_HOM': 0,
                            'BABA': 0, 'BABA_HOM': 0,
                            'BAAA': 0, 'BAAA_HOM': 0,
                            'ABAA': 0, 'ABAA_HOM': 0,
                        }
        # If there are overlapping positions...
        if len(overlap_list) > 0:
            # For every overlapping position...
            for pos in overlap_list:
                # For every quartet...
                for quartet in bootstrap_rep.keys():
                    # For every site pattern...
                    for key in bootstrap_rep[quartet].keys():
                        bootstrap_rep[quartet][key] += (site_patterns[pos][quartet][key] * pos_weight_dicc[pos])
        # For every quartet...
        for key in bootstrap_rep.keys():
            # Calculate numerators and denonimators for detection metrics.
            d_num = (bootstrap_rep[key]['ABBA'] - bootstrap_rep[key]['BABA'])
            d_den = (bootstrap_rep[key]['ABBA'] + bootstrap_rep[key]['BABA'])
            danc_num = (bootstrap_rep[key]['BAAA'] - bootstrap_rep[key]['ABAA'])
            danc_den = (bootstrap_rep[key]['BAAA'] + bootstrap_rep[key]['ABAA'])
            dplus_num = d_num + danc_num
            dplus_den = d_den + danc_den
            # Calculate numerators and denonimators for quantification metrics.
            fhom_num = d_num
            fhom_den = (bootstrap_rep[key]['ABBA_HOM'] - bootstrap_rep[key]['BABA_HOM'])
            fanc_num = danc_num
            fanc_den = (bootstrap_rep[key]['BAAA_HOM'] - bootstrap_rep[key]['ABAA_HOM'])
            fplus_num = dplus_num
            fplus_den = fhom_den + fanc_den
            # If D is undefined...
            if d_den == 0:
                # Set the value to nan.
                bootstrap_rep[key]['D'] = 'nan'
            # Else...
            else:
                # Calculate D.
                bootstrap_rep[key]['D'] = d_num / float(d_den)
            # If Danc is undefined...
            if danc_den == 0:
                # Set the value to nan.
                bootstrap_rep[key]['Danc'] = 'nan'
            # Else...
            else:
                # Calculate Danc.
                bootstrap_rep[key]['Danc'] = danc_num / float(danc_den)
            # If D+ is undefined...
            if dplus_den == 0:
                # Set the value to nan.
                bootstrap_rep[key]['D+'] = 'nan'
            # Else...
            else:
                # Calculate D+.
                bootstrap_rep[key]['D+'] = dplus_num / float(dplus_den)
            # If fhom is undefined...
            if fhom_den == 0:
                # Set the value to nan.
                bootstrap_rep[key]['fhom'] = 'nan'
            # Else...
            else:
                # Calculate fhom.
                bootstrap_rep[key]['fhom'] = fhom_num / float(fhom_den)
            # If fanc is undefined...
            if fanc_den == 0:
                # Set the value to nan.
                bootstrap_rep[key]['fanc'] = 'nan'
            # Else...
            else:
                # Calculate fanc.
                bootstrap_rep[key]['fanc'] = fanc_num / float(fanc_den)
            # If f+ is undefined...
            if fplus_den == 0:
                # Set the value to nan.
                bootstrap_rep[key]['f+'] = 'nan'
            # Else...
            else:
                # Calculate f+.
                bootstrap_rep[key]['f+'] = fplus_num / float(fplus_den)
            # Unpack the samples in the quartet.
            p1, p2, p3, p4 = key.split('-')
            # Intialize the results list.
            results_list = [
                p1, p2, p3, p4,
                str(bootstrap_rep[key]['ABBA']), str(bootstrap_rep[key]['BABA']),
                str(bootstrap_rep[key]['BAAA']), str(bootstrap_rep[key]['ABAA']),
                str(bootstrap_rep[key]['ABBA_HOM']), str(bootstrap_rep[key]['BABA_HOM']),
                str(bootstrap_rep[key]['BAAA_HOM']), str(bootstrap_rep[key]['ABAA_HOM']),
                str(bootstrap_rep[key]['D']), str(bootstrap_rep[key]['Danc']), str(bootstrap_rep[key]['D+']),
                str(bootstrap_rep[key]['fhom']), str(bootstrap_rep[key]['fanc']), str(bootstrap_rep[key]['f+']), str(rep),
            ]
            # Write the results list to the results file.
            out_file.write('\t'.join(results_list)+'\n')

        
# [4] Close output files.

# Close the results and log files.
out_file.close()
log_file.close()