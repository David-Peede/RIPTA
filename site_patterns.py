# Import packages.
import argparse
import gzip
import random
import sys
import itertools

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
    '-p', '--path', required=False,
    action='store', default='./',
    help='Path for results and qc file (default = current working directory).',
)

parser.add_argument("-GP",help="Uses GP when available", action='store_true', default = False)

args = parser.parse_args()

# [0] Intialize output files.

# Intailzie the output file prefix.
prefix = args.path+'{0}_{1}_{2}_{3}_'.format(args.p1_population, args.p2_population, args.p3_population, args.p4_population)
# Intialize the output files.
out_file = open(prefix+'observed_values.txt', 'w')
log_file = open(prefix+'observed_values_log.txt', 'w')


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
freq_flag = args.frequency == 'True' or args.GP
# If the file is gzipped...
if args.vcf_file.endswith('.gz'):
    # Intialize the gzipped vcf file.
    data = gzip.open(args.vcf_file, 'rt')
# Else...
else:
    # Intialize the vcf file.
    data = open(args.vcf_file, 'rt')

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
        # If site patterns are to be calculated from derived allele frequencies...
        if freq_flag:
            # Intialize a dictionary to store site patterns.
            site_patterns = {
                'ABBA': 0, 'ABBA_HOM': 0,
                'BABA': 0, 'BABA_HOM': 0,
                'BAAA': 0, 'BAAA_HOM': 0,
                'ABAA': 0, 'ABAA_HOM': 0,
            }
        # Else...
        else:
            # Intialize a dictionary to store site patterns.
            site_patterns = {}
            # For every quartet...
            for p1, p2, p3, p4 in itertools.product(pop_dicc[args.p1_population]['IND'],
                pop_dicc[args.p2_population]['IND'], pop_dicc[args.p3_population]['IND'],
                pop_dicc[args.p4_population]['IND']):
                # Intialize the quartet.
                quartet = '{0}-{1}-{2}-{3}'.format(p1, p2, p3, p4)
                # Fill the dictionary to store site patterns.
                site_patterns[quartet] = {
                    'ABBA': 0, 'ABBA_HOM': 0,
                    'BABA': 0, 'BABA_HOM': 0,
                    'BAAA': 0, 'BAAA_HOM': 0,
                    'ABAA': 0, 'ABAA_HOM': 0,
                }
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
            # Intialize a alternative allele frequency dictionary.
            freq_dicc = {}
            # For every focal population...
            for key in pop_dicc.keys():
                # Intialize an alternative allele counter.
                alt_allele_counter = 0
                # For every individual in the population...
                if not args.GP: #GP
                    for idx in pop_dicc[key]['IDX']:
		        # Count the number of alternative alleles.
                        alt_allele_counter += spline[idx][0:3].count('1')
                else:
                    posGP=-1
                    if 'GP' in spline[8].split(':'):
                    	posGP=spline[8].split(':').index('GP')
                    for idx in pop_dicc[key]['IDX']:
                        sspline = spline[idx].split(':')
                        if posGP!=-1 and sspline[posGP]!='.':
                            alt_allele_counter += float(sspline[posGP].split(',')[2])*2+float(sspline[posGP].split(',')[1])
                        else:
                            alt_allele_counter += float(spline[idx][0:3].count('1'))
		       	   
		    
			
		# Determine the alternative allele frequency.
                freq_dicc[key] = float(alt_allele_counter) / (len(pop_dicc[key]['IDX']) * 2)
            # If the ancestral allele is the refernce allele...
            if freq_dicc[args.p4_population] == 0.0:
                # Calculate site patterns.
                site_patterns['ABBA'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p2_population] * freq_dicc[args.p3_population]
                site_patterns['ABBA_HOM'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p3_population] * freq_dicc[args.p3_population]
                site_patterns['BABA'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p2_population]) * freq_dicc[args.p3_population]
                site_patterns['BABA_HOM'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p3_population]) * freq_dicc[args.p3_population]
                site_patterns['BAAA'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p2_population]) * (1 - freq_dicc[args.p3_population])
                site_patterns['BAAA_HOM'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p3_population]) * (1 - freq_dicc[args.p3_population])
                site_patterns['ABAA'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p2_population] * (1 - freq_dicc[args.p3_population])
                site_patterns['ABAA_HOM'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p3_population] * (1 - freq_dicc[args.p3_population])
            # Else-if the ancestral allele is the alternative allele...
            elif freq_dicc[args.p4_population] == 1.0:
                # Calculate site patterns.
                site_patterns['ABBA'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p2_population]) * (1 - freq_dicc[args.p3_population])
                site_patterns['ABBA_HOM'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p3_population]) * (1 - freq_dicc[args.p3_population])
                site_patterns['BABA'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p2_population] * (1 - freq_dicc[args.p3_population])
                site_patterns['BABA_HOM'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p3_population] * (1 - freq_dicc[args.p3_population])
                site_patterns['BAAA'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p2_population] * freq_dicc[args.p3_population]
                site_patterns['BAAA_HOM'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p3_population] * freq_dicc[args.p3_population]
                site_patterns['ABAA'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p2_population]) * freq_dicc[args.p3_population]
                site_patterns['ABAA_HOM'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p3_population]) * freq_dicc[args.p3_population]
            # Else...
            else:
                # Calculate site patterns.
                site_patterns['ABBA'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p2_population] * freq_dicc[args.p3_population] * (1 - freq_dicc[args.p4_population])
                site_patterns['ABBA_HOM'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p3_population] * freq_dicc[args.p3_population] * (1 - freq_dicc[args.p4_population])
                site_patterns['BABA'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p2_population]) * freq_dicc[args.p3_population] * (1 - freq_dicc[args.p4_population])
                site_patterns['BABA_HOM'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p3_population]) * freq_dicc[args.p3_population] * (1 - freq_dicc[args.p4_population])
                site_patterns['BAAA'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p2_population]) * (1 - freq_dicc[args.p3_population]) * (1 - freq_dicc[args.p4_population])
                site_patterns['BAAA_HOM'] += freq_dicc[args.p1_population] * (1 - freq_dicc[args.p3_population]) * (1 - freq_dicc[args.p3_population]) * (1 - freq_dicc[args.p4_population])
                site_patterns['ABAA'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p2_population] * (1 - freq_dicc[args.p3_population]) * (1 - freq_dicc[args.p4_population])
                site_patterns['ABAA_HOM'] += (1 - freq_dicc[args.p1_population]) * freq_dicc[args.p3_population] * (1 - freq_dicc[args.p3_population]) * (1 - freq_dicc[args.p4_population])
        # Else...
        else:
            # Intialize RNG values.
            rng_vals = [0, 2]
            # Randomly select a chromosome to sample.
            random_chr = random.choice(rng_vals)
            # For every index quartet...
            for p1_idx, p2_idx, p3_idx, p4_idx in itertools.product(range(len(pop_dicc[args.p1_population]['IND'])),
                range(len(pop_dicc[args.p2_population]['IND'])), range(len(pop_dicc[args.p3_population]['IND'])),
                range(len(pop_dicc[args.p4_population]['IND']))):
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
                    site_patterns[quartet]['ABBA'] += 1
                    site_patterns[quartet]['ABBA_HOM'] += 1
                elif ((p1 == '1') and (p2 == '0') and (p3 == '0') and (p4 == '1')):
                    site_patterns[quartet]['ABBA'] += 1
                    site_patterns[quartet]['ABBA_HOM'] += 1
                elif ((p1 == '0') and (p2 == '0') and (p3 == '1') and (p4 == '0')):
                    site_patterns[quartet]['ABBA_HOM'] += 1
                elif ((p1 == '1') and (p2 == '1') and (p3 == '0') and (p4 == '1')):
                    site_patterns[quartet]['ABBA_HOM'] += 1
                elif ((p1 == '1') and (p2 == '0') and (p3 == '1') and (p4 == '0')):
                    site_patterns[quartet]['BABA'] += 1
                elif ((p1 == '0') and (p2 == '1') and (p3 == '0') and (p4 == '1')):
                    site_patterns[quartet]['BABA'] += 1
                elif ((p1 == '1') and (p2 == '0') and (p3 == '0') and (p4 == '0')):
                    site_patterns[quartet]['BAAA'] += 1
                    site_patterns[quartet]['BAAA_HOM'] += 1
                elif ((p1 == '0') and (p2 == '1') and (p3 == '1') and (p4 == '1')):
                    site_patterns[quartet]['BAAA'] += 1
                    site_patterns[quartet]['BAAA_HOM'] += 1
                elif ((p1 == '1') and (p2 == '1') and (p3 == '0') and (p4 == '0')):
                    site_patterns[quartet]['BAAA_HOM'] += 1
                elif ((p1 == '0') and (p2 == '0') and (p3 == '1') and (p4 == '1')):
                    site_patterns[quartet]['BAAA_HOM'] += 1
                elif ((p1 == '0') and (p2 == '1') and (p3 == '0') and (p4 == '0')):
                    site_patterns[quartet]['ABAA'] += 1
                elif ((p1 == '1') and (p2 == '0') and (p3 == '1') and (p4 == '1')):
                    site_patterns[quartet]['ABAA'] += 1
                else:
                    continue

# Close VCF file
data.close()

# [3] Write the results file.

# Intialize a header list.
header_list = [
    'P1', 'P2', 'P3', 'P4',
    'ABBA', 'BABA', 'BAAA', 'ABAA',
    'ABBA_HOM', 'BABA_HOM', 'BAAA_HOM', 'ABAA_HOM',
    'D', 'Danc', 'D+', 'fhom', 'fanc', 'f+',
]
# Write the header list to the results file.
out_file.write('\t'.join(header_list)+'\n')
# If site patterns were calculated from derived allele frequencies...
if freq_flag:
    # Calculate numerators and denonimators for detection metrics.
    d_num = (site_patterns['ABBA'] - site_patterns['BABA'])
    d_den = (site_patterns['ABBA'] + site_patterns['BABA'])
    danc_num = (site_patterns['BAAA'] - site_patterns['ABAA'])
    danc_den = (site_patterns['BAAA'] + site_patterns['ABAA'])
    dplus_num = d_num + danc_num
    dplus_den = d_den + danc_den
    # Calculate numerators and denonimators for quantification metrics.
    fhom_num = d_num
    fhom_den = (site_patterns['ABBA_HOM'] - site_patterns['BABA_HOM'])
    fanc_num = danc_num
    fanc_den = (site_patterns['BAAA_HOM'] - site_patterns['ABAA_HOM'])
    fplus_num = dplus_num
    fplus_den = fhom_den + fanc_den
    # If D is undefined...
    if d_den == 0:
        # Set the value to nan.
        site_patterns['D'] = 'nan'
    # Else...
    else:
        # Calculate D.
        site_patterns['D'] = d_num / float(d_den)
    # If Danc is undefined...
    if danc_den == 0:
        # Set the value to nan.
        site_patterns['Danc'] = 'nan'
    # Else...
    else:
        # Calculate Danc.
        site_patterns['Danc'] = danc_num / float(danc_den)
    # If D+ is undefined...
    if dplus_den == 0:
        # Set the value to nan.
        site_patterns['D+'] = 'nan'
    # Else...
    else:
        # Calculate D+.
        site_patterns['D+'] = dplus_num / float(dplus_den)
    # If fhom is undefined...
    if fhom_den == 0:
        # Set the value to nan.
        site_patterns['fhom'] = 'nan'
    # Else...
    else:
        # Calculate fhom.
        site_patterns['fhom'] = fhom_num / float(fhom_den)
    # If fanc is undefined...
    if fanc_den == 0:
        # Set the value to nan.
        site_patterns['fanc'] = 'nan'
    # Else...
    else:
        # Calculate fanc.
        site_patterns['fanc'] = fanc_num / float(fanc_den)
    # If f+ is undefined...
    if fplus_den == 0:
        # Set the value to nan.
        site_patterns['f+'] = 'nan'
    # Else...
    else:
        # Calculate f+.
        site_patterns['f+'] = fplus_num / float(fplus_den)
    # Intialize the results list.
    results_list = [
        args.p1_population, args.p2_population, args.p3_population, args.p4_population,
        str(site_patterns['ABBA']), str(site_patterns['BABA']),
        str(site_patterns['BAAA']), str(site_patterns['ABAA']),
        str(site_patterns['ABBA_HOM']), str(site_patterns['BABA_HOM']),
        str(site_patterns['BAAA_HOM']), str(site_patterns['ABAA_HOM']),
        str(site_patterns['D']), str(site_patterns['Danc']), str(site_patterns['D+']),
        str(site_patterns['fhom']), str(site_patterns['fanc']), str(site_patterns['f+']),
    ]
    # Write the results list to the results file.
    out_file.write('\t'.join(results_list)+'\n')
# Else...
else:
    # For every quartet...
    for key in site_patterns.keys():
        # Calculate numerators and denonimators for detection metrics.
        d_num = (site_patterns[key]['ABBA'] - site_patterns[key]['BABA'])
        d_den = (site_patterns[key]['ABBA'] + site_patterns[key]['BABA'])
        danc_num = (site_patterns[key]['BAAA'] - site_patterns[key]['ABAA'])
        danc_den = (site_patterns[key]['BAAA'] + site_patterns[key]['ABAA'])
        dplus_num = d_num + danc_num
        dplus_den = d_den + danc_den
        # Calculate numerators and denonimators for quantification metrics.
        fhom_num = d_num
        fhom_den = (site_patterns[key]['ABBA_HOM'] - site_patterns[key]['BABA_HOM'])
        fanc_num = danc_num
        fanc_den = (site_patterns[key]['BAAA_HOM'] - site_patterns[key]['ABAA_HOM'])
        fplus_num = dplus_num
        fplus_den = fhom_den + fanc_den
        # If D is undefined...
        if d_den == 0:
            # Set the value to nan.
            site_patterns[key]['D'] = 'nan'
        # Else...
        else:
            # Calculate D.
            site_patterns[key]['D'] = d_num / float(d_den)
        # If Danc is undefined...
        if danc_den == 0:
            # Set the value to nan.
            site_patterns[key]['Danc'] = 'nan'
        # Else...
        else:
            # Calculate Danc.
            site_patterns[key]['Danc'] = danc_num / float(danc_den)
        # If D+ is undefined...
        if dplus_den == 0:
            # Set the value to nan.
            site_patterns[key]['D+'] = 'nan'
        # Else...
        else:
            # Calculate D+.
            site_patterns[key]['D+'] = dplus_num / float(dplus_den)
        # If fhom is undefined...
        if fhom_den == 0:
            # Set the value to nan.
            site_patterns[key]['fhom'] = 'nan'
        # Else...
        else:
            # Calculate fhom.
            site_patterns[key]['fhom'] = fhom_num / float(fhom_den)
        # If fanc is undefined...
        if fanc_den == 0:
            # Set the value to nan.
            site_patterns[key]['fanc'] = 'nan'
        # Else...
        else:
            # Calculate fanc.
            site_patterns[key]['fanc'] = fanc_num / float(fanc_den)
        # If f+ is undefined...
        if fplus_den == 0:
            # Set the value to nan.
            site_patterns[key]['f+'] = 'nan'
        # Else...
        else:
            # Calculate f+.
            site_patterns[key]['f+'] = fplus_num / float(fplus_den)
        # Unpack the samples in the quartet.
        p1, p2, p3, p4 = key.split('-')
        # Intialize the results list.
        results_list = [
            p1, p2, p3, p4,
            str(site_patterns[key]['ABBA']), str(site_patterns[key]['BABA']),
            str(site_patterns[key]['BAAA']), str(site_patterns[key]['ABAA']),
            str(site_patterns[key]['ABBA_HOM']), str(site_patterns[key]['BABA_HOM']),
            str(site_patterns[key]['BAAA_HOM']), str(site_patterns[key]['ABAA_HOM']),
            str(site_patterns[key]['D']), str(site_patterns[key]['Danc']), str(site_patterns[key]['D+']),
            str(site_patterns[key]['fhom']), str(site_patterns[key]['fanc']), str(site_patterns[key]['f+']),
        ]
        # Write the results list to the results file.
        out_file.write('\t'.join(results_list)+'\n')


# [4] Close output files.

# Close the results and log files.
out_file.close()
log_file.close()
