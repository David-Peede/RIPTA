![ripta_logo](./ripta_logo.png)

# Re-evaluating Introgression site Patterns Through Ancestral alleles

If you have any questions or feature requests please leave a detailed issue. If you use `RIPTA` in your work please cite doi: X.



## `observed_values.py`

```bash
usage: observed_values.py [-h] -v VCF_FILE -m META_DATA -f {True,False} -p1
                        P1_POPULATION -p2 P2_POPULATION -p3 P3_POPULATION -p4
                        P4_POPULATION [-p PATH]

optional arguments:
  -h, --help            show this help message and exit
  -v VCF_FILE, --vcf_file VCF_FILE
                        Path to the input vcf file.
  -m META_DATA, --meta_data META_DATA
                        Path to the tab delimited text file where the first
                        column contains the sample name as it appears in the
                        header of the input vcf file and the second column
                        contains the population the sample belongs to.
  -f {True,False}, --frequency {True,False}
                        If true calculate site patterns from derived allele
                        frequencies, if false calculate site patterns by
                        randomly sampling one chromosome.
  -p1 P1_POPULATION, --p1_population P1_POPULATION
                        Population in the meta data file to use as P1 for
                        comparisons (ie a potential recipient population).
  -p2 P2_POPULATION, --p2_population P2_POPULATION
                        Population in the meta data file to use as P2 for
                        comparisons (ie a potential recipient population).
  -p3 P3_POPULATION, --p3_population P3_POPULATION
                        Population in the meta data file to use as P3 for
                        comparisons (ie the source population).
  -p4 P4_POPULATION, --p4_population P4_POPULATION
                        Population in the meta data file to use as P3 for
                        comparisons (ie the outrgroup population used for
                        polarization).
  -p PATH, --path PATH  Path for results and qc file (default = current
                        working directory).
```



## `bootstrapping.py`

```bash
usage: bootstrapping.py [-h] -v VCF_FILE -m META_DATA -f {True,False} -p1
                        P1_POPULATION -p2 P2_POPULATION -p3 P3_POPULATION -p4
                        P4_POPULATION -cl CONTIG_LENGTH -bs BLOCK_SIZE -r
                        REPLICATES [-p PATH]

optional arguments:
  -h, --help            show this help message and exit
  -v VCF_FILE, --vcf_file VCF_FILE
                        Path to the input vcf file.
  -m META_DATA, --meta_data META_DATA
                        Path to the tab delimited text file where the first
                        column contains the sample name as it appears in the
                        header of the input vcf file and the second column
                        contains the population the sample belongs to.
  -f {True,False}, --frequency {True,False}
                        If true calculate site patterns from derived allele
                        frequencies, if false calculate site patterns by
                        randomly sampling one chromosome.
  -p1 P1_POPULATION, --p1_population P1_POPULATION
                        Population in the meta data file to use as P1 for
                        comparisons (ie a potential recipient population).
  -p2 P2_POPULATION, --p2_population P2_POPULATION
                        Population in the meta data file to use as P2 for
                        comparisons (ie a potential recipient population).
  -p3 P3_POPULATION, --p3_population P3_POPULATION
                        Population in the meta data file to use as P3 for
                        comparisons (ie the source population).
  -p4 P4_POPULATION, --p4_population P4_POPULATION
                        Population in the meta data file to use as P3 for
                        comparisons (ie the outrgroup population used for
                        polarization).
  -cl CONTIG_LENGTH, --contig_length CONTIG_LENGTH
                        Total length in base pairs of conting.
  -bs BLOCK_SIZE, --block_size BLOCK_SIZE
                        Block size to be sampled with replacemnt for builidng
                        bootstrapped contigs.
  -r REPLICATES, --replicates REPLICATES
                        Number of bootstrapped replicates to perform.
  -p PATH, --path PATH  Path for results and qc file (default = current
                        working directory).
```



## Example Usage

In the `tutorial` directory I have supplied a pre-filtered VCF file of chromosome 22 for the YRI, CEU, and Altai Neanderthal individuals along with the ancestral allele calls from the EPO pipeline encoded as an Ancestor individual. To calculate site patterns and all introgression metrics based on derived allele frequencies you can run:

```bash
python observed_values.py -v ./tutorial/tutorial_data.vcf.gz -m ./tutorial/tutorial_meta_data_freqs.txt -f True -p1 YRI -p2 CEU -p3 NEA -p4 ANC -p ./tutorial/freq_results/
```

and to calculate site patterns and all introgression metrics based on individual trios you can run:

```bash
python observed_values.py -v ./tutorial.tutorial_data.vcf.gz -m ./tutorial/tutorial_meta_data_trios.txt -f False -p1 YRI -p2 CEU -p3 NEA -p4 ANC -p ./tutorial/trio_results/
```

To assess significance by bootstrapping based on estimates from derived allele frequencies you can run:

```bash
python bootstrapping.py -v ./tutorial/tutorial_data.vcf.gz -m ./tutorial/tutorial_meta_data_freqs.txt --frequency -p1 YRI -p2 CEU -p3 NEA -p4 ANC -cl 51304566 -bs 10_000_000 -r 100 -p ./tutorial/freq_results/
```

and to perform bootstrapping on individual trios you can run:

```bash
python bootstrapping.py -v ./tutorial.tutorial_data.vcf.gz -m ./tutorial/tutorial_meta_data_trios.txt -p1 YRI -p2 CEU -p3 NEA -p4 ANC -cl 51304566 -bs 10_000_000 -r 100 -p ./tutorial/trio_results/
```

Both `observed_values.py` and `bootstrapping.py` expect a filtered VCF file and will output a results and a log file in the directory you specify with the `-p` option.
