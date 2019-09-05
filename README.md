# Pactyper
Snakemake pipeline for continious clone type prediction for WGS sequenced bacterial isolates based on their core genome.

## General information

All the code is in the Snakefile and is written in snakemake. Snakefile_non_computerome should be used by users which are not using [Computerome](https://www.computerome.dk/) to run the pipeline.

The pipeline takes (1) a fastq file of genome of interest and (2) a core genome fasta file as input and outputs:

Output file description | Output file location
------------ | ------------
Alignment quality statistics for the input sample | `[output_dir]/sample_alignments/alignment_statistics.txt`
Aligned fasta file with SNPs applied to the core genome sequence | `[output_dir]/sample_alignments/[sample_name].aligned.fa`
VCF file will all the called variants in the core genome | `[output_dir]/sample_alignments/[sample_name].vcf.gz`
Clone type preditcion for the input sample | `[output_dir]/sample_alignments/[sample_name]_clonetype_summary.txt`
SNP distance matrix which includes all the samples already present in the matrix | `[output_dir]/sample_alignments/[sample_name]_clonetype_snp_distances.txt`
SNP distance matrix visualization | IN PROGRESS, not implemented yet


## Required software

Before running the pipeline, make sure that the following programs are installed and added to the path:

[GNU parallel >=2013xxxx](https://www.gnu.org/software/parallel/) <br/>
[Perl>=5.12](https://www.perl.org/) <br/>
Perl Modules: Time::Piece (core with modern Perl) <br/>
[Bioperl >= 1.6](https://bioperl.org/) <br/>
[bwa mem>=0.7.12](http://bio-bwa.sourceforge.net/)<br/>
[readseq>=2.0](http://iubio.bio.indiana.edu/soft/molbio/readseq/java/)<br/>
[samclip>=0.2](https://github.com/tseemann/samclip)<br/>
[bedtools>2.0](https://bedtools.readthedocs.io/en/latest/)<br/>
[freebayes>=1.1](https://github.com/ekg/freebayes)<br/>
[vcflib>=1.0](https://github.com/vcflib/vcflib)<br/>
[vcftools>=0.1.16](http://vcftools.sourceforge.net/)<br/>
[snpeff>=4.3](http://snpeff.sourceforge.net/)<br/>
[minimap2>=2.6](https://github.com/lh3/minimap2)<br/>
[seqtk>=1.2](https://github.com/lh3/seqtk)<br/>
[snp-sites>=2.0](https://github.com/sanger-pathogens/snp-sites)<br/>
[snippy>=4.1.0](https://github.com/tseemann/snippy)<br/>
[vt>=0.5](https://genome.sph.umich.edu/wiki/Vt)<br/>
[samtools>=1.9](http://www.htslib.org/doc/samtools.html)<br/>
[seqkit>=0.7](https://bioinf.shenwei.me/seqkit/)<br/>
[snp-dists>=0.6.3](https://github.com/tseemann/snp-dists)
[datamash>=1.4](https://www.gnu.org/software/datamash/)

## Setting up the config.yaml file

In order for the pipeline to run, a configuration file is needed. A configuration file requires 6 fields to be present:

Field name | Description
------------ | ------------
input_dir | Directory in which all fastq files which will be analyzed are present
input_sample | The unique prefix of the sample for which the clone type has to be predicted
core_genome | Full path to the FASTA file containing all the core genome genes
output_dir | Directory where all the output files will be stored
include | If the input sample should be included to the final matrix with the predicted clone type and used in the future iterations
SNP_distance | Number of SNP difference in the core genome for isolates to be defined as of different clone type


Here is the example config.yaml file:

```
input_dir: "/home/project_name/fastqs/"
input_sample: "551_12062011-DK10-0"
core_genome: "/home/project_name/Pseudomonas_aeruginosa/Pseudomonas_aeruginosa_core_genome.fasta"
output_dir: "/home/project_name/output_files"
include: True
SNP_distance: 5000

```
### Overwritting configuration file in the command line


## Input files and input requirements

The code is written to (1) build and (2) apply and extend the clone type matrix over time.

1. The first input sample will be assigned the `[prefix]001` clone type. Include is automatically defined as "True" for the first input sample.
2. The second all other input samples will be compared to the existing samples and if the SNP distance is lower than defined in the configuration file, the existing clone type which passes the criteria is assigned. If none of the clone types are close enough to the input file, the new clone type is assigned to the input sample.
3. If it is stated in the configuration file that the input sample should be included in the final martix, the predicted clone type will be assigned to the input sample and it will be added to the final distance matrix.

### Disclaimer

No files should be deleted from the output_files directory or the code will fail during the next run. 

[Snippy4](https://github.com/tseemann/snippy) doesn't work with python3. Python3 should be disabled at that step and Python2 should be available. 

## Running the pipeline

In order to run the pipeline anaconda3 (version 4.0.0) has to be available. Snakemake is started from its directory:

```
-j option allows to choose the number of threads (1-28) used for the analysis (default:1)
--configfile option allows to chose the configuration file for the analysis
--config option allows to overwrite the config file
```

Here is the example code for running snakemake:

```
snakemake -j 10 --configfile config.yaml --config input_sample="test_sample"

```


## Rerunning the pipeline

Pipeline can be rerun when new samples are added by changing the sample name in the config.yaml file or by overwritting the config file in the command line.

## Author

Migle Gabrielaite | migle.gabrielaite@regionh.dk
