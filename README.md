# Pactyper
Snakemake pipeline for continiuos clone type prediction for WGS sequenced bacterial isolates based on their core genome.

## General information

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

## Setting up the config.yaml file

In order for the pipeline to run, a configuration file is needed. A configuration file requires 6 fields to be present:

Field name | Description
------------ | ------------
input_dir | Directory in which all fastq files which will be analyzed are present
input_sample | The name of the sample for which the clone type has to be predicted
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

## Input files and input requirements

No files should be deleted from the output directory or the code will fail during the next run. 

## Running the pipeline

## Rerunning the pipeline

Pipeline can be rerun when new samples are added by changing the sample name in the config.yaml file. 

## Author

Migle Gabrielaite | migle.gabrielaite@regionh.dk
