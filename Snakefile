import os
import sys
import argparse
import glob
from pathlib import Path
import re

if os.path.exists('config.yaml'):
	configfile: 'config.yaml'

# print function for debugging in snakemake
def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)


input_dir=config['input_dir']
input_sample=config['input_sample']
core_genome=config['core_genome']
output_dir=config['output_dir']
include_to_final=config['include']
SNP_distance=config['SNP_distance']
prefix=config['prefix']

# check if files are fastq.gz:
def check_fastq_end():
	fastqs = glob.glob(input_dir+'/'+input_sample+"*R1*"); fastqs += glob.glob(input_dir+'/'+input_sample+"*R2*")
	for fastq in fastqs:
		if not fastq.endswith('fastq.gz') and not fastq.endswith('fq.gz'):
			raise Exception("ERROR Input files must be fastq.gz or fq.gz !")
		return

check_fastq_end()

READS={}

READS['R1'] = glob.glob(input_dir+'/'+input_sample+"*R1*")[0]
READS['R2'] = glob.glob(input_dir+'/'+input_sample+"*R2*")[0]

# check if the file is already present in the matrix. Authomatically set include to False if the file is already in the matrix.
def check_if_exists_in_matrix():
	myfile_dist=Path(output_dir+'/output_files/complete_clonetype_distance_matrix.txt')
        if myfile_dist.is_file():
		with open(output_dir+'/output_files/complete_clonetype_distance_matrix.txt') as f:
			names = [line.split('\t')[0] for line in f][1:]
		trimmed_names = [name.split('::')[1] for name in names]
		if input_sample in trimmed_names:
			if include_to_final is True:
				config['include']=False
				eprint("The file already exists in the matrix. Parameter \"Include\" is set to False")

check_if_exists_in_matrix()

def number_of_existing_samples():
	myfile_aln=Path(output_dir+'/output_files/complete_clonetype_sorted_alignment.fa')
	if myfile_aln.is_file():
		sample_no=int(sum(1 for line in open(output_dir+'/output_files/complete_clonetype_sorted_alignment.fa') if line.rstrip())/2)
	else:
		sample_no=0
	return(sample_no)

sample_no=number_of_existing_samples()

# exclude all the names which are non-numeric
def follows_the_pattern(name):
	if(re.compile("[0-9][0-9][0-9]").match(name)):
		return True
	else:
		return False

def new_clonetype():
	myfile_dist=Path(output_dir+'/output_files/complete_clonetype_distance_matrix.txt')
	if myfile_dist.is_file():
		with open(output_dir+'/output_files/complete_clonetype_distance_matrix.txt') as f:
			names = [line.split('\t')[0] for line in f]
		num = [name.split('::')[0][-3:] for name in names]
		numbers=int(max(filter(follows_the_pattern, num)))+1
	elif sample_no==1:
		numbers=2
	else:
		numbers=1
	return("{:03d}".format(numbers))


new_ct=new_clonetype()
now_ct=float(new_ct)-1

eprint("Previously analyzed and added sample number:",sample_no)
eprint("Number of already existing clone types in the data set is:", int(now_ct))

rule all:
	input:
		vcf_gz=expand("{output_dir}/sample_alignments/{sample}/{sample}.vcf.gz", sample=config['input_sample'], output_dir=config['output_dir']),
		aligned=expand("{output_dir}/sample_alignments/{sample}/{sample}.aligned.fa", sample=config['input_sample'], output_dir=config['output_dir']),
		stats=expand("{output_dir}/sample_alignments/{sample}/alignment_statistics.txt", sample=config['input_sample'], output_dir=config['output_dir']),
		snp_dist=expand("{output_dir}/sample_alignments/{sample}/{sample}_clonetype_snp_distances.txt", output_dir=config['output_dir'], sample=config['input_sample']) if sample_no>=1 else [],
		initial_alignment=expand("{output_dir}/output_files/complete_clonetype_sorted_alignment.fa", output_dir=config['output_dir']) if sample_no==0 else [],
		clonetype_stats=expand("{output_dir}/sample_alignments/{sample}/{sample}_clonetype_summary.txt", output_dir=config['output_dir'], sample=config['input_sample']) if sample_no>=1 else []
	message:
		"Pipeline complete"

# Run snippy on the input sample using the core reference genome
rule run_snippy:
	input:
		R1=READS['R1'],
		R2=READS['R2']	
	output:
		vcf_gz=expand("{output_dir}/sample_alignments/{sample}/{sample}.vcf.gz", sample=config['input_sample'], output_dir=config['output_dir']),
		aligned=expand("{output_dir}/sample_alignments/{sample}/{sample}.aligned.fa", sample=config['input_sample'], output_dir=config['output_dir']),
		ref=temp(expand("{output_dir}/sample_alignments/{sample}/ref.fa", sample=config['input_sample'], output_dir=config['output_dir'])),
		ref_fai=temp(expand("{output_dir}/sample_alignments/{sample}/ref.fa.fai", sample=config['input_sample'], output_dir=config['output_dir'])),
		bam=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.bam", sample=config['input_sample'], output_dir=config['output_dir'])),
		bam_bai=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.bam.bai", sample=config['input_sample'], output_dir=config['output_dir'])),
		raw_vcf=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.raw.vcf", sample=config['input_sample'], output_dir=config['output_dir'])),
		filt_vcf=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.filt.vcf", sample=config['input_sample'], output_dir=config['output_dir'])),
		vcf=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.vcf", sample=config['input_sample'], output_dir=config['output_dir'])),
		tab=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.tab", sample=config['input_sample'], output_dir=config['output_dir'])),
		subs_vcf=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.subs.vcf", sample=config['input_sample'], output_dir=config['output_dir'])),
		consensus_fa=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.consensus.fa", sample=config['input_sample'], output_dir=config['output_dir'])),
		consensus_subs_fa=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.consensus.subs.fa", sample=config['input_sample'], output_dir=config['output_dir'])),
		log_file=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.log", sample=config['input_sample'], output_dir=config['output_dir'])),
		bed=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.bed", sample=config['input_sample'], output_dir=config['output_dir'])),
		gff=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.gff", sample=config['input_sample'], output_dir=config['output_dir'])),
		csv=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.csv", sample=config['input_sample'], output_dir=config['output_dir'])),
		html=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.html", sample=config['input_sample'], output_dir=config['output_dir'])),
		txt=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.txt", sample=config['input_sample'], output_dir=config['output_dir']))
		
	params:
		outdir=config['output_dir'],
		ref=config['core_genome'],
		prefix=config['input_sample'],
		cpus=10,
		mincov=10
	message:
		"Running snippy variant calling for {params.prefix}."
	log:
		expand("{output_dir}/logs/snippy_{sample}.log", output_dir=config["output_dir"], sample=config['input_sample'])
	shell:
		"""
		module unload anaconda3/4.0.0
		module load anaconda2/4.0.0
		module load tbl2asn/20190117
		module load parallel/20190122
		module load bwa/0.7.15
		module load samtools/1.9
		module load java/1.8.0
		module load jre/1.8.0
		module load perl/5.24.0
		module load freebayes/1.1.0-50-g61527c5
		module load vcflib/1.0.0-rc2
		module load vcftools/0.1.16
		module load snpeff/4.3r
		module load prokka/1.12
		module load minimap2/2.6
		module load seqtk/1.0-r82-dirty
		module load snp-sites/2.4.0
		module load emboss/6.6.0
		module load bcftools/1.9
		module load snippy/4.4.0
		module load vt/0.5772


		(snippy --outdir {params.outdir}/sample_alignments/{params.prefix} --ref {params.ref} --R1 {input.R1} --R2 {input.R2} --prefix {params.prefix} --cpus {params.cpus}  --force --mincov {params.mincov}) 2> {log}
		# how to add to a previous SNP distance? cp - if include is No, add if include is Yes?

		"""

def do_QC(fasta_file, output_file):
	with open(fasta_file, "r") as f:
		sequence = "".join(line.strip() for line in f if line[:1]!=">")
		seq_len = len(sequence)
		missing_len = sequence.count("N")+sequence.count("-")
		cons_positions = seq_len - missing_len
		missing_percent = (missing_len*100)/seq_len
		present_percent = 100 - missing_percent
	f.closed
	if missing_percent >5:
		raise Exception("ERROR: The sequence alignment is less than 95% of the reference sequence. Analysis can't be continued because of low coverage.")
	else:
		print("Sequence alignment covers 95% or more of the reference sequence")
		with open(output_file, "w") as outf:
			outf.write("Core genome consists of\t{total} positions.\n\
In total {considered} positions were covered by not less than 10 reads.\n\
{covered}% of postions were covered by not less than 10 reads.\n".format(total=seq_len, considered=cons_positions, covered=round(present_percent, 2)))
			
		
# Check if the alignment covers >=95% of the reference genome
rule do_QC:
	input:
		expand("{output_dir}/sample_alignments/{sample}/{sample}.aligned.fa", sample=config['input_sample'], output_dir=config['output_dir'])
	output:
		expand("{output_dir}/sample_alignments/{sample}/alignment_statistics.txt", sample=config['input_sample'], output_dir=config['output_dir'])
	message:
		"Performing quality check of the alignment coverage."
	run:
		do_QC(input[0], output[0])

if sample_no >= 2:
# SNP distance analysis on the provided samples
	rule snp_dists:
		input:
			aligned_file=expand("{output_dir}/sample_alignments/{sample}/{sample}.aligned.fa", sample=config['input_sample'], output_dir=config['output_dir']),
			initial_alignment=expand("{output_dir}/output_files/complete_clonetype_sorted_alignment.fa", output_dir=config['output_dir']),
			alignment_statistics=expand("{output_dir}/sample_alignments/{sample}/alignment_statistics.txt", sample=config['input_sample'], output_dir=config['output_dir'])
		output:
			consensus_sample=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.consensus.fa", sample=config['input_sample'], output_dir=config['output_dir'])),
			snp_dist=expand("{output_dir}/sample_alignments/{sample}/{sample}_clonetype_snp_distances.txt", output_dir=config['output_dir'], sample=config['input_sample']),
			temp_names=temp(expand("{output_dir}/output_files/temporary_names_snp_dist.txt", output_dir=config['output_dir'])),
			temp_dist=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}_temp_distances.txt", sample=config['input_sample'], output_dir=config['output_dir'])),
			temp_fasta=temp(expand("{output_dir}/sample_alignments/{sample}/temp_fasta.fa", output_dir=config['output_dir'], sample=config['input_sample'])),
		params:
			initial_distance=expand("{output_dir}/output_files/complete_clonetype_distance_matrix.txt", output_dir=config['output_dir']) if include_to_final == True else [],
			temp_consensus=temp(expand("{output_dir}/output_files/temp_complete_clonetype_sorted_alignment.fa", output_dir=config['output_dir'])),
			sample=config['input_sample'],
			incl=config['include']
		message:
			"Performing SNP distance analysis on the {params.sample}."
		shell:
			"""
			module load pigz/2.3.4
			module load seqkit/0.7.1	

			# creating a sorted fasta file with all the genes merged in one sequence
	
			seqkit sort --id-regexp "(.+)" {input.aligned_file} | grep -v "^>" | awk \'!/^>/ {{ printf "%s", $0; n = "\\n" }} /^>/ {{ print n $0; n = "" }} END {{ printf "%s", n }}\' |  sed "1s/^/>{params.sample}\\n/" > {output.consensus_sample}

			# performing snp-dist for each sample pair (new sample vs every older sample)
		
			cat {input.initial_alignment} | nl | grep ">" > {output.temp_names}	
			echo {params.sample} > {output.temp_dist}
		
			while read line name
			do
				# adding a 1 to the line number so both name and sequence would be included
				new_line=`echo "$line"+1 | bc`
				seq=`head -n $new_line {input.initial_alignment} | tail -2`
				new_seq=`cat {output.consensus_sample}`
				echo "$seq" > {output.temp_fasta}
				echo "$new_seq" >> {output.temp_fasta}
				dist=`snp-dists -q -b {output.temp_fasta} | tail -1  | cut -f 2`
				echo "$dist"
			done < {output.temp_names} >> {output.temp_dist}
	
			# merging the initial distance table with the new results
			paste {params.initial_distance} {output.temp_dist} > {output.snp_dist}	
		
			hline=`cat {output.temp_dist} | tr "\\n" "\\t" | sed "s/$/0\\n/"`
			echo "$hline" >> {output.snp_dist}	
			
		
			# if include is true output snps distance file and initial alignment file should be updated with the new file. Filename should be changed in the next rule.
			if [[ {params.incl} == "True" ]]
			then
				cat {input.initial_alignment} {output.consensus_sample} > {params.temp_consensus}
				mv {params.temp_consensus} {input.initial_alignment}
			
				cp {output.snp_dist} {params.initial_distance}
			fi
			"""
if sample_no == 1:
	rule snp_dists:
		input:
			aligned_file=expand("{output_dir}/sample_alignments/{sample}/{sample}.aligned.fa", sample=config['input_sample'], output_dir=config['output_dir']),
			initial_alignment=expand("{output_dir}/output_files/complete_clonetype_sorted_alignment.fa", output_dir=config['output_dir']),
			alignment_statistics=expand("{output_dir}/sample_alignments/{sample}/alignment_statistics.txt", sample=config['input_sample'], output_dir=config['output_dir'])
		output:
			consensus_sample=temp(expand("{output_dir}/sample_alignments/{sample}/{sample}.consensus.fa", sample=config['input_sample'], output_dir=config['output_dir'])),
			snp_dist=expand("{output_dir}/sample_alignments/{sample}/{sample}_clonetype_snp_distances.txt", output_dir=config['output_dir'], sample=config['input_sample']),
			temp_fasta=temp(expand("{output_dir}/sample_alignments/{sample}/temp_fasta.fa", output_dir=config['output_dir'], sample=config['input_sample'])),
			initial_distance=expand("{output_dir}/output_files/complete_clonetype_distance_matrix.txt", output_dir=config['output_dir'])
		params:
			temp_consensus=temp(expand("{output_dir}/output_files/temp_complete_clonetype_sorted_alignment.fa", output_dir=config['output_dir'])),
			sample=config['input_sample'],
			incl=config['include']
		message:
			"Performing SNP distance analysis on the {params.sample}."
		shell:
			"""
			module load pigz/2.3.4
			module load seqkit/0.7.1

			# creating a sorted fasta file with all the genes merged in one sequence

			seqkit sort --id-regexp "(.+)" {input.aligned_file} | grep -v "^>" | awk \'!/^>/ {{ printf "%s", $0; n = "\\n" }} /^>/ {{ print n $0; n = "" }} END {{ printf "%s", n }}\' |  sed "1s/^/>{params.sample}\\n/" > {output.consensus_sample}

			# performing snp-dists for the only two samples which are available

			cat {input.initial_alignment} {output.consensus_sample} > {output.temp_fasta}
			snp-dists -q -b {output.temp_fasta} > {output.snp_dist}
		

			# output snps distance file and initial alignment file should be updated with the new file. Filename should be changed in the next rule.
			
			cat {input.initial_alignment} {output.consensus_sample} > {params.temp_consensus}
			mv {params.temp_consensus} {input.initial_alignment}

			cp {output.snp_dist} {output.initial_distance}
			
			"""


# If there is only one alignment, the alignment should be added to the complete_alignment if include is YES file and no further action should be done.
rule create_complete_aligment:
	input:
		aligned_file=expand("{output_dir}/sample_alignments/{sample}/{sample}.aligned.fa", sample=config['input_sample'], output_dir=config['output_dir'])
	output:
		initial_alignment=expand("{output_dir}/output_files/complete_clonetype_sorted_alignment.fa", output_dir=config['output_dir']) if sample_no<1 else []
	params:
		prefix=config['prefix'],
		sample=config['input_sample']
	message:
		"Adding the sample to the multiple sequence FASTA file"
	shell:
		"""
		module load pigz/2.3.4
		module load seqkit/0.7.1

		# creating a sorted fasta file with all the genes merged in one sequence

		seqkit sort --id-regexp "(.+)" {input.aligned_file} |  grep -v "^>" | awk \'!/^>/ {{ printf "%s", $0; n = "\\n" }} /^>/ {{ print n $0; n = "" }} END {{ printf "%s", n }}\' |  sed "1s/^/>{params.prefix}001::{params.sample}\\n/" > {output.initial_alignment}
		"""


# Estimation of the closest clone types. If the SNP difference is <5000 SNPs, the isolate is given a clone type. Otherwise, it is considered a new clone type.
rule estimate_clonetype:
	input:
		snp_dist=expand("{output_dir}/sample_alignments/{sample}/{sample}_clonetype_snp_distances.txt", output_dir=config['output_dir'], sample=config['input_sample']),
		initial_alignment=expand("{output_dir}/output_files/complete_clonetype_sorted_alignment.fa", output_dir=config['output_dir']),
		initial_distance=expand("{output_dir}/output_files/complete_clonetype_distance_matrix.txt", output_dir=config['output_dir'])
	output:
		stats=expand("{output_dir}/sample_alignments/{sample}/{sample}_clonetype_summary.txt", output_dir=config['output_dir'], sample=config['input_sample']),
		temp_clonetype_list=temp(expand("{output_dir}/output_files/temp_{sample}_clonetype_list.txt", output_dir=config['output_dir'], sample=config['input_sample']))
	params:
		prefix=config['prefix'],
		distance=config['SNP_distance'],
		incl=config['include'],
		sample=config['input_sample'],
		new_clonetype=new_ct
	message:
		"Predicting the clone type for the input sample"
	shell:
		"""
		module load datamash/1.4

		close_clonetypes=`cat {input.snp_dist} | tail -1  | cut -f 2- | datamash transpose | head -n -1 | nl | awk '$2<{params.distance} {{print $1}}'`

		if [[ ! -z "$close_clonetypes" ]]
		then
			while read line
			do
				new_num=`echo "$line" +1 | bc`
				cat {input.initial_distance} | head -1 | cut -f "$new_num"| sed 's/::.*//'
			done <<< "$close_clonetypes" >> {output.temp_clonetype_list}

			clonetype_number=`cat {output.temp_clonetype_list} | sort | uniq -c | sort -k1,1n`
	
			line_number=`echo "$clonetype_number" | wc -l`
	
		else
			touch {output.temp_clonetype_list}
			line_number=0
		fi
	
		# if only one clone type is predicted:	
		if [[ "$line_number" == 1  ]]
		then
			number=`echo "$clonetype_number" | awk '{{print $1}}'`
			clonetype=`echo "$clonetype_number" | awk '{{print $2}}'`
			total_number=`cat {input.snp_dist} | cut -f 1 | grep -w "$clonetype" | wc -l`
			
			echo "The {params.sample} sample has less than {params.distance} SNP distance to the "$clonetype" clone type.
In total "$number" out of "$total_number" samples with this clone type fall below the SNP distance threshold." > {output.stats}
		elif [[ "$line_number" > 1 ]]
		then
			echo "WARNING! More than one clone type has smaller SNP distance than {params.distance} to the {params.sample}.
You might want to reconsider SNP distance threshold. If the sample is set to be included in the final SNP distance matrix, the most often occuring clone type will be used. It might cause problems in the future analysis. \The following clone types are reported to be close to the {params.sample} sample:" > {output.stats}
			while read clonetype_line
			do
				number=`echo "$clonetype_line" | awk '{{print $1}}'`
				clonetype=`echo "$clonetype_line" | awk '{{print $2}}'`
				total_number=`cat {input.snp_dist} | cut -f 1 | grep -w "$clonetype" | wc -l`
				echo "Clonetype "$clonetype" in "$number" out of "$total_number" samples with this clonetype which fall below the SNP distance threshold." >> {output.stats}
		
			done <<< "$clonetype_number"
	
		elif [[ "$line_number" == 0 ]]
		then
			echo "No already existing clone types are closer than the {params.distance} SNP distance threshold.
New clone type assigned: {params.new_clonetype}" > {output.stats}
			clonetype={params.prefix}{params.new_clonetype}
		fi

		# if include is true change the name of complete snp distance and sorted alignment files to the name of the assigned clone type
		if [[ {params.incl} == "True" ]]
		then
			sed -i "s/^>{params.sample}/>\"$clonetype\"::{params.sample}/" {input.initial_alignment}
			sed -i "s/{params.sample}/\"$clonetype\"::{params.sample}/" {input.initial_distance}
		fi
		"""


#rule visualize_snp_dists:
