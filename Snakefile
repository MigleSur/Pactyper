

rule all:
	input:
		vcf=expand("results/raw_vcf_calls/{{sample}}/{{sample}}.vcf.gz", sample=IDS)
	message:
		"Pipeline complete"

# Run snippy on all input samples using the core reference genome
rule run_snippy:
	input:
		R1=expand("{input_dir}/{{sample}}_R1_001.fastq.gz", input_dir=config["sample_dir"]),
                R2=expand("{input_dir}/{{sample}}_R2_001.fastq.gz", input_dir=config["sample_dir"])
	output:
		vcf=expand("results/raw_vcf_calls/{{sample}}/{{sample}}.vcf.gz"),
                bam=expand("results/raw_vcf_calls/{{sample}}/{{sample}}.bam"),
                ref=expand("results/raw_vcf_calls/{{sample}}/reference/ref.fa")
	params:
		outdir="results",
                ref=config["ref"],
                prefix="{sample}",
                cpus=10,
                mincov=1,
                minqual=50
	message:
		"Running snippy variant calling for {params.prefix}."
	log:
		expand("results/logs/snippy_{{sample}}.log",outdir=config["output_dir"])
	shell:
		"""
                (snippy --outdir {params.outdir}/raw_vcf_calls/{params.prefix} --ref {params.ref} --R1 {input.R1} --R2 {input.R2} --prefix {params.prefix} --cpus {params.cpus}  --force --mincov {params.mincov} --minqual {params.minqual}) 2> {log}
		"""


# SNP distance analysis on the provided samples
rule snp_dists:
	input:
	output:
	message:
		"Performing SNP distance analysis on the {params.sample}."
	shell:

# Estimation of the closest clone types. If the SNP difference is <5000 SNPs, the isolate is given a clone type. Otherwise, it is considered a new clone type.
rule estimate_clonetype:

rule visualize_snp_dists:
