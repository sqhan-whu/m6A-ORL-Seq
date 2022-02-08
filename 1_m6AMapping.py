#!/usr/bin/env runsnakemake
configfile: "config.yaml"

samples = config["samples"]
threads = config["threads"]
rRNA_database = config["rRNA_database"]
star_index = config["star_index"]
path = config["path"]

indir = path + "/fq"
outdir = path + "/results"
seqtk_fq1_e = 10
seqtk_fq2_b = 10
seqtk_fq2_e = 0

rule all:
	input:
		expand(outdir + "/bam/{sample}.filter.sorted.bam", sample = samples),
		expand(outdir + "/bam/{sample}.filter.sorted.bam.bai", sample = samples),
		expand(outdir + "/bam/{sample}.filter.sorted.bam.bw", sample = samples)

rule trim_galore:
	input:
		indir + "/{sample}_R1.fq.gz",
		indir + "/{sample}_R2.fq.gz"
	output:
		outdir + "/clean_fq/{sample}_R1_val_1.fq",
		outdir + "/clean_fq/{sample}_R2_val_2.fq"
	params:
		outdir + "/clean_fq"
	threads:
		threads
	shell:
		"trim_galore -j {threads} --quality 20 --phred33 --stringency 3 --length 30 --dont_gzip -o {params} --paired {input}"

rule echo_uniq:
	input:
		fq1 = outdir + "/clean_fq/{sample}_R1_val_1.fq",
		fq2 = outdir + "/clean_fq/{sample}_R2_val_2.fq"
	output:
		outdir + "/clean_fq/{sample}_uniq_list.txt"
	shell:
		"echo \"{input.fq1}\n{input.fq2}\" > {output}"
rule fastuniq:
	input:
		outdir + "/clean_fq/{sample}_uniq_list.txt"
	output:
		fq1 = temp(outdir + "/clean_fq/{sample}_uniq_1.fq"),
		fq2 = temp(outdir + "/clean_fq/{sample}_uniq_2.fq")
	shell:
		"fastuniq -i {input} -o {output.fq1} -p {output.fq2}"
rule seqtk_trimfq1:
	input:
		outdir + "/clean_fq/{sample}_uniq_1.fq"
	output:
		temp(outdir + "/clean_fq/{sample}_seqtk_1.fq")
	shell:
		"seqtk trimfq -e {seqtk_fq1_e} {input} > {output}"
rule seqtk_trimfq2:
	input:
		outdir + "/clean_fq/{sample}_uniq_2.fq"
	output:
		temp(outdir + "/clean_fq/{sample}_seqtk_2.fq")
	shell:
		"seqtk trimfq -b {seqtk_fq2_b} -e {seqtk_fq2_e} {input} > {output}"
rule rRNA_remove:
	input:
		fq1 = outdir + "/clean_fq/{sample}_seqtk_1.fq",
		fq2 = outdir + "/clean_fq/{sample}_seqtk_2.fq"
	output:
		fq1 = outdir + "/clean_fq/{sample}_rmrRNA_1.fq",
		fq2 = outdir + "/clean_fq/{sample}_rmrRNA_2.fq",
		sam = temp(outdir + "/clean_fq/{sample}_rmrRNA.sam")
	params:
		outdir + "/clean_fq/{sample}_rmrRNA.fq"
	threads:
		threads
	shell:
		"bowtie -p {threads} {rRNA_database} -1 {input.fq1} -2 {input.fq2} --un {params} -S {output.sam}"

rule star_mapping:
	input:
		outdir + "/clean_fq/{sample}_rmrRNA_1.fq",
		outdir + "/clean_fq/{sample}_rmrRNA_2.fq"
	output:
		outdir + "/bam/{sample}.Aligned.sortedByCoord.out.bam"
	params:
		outdir + "/bam/{sample}."
	threads:
		threads
	shell:
		"STAR --genomeDir {star_index} --readFilesIn {input} --runThreadN {threads} --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --readFilesCommand cat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params}"

rule filter_bam:
	input:
		outdir + "/bam/{sample}.Aligned.sortedByCoord.out.bam"
	output:
		temp(outdir + "/bam/{sample}.filter.bam")
	threads:
		threads
	shell:
		#"bamtools filter -in {input} -out {output} -isPrimaryAlignment true -isProperPair true -isMapped true -mapQuality '>1'"
		"samtools view -@ {threads} -bS -f 3 -F 4 -F 8 -F 1024 -q 5 -o {output} {input}"

rule bam_sort:
	input:
		outdir + "/bam/{sample}.filter.bam"
	output:
		outdir + "/bam/{sample}.filter.sorted.bam"
	threads:
		threads
	shell:
		"samtools sort -@ {threads} -o {output} {input}"

rule bam_index:
	input:
		outdir + "/bam/{sample}.filter.sorted.bam"
	output:
		outdir + "/bam/{sample}.filter.sorted.bam.bai"
	threads:
		threads
	shell:
		"samtools index -@ {threads} {input}"

rule bam_bw:
	input:
		outdir + "/bam/{sample}.filter.sorted.bam"
	output:
		outdir + "/bam/{sample}.filter.sorted.bam.bw"
	shell:
		"bamCoverage --bam {input} -o {output} --binSize 1 --normalizeUsing CPM"
