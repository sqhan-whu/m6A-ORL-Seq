#!/usr/bin/env runsnakemake
configfile: "config.yaml"
indir = "results"
outdir = "results/final_sites"
script_dir = "~/bin_python/script"
GENOME_fa = "~/hg38/bowtie2_index/hg38.fa"
rule all:
	input:
		indir + "/sites/IN/input.sites",
		indir + "/sites/ip.sites",
		outdir + "/m6A.sites",
		outdir + "/m6A.sites.ACA"


rule ip_sites_pos:
	input:
		indir + "/sites/sample.pos.filter.bam"
	output:
		indir + "/sites/ip.sites.pos"
	params:
		script_dir + "/get_bam_stop_rate.py"
	shell:
		"python {params} {input} {GENOME_fa} > {output}"

rule ip_sites_neg:
	input:
		indir + "/sites/sample.neg.filter.bam"
	output:
		indir + "/sites/ip.sites.neg"
	params:
		script_dir + "/get_bam_stop_rate_neg.py"
	shell:
		"python {params} {input} {GENOME_fa} > {output}"

rule ip_sites:
	input:
		indir + "/sites/ip.sites.pos",
		indir + "/sites/ip.sites.neg"
	output:
		indir + "/sites/ip.sites"
	shell:
		"cat {input} |sort -k1,1 -k2,2n > {output}"

rule input_sites_pos:
	input:
		indir + "/sites/IN/sample.pos.filter.bam"
	output:
		indir + "/sites/IN/input.sites.pos"
	params:
		script_dir + "/get_bam_stop_rate.py"
	shell:
		"python {params} {input} {GENOME_fa} > {output}"

rule input_sites_neg:
	input:
		indir + "/sites/IN/sample.neg.filter.bam"
	output:
		indir + "/sites/IN/input.sites.neg"
	params:
		script_dir + "/get_bam_stop_rate_neg.py"
	shell:
		"python {params} {input} {GENOME_fa} > {output}"

rule input_sites:
	input:
		indir + "/sites/IN/input.sites.pos",
		indir + "/sites/IN/input.sites.neg"
	output:
		indir + "/sites/IN/input.sites"
	shell:
		"cat {input} |sort -k1,1 -k2,2n > {output}"

rule overlap_sites:
	input:
		f1 = indir + "/sites/ip.sites",
		f2 = indir + "/sites/IN/input.sites",
	output:
		outdir + "/m6A.sites"
	params:
		script = script_dir + "/overlap.py",
		uniq_f2 = indir + "/sites/IN/input.sites",
		overlap = indir + "/sites/IN/input.sites.overlap.ip.txt"
	shell:
		"python {params.script} {input.f1} {input.f2} {output} {params.uniq_f2} {params.overlap}"

rule get_ACA_sites:
	input:
		outdir + "/m6A.sites"
	output:
		outdir + "/m6A.sites.ACA"
	shell:
		"grep -P \".*\\t..ACA\" {input} > {output}"
