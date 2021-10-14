## /scratch/sqhan/11.xyl/20210702/results/sites
## time: 2021.10.14
## usage: python get_bam_stop_rate.py total.pos.filter.bam ~/project/project/00.DATABASE/hg38/bowtie2_index/hg38.fa

import sys
import pysam
min_value = 5
rate = 0.2
motif = {'GGACC': 6, 'CCACC': 6, 'CGACC': 6, 'GCACC': 6, 
                'GCACT': 6, 'CCACT': 6, 'GGACT': 6,'CGACT': 6,
                'GCACA': 6, 'CCACA': 6, 'GGACA': 6,'CGACA': 6}

def complement(s):
	basecomplement = {
         "A":"T",
          "T":"A",
          "G":"C",
          "C":"G",
          "a":"t",
          "t":"a",
          "g":"c",
          "c":"g",
			"N":"N",}
	letters = list(s)
	letters = [basecomplement[base] for base in letters]
	letters = letters[::-1]

	return ''.join(letters)


def get_bam_startsites(bam_file, genome):
	out_dict = []
	fasta = pysam.FastaFile(genome)
	n = 0
	dict_start = dict()
	bf = pysam.AlignmentFile(bam_file, 'rb')
	for r in bf:
		k = 0
		for mark, num in r.cigartuples:
			if mark == 0:
				k +=num
			if mark == 1:
				k -=num
			if mark == 2:
				k +=num
		content = str(r.reference_name)+'\t'+str(r.pos+1+k)
		
		dict_start.setdefault(content,0)
		dict_start[content] +=1
	
	for k, v in dict_start.items():
		sequence = ''
		k = k.split('\t')
		if v >= min_value:
			cover = bf.count(contig=k[0], start=int(k[1])-2, stop=int(k[1]))
			if v/cover > rate:
				sequence = fasta.fetch(str(k[0]),int(k[1])-3,int(k[1])+2)
				sequence = complement(sequence)
				if sequence in motif.keys():
					print(str(k[0])+'\t'+str(k[1])+"\t"+str(v)+"\t"+str(cover)+'\t'+str(v/cover)+'\t'+sequence)
bam_file,genome = sys.argv[1:]

get_bam_startsites(bam_file,genome)
