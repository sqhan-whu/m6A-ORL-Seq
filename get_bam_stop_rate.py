## author: Han Shaoqing
## time: 2021.10.14
## usage: python get_bam_stop_rate.py total.pos.filter.bam ~/project/project/00.DATABASE/hg38/bowtie2_index/hg38.fa

import sys
import pysam
min_value = 5
rate = 0.2
motif = {'GGACC': 6, 'CCACC': 6, 'CGACC': 6, 'GCACC': 6, 
                'GCACT': 6, 'CCACT': 6, 'GGACT': 6,'CGACT': 6,
                'GCACA': 6, 'CCACA': 6, 'GGACA': 6,'CGACA': 6}

def get_bam_startsites(bam_file, genome):
	out_dict = []
	fasta = pysam.FastaFile(genome)
	n = 0
	dict_start = dict()
	bf = pysam.AlignmentFile(bam_file, 'rb')
	for r in bf:
		content = str(r.reference_name)+'\t'+str(r.pos)
		dict_start.setdefault(content,0)
		dict_start[content] +=1
	
	for k, v in dict_start.items():
		sequence = ''
		k = k.split('\t')
		if v >= min_value:

			cover = bf.count(contig=k[0], start=int(k[1]), stop=int(k[1])+1)
			if v/cover > rate:
				sequence = fasta.fetch(str(k[0]),int(k[1])-3,int(k[1])+2)

				if sequence in motif.keys():
					print(str(k[0])+'\t'+str(k[1])+"\t"+str(v)+"\t"+str(cover)+'\t'+str(v/cover)+'\t'+sequence)
          
bam_file,genome = sys.argv[1:]
get_bam_startsites(bam_file,genome)
