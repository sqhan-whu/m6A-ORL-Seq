## m6A-ORL-Seq
### Content
m6A-ORL-Seq includes following steps:

1_m6AMapping.py: fastuniq and seqtk trimfq were used before map the reads to the genome.

2_callpeaks_exomePeak.R: m6A peaks calling method were implemented through the R package exomepeak.

3_SnakeSites.py: Distinguish between positive and negative chains, it is used to distinguish between positive and negative chains, and obtain the Stop sites on the positive and negative chains respectively.

HOMERMotif.sh:findMotifsGenome.pl is used to calculate m6A motifs.

get_bam_stop_rate.py and get_bam_stop_rate_neg.py: screen RT stop sites from filtered bam files.

### Raw and processed data files for this study:

Raw data files are available at NCBI Gene Expression Omnibus (GEO) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185753).
