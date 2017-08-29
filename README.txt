circRNA_eCLIP

The purpose of this repo is to find reads mapping to circular RNAs within an eCLIP dataset. 

Backsplicing Junctions Database

A customized "backsplicing junctions" database was created, which contains the sequences of all possible backsplicing junctions in the human genome ( hg19 ). A backsplicing junction was formed by merging 30bp from the 3' end of the downstream exon and 5` end of the upstream exon. This was done for all exons within a gene, excluding the first and last exons of a gene.
You can create your own backsplicing junctions database from the FASTA file using the 'create_backsplicing_db_from_annotation.py' script.Currently it takes roughly 3 - 4 hours to create the FASTA file for backsplicing junctions
Index the junction database using 'bwa index -p <prefix> <in.fasta>

BWA Alignment

The fastq reads from the eCLIP experiment first need to be aligned to the backsplicing junction database using bwa. The following
command can be used:
bwa aln <db_prefix> <in.fq>
The previous command will product gapped / ungapped alignments for each read. The 'bwa samse' command will convert the previous 
output into a sam file.
bwa samse -f <out.sam> <db_prefix> <in.sai> <in.fq>

Identify Reads Mapping To Backsplicing Junctions

The final step is to execute the "find_possible_circRNA_reads_from_bwa_alignment.py" script on the resulting sam file. The script will extract reads that map to backsplicing junctions and also quatify the number of reads mapping to each junction in the database. Junctions with at least 10 reads are reported.
If an Input sample is provided, the script withh also perform input normalization on the junctions found in the IP sample.
