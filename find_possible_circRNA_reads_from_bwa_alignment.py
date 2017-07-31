# This script accepts a sam file for reads aligned to the "backsplicing database" using BWA. The script will search for
# reads that map acroos any of the backsplicing junctions. The backsplicing junction database was generated for all
# combinations of backsplicing exons. The database currently does not contain intron information, only exon-exon
# backsplicing junctions are currently available.

# TODO: Analyze soft/hard clipped reads. For now excluding them from analysis.
# TODO: Add analysis for reads showing insertions and deletions.
# TODO: Add analysis for secondary alignments. For now, only considering primary alignments.

# Import required libraries
import pandas
import numpy
import pybedtools
import os
from argparse import ArgumentParser


def parse_sam_file (samfile):
    samfile = pandas.read_csv(samfile, sep='\t', comment='@', names=['rname','flag','junc_name','pos','mapq','cigar',
                                                                     'rnext','pnext','tlen','seq','qual','AS','XS'])
    # Sam file flag 4 reads are unmapped. So, excluding such reads.
    samfile = samfile[samfile['flag'] != 4]
    # Annotating transcript across which each read was mapped.
    samfile['tx'] = map(lambda x: x.split('__')[0], samfile['junc_name'])
    # Generate read and transcript key.
    samfile['rname-tx'] = samfile['rname'] +'-'+ samfile['tx']
    # Extract reads that show soft or hard splicing or insertions or deletions. For now not considering such reads. Need
    #  to analyze them further
    X = samfile[samfile['cigar'].str.contains('S') | samfile['cigar'].str.contains('H') | samfile['cigar'].str.
        contains('I') | samfile['cigar'].str.contains('D')]
    samfile = samfile.loc[set(samfile.index).difference(set(X.index))]
    # Keep only primary alignments. Exclusing secondary alignmens for now
    samfile = samfile[(samfile['flag'] == 0) | (samfile['flag'] == 16)]
    # Search for reads that start mapping after position '25' or those that finish mapping before position '25'. Such
    # reads essentially map to only one exon of the junction, so excluding them from analysis.
    samfile['match_len'] = samfile['cigar'].str.extract('(\d\d)', expand=True)
    samfile['match_end_pos'] = samfile['pos'].astype(int) + samfile['match_len'].astype(int)
    samfile = samfile[(samfile['match_end_pos'] > 25) & (samfile['pos'] < 25)]
    return samfile


def count_backsplicing_junc_coverage (samfile):
    return samfile['junc_name'].value_counts()


def main():
    parser = ArgumentParser(description="Find reads that potentially map across backsplicing junctions.")
    parser.add_argument("--alignment_file", help="Full path and name of sam file")
    parser.add_argument("--output_dir", help="Full path to output directory")
    args = parser.parse_args()

    if os.path.isfile(args.alignment_file):
        print("Found Sam File...\n")
        print("Reading bed file...\n")
        samfile_filtered = parse_sam_file(args.alignment_file)
        backsplicing_junc_cov = count_backsplicing_junc_coverage(samfile_filtered)
        backsplicing_junc_cov.columns = ['Coverage']
        print("Found %d possible reads that map to backspicing junctions."%len(set(samfile_filtered.rname)))
        print("Found %d possible backspicing junctions."%len(backsplicing_junc_cov.index))
        output_file = args.output_dir+'/Backsplicing_Junction_Coverages.txt'
        backsplicing_junc_cov.to_csv(output_file, sep='\t')
        print("Coverage for each backsplicing junction found written to following file: %s"%output_file)


if __name__ == '__main__':
    main()

