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
    """
    This function will parse the input sam file and search for reads mapping to backsplicing junctions. Unmapped reads
    are exluded from analysis. For now, only primary alignments are being considered. Also, only those reads that span
    across both exons of a backsplicing junctions are considered. The rest are excluded from all analysis.
    :param samfile: Input BWA alignment sam file. This has to be the alignment file for the backsplicing database.
    :return: Filtered sam file that contains only those reads that pass all parameters.
    """
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

    # Keep only primary alignments. Excluding secondary alignmens for now
    samfile = samfile[(samfile['flag'] == 0) | (samfile['flag'] == 16)]

    # Search for reads that start mapping after position '36' or those that finish mapping before position '36'. Such
    # reads essentially map to only one exon of the junction, so excluding them from analysis.
    samfile['match_len'] = samfile['cigar'].str.extract('(\d\d)', expand=True)
    samfile['match_end_pos'] = samfile['pos'].astype(int) + samfile['match_len'].astype(int)
    samfile = samfile[(samfile['match_end_pos'] >= 40) & (samfile['pos'] <= 31)]

    # Remove reads that have a mapping quality less than 20 ( -log10(0.01) ). Also remove reads that have a mapping
    # quality of 255. 255 means that a mapping quality score is not available for the read.
    samfile = samfile[(samfile['mapq'] >= 20) & (samfile['mapq'] != 255)]

    # Remove duplicate reads. Reads are considered duplicates if they have the same start and end position on the same
    # junction
    samfile = check_for_duplicate_reads(samfile)

    coverage = count_backsplicing_junc_coverage(samfile)
    return samfile, coverage


def check_for_duplicate_reads(samfile):
    samfile['key'] = samfile['junc_name']+'_'+samfile['flag'].astype(str)+'_'+samfile['pos'].astype(str)+'_'+\
                     samfile['cigar'].astype(str)+'_'+samfile['seq']
    samfile = samfile.drop_duplicates('key')
    return samfile.drop('key',1)


def count_backsplicing_junc_coverage (samfile):
    """
    This function counts the number of reads mapped across each backsplicing junction in the alignment file.
    :param samfile: The filtered sam file from the 'parse_sam_file' function.
    :return: Dataframe for coverge across each junction
    """
    return samfile['junc_name'].value_counts()


def main():
    parser = ArgumentParser(description="Find reads that potentially map across backsplicing junctions.")
    parser.add_argument("--alignment_file", help="Full path and name of sam file")
    parser.add_argument("--output_dir", help="Full path to output directory")
    args = parser.parse_args()

    if os.path.isfile(args.alignment_file):
        print("Found Sam File...\n")
        samfile_filtered, backsplicing_junc_cov = parse_sam_file(args.alignment_file)
        output_file = args.output_dir+'/Filtered_Sam_File.txt'
        samfile_filtered.to_csv(output_file, sep='\t')
        print("Filter sam file written to following file: %s"%output_file)

        backsplicing_junc_cov.columns = ['Coverage']
        print("Found %d possible reads that map to backspicing junctions."%len(set(samfile_filtered.rname)))
        print("Found %d possible backspicing junctions."%len(backsplicing_junc_cov.index))
        output_file = args.output_dir+'/Backsplicing_Junction_Coverages.txt'
        backsplicing_junc_cov.to_csv(output_file, sep='\t')
        print("Coverage for each backsplicing junction found written to following file: %s"%output_file)


if __name__ == '__main__':
    main()

