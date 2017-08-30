# This script accepts a sam file for reads aligned to the "backsplicing database" using BWA. The script will search for
# reads that map acroos any of the backsplicing junctions. The backsplicing junction database was generated for all
# combinations of backsplicing exons. The database currently does not contain intron information, only exon-exon
# backsplicing junctions are currently available.

# TODO: Analyze soft/hard clipped reads. For now excluding them from analysis.
# TODO: Add analysis for reads showing insertions and deletions.
# TODO: Add analysis for secondary alignments. For now, only considering primary alignments.
# TODO: Add function to check against known circular RNAs. Bed file has been downloaded. Use bedtools to intersect.

# Import required libraries
import pandas
import numpy
import pybedtools
import os
from pybedtools import *
from argparse import ArgumentParser


def parse_sam_file(samfile, cut_off, sample):
    """
    This function will parse the input sam file and search for reads mapping to backsplicing junctions. Unmapped reads
    are exluded from analysis. For now, only primary alignments are being considered. Also, only those reads that span
    across both exons of a backsplicing junctions are considered. The rest are excluded from all analysis.
    :param samfile: BWA alignment sam file. This has to be the alignment file for the backsplicing database.
    :param cut_off: Read cut-off for filtering junction. Unless given, cut-off is set at 10 reads.
    :param sample: Type of sample. IP or INPUT. If INPUT, no filtering of junctions is performed.
    :return: Filtered sam file that contains only those reads that pass all parameters.
    """
    # TODO: fix this, different SAM files might have different columns. Look at pysam to easily parse this
    samfile = pandas.read_table(samfile, sep='\t', comment='@', names=[
        'rname', 'flag', 'junc_name', 'pos', 'mapq', 'cigar',
        'rnext', 'pnext', 'tlen', 'seq', 'qual', 'XT', 'NM',
        'X0', 'X1', 'XM', 'XO', 'XG', 'MDZ'
    ])
    # Sam file flag 4 reads are unmapped. So, excluding such reads.
    samfile = samfile[samfile['flag'] != 4]
    # Annotating transcript across which each read was mapped.
    samfile['tx'] = map(lambda x: x.split('__')[0], samfile['junc_name'])

    # Generate read and transcript key.
    samfile['rname-tx'] = samfile['rname'] + '-' + samfile['tx']

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
    samfile = samfile[(samfile['match_end_pos'] > 35) & (samfile['pos'] <= 25)]

    # Remove reads that have a mapping quality less than 20 ( -log10(0.01) ). Also remove reads that have a mapping
    # quality of 255. 255 means that a mapping quality score is not available for the read.
    samfile = samfile[(samfile['mapq'] >= 20) & (samfile['mapq'] != 255)]

    # Remove duplicate reads. Reads are considered duplicates if they have the same start and end position on the same
    # junction
    samfile = check_for_duplicate_reads(samfile)

    # Count coverage across all backsplicing junctions. Coverage == # of reads spanning that junction
    coverage = count_backsplicing_junc_coverage(samfile)

    # Filter junctions based on read cut-off. This is done only for IP sample.
    if sample == 'ip':
        coverage = coverage[coverage >= cut_off]
        samfile = samfile[samfile['junc_name'].isin(list(coverage.index))]

    if len(samfile) == 0:
        return samfile, coverage, 0
    else:
        # Create bed file from all filtered junctions
        bed_file = create_bed_for_junctions(set(samfile.junc_name))
        # Intersect bed file with known circular RNA taken from circBase (http://www.circbase.org/)
        known = bed_file.intersect('/home/shsathe/backsplicing_exons_db/hsa_hg19_circRNA.bed', u=True)
        return samfile, coverage, known


def create_bed_for_junctions(juncs):
    """
    This function will create a pybedtools bed object. Each record in the bed object will be a unique exon / intron
    involved in backsplicing.
    :param juncs: A set of backsplicing junctions. Each junction contains the coordinates for the start / stop regions
    of the exons involved in forming the junction
    :return: pybedtools bed object for the exons involved in forming backsplicing junctions
    """
    df = pandas.DataFrame(columns=['chr', 'start', 'stop', 'gene', 'num', 'strand'])
    for i in juncs:
        for k in i.split('__and__'):
            j = k.split('__')
            df.loc[i, 'chr'] = j[0]
            df.loc[i, 'start'] = j[1]
            df.loc[i, 'stop'] = j[2]
            df.loc[i, 'gene'] = j[3]
            df.loc[i, 'num'] = j[4]
            df.loc[i, 'strand'] = j[5]
    return BedTool.from_dataframe(df)


def check_for_duplicate_reads(samfile):
    """
    This function will check for duplicate reads and remove them from the samfile. The function will return a samfile
    that does not contain any duplicate reads. Duplicate reads have the same start, strand, same cigar string, same
    sequence. When duplicate reads are encountered, the function will keep one copy of the reads, while discard the
    others.
    :param samfile: The filtered sam file from the 'parse_sam_file' function.
    :return: samfile without duplicate reads
    """
    samfile['key'] = samfile['junc_name'] + '_' + samfile['flag'].astype(str) + '_' + samfile['pos'].astype(str) + '_' + \
                     samfile['cigar'].astype(str)
    samfile = samfile.drop_duplicates('key')
    return samfile.drop('key', 1)


def count_backsplicing_junc_coverage(samfile):
    """
    This function counts the number of reads mapped across each backsplicing junction in the alignment file.
    :param samfile: The filtered sam file from the 'parse_sam_file' function.
    :return: Dataframe for coverge across each junction
    """
    return samfile['junc_name'].value_counts()


def normalize_by_input(ip, inpt):
    """
    This function accepts junction coverage files for IP and Input and normalizes each junction read count in IP file
    with the corresponding read count in the Input sample
    :param ip: IP junction coverage file
    :param inpt: Input junction coverage file
    :return: Dataframe with number of IP and Input reads for each junction and the normalized IP read counts.
    """
    df = pandas.DataFrame(index=list(ip.index), columns=['IP', 'INPUT', 'IP/INPUT'])
    df['IP'] = ip
    df.loc[set(ip.index).intersection(set(inpt.index)), 'INPUT'] = inpt.loc[set(ip.index).intersection(set(inpt.index))]
    df['IP/INPUT'] = df['IP'].div(df['INPUT'])
    df['INPUT'] = df['INPUT'].fillna(0)
    return df


def main():
    parser = ArgumentParser(description="Find reads that potentially map across backsplicing junctions.")
    parser.add_argument("--ip", help="Full path and name of sam file for IP sample (Required)", default='')
    parser.add_argument("--input", help="Full path and name of sam file INPUT sample")
    parser.add_argument("--output_dir", help="Full path to output directory. Default is same directory as IP sample "
                                             "(Required)")
    parser.add_argument("--cut_off", help="Cut-Off for filtering junctions. Default set to 10 reads", default=10)
    parser.add_argument("--type", help="Type of junctions. exon-exon or intron-intron")
    parser.add_argument("--known_circRNA", help="bed file of known circRNA positions")
    args = parser.parse_args()

    # Check if IP sam file exists
    if os.path.isfile(args.ip):
        print("Found IP Sam File...\n")
        ip_filtered, ip_backsplicing_junc_cov, ip_backsplicing_junc_bed_file = parse_sam_file(args.ip,
                                                                                              args.cut_off, 'ip')

        # Check if any junctions passed the 10 reads cut-off.
        if ip_backsplicing_junc_bed_file == 0:
            print("No backsplicing %s junctions were found in the IP sample." % args.type)
        else:
            if args.input == '':
                print('INPUT SAM file not provided. Junction read counts will not be input normalized.')
            # If input sam file is given, parse input sam file and normalized ip with input
            elif os.path.isfile(args.input):
                print("Found INPUT Sam File. Input normalization of IP sample will be performed")
                input_samfile_filtered, input_backsplicing_junc_cov, input_backsplicing_junc_bed_file = parse_sam_file(
                    args.input, args.cut_off, 'input')
                ip_normalized_juncs = normalize_by_input(ip_backsplicing_junc_cov, input_backsplicing_junc_cov)
                output_file = args.output_dir + '/IP_Backsplicing_Junctions_Coverage_Input_Normalized.txt'
                ip_normalized_juncs.to_csv(output_file, sep='\t')
                print("INPUT Normalized Junctions written to following file: %s" % output_file)
            else:
                parser.error("INPUT SAM file not found. Enter correct path and name")

            output_file = args.output_dir + '/Filtered_Sam_File.txt'
            ip_filtered.to_csv(output_file, sep='\t')
            print("Filter sam file written to following file: %s" % output_file)

            output_file = args.output_dir + '/Known_Backsplicing_Junctions.bed'
            ip_backsplicing_junc_bed_file.saveas(output_file)
            print("Bed file for known backsplicing %ss written to following file: %s" % (args.type.split('-')[0],
                                                                                         output_file))

            print("Found %d possible reads that map to backspicing junctions." % len(set(ip_filtered.rname)))
            print("Found %d possible backspicing junctions." % len(ip_backsplicing_junc_cov.index))
            output_file = args.output_dir + '/Backsplicing_Junction_Coverages.txt'
            ip_backsplicing_junc_cov.to_csv(output_file, sep='\t')
            print("Coverage for each backsplicing junction found written to following file: %s" % output_file)

    else:
        parser.error("IP SAM file not found. Enter correct path and name.")


if __name__ == '__main__':
    main()
