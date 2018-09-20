# This script accepts a sam file for reads aligned to the "backsplicing database" using BWA. The script will search for
# reads that map across any of the backsplicing junctions. The backsplicing junction database was generated for all
# combinations of backsplicing exons. The database currently does not contain intron information, only exon-exon
# backsplicing junctions are currently available.

# TODO: Analyze soft/hard clipped reads. For now excluding them from analysis.
# TODO: Add analysis for reads showing insertions and deletions.
# TODO: Add analysis for secondary alignments. For now, only considering primary alignments.


# Import required libraries
from __future__ import division
import pandas
import numpy
import os
from pybedtools import *
from argparse import ArgumentParser
from ntpath import basename
import json
from Bio import SeqIO


def parse_sam_file(s, cut_off, sample):
    """
    This function will parse the input sam file and search for reads mapping to backsplicing junctions. Unmapped reads
    are excluded from analysis. For now, only primary alignments are being considered. Also, only those reads that span
    across both exons of a backsplicing junctions are considered. The rest are excluded from all analysis.
    :param samfile: BWA alignment sam file. This has to be the alignment file for the backsplicing database.
    :param cut_off: Read cut-off for filtering junction. Unless given, cut-off is set at 10 reads.
    :param sample: Type of sample. IP or INPUT. If INPUT, no filtering of junctions is performed.
    :return: Filtered sam file that contains only those reads that pass all parameters.
    """
    # TODO: fix this, different SAM files might have different columns. Look at pysam to easily parse this
    samfile = pandas.read_csv(s, sep='\t', comment='@', index_col=False, header=None, engine='python')

    # Sam file flag 4 reads are unmapped. So, excluding such reads.
    #samfile = samfile[samfile['flag'] != 4]

    # Annotating transcript across which each read was mapped.
    samfile['tx'] = map(lambda x: x.split('__')[3], samfile[2])

    # Generate read and transcript key.
    samfile['rname-tx'] = samfile[0] + '-' + samfile['tx']

    # Extract reads that show soft or hard splicing or insertions or deletions. For now not considering such reads. Need
    #  to analyze them further
    X = samfile[samfile[5].str.contains('S') | samfile[5].str.contains('H') | samfile[5].str.contains('I') |
                samfile[5].str.contains('D')]
    samfile = samfile.loc[set(samfile.index).difference(set(X.index))]

    # TODO: remove this step. probably not required with the new processing pipeline. For now, commented it out.
    # Keep only primary alignments. Excluding secondary alignmens for now
    #print set(samfile['flag'])
    #print samfile.head()
    #samfile = samfile[(samfile['flag'] == 0) | (samfile['flag'] == 16)]

    # Search for reads that start mapping after position '36' or those that finish mapping before position '36'. Such
    # reads essentially map to only one exon of the junction, so excluding them from analysis.
    samfile['match_len'] = samfile[5].str.extract('(\d\d)', expand=True)
    samfile['match_end_pos'] = samfile[3].astype(int) + samfile['match_len'].astype(int)
    samfile = samfile[(samfile['match_end_pos'] > 35) & (samfile[3] <= 25)]

    # Remove reads that have a mapping quality less than 20 ( -log10(0.01) ). Also remove reads that have a mapping
    # quality of 255. 255 means that a mapping quality score is not available for the read.
    samfile = samfile[(samfile[4] >= 20) & (samfile[4] != 255)]

    # Remove duplicate reads. Reads are considered duplicates if they have the same start and end position on the same
    # junction
    #samfile = check_for_duplicate_reads(samfile)

    # Count coverage across all backsplicing junctions. Coverage == # of reads spanning that junction
    coverage = count_backsplicing_junc_coverage(samfile)

    # Filter junctions based on read cut-off. This is done only for IP backsplicing sample.
    if sample == 'ip_back':
        coverage = coverage[coverage >= cut_off]
        samfile = samfile[samfile[2].isin(list(coverage.index))]

    return samfile, coverage


def intersect_with_known_juncs(bed):
    # Intersect bed file with known circular RNA taken from circBase (http://www.circbase.org/)
    return bed_file.intersect('/home/shsathe/backsplicing_exons_db/hsa_hg19_circRNA.bed', u=True)


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
    samfile['key'] = samfile[2] + '_' + samfile[1].astype(str) + '_' + samfile[3].astype(str) + '_' + \
                     samfile[5].astype(str)
    samfile = samfile.drop_duplicates('key')
    return samfile.drop('key', 1)


def count_backsplicing_junc_coverage(samfile):
    """
    This function counts the number of reads mapped across each backsplicing junction in the alignment file.
    :param samfile: The filtered sam file from the 'parse_sam_file' function.
    :return: Dataframe for coverge across each junction
    """
    return samfile[2].value_counts()


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


def diff_exp_of_back_junc(back_cov, linear_cov):
    """
    This function accepts junctions coverage dataframes for back and linear splicing junctions. The function will
    iterate over each back-splicing junction to find a corresponding linear junctions. If the coverage for the
    corresponding linear junction is available, ratio of back / linear junction counts will be calculated as a measure
    of differential expression of the back-splicing junction
    :param back_cov: Coverage dataframe for back-splicing junctions
    :param linear_cov: Coverage dataframe for linear-splicing junctions
    :return df: DataFrame containing differential expression ratios of back-splicing junctions
    :return d: A dictionary specifying corresponding linear junctions for each backsplicing junction
    """
    # Parse linear junctions to extract transcript and exon number information
    linear_cov_df = pandas.DataFrame(index=linear_cov.index, columns=['cov','tx','exon_1','exon_2','junc'])
    linear_cov_df['cov'] = linear_cov
    linear_cov_df['tx'] = map(lambda x: x.split('__')[3], linear_cov_df.index)
    linear_cov_df['exon_1'] = map(lambda x: x.split('__')[4], linear_cov_df.index)
    linear_cov_df['exon_2'] = map(lambda x: x.split('__')[10], linear_cov_df.index)
    linear_cov_df['junc'] = linear_cov_df['tx'] + '__' + linear_cov_df['exon_1'].astype(str) + '__' + linear_cov_df['tx'] \
                         + '__' + linear_cov_df['exon_2'].astype(str)
    linear_cov_df = linear_cov_df.set_index('junc')

    # Create empty dataframe to fill ratios into
    df = pandas.DataFrame(index=back_cov.index, columns=['IP_Back', 'IP_Linear', 'Back/Linear'])
    d = dict()
    # Iterate over each back-splicing junction in a for loop
    for i in back_cov.index:
        df.loc[i, 'IP_Back'] = back_cov.loc[i]
        x = i.split('__')
        # Find corresponding linear-splicing junction depending on the strand.
        if '+' in i:
            corresponding_linear_junc = x[3] + '__' + x[11] + '__' + x[3] + '__' + str(int(x[11]) + 1)
        elif '-' in i:
            corresponding_linear_junc = x[3] + '__' + str(int(x[4]) + 1) + '__' + x[3] + '__' + x[4]
        else:
            print(i)
        d[i] = corresponding_linear_junc
        # Calculate back / linear ratio if the coverage for corresponding linear junction is available
        if corresponding_linear_junc in list(linear_cov_df.index):
            df.loc[i,'IP_Linear'] = linear_cov_df.loc[corresponding_linear_junc, 'cov']
            df.loc[i, 'Back/Linear'] = back_cov.loc[i] / linear_cov_df.loc[corresponding_linear_junc, 'cov']
            print('%s -- %s' % (i, corresponding_linear_junc))
        else:
            df.loc[i, 'IP_Linear'] = 0
            df.loc[i, 'Back/Linear'] = numpy.nan
    return df, d


def generate_sam_index (inds, out_file):
    f = open(out_file, 'w')
    f.write('@HD\tVN:1.5\tSO:queryname\n')
    for i in inds:
        f.write('@SQ\tSN:%s\tLN:61\n'%i)
    f.close()


def fasta_file_for_mapped_junctions (inds, out_dir):
    fasta_sequences = SeqIO.parse(open('/home/shsathe/backsplicing_exons_db/v4/exon_exon/backsplicing_exons_seqs.fasta'), 'fasta')
    f = open(out_dir + '/mapped_junctions_sequences.fa', 'w')
    for fasta in fasta_sequences:
        if fasta.id in inds:
            f.write('>%s\n' %(fasta.id))
            f.write('%s\n' % (str(fasta.seq)))
    f.close()


def main():
    parser = ArgumentParser(description="Find reads that potentially map across backsplicing junctions.")
    parser.add_argument("--ip", help="Full path and name of sam file for IP sample (Required)")
    parser.add_argument("--ip_linear", help="Full path and name of sam file for IP sample mapped to linear junctions",
                        default='')
    parser.add_argument("--input", help="Full path and name of sam file INPUT sample", default='')
    parser.add_argument("--input_linear", help="Full path and name of sam file INPUT sample mapped to linear junctions",
                        default='')
    parser.add_argument("--output_dir", help="Full path to output directory. Default is same directory as IP sample "
                                             "(Required)")
    parser.add_argument("--cut_off", help="Cut-Off for filtering junctions. Default set to 5 reads", default=5, type=int)
    parser.add_argument("--type", help="Type of junctions. exon-exon or intron-intron", default='exon-exon')
    args = parser.parse_args()

    # Check if IP sam file exists
    if os.path.isfile(args.ip):
        print("Found IP Sam File...\n")
        fname = basename(args.ip).split('.')[0]
        ip_back_filtered, ip_back_coverage = parse_sam_file(args.ip, args.cut_off, 'ip_back')

        # Check if any junctions passed the 10 reads cut-off.
        if len(ip_back_filtered) == 0:
            print("No backsplicing %s junctions were found in the IP sample." % args.type)
        else:
            # Intersect junction passing read cut-off with known back-splicing exons. Outputn as bed file
            ip_back_known_bed_file = create_bed_for_junctions(set(ip_back_filtered[2]))
            output_file = args.output_dir + '/' + fname + '_Known_Backsplicing_Junctions.bed'
            ip_back_known_bed_file.saveas(output_file)
            print("Bed file for known backsplicing %ss written to following file: %s" % (args.type.split('-')[0],
                                                                                         output_file))

            if args.ip_linear == '':
                print('IP SAM file for linear junctions not provided. Differential analysis of back-splicing junctions '
                      'over linear junctions will not be performed')
            elif os.path.isfile(args.ip_linear):
                # Parse IP linear junctions samfile
                ip_linear_samfile, ip_linear_junc_cov = parse_sam_file(args.ip_linear, args.cut_off, 'ip_linear')

                # Calculate differential expression of IP back-splicing junctions compared to IP linear junctions
                ip_back_junc_diff_exp, d = diff_exp_of_back_junc(ip_back_coverage, ip_linear_junc_cov)
                with open(args.output_dir + '/' + fname + '_Corresponding_Linear_Junctions.txt', 'w') as file:
                    file.write(json.dumps(d))

                output_file = args.output_dir + '/' + fname + '_Linear_Juncs_Filtered_Sam_File.txt '
                ip_linear_samfile.to_csv(output_file, sep='\t')
                print("Filter sam file written to following file: %s" % output_file)

                print("Found %d possible reads that map to backspicing junctions." % len(set(ip_linear_samfile.rname)))
                print("Found %d possible backspicing junctions." % len(ip_linear_junc_cov.index))
                output_file = args.output_dir + '/' + fname + '_LinearSplicing_Junction_Coverages.txt'
                ip_linear_junc_cov.to_csv(output_file, sep='\t')
                print("Coverage for each linear-splicing junction found written to following file: %s" % output_file)

                # Output IP back-splicing junction differential expression results
                output_file = args.output_dir + \
                              '/' + fname + '_IP_Backsplicing_Junctions_Differential_Expression_Over_Linear_Junctions.txt'
                ip_back_junc_diff_exp.to_csv(output_file, sep='\t')
                print("Differential Expression analysis of back-splicing junctions written to following file: %s" %
                      output_file)

            else:
                parser.error('IP SAM for linear junctions not found. Enter correct path and name to SAM file')

            if args.input == '':
                print('INPUT SAM file not provided. Junction read counts will not be input normalized.')
            # If input sam file is given, parse input sam file and normalized ip with input
            elif os.path.isfile(args.input):
                print("Found INPUT Sam File. Input normalization of IP sample will be performed")
                input_samfile, input_backsplicing_junc_cov = parse_sam_file(args.input, args.cut_off, 'input')

                # Normalize IP back-splicing junction counts with Input counts
                ip_normalized_juncs = normalize_by_input(ip_back_coverage, input_backsplicing_junc_cov)

                # Output input normalized counts
                output_file = args.output_dir + '/' + fname + '_IP_Backsplicing_Junctions_Coverage_Input_Normalized.txt'
                ip_normalized_juncs.to_csv(output_file, sep='\t')
                print("INPUT Normalized Junctions written to following file: %s" % output_file)
            else:
                parser.error("INPUT SAM file not found. Enter correct path and name")

            output_file = args.output_dir + '/' + fname + '_Filtered_Sam_File.sam'

            # Generate index for output sam file
            generate_sam_index(set(ip_back_filtered[2]), output_file)

            # Generate fasta file for mapped junctions
            fasta_file_for_mapped_junctions(set(ip_back_filtered[2]), args.output_dir)

            ip_back_filtered.drop([u'tx', u'rname-tx', u'match_len', u'match_end_pos'], 1).to_csv(output_file, sep='\t', header=False, index=False, mode='a')
            print("Filter sam file written to following file: %s" % output_file)

            print("Found %d possible reads that map to backspicing junctions." % len(set(ip_back_filtered[0])))
            print("Found %d possible backspicing junctions." % len(ip_back_coverage.index))
            output_file = args.output_dir + '/' + fname + '_Backsplicing_Junction_Coverages.txt'
            ip_back_coverage.to_csv(output_file, sep='\t')
            print("Coverage for each backsplicing junction found written to following file: %s" % output_file)

    else:
        parser.error("IP SAM file not found. Enter correct path and name.")


if __name__ == '__main__':
    main()
