# This script accepts a bed file for aligned reads. The bed file can be made by using "bedtools bamtobed". Convert the
# bam file from aligner to a bed file. Other annotation files used in the script include bed files for exon and intron
# junction starts and stops and bed files for entire exon and intron regions. The script will search for reads that
# span one of the possible combinations:
# 1) Exon - Exon backsplicing ( junctions only )
# 2) Exon - Intron backsplicing ( junctions only )
# 3) Intron - Intron backsplicing ( junctions only )
# 4) Exon - Exon backsplicing ( junction of one exon and internal of another exon -- not sure how a situation like this
# would arise )
# 5) Exon - Intron backsplicing (junction of one exon/intron and internal of another intron/exon -- could be an intron
# lariat ?)
# 6) Intron - Intron backsplicing ( junction of one intron and internal to another intron -- Again not sure how this
# would arise )

import pandas
import numpy
import pybedtools
import os
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter


def intersect_bed_files(f1, f2, col):
    """
    This functions accepts two bed files: alignment bed file and an annotation bed file. The function will perform an
    intersect, convert output to dataframe format, search for records that show overlap of at least 1bp, drop duplicate
    reads and then return the remaining dataframe.
    :param col:
    :param f1: Alignment bed file
    :param f2: Annotation bed file
    :return: Intersect dataframe
    """
    f_intersect = f1.intersect(f2, wao=True, s=True)
    f_intersect_df = f_intersect.to_dataframe(
        names=['chr', 'start', 'stop', 'rname', 'score', 'strand', col + '_chr', col + '_start',
               col + '_stop', col + '_txID', col + '_num', col + '_strand', col + '_overlap'])
    f_intersect_df = f_intersect_df[f_intersect_df[col + '_overlap'] > 0]
    f_intersect_df['rname-key'] = f_intersect_df['rname'] + '-' + f_intersect_df['chr'] + '-' + f_intersect_df['start'] \
        .astype(str) + '-' + f_intersect_df['stop'].astype(str)
    return f_intersect_df.drop_duplicates('rname-key')


def keep_multi_jxn_spanning_reads(df, col1, col2):
    """
    This function searches for reads that map to multiple exons or introns within a transcript. The function does not
    search for reads mapping to multiple exons and introns at the same time. It searches for exons or introns
    individually.
    Once multi-mapping reads are found, their positions are annotated as 'start', 'stop' or 'internal' to an exon or
    intron. (annotate_mapping_reads)
    The function then searches for reads mapping to consecutive exons or introns. (find_circular_rna_junction_reads)
    :param df: Input dataframe returned by the 'intersect_bed_files' function.
    :param col1: column to form read name and transcript key from ( can be exon_txID or intron_txID )
    :param col2: column name that secifies the exon or intron number in the transcript ( can be exon_num or intron_num )
    :return: Returns 1) dataframe containing only reads that map to multiple exons or introns within a given transcript,
     2) list of potential junctions with each value in the list containing the read name, txID, and exon/intron numbers,
     3) list of potential junctions with each value in the list containing the read name and txID
    """
    # Form read name and mapped transcript ID key
    df['rname-tx'] = df['rname'] + '-' + df[col1]
    # Find reads that map to multiple locations ( not necessarily multiple exons or introns )within a transcript
    df = df[df.duplicated('rname-tx', keep=False)]
    keep_reads = list()
    for name, grp in df.groupby('rname-tx'):
        s_temp = grp[col2].value_counts()
        # For a given read mapping to a transcript, check if the read maps to multiple exons or introns. If the reads
        # maps to a single exon or intron, read-trasncript will be discarded from further analysis.
        if len(s_temp) > 1:
            keep_reads.append(name)
    if len(keep_reads) > 0:
        df = df[df['rname-tx'].isin(keep_reads)]
        df = annotate_mapping_positions(df, col1.split('_')[0])
        potential_circrna_juncs, potential_circrna_rname_tx = find_circular_rna_junction_reads(df.sort_values(col2),
                                                                                               col1.split('_')[0])
        return df, potential_circrna_juncs, potential_circrna_rname_tx
    else:
        return df, [], []


def annotate_mapping_positions(df, col):
    """
    This function will annotate the mapping position of read. Mapping positions are annotated as "exon/intron_start",
    "exon/intron_stop" or "exon/intron_internal". Mapping positions are annotated as 'start' or stop' if they fall
    within 25bp of the start or stop of the exon/intron respectively.
    :param df: Dataframe returned from the 'keep_multimapping_reads'
    :param col: Keyword specifying whether searching for exon or intron
    :return: Returns dataframe with the annotation for each mapping position
    """
    # Define custom annotation ( exon_start, exon_stop OR intron_start, intron_stop)
    strt = col + '_start'
    stp = col + '_stop'

    # Get distance between mapping position and start or stop of an exon/intron
    df['start_dist'] = numpy.abs(df[strt] - df['start'])
    df['stop_dist'] = numpy.abs(df[stp] - df['stop'])

    # All mapping positions are defined initially as exon/intron_internal. Annotations will be changed as we encounter
    # reads mapping close to start or stop positions.
    df['loc'] = col + '_internal'

    # Find reads within 25 bps of start of exon/intron.
    inds = df[df['start_dist'] <= 25].index
    df.loc[inds, 'loc'] = col + '_start'

    # Find reads within 25 bps of stop of exon/intron.
    inds = df[df['stop_dist'] <= 25].index
    df.loc[inds, 'loc'] = col + '_stop'

    return df.drop(['rname-key', 'start_dist', 'stop_dist'], 1)


def find_circular_rna_junction_reads(df, col):
    """
    Find reads that potentially map to circular RNA junctions ( exon-exon or intron-intron junctions only for now ).
    The functions currently handles only those circular RNAs that form on consecutive exons or introns.
    :param df: Dataframe returned from the 'annotation_mapping_positions' function. The dataframe is sorted by the
    exon/intron numbers within a transcript
    :param col: Keyword specifying whether to search for exons or introns
    :return: Returns 1) list of potential reads that map to circRNA junctions. Each value in the list contains the read
    name, tx ID and the exon/intron numbers on which the circular RNA is being formed. 2) List containing the read
    names and tx ID only.
    """
    # Create groups to search by unique combinations of a read and transcript
    print(col)
    rname_tx_grouped = df.groupby('rname-tx')
    # Create empty lists. The lists will be filled as circRNA events are found
    potential_juncs = list()
    potential_rname_tx = list()
    for name, grp in rname_tx_grouped:
        # Extract exon/intron numbers and annotations for each combination of read and transcript
        nums = list(grp[col + '_num'])
        locs = list(grp['loc'])
        # Treat '-' and '+' strand differently
        if '-' in list(grp['strand']):
            # Loop over each read and the next read in the group
            for i in range(len(nums) - 1):
                for j in range(i + 1, len(nums)):
                    # Check if the two reads map to consecutive exons. The second read had to map to the next immediate
                    # exon/intron of the previous read.
                    if (nums[j] - nums[i]) == 1:
                        # Check if the read maps to the stop of the first exon and the start of the next exon. This
                        # trend will be reversed on the positive strand
                        if ('start' in locs[j]) & ('stop' in locs[i]):
                            potential_rname_tx.append(name)
                            potential_juncs.append(name.split('-')[2] + '-' + str(nums[i]) + '-' + str(nums[j]))
        else:
            # Loop over each read and the next read in the group
            for i in range(len(nums) - 1):
                for j in range(i + 1, len(nums)):
                    # Check if the two reads map to consecutive exons. The second read had to map to the next immediate
                    # exon/intron of the previous read.
                    if (nums[j] - nums[i]) == 1:
                        # Check if the read maps to the start of the first exon and the stop of the next exon. This
                        # trend will be reversed on the positive strand
                        if ('stop' in locs[j]) & ('start' in locs[i]):
                            potential_rname_tx.append(name)
                            potential_juncs.append(name.split('-')[2] + '-' + str(nums[i]) + '-' + str(nums[j]))
    return potential_juncs, potential_rname_tx


def find_confident_circular_rna_junctions(df, reads, col):
    """
    This function checks whether the reads mapping to circRNA junctions are unique to those junctions. If the read maps
    to any other location ( could be another exon/intron) within the same transcript, the read and circRNA junction
    will be discarded. We can be more confident of the remaining reads since they only map to the circRNA junction.
    :param df: Dataframe returned from 'annotate_mapping_positions' function
    :param reads: List of reads that map exclusively to exons or exclusively to introns within a transcript
    :param col: 'exon' or 'intron'
    :return: Returns a list of read name and transcript combinations that show evidence of circular RNAs
    """
    df = df[df['rname-tx'].isin(reads)]
    grouped = df.groupby('rname-tx')
    confident_reads = list()
    for name, grp in grouped:
        if len(grp[col].unique()) == 2:
            confident_reads.append(name)
    return confident_reads


def read_counts_for_circular_rna_junctions (confident_juncs, potential_reads):
    juncs_counts = dict()
    for i in confident_juncs:
        for j in potential_reads:
            if i in j:
                candidate = j.split('-')[2] +'-'+ j.split('-')[3] +'-'+ j.split('-')[4]
                if candidate in juncs_counts.keys():
                    juncs_counts[candidate] = juncs_counts[candidate] + 1
                else:
                    juncs_counts[candidate] = 1
    return juncs_counts


def main():
    parser = ArgumentParser(description="Find reads that potentially map across circularRNAs.")
    parser.add_argument("--alignment_bed", help="Full path and name of bed file containing alignments")
    parser.add_argument("--output_dir", help="Full path to output directory")
    parser.add_argument(
        "--exons",
        default='/home/shsathe/annotation_files/hg19_exons_numbers.bed',
        help="Full path to bed file containing exons"
    )
    parser.add_argument(
        "--introns",
        default='/home/shsathe/annotation_files/hg19_introns_numbers.bed',
        help="Full path to bed file containing introns"
    )

    args = parser.parse_args()

    if os.path.isfile(args.alignment_bed):
        print("Found Bed File...\n")
        print("Reading bed file...\n")
        bedfile = pybedtools.BedTool(args.alignment_bed)

        print("Intersecting mapped reads with exon and intron annotations...")
        exon_bed = pybedtools.BedTool(args.exons)
        reads_and_exons_df = intersect_bed_files(bedfile, exon_bed, 'exon')

        intron_bed = pybedtools.BedTool(args.introns)
        reads_and_introns_df = intersect_bed_files(bedfile, intron_bed, 'intron')

        ### If aligned reads do not overlap with either exons or introns, pass
        if (len(reads_and_exons_df) == 0) & (len(reads_and_introns_df) == 0):
            parser.error("No reads found that map to introns or exons ! Check input bed file...")

        ### If aligned reads overlap with ONLY exons
        elif (len(reads_and_exons_df) != 0) & (len(reads_and_introns_df) == 0):
            print("Found exon-mapping reads, but no intron-mapping reads found. Proceeding to analyze exon-mapping "
                  "reads...")
            reads_multi_mapping_to_exons, potential_circrna_exon_exon_juncs, potential_circrna_exon_exon_rname_tx = \
                keep_multi_jxn_spanning_reads(reads_and_exons_df, 'exon_txID', 'exon_num')
            if len(reads_multi_mapping_to_exons) == 0:
                parser.error("No reads found that map multiple exons within the same transcript. Dataset likely does "
                             "not contain reads that map to circular RNAs.")
            else:
                print("Found reads that map to multiple exons within the same transcript..")
                confident_exon_exon_juncs = find_confident_circular_rna_junctions(
                    reads_multi_mapping_to_exons, set(potential_circrna_exon_exon_rname_tx), 'exon_num')
                reads_multi_mapping_to_exons[reads_multi_mapping_to_exons['rname-tx'].isin(
                    confident_exon_exon_juncs)].to_csv(args.output_dir + '/CircRNA_Exon_Exon_Junction_Events.txt',
                                                       sep='\t')
                exon_exon_juncs_counts = read_counts_for_circular_rna_junctions(confident_exon_exon_juncs,
                                                                                potential_circrna_exon_exon_juncs)
                if len(exon_exon_juncs_counts) != 0:
                    exon_exon_juncs_counts_df = pandas.DataFrame.from_dict(exon_exon_juncs_counts, orient='index')
                    exon_exon_juncs_counts_df.columns = ['Read Count']
                    exon_exon_juncs_counts_df.to_csv(args.output_dir + '/CircRNA_Exon_Exon_Junction_Read_Counts.txt',
                                                     sep='\t')
        ### If aligned reads overlap with ONLY introns
        elif (len(reads_and_exons_df) == 0) & (len(reads_and_introns_df) != 0):
            print("Found intron-mapping reads, but no exon-mapping reads found. Proceeding to analyze intron-mapping "
                  "reads...")
            reads_multi_mapping_to_introns, potential_circrna_intron_intron_juncs, \
            potential_circrna_intron_intron_rname_tx = keep_multi_jxn_spanning_reads(reads_and_introns_df, 'intron_txID',
                                                                               'intron_num')
            if (reads_multi_mapping_to_introns) == 0:
                parser.error("No reads found that map multiple introns within the same transcript. Dataset likely does "
                             "not contain reads that map to circular RNAs.")
            else:
                print("Found reads that map to multiple introns within the same transcript..")
                confident_intron_intron_juncs = find_confident_circular_rna_junctions(
                    reads_multi_mapping_to_introns, set(potential_circrna_intron_intron_rname_tx), 'intron_num')
                reads_multi_mapping_to_introns[reads_multi_mapping_to_introns['rname-tx'].isin(
                    confident_intron_intron_juncs)].to_csv(args.output_dir +
                                                           '/CircRNA_Intron_Intron_Junction_Events.txt', sep='\t')
                intron_intron_juncs_counts = read_counts_for_circular_rna_junctions(
                    confident_intron_intron_juncs, potential_circrna_intron_intron_juncs)
                if len(intron_intron_juncs_counts) != 0:
                    intron_intron_juncs_counts_df = pandas.DataFrame.from_dict(intron_intron_juncs_counts, orient='index')
                    intron_intron_juncs_counts_df.columns = ['Read Count']
                    intron_intron_juncs_counts_df.to_csv(args.output_dir +
                                                         '/CircRNA_Intron_Intron_Junction_Read_Counts.txt', sep='\t')

        else:
            print("Found exon-mapping reads and intron-mapping reads found. Proceeding with analysis...")
            reads_multi_mapping_to_exons, potential_circrna_exon_exon_juncs, \
            potential_circrna_exon_exon_rname_tx = keep_multi_jxn_spanning_reads(
                reads_and_exons_df, 'exon_txID', 'exon_num'
            )
            reads_multi_mapping_to_introns, potential_circrna_intron_intron_juncs, \
            potential_circrna_intron_intron_rname_tx = keep_multi_jxn_spanning_reads(
                reads_and_introns_df, 'intron_txID', 'intron_num'
            )

            if (len(reads_multi_mapping_to_exons) == 0) & (len(reads_multi_mapping_to_introns) == 0):
                parser.error("No reads found that map to multiple exons or introns within the same transcript. Dataset "
                             "likely does not contain circular RNAs.")

            elif (len(reads_multi_mapping_to_exons) != 0) & (len(reads_multi_mapping_to_introns) == 0):
                print("Found reads mapping to multiple exons in the same transcript. No such intronic reads found. "
                      "Proceeding with analysis...")
                confident_exon_exon_juncs = find_confident_circular_rna_junctions(
                    reads_multi_mapping_to_exons, set(potential_circrna_exon_exon_rname_tx), 'exon_num')
                reads_multi_mapping_to_exons[reads_multi_mapping_to_exons['rname-tx'].isin(
                    confident_exon_exon_juncs)].to_csv(args.output_dir + '/CircRNA_Exon_Exon_Junction_Events.txt',
                                                       sep='\t')
                exon_exon_juncs_counts = read_counts_for_circular_rna_junctions(confident_exon_exon_juncs,
                                                                                potential_circrna_exon_exon_juncs)
                if len(exon_exon_juncs_counts) != 0:
                    exon_exon_juncs_counts_df = pandas.DataFrame.from_dict(exon_exon_juncs_counts, orient='index')
                    exon_exon_juncs_counts_df.columns = ['Read Count']
                    exon_exon_juncs_counts_df.to_csv(args.output_dir + '/CircRNA_Exon_Exon_Junction_Read_Counts.txt',
                                                     sep='\t')

            elif (len(reads_multi_mapping_to_exons) == 0) & (len(reads_multi_mapping_to_introns) != 0):
                print("Found reads mapping to multiple introns in the same transcript. No such exonic reads found. "
                      "Proceeding with analysis...")
                confident_intron_intron_juncs = find_confident_circular_rna_junctions(
                    reads_multi_mapping_to_introns, set(potential_circrna_intron_intron_rname_tx), 'intron_num')
                reads_multi_mapping_to_introns[reads_multi_mapping_to_introns['rname-tx'].isin(
                    confident_intron_intron_juncs)].to_csv(args.output_dir +
                                                           '/CircRNA_Intron_Intron_Junction_Events.txt', sep='\t')
                intron_intron_juncs_counts = read_counts_for_circular_rna_junctions(
                    confident_intron_intron_juncs, potential_circrna_intron_intron_juncs)
                if len(intron_intron_juncs_counts) != 0:
                    intron_intron_juncs_counts_df = pandas.DataFrame.from_dict(intron_intron_juncs_counts,
                                                                               orient='index')
                    intron_intron_juncs_counts_df.columns = ['Read Count']
                    intron_intron_juncs_counts_df.to_csv(args.output_dir +
                                                         '/CircRNA_Intron_Intron_Junction_Read_Counts.txt', sep='\t')

            else:
                print("Found reads that map to multiple exons and introns within the same transcript..")
                confident_exon_exon_juncs = find_confident_circular_rna_junctions(
                    reads_multi_mapping_to_exons, set(potential_circrna_exon_exon_rname_tx), 'exon_num')
                reads_multi_mapping_to_exons[reads_multi_mapping_to_exons['rname-tx'].isin(
                    confident_exon_exon_juncs)].to_csv(args.output_dir + '/CircRNA_Exon_Exon_Junction_Events.txt',
                                                       sep='\t')
                exon_exon_juncs_counts = read_counts_for_circular_rna_junctions(confident_exon_exon_juncs,
                                                                                potential_circrna_exon_exon_juncs)
                if len(exon_exon_juncs_counts) != 0:
                    exon_exon_juncs_counts_df = pandas.DataFrame.from_dict(exon_exon_juncs_counts, orient='index')
                    exon_exon_juncs_counts_df.columns = ['Read Count']
                    exon_exon_juncs_counts_df.to_csv(args.output_dir + '/CircRNA_Exon_Exon_Junction_Read_Counts.txt',
                                                     sep='\t')

                confident_intron_intron_juncs = find_confident_circular_rna_junctions(
                    reads_multi_mapping_to_introns, set(potential_circrna_intron_intron_rname_tx), 'intron_num')
                reads_multi_mapping_to_introns[reads_multi_mapping_to_introns['rname-tx'].isin(
                    confident_intron_intron_juncs)].to_csv(args.output_dir +
                                                           '/CircRNA_Intron_Intron_Junction_Events.txt', sep='\t')
                intron_intron_juncs_counts = read_counts_for_circular_rna_junctions(
                    confident_intron_intron_juncs, potential_circrna_intron_intron_juncs)
                if len(intron_intron_juncs_counts) != 0:
                    intron_intron_juncs_counts_df = pandas.DataFrame.from_dict(intron_intron_juncs_counts,
                                                                               orient='index')
                    intron_intron_juncs_counts_df.columns = ['Read Count']
                    intron_intron_juncs_counts_df.to_csv(args.output_dir +
                                                         '/CircRNA_Intron_Intron_Junction_Read_Counts.txt', sep='\t')


if __name__ == '__main__':
    main()
