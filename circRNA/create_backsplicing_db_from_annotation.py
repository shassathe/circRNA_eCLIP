# Import required libraries
import pandas
import numpy
import os
import pybedtools
from argparse import ArgumentParser
from tqdm import trange

def parse_exons_bed_file(f, min_exon_length=35):
    """
    This function reads and parses the initial exon coordinates bed file provided by the user. The function will read
    the file and assign a sequential 'exon number' to each exon within the group. Group is determined from column 4 of
    the bed file.
    :param f: Path to exon BED file
    :return: pandas dataframe with bed file columns
    """
    # TODO:The whole function works on pandas, even for parsing a BED file. Can be probably implemented in a better way
    # TODO:using pybedtools.

    # Read bed file using pandas
    exons = sort_bed_file(f)
    # Assign exon number to each exon
    exons['exon_num'] = exons.groupby('name').cumcount()
    exons['key'] = exons['chrom'] + '__' + \
                   exons['start'].astype(str) + '__' + \
                   exons['end'].astype(str) + '__' + \
                   exons['name'] + '__' + \
                   exons['exon_num'].astype(str) + '__' + \
                   exons['strand']
    #
    exons['len'] = exons['end'] - exons['start']
    exons = filter_min_exon_length(exons, min_exon_length)
    return exons[['chrom', 'start', 'end', 'key', 'exon_num', 'strand']]

def sort_bed_file(f):
    """
    sorts a bed file and returns as a dataframe.
    :param f: string
    :return: pandas.DataFrame()
    """
    bed_tool = pybedtools.BedTool(f)
    bed_tool = bed_tool.sort()
    return bed_tool.to_dataframe()

def filter_min_exon_length(df, min_length):
    """
    Removes any exon whose length is less than min_length
    :param df: pandas.DataFrame()
    :param min_length: int
    :return:
    """
    return df[df['len'] >= min_length]

def create_exon_starts_bed_file(df):
    """
    Creates a bed file that contains coordinates for the first 35 bp of each exon
    :param df: Exon coordinates dataframe
    :return: Coordinates for the start of each exon
    """
    df['start'] = df['start'] - 1
    df['end'] = df['start'] + 35

    return pybedtools.BedTool.from_dataframe(df)


def create_exon_stops_bed_file(df):
    """
    Creates a bed file that contains coordinates for the last 35 bp of each exon
    :param df: Exon coordinates dataframe
    :return: Coordinates for the end of each exon
    """
    df['start'] = df['end'] - 36
    return pybedtools.BedTool.from_dataframe(df)


def extract_genome_seqs_for_regions(bed, genome_seqs):
    """
    This function uses the pybedtools wrapper for 'bedtools getfasta' to extract sequences for the regions within a
    given bed file. Here, the bed file is either the starting or ending 35 bp of an exon
    :param bed: Start / End coordinates bed file
    :param genome_seqs: Genome fasta file to extract the sequences from
    :return: Dataframe containing the sequences for regions in the passed bed file
    """
    region_seqs = bed.sequence(fi=genome_seqs, name=True, tab=True)
    region_seqs = open(region_seqs.seqfn).read()
    l = [map(str, x.split("\t")) for x in region_seqs.split('\n')]
    region_seqs = pandas.DataFrame(data=l, columns=['key','seq'])
    region_seqs = region_seqs[region_seqs['key']!=""]
    region_seqs = add_keys(region_seqs)
    return region_seqs

def add_keys(region_seqs):
    """
    Parses the 'key' field in region_seqs and appends columns to row
    :param region_seqs:
    :return:
    """
    region_seqs.index = map(lambda x: x.split('__')[3]+'__'+x.split('__')[4], region_seqs['key'])
    region_seqs['gene'] = map(lambda x: x.split('__')[3], region_seqs['key'])
    region_seqs['exon_num'] = map(lambda x: x.split('__')[4], region_seqs['key'])
    region_seqs['strand'] = map(lambda x: x.split('__')[5], region_seqs['key'])
    return region_seqs


def generate_backsplicing_junction_seqs(starts, stops):
    """
    This function concatenates the sequences from two backsplicing exons to generate a single backsplicing junction
    sequence. This will be done for all possible combinations of backsplicing exons within a region, eg: gene or
    transcript.
    :param starts: Dataframe for exon starts.
    :param ends: Dataframe for exon ends.
    :return dd: pandas Dataframe with the backsplicing junction coordinate and name as key and the corresponding
    sequence as column values
    """
    # Group dataframe by region.
    progress = trange(len(set(starts['gene'])))
    start_grouped = starts.groupby('gene')
    d = dict()

    for name, grp in start_grouped:
        # Since we are not considering terminal exons, only consider those regions that have more than 3 exons. If a
        # transcript has 3 or less exons, skip that transcript
        if list(grp['exon_num'])[-1] > 3:
            df = stops[stops['gene'] == name]
            # Loop over an exon and every exon succeeding in
            if '+' in list(grp['strand']):
                for i in range(1, len(grp)-2):
                    for j in range(i+1, len(grp)-1):
                        seq = list(df['seq'])[j]+list(grp['seq'])[i]
                        ind = list(df.key)[j]+'__and__'+list(grp.key)[i]
                        d[ind] = seq

            if '-' in list(grp['strand']):
                for i in range(1, len(grp)-2):
                    for j in range(i+1, len(grp)-1):
                        seq = list(df['seq'])[i]+list(grp['seq'])[j]
                        ind = list(df.key)[i]+'__and__'+list(grp.key)[j]
                        d[ind] = seq
        progress.update(1)
    df = pandas.DataFrame.from_dict(d, orient='index')
    df.columns = ['Seq']
    return df


def main():
    parser = ArgumentParser(
        description="Generate backsplicing junction index."
    )
    parser.add_argument(
        '--exons',
        help='Full path and name of bed file containing '
             'exon coordinates (Required)',
        required=True,
    )
    parser.add_argument(
        '--genome',
        help='Full path and name of fasta file containing '
             'all genome chromosome sequences '
             '(Required)',
        required=True,
    )
    parser.add_argument(
        '--output_dir',
        help='Full path to output directory. '
             'Default = current dir',
        default='./',
        required=False
    )
    args = parser.parse_args()

    # Parse exons bed file
    exons = parse_exons_bed_file(args.exons)

    # Create exon starts and stops bed files
    starts = create_exon_starts_bed_file(exons)
    stops = create_exon_stops_bed_file(exons)

    # Create dataframe containing sequences for starts and stops of each exon
    starts = extract_genome_seqs_for_regions(starts, args.genome)
    stops = extract_genome_seqs_for_regions(stops, args.genome)

    # Generate backsplicing junction sequences
    d = generate_backsplicing_junction_seqs(starts, stops)

    f = open(args.output_dir + 'exon_exon/backsplicing_exons_seqs.fasta', 'w')
    for i in d.keys():
        f.write('>%s\n'%i)
        f.write('%s\n'%d[i])
    f.close()


if __name__ == '__main__':
    main()
