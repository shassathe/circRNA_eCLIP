# Import required libraries
import pandas
import numpy
import os
import pybedtools
from argparse import ArgumentParser


def parse_exons_bed_file(f):
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
    exons = pandas.read_csv(f, sep='\t', names=['chr', 'start', 'stop', 'gene', 'x', 'strand'])
    # Sort bed file by group and start coordinate
    exons = exons.sort_values(['gene', 'start'])
    # Assign exon number to each exon
    exons['exon_num'] = exons.groupby('gene').cumcount()
    exons['key'] = exons['chr'] + '__' + exons['start'].astype(str) + '__' + exons['stop'].astype(str) + '__' + \
                   exons['gene'] + '__' + exons['exon_num'].astype(str) + '__' + exons['strand']
    return exons[['chr', 'start', 'stop', 'key', 'exon_num', 'strand']]


def create_exon_starts_bed_file(df):
    """
    Creates a bed file that contains coordinates for the first 35 bp of each exon
    :param df: Exon coordinates dataframe
    :return: Coordinates for the start of each exon
    """
    df['start'] = df['start'] - 1
    df['stop'] = df['start'] + 35
    return BedTool.from_dataframe(df)


def create_exon_stops_bed_file(df):
    """
    Creates a bed file that contains coordinates for the last 35 bp of each exon
    :param df: Exon coordinates dataframe
    :return: Coordinates for the end of each exon
    """
    df['start'] = df['stop'] - 36
    return BedTool.from_dataframe(df)


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
    region_seqs.index = map(lambda x: x.split('__')[3]+'__'+x.split('__')[4], region_seqs['key'])
    region_seqs['gene'] = map(lambda x: x.split('__')[3], region_seqs['key'])
    region_seqs['exon_num'] = map(lambda x: x.split('__')[4], region_seqs['key'])
    region_seqs['strand'] = map(lambda x: x.split('__')[5], region_seqs['key'])
    return region_seqs


def parse_seq_file(f):
    """
    This function parses relevant information from the 'key' column in the input dataframe.
    :param f: Full path of the input dataframe as given in args.exon_start and args.exon_ends
    :return df: Parsed dataframe
    """
    # Read input dataframe.
    df = pandas.read_csv(f, sep='\t', names=['key', 'seq'])
    # Split 'key' column to extract relevant information.
    df.index = map(lambda x: x.split('__')[3] + '__' + x.split('__')[4], df['key'])
    df['tx'] = map(lambda x: x.split('__')[3], df['key'])
    df['exon_num'] = map(lambda x: x.split('__')[4], df['key'])
    df['strand'] = map(lambda x: x.split('__')[5], df['key'])
    return df


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

    df = pandas.DataFrame.from_dict(d, orient='index')
    df.columns = ['Seq']
    return df


def generate_linear_splicing_junction_seqs(starts, stops):
    start_grouped = starts.groupby('gene')
    d = dict()
    for name, grp in start_grouped:
        if 2 in list(starts.exon_num):
            df = stops[stops['gene'] == name]
            if '+' in list(df['strand']):
                for i in range(1, max(df['exon_num'])):
                    d[df[df['exon_num'] == i].key[0] +'__'+ grp[grp['exon_num'] == i+1].key[0]] = \
                        df[df['exon_num'] == i]['seq'][0] + grp[grp['exon_num'] == i+1]['seq'][0]
            if '-' in list(df['strand']):
                for i in range(1, max(df['exon_num'])):
                    d[df[df['exon_num'] == i+1].key[0] +'__'+ grp[grp['exon_num'] == i].key[0]] = \
                        df[df['exon_num'] == i+1]['seq'][0] + grp[grp['exon_num'] == i]['seq'][0]

    df = pandas.DataFrame.from_dict(d, orient='index')
    df.columns = ['Seq']
    return df


def main():
    parser = ArgumentParser(description="Generate backsplicing junction index.")
    parser.add_argument('--exons', help='Full path and name of bed file containing exon coordinates (Required)')
    parser.add_argument('--genome', help='Full path and name of fasta file containing all genome chromosome sequences '
                                         '(Required)')
    parser.add_argument('--output_dir', help='Full path to output directory. Default = current dir', default='./')
    args = parser.parse_args()

    # Check if exons bed file exists
    if os.path.isfile(args.exons):
        # Check if genome sequences fasta file exists
        if os.path.isfile(args.genome):
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

            # Generate linear splicing junction sequences
            d = generate_linear_splicing_junction_seqs(starts, stops)
            f = open(args.output_dir + 'exon_exon/linear_splicing_exons_seqs.fasta', 'w')
            for i in d.keys():
                f.write('>%s\n'%i)
                f.write('%s\n'%d[i])
            f.close()
        else:
            parser.error('Genome sequences Fasta file not found. Enter full path and name of the Fasta file.')
    else:
        parser.error('Exons BED file not found. Enter full path and name of the BED file.')


if __name__ == '__main__':
    main()
