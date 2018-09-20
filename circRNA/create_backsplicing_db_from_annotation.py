# Import required libraries
import pandas
import numpy
import os
import pybedtools
from argparse import ArgumentParser


def parse_exons_bed_file(f):
    """
    This function reads and parses the initial exon coordinates bed file provided by the user. The function will read 
    the file and assign a sequential 'exon number' to each exon within the group. Group is determined from column 4 
    of the bed file. If column 4 contains gene ids, exon will be grouped according to the corresponding gene. Use 
    transcript ids to group by transcript isoforms 
    :param f: Path to exon BED file 
    :return: pandas dataframe with bed file columns 
    """
    # TODO:The whole function works on pandas, even for parsing a BED file. Can be probably implemented in a better way
    # TODO:using pybedtools.

    # Read bed file using pandas
    df = pandas.read_csv(f, sep='\t', names=['chr', 'start', 'stop', 'gene', 'x', 'strand'])
    # Assign exon number to each exon
    df['exon_num'] = df.groupby('gene').cumcount()
    # Generate a key containing the coordinates and tx information for each exon. This key will help get rid of
    # duplicates records from the bed file, if any are present
    df['key'] = df['chr'] + '__' + df['start'].astype(str) + '__' + df['stop'].astype(str) + '__' + \
                   df['gene'] + '__' + df['exon_num'].astype(str) + '__' + df['strand']
    df = df[['chr', 'start', 'stop', 'key', 'exon_num', 'strand']]
    df['stop2'] = df['start'] + 30
    df['start2'] = df['stop'] - 30
    return df[['chr','start','stop2','key','exon_num','strand']], df[['chr','start2','stop','key','exon_num','strand']]


def extract_genome_seqs_for_regions(b, genome_seqs):
    """
    This function uses the pybedtools wrapper for 'bedtools getfasta' to extract sequences for the regions within a
    given bed file. Here, the bed file is either the starting or ending 30 bp of an exon
    :param bed: Start / End coordinates BedTools object
    :param genome_seqs: Genome fasta file to extract the sequences from
    :return: Dataframe containing the sequences for regions in the passed bed file
    """
    bed = pybedtools.BedTool.from_dataframe(b)
    region_seqs = bed.sequence(fi=genome_seqs, name=True, tab=True)
    region_seqs = open(region_seqs.seqfn).read()
    l = [map(str, x.split("\t")) for x in region_seqs.split('\n')]
    region_seqs = pandas.DataFrame(data=l, columns=['key','seq']).dropna()
    region_seqs.index = map(lambda x: x.split('__')[3]+'__'+x.split('__')[4], region_seqs['key'])
    region_seqs['gene'] = map(lambda x: x.split('__')[3], region_seqs['key'])
    region_seqs['exon_num'] = map(lambda x: x.split('__')[4], region_seqs['key'])
    region_seqs['exon_num'] = region_seqs['exon_num'].astype(int)
    region_seqs['strand'] = map(lambda x: x.split('__')[5], region_seqs['key'])
    return region_seqs


def generate_backsplicing_junction_seqs(lefts, rights):
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
    d = dict()
    for name, grp in lefts.groupby('gene'):
        # Since we are not considering terminal exons, only consider those regions that have more than 3 exons. If a
        # transcript has 3 or less exons, skip that transcript
        # grp = grp.sort_values('exon_num')
        if 4 in list(grp['exon_num']):
            l1 = list(grp['seq'])
            k1 = list(grp['key'])

            l2 = list(rights[rights['gene'] == name].sort_values('exon_num')['seq'])
            k2 = list(rights[rights['gene'] == name].sort_values('exon_num')['key'])
            # Loop over an exon and all possible successive exon, starting with itself
            for i in range(0, len(l1) - 2):
                seq = l2[i] + l1[i]
                ind = k2[i] + '__and__' + k1[i]
                d[ind] = seq

                for j in range(i + 1, len(l1) - 1):
                    seq = l2[j] + l1[i]
                    ind = k2[j] + '__and__' + k1[i]
                    d[ind] = seq
        else:
            continue
    return d


def generate_linear_splicing_junction_seqs(lefts, rights):
    start_grouped = lefts.groupby('gene')
    d = dict()
    for name, grp in start_grouped:
        if 2 in list(lefts.exon_num):
            df = rights[rights['gene'] == name]
            if '+' in list(df['strand']):
                for i in range(1, max(df['exon_num'])):
                    d[df[df['exon_num'] == i].key[0] + grp[grp['exon_num'] == i+1].key[0]] = \
                        df[df['exon_num'] == i]['seq'][0] + grp[grp['exon_num'] == i+1]['seq'][0]
            if '-' in list(df['strand']):
                for i in range(1, max(df['exon_num'])):
                    d[df[df['exon_num'] == i+1].key[0] + grp[grp['exon_num'] == i].key[0]] = \
                        df[df['exon_num'] == i+1]['seq'][0] + grp[grp['exon_num'] == i]['seq'][0]

    df = pandas.DataFrame.from_dict(d, orient='index')
    print df.head()
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
            lefts, rights = parse_exons_bed_file(args.exons)
            lefts['start'] = lefts['start'] - 1

            # Create dataframe containing sequences for lefts and rights of each exon
            lefts = extract_genome_seqs_for_regions(lefts, args.genome)
            rights = extract_genome_seqs_for_regions(rights, args.genome)

            # Generate backsplicing junction sequences
            d = generate_backsplicing_junction_seqs(lefts, rights)
            f = open(args.output_dir + 'exon_exon/backsplicing_exons_seqs.fasta', 'w')
            for i in d.keys():
                f.write('>%s\n'%i)
                f.write('%s\n'%d[i])
            f.close()

            # Generate linear splicing junction sequences
            #d = generate_linear_splicing_junction_seqs(lefts, rights)
            #f = open(args.output_dir + 'exon_exon/linear_splicing_exons_seqs.fasta', 'w')
            #for i in d.keys():
            #    f.write('>%s\n'%i)
            #    f.write('%s\n'%d[i])
            #f.close()
        else:
            parser.error('Genome sequences Fasta file not found. Enter full path and name of the Fasta file.')
    else:
        parser.error('Exons BED file not found. Enter full path and name of the BED file.')


if __name__ == '__main__':
    main()
