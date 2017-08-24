from circRNA import find_possible_circRNA_reads_from_bwa_alignment as fpc
import pandas as pd

def duped_sam_1():
    samfile = 'test/data/226_igf2bp2.sai.sorted.read.sam'
    df = pd.read_table(
        samfile, sep='\t', comment='@', names=[
            'rname', 'flag', 'junc_name', 'pos', 'mapq', 'cigar',
            'rnext', 'pnext', 'tlen', 'seq', 'qual', 'XT', 'NM',
            'X0', 'X1', 'XM', 'XO', 'XG', 'MDZ'
        ]
    )
    return df

def test_check_for_duplicate_reads_1():
    print("assume function keeps only unique reads")
    df = duped_sam_1()
    deduped = fpc.check_for_duplicate_reads(df)
    assert deduped.shape[0] == 1

def test_check_for_duplicate_reads_2():
    print("assume function removes both duplicated reads")
    df = duped_sam_1()
    deduped = fpc.check_for_duplicate_reads(df)
    assert deduped.shape[0] == 0