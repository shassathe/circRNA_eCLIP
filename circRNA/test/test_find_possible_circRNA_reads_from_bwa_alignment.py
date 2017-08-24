from circRNA import find_possible_circRNA_reads_from_bwa_alignment as fpc
import pandas as pd

def duped_sam_1():
    samfile = 'test/data/226_igf2bp2.sai.sorted.read.sam'
    return fpc.parse_sam_file(samfile)

def test_check_for_duplicate_reads_1():
    print("assume function keeps only unique reads")
    df, _ = duped_sam_1()
    deduped = fpc.check_for_duplicate_reads(df)
    assert deduped.shape[0] == 1

def test_check_for_duplicate_reads_2():
    print("assume function removes both duplicated reads")
    df, _ = duped_sam_1()
    deduped = fpc.check_for_duplicate_reads(df)
    assert deduped.shape[0] == 0