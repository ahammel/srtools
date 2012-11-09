from srtools import seq
import pytest

TEST_FOLDER = "/home/alex/repos/py-srtools/srtools/test"

class SequenceTestSetup(object):
    fasta_reads = seq.read_fasta(TEST_FOLDER + "/test_data/fasta.fa")

class TestSequenceFunctions(SequenceTestSetup):
    def test_reverse_complement(self):
        assert seq.reverse_complement("") == ""
        assert seq.reverse_complement("N") == "N"
        assert seq.reverse_complement("ACGT") == "ACGT"
        assert seq.reverse_complement("GCCAT") == "ATGGC"

    def test_read_fasta(self):
        assert len(self.fasta_reads) == 2
        assert self.fasta_reads["read1"] == "AAAAAAAAAAAAAA"
        assert self.fasta_reads["read2"] == "ATGAATGAATAAAAAAAAATAATTATTTCAT" 

    def test_gc_content(self):
        pytest.raises(seq.NullSequenceError, seq.gc_content, "")
        pytest.raises(seq.NullSequenceError, seq.gc_content, "NNNNN")
        assert seq.gc_content("AAATTT") == 0.0
        assert seq.gc_content("GGGGCC") == 1.0
        assert seq.gc_content("ACGTACGT") == 0.5
        assert seq.gc_content("ACCTACGT") == 0.5
    
    def test_random_sequence(self):
        sequence = seq.random_sequence(10)
        assert len(sequence) == 10
        for letter in sequence:
            assert letter in "ACGT"

    def test_randomize_sequence(self):
        sequence = "AACGNNNNAACG"
        rseq = seq.randomize_sequence(sequence)
        assert rseq[4:8] == "NNNN"
        assert rseq != seq

    def test_reading_frames(self):
        assert seq.reading_frames(self.fasta_reads["read1"]) == \
            [["AAA", "AAA", "AAA", "AAA", "AA"],
             ["TTT", "TTT", "TTT", "TTT", "TT"],
             ["A", "AAA", "AAA", "AAA", "AAA", "A"],
             ["T", "TTT", "TTT", "TTT", "TTT", "T"],
             ["AA", "AAA", "AAA", "AAA", "AAA"],
             ["TT", "TTT", "TTT", "TTT", "TTT"]]

        assert seq.reading_frames(self.fasta_reads["read2"]) == \
            [["ATG", "AAT", "GAA", "TAA", "AAA", "AAA", "ATA", "ATT", "ATT", "TCA", "T"],
             ["ATG", "AAA", "TAA", "TTA", "TTT", "TTT", "TTT", "ATT", "CAT", "TCA", "T"],
             ["A", "TGA", "ATG", "AAT", "AAA", "AAA", "AAA", "TAA", "TTA", "TTT", "CAT"],
             ["A", "TGA", "AAT", "AAT", "TAT", "TTT", "TTT", "TTA", "TTC", "ATT", "CAT"],
             ["AT", "GAA", "TGA", "ATA", "AAA", "AAA", "AAT", "AAT", "TAT", "TTC", "AT"],
             ["AT", "GAA", "ATA", "ATT", "ATT", "TTT", "TTT", "TAT", "TCA", "TTC", "AT"]]

    def test_open_reading_frames(self):
        assert seq.open_reading_frames(self.fasta_reads["read1"]) == []
        assert seq.open_reading_frames(self.fasta_reads["read2"]) == \
            ["ATGAATGAA", "ATGAAA", "ATGAATAAAAAAAAA"]


