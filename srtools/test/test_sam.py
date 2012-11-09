try:
    import pytest
except ImportError:
    import nose as pytest
from srtools import sam, stats, seq, gff
import os
import glob
from srtools.test import TEST_FOLDER


class ReadTestSetup(object):
    """Does the setup for the tests of the Read classes and associated functions."""
    sam_string = ("SRR360147.1\t77\t*\t0\t0\t*\t*\t0\t0\tAGATGACGAAGAAGCTTGATC"
                  "TCACGAANNNNNNNNTTNNCATCCNNNTNNTNNNNNNNNNNNNNNNNNNNNNNN\tHHH"
                  "HHHHHHHHHHHHHHHHHHHHHH==@##################################"
                  "#############\tYT:Z:UP\tYF:Z:NS")

    single_read = sam.Read(sam_string)

    head_string = ("@HD\tVN:1.0\tSO:unsorted\n"
                   "@SQ\tSN:1\tLN:30427671\n"
                   "@SQ\tSN:2\tLN:19698289\n"
                   "@SQ\tSN:3\tLN:23459830\n"
                   "@SQ\tSN:4\tLN:18585056\n"
                   "@SQ\tSN:5\tLN:26975502\n"
                   "@SQ\tSN:Mt\tLN:366924\n"
                   "@SQ\tSN:Pt\tLN:154478\n"
                   "@PG\tID:bowtie2\tPN:bowtie2\tVN:2.0.0-beta5\n")

    def test_make_single_read(self):
        with open("single_read_test.sam", "w") as f:
            print(self.head_string, end="", file=f)
            print(self.sam_string, file=f)
            assert True

    single_read_alignment = sam.SamAlignment("single_read_test.sam") 

    indel_algn = sam.SamAlignment(TEST_FOLDER + "/test_data/test_indel.sam")

    reverse_complement_align = sam.SamAlignment(TEST_FOLDER + "/test_data/test_rc_consensus.sam")
    
    reverse_complement_align_copy = sam.SamAlignment(TEST_FOLDER + "/test_data/test_rc_consensus.sam")

    expressed_locus_alignment = sam.SamAlignment(TEST_FOLDER + "/test_data/test_expressed_locus.sam")

    mate_pair_align = sam.SamAlignment(TEST_FOLDER + "/test_data/test_mate_pair.sam")

    mate_pair_consensus = ("AATGCGAGTGGAAAATGGCTGAAGACTCGATCAAGATACCTCCATCCAC"
                           "CAACACGGTGAAGCAGAGCTGGATTG" + "N" * 36 + "AGAATGCA"
                           "CAGTGGAACTATCTCAAGAACATGATCATTGGTGTCTTGTTGTTCATCT"
                           "CCGTCATTAGTTGGATCATA")

    
class TestReadMethods(ReadTestSetup):
    """Tests the Read class methods."""
    def test_init(self):
        assert self.single_read.qname == "SRR360147.1"
        assert self.single_read.flag == 77
        assert self.single_read.rname == '*'
        assert self.single_read.pos == 0
        assert self.single_read.mapq == 0
        assert self.single_read.cigar.elements == []
        assert self.single_read.rnext == '*'
        assert self.single_read.pnext == 0
        assert self.single_read.tlen == 0
        assert self.single_read.seq == ("AGATGACGAAGAAGCTTGATCTCACGAANNNNNNNNTTNNCA"
                                        "TCCNNNTNNTNNNNNNNNNNNNNNNNNNNNNNN")
        assert self.single_read.qual == ("HHHHHHHHHHHHHHHHHHHHHHHHH==@############"
                                         "###################################")
        assert self.single_read.tags == ["YT:Z:UP", "YF:Z:NS"]

    def test_eq(self):
        read1, read2, read3 = list(self.indel_algn)
        assert read1 != read2
        assert read1 != read3
        assert read2 != read3

    def test_str(self):
        assert str(self.single_read) == self.sam_string

    def test_get_covered_range(self):
        self.expressed_locus_alignment.rewind()
        test_reads = list(self.expressed_locus_alignment) 
        assert test_reads[0].get_covered_range() == (1, 5)
        assert test_reads[1].get_covered_range() == (3, 7)
        assert test_reads[2].get_covered_range() == (13, 17)
        assert test_reads[3].get_covered_range() == (14, 18)

    def test_has_mate_pair(self):
        assert not self.single_read.has_mate_pair()
        for read in self.mate_pair_align:
            assert read.has_mate_pair()


class TestReadFunctions(ReadTestSetup):
    """Tests of the functions associated with the Read class."""

    def test_consensus(self):
        self.indel_algn.rewind()
        with pytest.raises(sam.UnmappedReadError):
            sam.consensus([self.single_read])
        assert sam.consensus(list(self.indel_algn)) ==\
             "AGATGACGGAAGCTTGATCTCACGAANNNNNNNNTTNNCATCCNNNTNNT"

        assert sam.consensus(self.reverse_complement_align) == \
                                                                "AAAGGGAAAA"

        self.mate_pair_align.rewind()
        mate_pair_loci = sam.expressed_loci(self.mate_pair_align)
        con = sam.consensus(next(mate_pair_loci))
        assert con == self.mate_pair_consensus

    def test_dot_indels(self):
        self.indel_algn.rewind()
        assert list(sam.dot_indels(self.indel_algn)) == \
            [("AGATGACG..GAAGCTTGATCTCACGAA..NNNNNNNNTTNNCATCCNNNTNNT",
              sam.Cigar("8M2D65M"),
              1),
             ("AGATGACGAAGAAGCTTGATCTCACGAA..NNNNNNNNTTNNCATCCNNNTNNA",
              sam.Cigar("77M"),
              1),
             ("AGATGACG..GAAGCTTGATCTCACGAATTNNNNNNNNTTNNCATCCNNNTNNT",
              sam.Cigar("8M2D18M2I24M"),
              1)]

    def test_overlaps(self):
        read1, read2, read3, read4 = \
            list(sam.SamAlignment(TEST_FOLDER + "/test_data/test_overlap.sam"))
        assert not sam.overlaps(read1, read2)
        assert sam.overlaps(read1, read3)
        assert not sam.overlaps(read1, read4)
        assert sam.overlaps(read2, read3)
        assert sam.overlaps(read2, read4)
        assert sam.overlaps(read3, read4)

        self.expressed_locus_alignment.rewind()
        read1, read2, read3, read4 = list(self.expressed_locus_alignment)
        assert sam.overlaps(read1, read2)
        assert not sam.overlaps(read1, read3)
        assert not sam.overlaps(read1, read4)
        assert not sam.overlaps(read2, read3)
        assert not sam.overlaps(read2, read4)
        assert sam.overlaps(read3, read4)

    def test_coverage(self):
        self.expressed_locus_alignment.rewind()
        test_reads = list(self.expressed_locus_alignment)
        assert sam.coverage(test_reads[:1]) == (1,5)
        assert sam.coverage(test_reads[:2]) == (1,7)
        assert sam.coverage(test_reads[:3]) == (1,17)
        assert sam.coverage(test_reads[:4]) == (1,18)
        assert sam.coverage(test_reads[2:]) == (13,18)
        assert sam.coverage(test_reads[1:]) == (3,18)

    def test_convert_indecies(self):
        cigar = [(1, "M"), (2, "M"), (3, "M")]
        c_cigar = sam.convert_indecies(cigar)
        assert c_cigar[0] == (0, 1, "M")
        assert c_cigar[1] == (1, 2, "M")
        assert c_cigar[2] == (3, 3, "M")

    def test_make_dot_queue(self):
        trivial_read = ("AAAA", [(4, "M")], 1)
        trivial_list = [("AAAA", [(4, "M")], 1)]
        assert sam.make_dot_queue(trivial_read, trivial_list) == []
        hard_read = ("AAAA", [(3, "M"), (2, "D"), (1, "M")], 1)
        hard_list = [("AAAA", [(3, "M"), (2, "D"), (1, "M")], 1),
                     ("AAATTA", [(6, "M")], 1),
                     ("CAATTA", [(1, "I"), (4, "M")], 2)]
        assert sam.make_dot_queue(hard_read, hard_list) == \
                [(3, 2), (1, 1)]


class CigarTestSetup(object):
    pass


class TestCigarMethods(CigarTestSetup):
    def test_string_parsing(self):
        assert sam.Cigar('*').elements == []
        assert sam.Cigar('40M').elements == [(40, 'M')]
        assert sam.Cigar('18M2D20M').elements == \
                            [(18, 'M'), (2, 'D'), (20, 'M')]
    
    def test_read_length(self):
        assert sam.Cigar("*").read_length() == 0
        assert sam.Cigar('18M2I20M').read_length() == 40

    def test_aligned_length(self):
        assert sam.Cigar("*").aligned_length() == 0
        assert sam.Cigar('18M2I20M').aligned_length() == 38


class AlignmentTestSetup(ReadTestSetup):
    pass


class TestAlignmentMethods(AlignmentTestSetup):
    def test_str(self):
        for test_file in glob.glob(TEST_FOLDER + "/test_data/*.sam"):
            test_algn = sam.SamAlignment(test_file)
            with open("tmp.sam", "w") as f:
                print(test_algn, file=f)
            test_algn = sam.SamAlignment(test_file)
            assert str(test_algn) == str(sam.SamAlignment("tmp.sam"))
        os.remove("tmp.sam")

    def test_iter(self):
        self.expressed_locus_alignment.rewind()
        for read in self.expressed_locus_alignment:
            assert isinstance(read, sam.Read)

    def test_filter_reads(self):
        self.expressed_locus_alignment.rewind()
        test_reads = self.expressed_locus_alignment
        filtered_reads = test_reads.filter_reads(lambda x: 2 <= x.pos <= 15)
        assert next(filtered_reads).qname == "SRR360147.3"
        assert next(filtered_reads).qname == "SRR360147.2"

    def test_collect_reads(self):
        self.expressed_locus_alignment.rewind()
        test_reads = self.expressed_locus_alignment
        collected_reads = test_reads.collect_reads(lambda x: x.pos < 5)
        first_batch = next(collected_reads)
        second_batch = next(collected_reads)
        first_batch_0 = next(first_batch)
        first_batch_1 = next(first_batch)
        second_batch_0 = next(second_batch)
        second_batch_1 = next(second_batch)
        assert first_batch_0.qname == "SRR360147.1"
        assert first_batch_1.qname == "SRR360147.3"
        assert second_batch_0.qname == "SRR360147.2"
        assert second_batch_1.qname == "SRR360147.4"

    def test_expressed_loci(self):
        self.expressed_locus_alignment.rewind()
        alignment = self.expressed_locus_alignment
        locus_gen = sam.expressed_loci(alignment)
        group1 = next(locus_gen)
        group2 = next(locus_gen)
        assert group1[0].qname == "SRR360147.1"
        assert group1[1].qname == "SRR360147.3"
        assert group2[0].qname == "SRR360147.2"
        assert group2[1].qname == "SRR360147.4"

    def test_mate_pairs(self):
        self.mate_pair_align.rewind()
        first = next(self.mate_pair_align)
        second = next(self.mate_pair_align)

        self.mate_pair_align.rewind()
        pairs = self.mate_pair_align.mate_pairs()

        pair = next(pairs)
        assert pair[0] == first
        assert pair[1] == second

        with pytest.raises(StopIteration):
            next(pairs)

        self.single_read_alignment.rewind()
        pairs = self.single_read_alignment.mate_pairs()

        with pytest.raises(StopIteration):
            next(pairs)


def test_speed_test():
   align = sam.SamAlignment(TEST_FOLDER + "/test_data/speed_test.sam")
   annotation = gff.read_gff(TEST_FOLDER + "/test_data/speed_test.gff")
   for locus in sam.expressed_loci(align):
        gff.in_features(locus, annotation.features)
        sam.consensus(locus)


def test_tear_down():
    os.remove("single_read_test.sam")


if __name__ == '__main__':
    test_speed_test()
