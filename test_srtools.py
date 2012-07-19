import pytest
import srstats
import srtools
import os
import glob


class ReadTestSetup(object):
    """Does the setup for the tests of the Read classes and associated functions."""
    sam_string = ("SRR360147.1\t77\t*\t0\t0\t*\t*\t0\t0\tAGATGACGAAGAAGCTTGATC"
                  "TCACGAANNNNNNNNTTNNCATCCNNNTNNTNNNNNNNNNNNNNNNNNNNNNNN\tHHH"
                  "HHHHHHHHHHHHHHHHHHHHHH==@##################################"
                  "#############\tYT:Z:UP\tYF:Z:NS")

    single_read = srtools.Read(
            qname="SRR360147.1",
            flag=77,
            rname='*',
            pos=0,
            mapq=0,
            cigar='*',
            rnext='*',
            pnext=0,
            tlen=0,
            seq=("AGATGACGAAGAAGCTTGATCTCACGAANNNNNNNNTTNNCATCCNNNT"
                 "NNTNNNNNNNNNNNNNNNNNNNNNNN"),
            qual=("HHHHHHHHHHHHHHHHHHHHHHHHH==@####################"
                  "###########################"),
            tags=["YT:Z:UP", "YF:Z:NS"]
            )

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

    single_read_alignment = srtools.SamAlignment("single_read_test.sam") 

    indel_algn = srtools.SamAlignment("./test/test_indel.sam")

    #indel_algn.reads = list(indel_algn.reads)

    reverse_complement_align = srtools.SamAlignment("test/test_rc_consensus.sam")

    expressed_locus_alignment = srtools.SamAlignment("test/test_expressed_locus.sam")

    
class TestReadMethods(ReadTestSetup):
    """Tests the Read class methods."""
    def test_init(self):
        assert self.single_read.qname == "SRR360147.1"
        assert self.single_read.flag == 77
        assert self.single_read.rname == '*'
        assert self.single_read.pos == 0
        assert self.single_read.mapq == 0
        assert self.single_read.cigar == []
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
        test_reads = list(self.expressed_locus_alignment) 
        assert test_reads[0].get_covered_range() == (1, 5)
        assert test_reads[1].get_covered_range() == (3, 7)
        assert test_reads[2].get_covered_range() == (13, 17)
        assert test_reads[3].get_covered_range() == (14, 18)


class TestReadFunctions(ReadTestSetup):
    """Tests of the functions associated with the Read class."""

    def test_read_cigar(self):
        assert srtools.read_cigar('*') == []
        assert srtools.read_cigar('40M') == [(40, 'M')]
        assert srtools.read_cigar('18M2D20M') == \
                    [(18, 'M'), (2, 'D'), (20, 'M')]

    def test_parse_sam_read(self):
        test_read = srtools.parse_sam_read(self.sam_string)
        assert test_read == self.single_read

    def test_consensus(self):
        self.indel_algn.rewind()
        with pytest.raises(srtools.UnmappedReadError):
            srtools.consensus([self.single_read])
        assert srtools.consensus(list(self.indel_algn)) ==\
             "AGATGACGGAAGCTTGATCTCACGAANNNNNNNNTTNNCATCCNNNTNNT"

        assert srtools.consensus(self.reverse_complement_align) == \
                                                                "AAAGGGAAAA"

    def test_dot_indels(self):
        self.indel_algn.rewind()
        assert list(srtools.dot_indels(self.indel_algn)) == \
            [("AGATGACG..GAAGCTTGATCTCACGAA..NNNNNNNNTTNNCATCCNNNTNNT",
              [(8, 'M'), (2, 'D'), (65, 'M')],
              1),
             ("AGATGACGAAGAAGCTTGATCTCACGAA..NNNNNNNNTTNNCATCCNNNTNNA",
              [(77, 'M')],
              1),
             ("AGATGACG..GAAGCTTGATCTCACGAATTNNNNNNNNTTNNCATCCNNNTNNT",
              [(8, 'M'), (2, 'D'), (18, 'M'), (2, 'I'), (24, 'M')],
              1)]

    def test_overlaps(self):
        read1, read2, read3, read4 = list(srtools.SamAlignment("test/test_overlap.sam"))
        assert not srtools.overlaps(read1, read2)
        assert srtools.overlaps(read1, read3)
        assert not srtools.overlaps(read1, read4)
        assert srtools.overlaps(read2, read3)
        assert srtools.overlaps(read2, read4)
        assert srtools.overlaps(read3, read4)

        self.expressed_locus_alignment.rewind()
        read1, read2, read3, read4 = list(self.expressed_locus_alignment)
        assert srtools.overlaps(read1, read2)
        assert not srtools.overlaps(read1, read3)
        assert not srtools.overlaps(read1, read4)
        assert not srtools.overlaps(read2, read3)
        assert not srtools.overlaps(read2, read4)
        assert srtools.overlaps(read3, read4)

    def test_coverage(self):
        self.expressed_locus_alignment.rewind()
        test_reads = list(self.expressed_locus_alignment)
        assert srtools.coverage(test_reads[:1]) == (1,5)
        assert srtools.coverage(test_reads[:2]) == (1,7)
        assert srtools.coverage(test_reads[:3]) == (1,17)
        assert srtools.coverage(test_reads[:4]) == (1,18)
        assert srtools.coverage(test_reads[2:]) == (13,18)
        assert srtools.coverage(test_reads[1:]) == (3,18)

    def test_in_features(self):
        self.expressed_locus_alignment.rewind()
        alignment = self.expressed_locus_alignment
        annotation = srtools.read_gff("test/test.gff")
        genes = annotation.filter_features(lambda x: x.f_type == "gene")
        locus_gen = srtools.expressed_loci(alignment)
        group1 = next(locus_gen)
        group2 = next(locus_gen)
        assert not srtools.in_features(group1, genes)
        assert srtools.in_features(group2, genes)


class SequenceTestSetup(object):
    fasta_reads = srtools.read_fasta("test/fasta.fa")

class TestSequenceFunctions(SequenceTestSetup):
    def test_reverse_complement(self):
        assert srtools.reverse_complement("") == ""
        assert srtools.reverse_complement("N") == "N"
        assert srtools.reverse_complement("ACGT") == "ACGT"
        assert srtools.reverse_complement("GCCAT") == "ATGGC"

    def test_read_fasta(self):
        assert len(self.fasta_reads) == 2
        assert self.fasta_reads["read1"] == "AAAAAAAAAAAAAA"
        assert self.fasta_reads["read2"] == "ATGAATGAATAAAAAAAAATAATTATTTCAT" 

    def test_convert_indecies(self):
        cigar = [(1, "M"), (2, "M"), (3, "M")]
        c_cigar = srtools.convert_indecies(cigar)
        assert c_cigar[0] == (0, 1, "M")
        assert c_cigar[1] == (1, 2, "M")
        assert c_cigar[2] == (3, 3, "M")

    def test_make_dot_queue(self):
        trivial_read = ("AAAA", [(4, "M")], 1)
        trivial_list = [("AAAA", [(4, "M")], 1)]
        assert srtools.make_dot_queue(trivial_read, trivial_list) == []
        hard_read = ("AAAA", [(3, "M"), (2, "D"), (1, "M")], 1)
        hard_list = [("AAAA", [(3, "M"), (2, "D"), (1, "M")], 1),
                     ("AAATTA", [(6, "M")], 1),
                     ("CAATTA", [(1, "I"), (4, "M")], 2)]
        assert srtools.make_dot_queue(hard_read, hard_list) == \
                [(3, 2), (1, 1)]

    def test_gc_content(self):
        pytest.raises(srtools.NullSequenceError, srtools.gc_content, "")
        pytest.raises(srtools.NullSequenceError, srtools.gc_content, "NNNNN")
        assert srtools.gc_content("AAATTT") == 0.0
        assert srtools.gc_content("GGGGCC") == 1.0
        assert srtools.gc_content("ACGTACGT") == 0.5
        assert srtools.gc_content("ACCTACGT") == 0.5
    
    def test_random_sequence(self):
        seq = srtools.random_sequence(10)
        assert len(seq) == 10
        for letter in seq:
            assert letter in "ACGT"

    def test_randomize_sequence(self):
        seq = "AACGNNNNAACG"
        rseq = srtools.randomize_sequence(seq)
        assert rseq[4:8] == "NNNN"
        assert rseq != seq

    def test_reading_frames(self):
        assert srtools.reading_frames(self.fasta_reads["read1"]) == \
            [["AAA", "AAA", "AAA", "AAA", "AA"],
             ["TTT", "TTT", "TTT", "TTT", "TT"],
             ["A", "AAA", "AAA", "AAA", "AAA", "A"],
             ["T", "TTT", "TTT", "TTT", "TTT", "T"],
             ["AA", "AAA", "AAA", "AAA", "AAA"],
             ["TT", "TTT", "TTT", "TTT", "TTT"]]

        assert srtools.reading_frames(self.fasta_reads["read2"]) == \
            [["ATG", "AAT", "GAA", "TAA", "AAA", "AAA", "ATA", "ATT", "ATT", "TCA", "T"],
             ["ATG", "AAA", "TAA", "TTA", "TTT", "TTT", "TTT", "ATT", "CAT", "TCA", "T"],
             ["A", "TGA", "ATG", "AAT", "AAA", "AAA", "AAA", "TAA", "TTA", "TTT", "CAT"],
             ["A", "TGA", "AAT", "AAT", "TAT", "TTT", "TTT", "TTA", "TTC", "ATT", "CAT"],
             ["AT", "GAA", "TGA", "ATA", "AAA", "AAA", "AAT", "AAT", "TAT", "TTC", "AT"],
             ["AT", "GAA", "ATA", "ATT", "ATT", "TTT", "TTT", "TAT", "TCA", "TTC", "AT"]]

    def test_open_reading_frames(self):
        assert srtools.open_reading_frames(self.fasta_reads["read1"]) == []
        assert srtools.open_reading_frames(self.fasta_reads["read2"]) == \
            ["ATGAATGAA", "ATGAAA", "ATGAATAAAAAAAAA"]


class AlignmentTestSetup(ReadTestSetup):
    pass


class TestAlignmentMethods(AlignmentTestSetup):
    def test_str(self):
        for test_file in glob.glob("test/*.sam"):
            test_algn = srtools.SamAlignment(test_file)
            with open("tmp.sam", "w") as f:
                print(test_algn, file=f)
            test_algn = srtools.SamAlignment(test_file)
            assert str(test_algn) == str(srtools.SamAlignment("tmp.sam"))
        os.remove("tmp.sam")

    def test_iter(self):
        self.expressed_locus_alignment.rewind()
        for read in self.expressed_locus_alignment:
            assert isinstance(read, srtools.Read)

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
        locus_gen = srtools.expressed_loci(alignment)
        group1 = next(locus_gen)
        group2 = next(locus_gen)
        assert group1[0].qname == "SRR360147.1"
        assert group1[1].qname == "SRR360147.3"
        assert group2[0].qname == "SRR360147.2"
        assert group2[1].qname == "SRR360147.4"


class FeatureTestSetup(object):
    test_feature = srtools.parse_gff_feature("Chr1\tTAIR9\tchromosome\t1\t"
                                             "30427671\t.\t.\t.\tID=Chr1;"
                                             "Name=Chr1")

    test_annotation = srtools.read_gff("test/test.gff")


class TestFeatureFunctions(FeatureTestSetup):
    def test_parse_gff_feature(self):
        assert self.test_feature.sequence == "Chr1"
        assert self.test_feature.source == "TAIR9"
        assert self.test_feature.f_type == "chromosome"
        assert self.test_feature.start == 1
        assert self.test_feature.end == 30427671
        assert self.test_feature.score == None
        assert self.test_feature.strand == None
        assert self.test_feature.frame == None
        assert self.test_feature.attribute == "ID=Chr1;Name=Chr1"


class TestFeatureMethods(FeatureTestSetup):
    pass


class TestGenomeAnnotationMethods(FeatureTestSetup):
    def test_collect_features(self):
        collected_features =\
            self.test_annotation.collect_features(lambda x: "1" in x.sequence)
        first_batch = next(collected_features)
        second_batch = next(collected_features)
        assert first_batch[0].attribute == "f1"
        assert first_batch[1].attribute == "f2"
        assert second_batch[0].attribute == "f3"


class TestGenomeAnnotationFunctions(FeatureTestSetup):
    def test_read_gff(self):
        assert self.test_annotation.head == "## Header\n"
        assert self.test_annotation.features[0].sequence == "Chr1"
        assert self.test_annotation.features[0].source == "TAIR9"
        assert self.test_annotation.features[0].f_type == "chromosome"
        assert self.test_annotation.features[0].start == 1
        assert self.test_annotation.features[0].end == 30427671
        assert self.test_annotation.features[0].score == None
        assert self.test_annotation.features[0].strand == None
        assert self.test_annotation.features[0].frame == None
        assert self.test_annotation.features[0].attribute == "f1"


class TestORFFunctions():
    reads = srtools.read_fasta("test/fasta.fa")


def test_summary_statistics():
    srstats.print_summary_statistics("test/speed_test.sam",
                                     output_file="test_summary.txt")
    test_file = open("test_summary.txt")
    known_file = open("test/david_summary.txt")
    test_lines = test_file.readlines()
    known_lines = known_file.readlines()
    for line in test_lines:
        assert line in known_lines
    for line in known_lines:
        assert line in test_lines
    test_file.close()
    known_file.close()
    os.remove("test_summary.txt")


def test_speed_test():
   align = srtools.read_sam("test/speed_test.sam")
   annotation = srtools.read_gff("test/speed_test.gff")
   for locus in srtools.expressed_loci(align.reads):
        srtools.in_features(locus, annotation.features)
        srtools.consensus(locus)

def test_tear_down():
    os.remove("single_read_test.sam")
