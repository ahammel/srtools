import pytest
import srtools
import os
import glob

test_1read = "./test/test_1read.sam"

sam_str = """\
SRR360147.1\t77\t*\t0\t0\t*\t*\t0\t0\
\tAGATGACGAAGAAGCTTGATCTCACGAANNNNNNNNTTNNCATCCNNNTNNTNNNNNNNNNNNNNNNNNNNNNNN\
\tHHHHHHHHHHHHHHHHHHHHHHHHH==@###############################################\
\tYT:Z:UP\tYF:Z:NS"""

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

headstr = """\
@HD\tVN:1.0\tSO:unsorted
@SQ\tSN:1\tLN:30427671
@SQ\tSN:2\tLN:19698289
@SQ\tSN:3\tLN:23459830
@SQ\tSN:4\tLN:18585056
@SQ\tSN:5\tLN:26975502
@SQ\tSN:Mt\tLN:366924
@SQ\tSN:Pt\tLN:154478
@PG\tID:bowtie2\tPN:bowtie2\tVN:2.0.0-beta5
"""

algn = srtools.Alignment(head=headstr, reads=[single_read])

indel_algn = srtools.read_sam("./test/test_indel.sam")
rc_align = srtools.read_sam("test/test_rc_consensus.sam")

class TestRead:
    def test_read_cigar(self):
        """Note: there should never be a need to use Read.read_cigar like
        this.

        """
        assert srtools.Read.read_cigar('*') == []
        assert srtools.Read.read_cigar('40M') == [(40, 'M')]
        assert srtools.Read.read_cigar('18M2D20M') == \
                    [(18, 'M'), (2, 'D'), (20, 'M')]

    def test_init(self):
        assert single_read.qname == "SRR360147.1"
        assert single_read.flag == 77
        assert single_read.rname == '*'
        assert single_read.pos == 0
        assert single_read.mapq == 0
        assert single_read.cigar == []
        assert single_read.rnext == '*'
        assert single_read.pnext == 0
        assert single_read.tlen == 0
        assert single_read.seq == ("AGATGACGAAGAAGCTTGATCTCACGAANNNNNNNNTTNNCA"
                                   "TCCNNNTNNTNNNNNNNNNNNNNNNNNNNNNNN")
        assert single_read.qual == ("HHHHHHHHHHHHHHHHHHHHHHHHH==@############"
                                    "###################################")
        assert single_read.tags == ["YT:Z:UP", "YF:Z:NS"]
        assert single_read == single_read

    def test_eq(self):
        read1, read2, read3 = indel_algn.reads
        assert read1 != read2
        assert read1 != read3
        assert read2 != read3

    def test_str(self):
        assert str(single_read) == sam_str

    def test_reverse(self):
        test_rc = srtools.read_sam("test/test_rc.sam")
        seq, rc = test_rc.reads
        print("Seq: ", seq)
        print("RRC: ", rc.reverse())
        print("RC: ", rc)
        print("RSeq: ", seq.reverse())
        assert seq.reverse() == rc
        assert rc.reverse() == seq

class TestAlignment:
    def test_str(self):
        for test_file in glob.glob("test/*.sam"):
            test_algn = srtools.read_sam(test_file)
            with open("tmp.sam", "w") as f:
                print(test_algn, file=f)
            assert str(test_algn) == str(srtools.read_sam("tmp.sam"))
        os.remove("tmp.sam")


def TestFeature():
    pass


def test_reverse_complement():
    assert srtools.reverse_complement("") == ""
    assert srtools.reverse_complement("N") == "N"
    assert srtools.reverse_complement("ACGT") == "ACGT"
    assert srtools.reverse_complement("GCCAT") == "ATGGC"

def test_parse_sam_read():
    test_read = srtools.parse_sam_read(sam_str)
    assert test_read == single_read


def test_read_sam():
    test_algn = srtools.read_sam(test_1read)
    assert test_algn.head == algn.head
    for test_read in test_algn.reads:
        assert single_read == test_read


def test_parse_gff_feature():
    feature = srtools.parse_gff_feature("Chr1\tTAIR9\tchromosome\t1\t30427671"
                                        "\t.\t.\t.\tID=Chr1;Name=Chr1")
    assert feature.sequence == "Chr1"
    assert feature.source == "TAIR9"
    assert feature.f_type == "chromosome"
    assert feature.start == 1
    assert feature.end == 30427671
    assert feature.score == None
    assert feature.strand == None
    assert feature.frame == None
    assert feature.attribute == "ID=Chr1;Name=Chr1"


def test_read_gff():
    features = srtools.read_gff("test/test.gff")
    assert features.head == "## Header\n"
    assert features.features[0].sequence == "Chr1"
    assert features.features[0].source == "TAIR9"
    assert features.features[0].f_type == "chromosome"
    assert features.features[0].start == 1
    assert features.features[0].end == 30427671
    assert features.features[0].score == None
    assert features.features[0].strand == None
    assert features.features[0].frame == None
    assert features.features[0].attribute == "ID=Chr1;Name=Chr1"


def test_convert_indecies():
    cigar = [(1, "M"), (2, "M"), (3, "M")]
    c_cigar = srtools.convert_indecies(cigar)
    assert c_cigar[0] == (0, 1, "M")
    assert c_cigar[1] == (1, 2, "M")
    assert c_cigar[2] == (3, 3, "M")


def test_make_dot_queue():
    trivial_read = ("AAAA", [(4, "M")], 1)
    trivial_list = [("AAAA", [(4, "M")], 1)]
    assert srtools.make_dot_queue(trivial_read, trivial_list) == []
    hard_read = ("AAAA", [(3, "M"), (2, "D"), (1, "M")], 1)
    hard_list = [("AAAA", [(3, "M"), (2, "D"), (1, "M")], 1),
                 ("AAATTA", [(6, "M")], 1),
                 ("CAATTA", [(1, "I"), (4, "M")], 2)]
    assert srtools.make_dot_queue(hard_read, hard_list) == \
            [(3, 2), (1, 1)]


def test_dot_indels():
    assert list(srtools.dot_indels(indel_algn.reads)) == \
        [("AGATGACG..GAAGCTTGATCTCACGAA..NNNNNNNNTTNNCATCCNNNTNNT",
          [(8, 'M'), (2, 'D'), (65, 'M')],
          1),
         ("AGATGACGAAGAAGCTTGATCTCACGAA..NNNNNNNNTTNNCATCCNNNTNNA",
          [(77, 'M')],
          1),
         ("AGATGACG..GAAGCTTGATCTCACGAATTNNNNNNNNTTNNCATCCNNNTNNT",
          [(8, 'M'), (2, 'D'), (18, 'M'), (2, 'I'), (24, 'M')],
          1)]


def test_consensus():
    with pytest.raises(srtools.UnmappedReadError):
        srtools.consensus([single_read])
    assert srtools.consensus(indel_algn.reads) ==\
         "AGATGACGGAAGCTTGATCTCACGAANNNNNNNNTTNNCATCCNNNTNNT"

    assert srtools.consensus(rc_align.reads) == "AAACCCTTTT"


def test_overlaps():
    read1, read2, read3, read4 = srtools.read_sam("test/test_overlap.sam").reads
    assert not srtools.overlaps(read1, read2)
    assert srtools.overlaps(read1, read3)
    assert not srtools.overlaps(read1, read4)
    assert srtools.overlaps(read2, read3)
    assert srtools.overlaps(read2, read4)
    assert srtools.overlaps(read3, read4)
