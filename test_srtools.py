import pytest
import srtools

test_1read = "./test/test_1read.sam"

sam_str = """
SRR360147.1     77      *       0       0       *       *       0       0\
 AGATGACGAAGAAGCTTGATCTCACGAANNNNNNNNTTNNCATCCNNNTNNTNNNNNNNNNNNNNNNNNNNNNNN\
 HHHHHHHHHHHHHHHHHHHHHHHHH==@###############################################\
 YT:Z:UP YF:Z:NS"""

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


class TestRead:
    def test_parse_cigar(self):
        """Note: there should never be a need to use Read.parse_cigar like
        this.

        """
        assert srtools.Read.parse_cigar('*') == []
        assert srtools.Read.parse_cigar('40M') == [(40, 'M')]
        assert srtools.Read.parse_cigar('18M2D20M') == \
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


def test_parse_sam_read():
    test_read = srtools.parse_sam_read(sam_str)
    assert test_read.qname == single_read.qname
    assert single_read.flag == test_read.flag
    assert single_read.rname == test_read.rname
    assert single_read.pos == test_read.pos
    assert single_read.mapq == test_read.mapq
    assert single_read.cigar == test_read.cigar
    assert single_read.rnext == test_read.rnext
    assert single_read.pnext == test_read.pnext
    assert single_read.tlen == test_read.tlen
    assert single_read.seq == test_read.seq
    assert single_read.qual == test_read.qual
    assert single_read.tags == test_read.tags
    assert single_read == test_read


def test_read_sam():
    test_algn = srtools.read_sam(test_1read)
    assert test_algn.head == algn.head
    for test_read in test_algn.reads:
        assert single_read.qname == test_read.qname
        assert single_read.flag == test_read.flag
        assert single_read.rname == test_read.rname
        assert single_read.pos == test_read.pos
        assert single_read.mapq == test_read.mapq
        assert single_read.cigar == test_read.cigar
        assert single_read.rnext == test_read.rnext
        assert single_read.pnext == test_read.pnext
        assert single_read.tlen == test_read.tlen
        assert single_read.seq == test_read.seq
        assert single_read.qual == test_read.qual
        assert single_read.tags == test_read.tags


def test_dot_deletions():
    assert srtools.dot_deletions(indel_algn.reads[0]) == \
        "AGATGACG..GAAGCTTGATCTCACGAANNNNNNNNTTNNCATCCNNNTNNT"
    assert srtools.dot_deletions(single_read) == single_read.seq


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
          [(8, 'M'), (2, 'D'), (65,'M')],
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
