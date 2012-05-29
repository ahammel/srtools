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

two_read_algn = srtools.read_sam("./test/test_2reads.sam")

def test_Read_class():
    """Tests whether the Read.__init__ method is working correctly."""
    assert single_read.qname == "SRR360147.1"
    assert single_read.flag == 77
    assert single_read.rname == '*'
    assert single_read.pos == 0
    assert single_read.mapq == 0
    assert single_read.cigar == '*'
    assert single_read.rnext == '*'
    assert single_read.pnext == 0
    assert single_read.tlen == 0
    assert single_read.seq == ("AGATGACGAAGAAGCTTGATCTCACGAANNNNNNNNTTNNCATCCNN"
                               "NTNNTNNNNNNNNNNNNNNNNNNNNNNN")
    assert single_read.qual == ("HHHHHHHHHHHHHHHHHHHHHHHHH==@#################"
                                "##############################")
    assert single_read.tags == ["YT:Z:UP", "YF:Z:NS"]


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


def test_read_sam():
    test_algn = srtools.read_sam(test_1read)
    assert test_algn.head == algn.head
    for test_read in test_algn.reads:
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

def test_consensus():
    with pytest.raises(srtools.UnmappedReadError):
        srtools.consensus([single_read])
    assert srtools.consensus(two_read_algn.reads) == single_read.seq
