import pytest
import srtools

testfile = "./marq_aln.test.sam"

sam_str =  """
SRR360147.1     77      *       0       0       *       *       0       0\
 AGATGACGAAGAAGCTTGATCTCACGAANNNNNNNNTTNNCATCCNNNTNNTNNNNNNNNNNNNNNNNNNNNNNN\
 HHHHHHHHHHHHHHHHHHHHHHHHH==@###############################################\
 YT:Z:UP YF:Z:NS"""

read = srtools.Read(
        qname = "SRR360147.1",
        flag = 77,
        rname = '*',
        pos = 0,
        mapq = 0,
        cigar = '*',
        rnext = '*',
        pnext = 0,
        tlen = 0,
        seq = ("AGATGACGAAGAAGCTTGATCTCACGAANNNNNNNNTTNNCATCCNNNT"
               "NNTNNNNNNNNNNNNNNNNNNNNNNN"),
        qual = ("HHHHHHHHHHHHHHHHHHHHHHHHH==@####################"
                "###########################"),
        tags = ["YT:Z:UP","YF:Z:NS"]
        )

headstr="""\
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

algn = srtools.Alignment(head = headstr, reads = [read])

def test_Read_class():
    """Tests whether the Read.__init__ method is working correctly."""
    assert read.qname == "SRR360147.1"
    assert read.flag == 77
    assert read.rname == '*'
    assert read.pos == 0
    assert read.mapq == 0
    assert read.cigar == '*'
    assert read.rnext == '*'
    assert read.pnext == 0
    assert read.tlen == 0
    assert read.seq == ("AGATGACGAAGAAGCTTGATCTCACGAANNNNNNNNTTNNCATCCNNNT"
                        "NNTNNNNNNNNNNNNNNNNNNNNNNN")
    assert read.qual == ("HHHHHHHHHHHHHHHHHHHHHHHHH==@####################"
                         "###########################")
    assert read.tags == ["YT:Z:UP","YF:Z:NS"]

def test_parse_sam_read():
    test_read = srtools.parse_sam_read(sam_str)
    assert test_read.qname == read.qname
    assert read.flag == test_read.flag 
    assert read.rname == test_read.rname 
    assert read.pos == test_read.pos 
    assert read.mapq == test_read.mapq 
    assert read.cigar == test_read.cigar 
    assert read.rnext == test_read.rnext 
    assert read.pnext == test_read.pnext 
    assert read.tlen == test_read.tlen 
    assert read.seq == test_read.seq 
    assert read.qual == test_read.qual 
    assert read.tags == test_read.tags
    

def test_read_sam():
    test_algn = srtools.read_sam(testfile)
    assert test_algn.head == algn.head
    for test_read in test_algn.reads:
        assert test_read.qname == read.qname
        assert read.flag == test_read.flag 
        assert read.rname == test_read.rname 
        assert read.pos == test_read.pos 
        assert read.mapq == test_read.mapq 
        assert read.cigar == test_read.cigar 
        assert read.rnext == test_read.rnext 
        assert read.pnext == test_read.pnext 
        assert read.tlen == test_read.tlen 
        assert read.seq == test_read.seq 
        assert read.qual == test_read.qual 
        assert read.tags == test_read.tags
