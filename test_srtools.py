import pytest
import srtools


def test_Read_class():
    """Tests whether the Read.__init__ method is working correctly.

    """
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
            seq = ("AGATGACGAAGAAGCTTGATCTCACGAANNNNNNNNTTNNCATCCNNNT",
                "NNTNNNNNNNNNNNNNNNNNNNNNNN"),
            qual = ("HHHHHHHHHHHHHHHHHHHHHHHHH==@####################",
                "###########################")
            )
    assert read.qname == "SRR360147.1"
    assert read.flag == 77
    assert read.rname == '*'
    assert read.pos == 0
    assert read.mapq == 0
    assert read.cigar == '*'
    assert read.rnext == '*'
    assert read.pnext == 0
    assert read.tlen == 0
    assert read.seq == ("AGATGACGAAGAAGCTTGATCTCACGAANNNNNNNNTTNNCATCCNNNT",
                       "NNTNNNNNNNNNNNNNNNNNNNNNNN")
    assert read.qual == ("HHHHHHHHHHHHHHHHHHHHHHHHH==@####################",
                         "###########################")
