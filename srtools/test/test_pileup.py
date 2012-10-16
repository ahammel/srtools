from srtools import pileup
try:
    import pytest
except ImportError:
    import nose as pytest

TEST_FOLDER = "/home/alex/repos/py-srtools/srtools/test"

class PileupTestSetup(object):
    test_pileup = pileup.Pileup(TEST_FOLDER + "/test_data/test.pileup")

    test_pileup_read = pileup.PileupRead("Chr1    3627    N       "
                                         "1       ^SG     B")

class TestPileupReadMethods(PileupTestSetup):
    def test_init(self):
        assert self.test_pileup_read.chromosome_name == "Chr1"
        assert self.test_pileup_read.coordinate == 3627
        assert self.test_pileup_read.reference_base == "N"
        assert self.test_pileup_read.read_bases == "1"
        assert self.test_pileup_read.read_bases == "1"
        assert self.test_pileup_read.read_qualities == "^SG"
        assert self.test_pileup_read.alignment_mapping_qualities == "B"

    def test_bad_data(self):
        with pytest.raises(ValueError):
            pileup.PileupRead("")
        with pytest.raises(ValueError):
            pileup.PileupRead("Banana")
        with pytest.raises(ValueError):
            pileup.PileupRead(42)


class TestPileupMethods(PileupTestSetup):
    def test_init(self):
        reads = list(self.test_pileup)
        assert reads[0] == self.test_pileup_read
        assert len(reads) == 10
