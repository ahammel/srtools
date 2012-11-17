from srtools import gff, sam
from srtools.test import TEST_FOLDER
from srtools.test.test_sam import AlignmentTestSetup

class FeatureTestSetup(AlignmentTestSetup):
    test_feature = gff.Feature("Chr1\tTAIR9\tchromosome\t1\t"
                                            "30427671\t.\t.\t.\tID=Chr1;"
                                            "Name=Chr1")

    test_annotation = gff.read_gff(TEST_FOLDER + "/test_data/test.gff")


class TestFeatureMethods(FeatureTestSetup):
    def test_init(self):
        assert self.test_feature.sequence == "Chr1"
        assert self.test_feature.source == "TAIR9"
        assert self.test_feature.f_type == "chromosome"
        assert self.test_feature.start == 1
        assert self.test_feature.end == 30427671
        assert self.test_feature.score == None
        assert self.test_feature.strand == None
        assert self.test_feature.frame == None
        assert self.test_feature.attribute == "ID=Chr1;Name=Chr1"


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


    def test_in_features(self):
        self.expressed_locus_alignment.rewind()
        alignment = self.expressed_locus_alignment
        annotation = gff.read_gff(TEST_FOLDER + "/test_data/test.gff")
        genes = annotation.filter_features(lambda x: x.f_type == "gene")
        locus_gen = sam.expressed_loci(alignment)
        group1 = next(locus_gen)
        group2 = next(locus_gen)
        assert not gff.in_features(group1, genes)
        assert gff.in_features(group2, genes)
