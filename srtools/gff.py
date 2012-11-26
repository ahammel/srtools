"""Tools for reading and manipulating General Feature Format files and for 
comparing SAM reads aligned to the relevant reference genome.

"""
from srtools import sam

class Feature(object):
    """A GFF genomic feature."""
    def __init__(self, feature_string):
        fields = feature_string.strip().split("\t")

        while True:
            try:
                fields[fields.index('.')] = None
            except ValueError:
                break

        self.sequence = fields[0]
        self.source = fields[1]
        self.f_type = fields[2]
        self.start = int(fields[3])
        self.end = int(fields[4])
        if fields[5]:
            self.score = float(fields[5])
        else:
            self.score = fields[5]
        self.strand = fields[6]
        self.frame = fields[7]
        self.attribute = fields[8]


class GenomeAnnotation(object):
    """A genome-spanning collection of Features"""
    def __init__(self, gff_file):
        self.data_file = gff_file
        self.features = self.read_features()

    def __iter__(self):
        return iter(self.features)

    def read_features(self):
        with open(self.data_file) as f:
            return [Feature(line) for line in f if not line.startswith("#")]

    def head(self):
        with open(self.data_file) as f:
            return "".join(line for line in f if line.startswith("##"))

    def filter_features(self, function):
        """Returns a list of features where function(feature) reutrns a truthy
        value.

        """
        return [f for f in self.features if function(f)]

    def collect_features(self, function):
        """Returns a generator which yields lists of consecutive features which
        all return the same value when the specified function is applied.

        """
        collection = [self.features[0]]
        for feature in self.features[1:]:
            if function(feature) == function(collection[-1]):
                collection.append(feature)
            else:
                yield collection
                collection = [feature]
        yield collection


def in_features(reads, features):
    """Returns a boolean indicating whether any of the reads in the first
    argument overlap with any of the features in the second.

    """
    overlap = False
    r0, r1 = sam.coverage(reads)
    for f in features:
        if f.start <= r0 <= f.end or f.start <= r1 <= f.end:
            overlap = True
            break
        elif f.start > r1:
            break
    return overlap
