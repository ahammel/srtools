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
    def __init__(self, head, features):
        self.head = head
        self.features = features

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


def read_gff(gff_file):
    """Creates a GenomeAnnotation object from a GFF file"""
    headlines = []
    features = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("##"):
                headlines.append(line)
            elif line.startswith("#"):
                pass
            else:
                features.append(Feature(line))
    return GenomeAnnotation(head="".join(headlines), features=features)


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
