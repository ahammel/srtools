from srtools import sr_sam

class Feature(object):
    """A GFF genomic feature."""
    def __init__(self, sequence, source, f_type, start, end, score,
                 strand, frame, attribute):
        self.sequence = str(sequence)
        self.source = str(source)
        self.f_type = str(f_type)
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attribute = attribute


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


def parse_gff_feature(feature_string):
    """Creates a Feature object from a GFF feature string."""
    fields = feature_string.strip().split("\t")

    while True:
        try:
            fields[fields.index('.')] = None
        except ValueError:
            break

    sequence = fields[0]
    source = fields[1]
    f_type = fields[2]
    start = int(fields[3])
    end = int(fields[4])
    if fields[5]:
        score = float(fields[5])
    else:
        score = fields[5]
    strand = fields[6]
    frame = fields[7]
    attribute = fields[8]

    return Feature(sequence, source, f_type, start, end, score, strand,
                   frame, attribute)


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
                features.append(parse_gff_feature(line))
    return GenomeAnnotation(head="".join(headlines), features=features)


def in_features(reads, features):
    """Returns a boolean indicating whether any of the reads in the first
    argument overlap with any of the features in the second.

    """
    overlap = False
    r0, r1 = sr_sam.coverage(reads)
    for f in features:
        if f.start <= r0 <= f.end or f.start <= r1 <= f.end:
            overlap = True
            break
        elif f.start > r1:
            break
    return overlap
