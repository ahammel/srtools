import re


COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


class UnmappedReadError(ValueError):
    """The exception raised when attempting an illegal operation on an unmapped
    read. A consensus sequence cannot be derived from an unmapped read, for
    example.

    """
    pass


class Read():
    """A sam-format sequence read."""

    def read_cigar(cigar):
        """Takes a cigar string, and returns a list of 2-tuples consisting
        of the index (int) and the operation (one-character str).

        """
        return [(int(a), b) for (a, b) in re.findall(r'(\d+)(\D)', cigar)]

    def print_cigar(cigar):
        """Prints a cigar in the standard SAM format"""
        if not cigar:
            cigar_str = '*'
        else:
            cigar_str =\
                ''.join([str(char) for substr in cigar for char in substr])
                # Strings and flattens the list of tuples
        return cigar_str

    def __init__(self, qname, flag, rname, pos, mapq, cigar, rnext, pnext,
                 tlen, seq, qual, tags=[]):
        self.qname = str(qname)
        self.flag = int(flag)
        self.rname = str(rname)
        self.pos = int(pos)
        self.mapq = int(mapq)
        self.cigar = Read.read_cigar(cigar)
        self.rnext = str(rnext)
        self.pnext = int(pnext)
        self.tlen = int(tlen)
        self.seq = str(seq)
        self.qual = str(qual)
        self.tags = [str(x) for x in tags]

    def __eq__(self, other):
        tests = [self.qname == other.qname,
                 self.flag == other.flag,
                 self.rname == other.rname,
                 self.pos == other.pos,
                 self.mapq == other.mapq,
                 self.cigar == other.cigar,
                 self.rnext == other.rnext,
                 self.pnext == other.pnext,
                 self.tlen == other.tlen,
                 self.seq == other.seq,
                 self.qual == other.qual,
                 self.tags == other.tags]

        return all(result == True for result in tests)

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        attrs = [self.qname, self.flag, self.rname, self.pos, self.mapq,
                 Read.print_cigar(self.cigar), self.rnext, self.pnext,
                 self.tlen, self.seq, self.qual] + self.tags
        return "\t".join([str(x) for x in attrs])

    def get_covered_range(self):
        """Returns a tuple consisiting of the first and last position covered
        by the read.

        """
        first_base = self.pos
        last_base = self.pos + sum([i for i, o in self.cigar if o == "M"]) - 1
        return (first_base, last_base)


class Alignment():
    """A sam-format sequence alignment"""
    def __init__(self, head, reads):
        self.head = head
        self.reads = reads

    def __str__(self):
        headstr = self.head
        readstr = "\n".join([str(read) for read in self.reads])
        return headstr + readstr

    def filter_reads(self, function):
        """Returns a list of reads where function(read) returns a truthy value.

        """
        return [r for r in self.reads if function(r)]

    def collect_reads(self, function):
        """Returns a generator which yields lists of consecutive reads which
        all return the same value when the specified function is applied.

        """
        collection = [next(self.reads)]
        for read in self.reads:
            if function(read) == function(collection[-1]):
                collection.append(read)
            else:
                yield collection
                collection = [read]
        yield collection


class Feature():
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


class GenomeAnnotation():
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
        """REturns a generator which yields lists of consecutive features which
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

def reverse_complement(sequence):
    rc = []
    seq = list(sequence)
    while seq:
        rc.append(COMPLEMENT[seq.pop()])
    return "".join(rc)


def parse_sam_read(string):
    """Takes a string in SAMfile format and returns a Read object."""
    fields = string.strip().split()
    return Read(fields[0], fields[1], fields[2], fields[3], fields[4],
                fields[5], fields[6], fields[7], fields[8], fields[9],
                fields[10], tags=fields[11:])


def read_sam(samfile):
    """Creates an Alignment object from a correctly formatted SAM file.
    
    Note: the Alignment.reads object is a generator expression. This (hopefully)
    reduces space complexity to the point where it's possible to actually read
    a real SAM file (which can easily be 50 Gb) in reasonable memory, but there
    is no backtracking. Once you've read a read, it gets garbage collected.

    """
    headlines = []
    with open(samfile) as f:
        while True:
            line = f.readline()
            if line.startswith("@"):
                headlines.append(line)
            else:
                reads = ((parse_sam_read(read)\
                          for read in [line] + f.readlines()\
                          if read))
                break
    return Alignment(head="".join(headlines), reads=reads)


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



def convert_indecies(cigar):
    """Converts a cigar from (n, operator) format to (index, n, operator).
    The index is the zero-based position of the operator, and n is its length.

    """
    index = 0
    c_cigar = []  # c for converted
    for i, o in cigar:
        c_cigar.append((index, i, o))
        index += i
    return c_cigar


def make_dot_queue(stripped_read, stripped_read_list):
    """Returns a queue of positions at which the stripped_read should be dotted
    to indicate an indel. A stripped_read is a 3-tuple of the sequence, the
    cigar and the position. Helper function for dot_indels.

    """
    queue = []
    s_seq, s_cigar, s_pos = stripped_read
    for read in stripped_read_list:
        seq, cigar, pos = read
        c_cigar = convert_indecies(cigar)
        offset = pos - s_pos
        if seq == s_seq and pos == s_pos:
            queue.extend([(i, n) for i, n, o in c_cigar if o in "ND"])
        else:
            queue.extend([(i + offset, n) for i, n, o in c_cigar if o == "I"])
    return queue


def dot_from_queue(stripped_read, queue):
    """Returns a short-read sequence with dots indicating the positions of
    indels. A stripped_read is a 3-tuple of the sequence, the cigar and the
    position.Helper function for dot_indels.

    """
    seq, cigar, pos = stripped_read
    for i, n in queue:
        seq = seq[:i] + ('.' * n) + seq[i:]
    return seq


def dot_indels(reads):
    """Given an iterable of reads, adds dots to the sequences where there are
    indels. Returns a list of 3-tuples of the dotted sequence, the cigar and
    the position.

    """
    stripped_reads = [(read.seq, read.cigar, read.pos) for read in reads]
    dotted_reads = []
    for sr in stripped_reads:
        queue = make_dot_queue(sr, stripped_reads)
        normalized_read = dot_from_queue(sr, queue)
        dotted_reads.append(normalized_read)

    cigars = [c for s, c, p in stripped_reads]
    positions = [p for s, c, p in stripped_reads]
    return zip(dotted_reads, cigars, positions)


def majority(nucleotides, cutoff=0.5):
    """Given a collection of strings, returns the majority rule consensus among
    them. If there is no majority above the cutoff fraction, returns "N".

    """
    for i in nucleotides:
        if len([x for x in nucleotides if x == i]) / len(nucleotides) > cutoff:
            consensus = i
            break
    else:
        consensus = "N"
    return consensus

    
def consensus(reads, cutoff=0.5):
    """Returns the consensus sequence of a collection of reads."""
    all_nucleotides = {}
    for read in dot_indels(reads):
        seq, cigar, pos = read
        if pos == 0:
            raise UnmappedReadError
        else:
            index = pos
            for nuc in seq:
                all_nucleotides.setdefault(index, [])
                all_nucleotides[index].append(nuc)
                index += 1

    consensus_sequence = ""
    for position in range(min(all_nucleotides), max(all_nucleotides) + 1):
        try:
            consensus_sequence += \
                majority(all_nucleotides[position], cutoff=cutoff)
        except KeyError:
            consensus_sequence += 'N'
    return consensus_sequence.replace('.', '')


def overlaps(read1, read2):
    """Returns True if the two reads cover at least one base in common and
    False otherwise.

    """
    r1 = read1.get_covered_range()
    r2 = read2.get_covered_range()

    r_overlap = r1[1] in set(range(r2[0], r2[1]+1))
    l_overlap = r2[1] in set(range(r1[0], r1[1]+1))
    return r_overlap or l_overlap


def expressed_loci(reads):
    """Returns a generator object which yields lists of overlapping reads"""
    locus = []
    for read in reads:
        if not locus or overlaps(locus[-1], read):
            locus.append(read)
        else:
            yield locus
            locus = [read]
    yield locus


def coverage(reads):
    """Returns a tuple consisting of the positions of the first and last base
    covered by the list of reads.

    """
    if not reads:
        return (0,0)

    first, last = reads[0].get_covered_range() 
    for read in reads[1:]:
        x0, x1 = read.get_covered_range()
        if x0 < first:
            first = x0
        if x1 > last:
            last = x1
    return (first, last)


def in_features(reads, features):
    """Returns a boolean indicating whether any of the reads in the first
    argument overlap with any of the features in the second.

    """
    overlap = False
    r0, r1 = coverage(reads)
    for f in features:
        if f.start <= r0 <= f.end or f.start <= r1 <= f.end:
            overlap = True
            break
        elif f.start > r1:
            break
    return overlap
