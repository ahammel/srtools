import re


class UnmappedReadError(ValueError):
    """The exception raised when attempting an illegal operation on an unmapped
    read. A consensus sequence cannot be derived from an unmapped read, for
    example.

    """
    pass


class UnpairedReadError(ValueError):
    """The exception raised when performing a paired-read-specific operation
    on an unpaired read.

    """
    pass


class Read(object):
    """A sam-format sequence read."""

    def __init__(self, qname, flag, rname, pos, mapq, cigar, rnext, pnext,
                 tlen, seq, qual, tags=[]):
        self.qname = str(qname)
        self.flag = int(flag)
        self.rname = str(rname)
        self.pos = int(pos)
        self.mapq = int(mapq)
        self.cigar = Cigar(cigar)
        if rnext == "=":
            self.rnext = self.rname
        else:
            self.rnext = str(rnext)
        self.pnext = int(pnext)
        self.tlen = int(tlen)
        self.seq = str(seq)
        self.qual = str(qual)
        self.tags = [str(x) for x in tags]

    def __eq__(self, other):
        return str(self) == str(other)

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        attrs = [self.qname, self.flag, self.rname, self.pos, self.mapq,
                 str(self.cigar), self.rnext, self.pnext,
                 self.tlen, self.seq, self.qual] + self.tags
        return "\t".join([str(x) for x in attrs])

    def get_covered_range(self):
        """Returns a tuple consisiting of the first and last position covered
        by the read.

        """
        first_base = self.pos
        last_base = self.pos + sum([i for i, o in self.cigar if o == "M"]) - 1
        return (first_base, last_base)

    def has_mate_pair(read):
        """Returns true if the read has a mate pair in the alignment according
        to the bitflag, rnext, and pnext fields.

        """
        flag_set = read.flag & 3 == 3       # Returns True if first two 
                                            # bitflags are set (i.e., if the
                                            # read is paired and mapped in its
                                            # proper pair).
        pnext_set = read.pnext != 0
        rnext_set = read.rnext != '*'

        return flag_set and pnext_set and rnext_set


class Cigar(object):
    """A cigar, as used in SAM-format short reads."""
    def __init__(self, cigar_string):
        self.elements = [(int(a), b) for (a, b) in
                         re.findall(r'(\d+)(\D)', cigar_string)]

    def __iter__(self):
        return iter(self.elements)

    def __next__(self):
        return next(iter(self))

    def __str__(self):
        if not self.elements:
            cigar_str = '*'
        else:
            cigar_str = "".join([str(char) for substr in self.elements
                                 for char in substr])
        return cigar_str

    def __eq__(self, other):
        return self.elements == other.elements

    def __ne__(self, other):
        return self.elements != other.elements


class Alignment(object):
    """A sam-format sequence alignment"""
    def __init__(self, data_file):
        self.data_file = data_file
        self.stream = self.read_generator()

    def __next__(self):
        return next(self.stream)

    def __iter__(self):
        return self

    def filter_reads(self, function):
        """Returns a generator of reads where function(read) returns a truthy
        value.

        """
        for r in self:
            if function(r):
                yield r

    def filter_consecutive_reads(self, function):
        """Returns a generator of consecutive reads where function(read)
        returns a truthy value. Helper method to Alignment.collect_reads.

        """
        reads = iter(self)
        first_read = next(reads)
        test_value = function(first_read)
        yield first_read
        for read in self:
            if function(read) == test_value:
                yield read

    def collect_reads(self, function):
        """Returns a generator which yields generators of consecutive reads
        which all return the same value when the specified function is applied.

        """
        while True:
            yield self.filter_consecutive_reads(function)

    def rewind(self):
        """Calls the read_generator method, thereby reseting the stream of
        reads.

        """
        self.stream = self.read_generator()

    def read_generator(self):
        raise NotImplementedError("Child class must provide read_generator!")


class SamAlignment(Alignment):
    """A stream of reads from a sam file."""
    def __str__(self):
        headstr = self.head()
        readlines = []
        for read in self:
            readlines.append(str(read))
        return headstr + "\n".join(readlines)

    def head(self):
        headlines = []
        with open(self.data_file) as f:
            for line in f:
                if line and line.startswith("@"):
                    headlines.append(line)
                else:
                    break
        return "".join(headlines)

    def read_generator(self):
        with open(self.data_file) as f:
            for line in f:
                if line and not line.startswith("@"):
                    yield parse_sam_read(line)

    def mate_pairs(self):
        """Returns a mate pair generator, which yields mated pairs of reads.
        Calling this method on an unpaired alignment will return an empty
        generator.

        """
        unpaired_reads = {}
        
        for read in self:
            try:
                mate_pair = unpaired_reads[(read.rname, read.pos)]
                yield (mate_pair, read)
                del unpaired_reads[(read.rname, read.pos)]
            except KeyError:
                if read.has_mate_pair():
                    unpaired_reads[(read.rnext, read.pnext)] = read


def parse_sam_read(string):
    """Takes a string in SAMfile format and returns a Read object."""
    fields = string.strip().split()
    return Read(fields[0], fields[1], fields[2], fields[3], fields[4],
                fields[5], fields[6], fields[7], fields[8], fields[9],
                fields[10], tags=fields[11:])


def convert_indecies(cigar):
    """Converts a cigar from (n, operator) format to (index, n, operator).
    The index is the zero-based position of the operator, and n is its length.

    """
    index = 0
    converted_cigar = []
    for i, o in cigar:
        converted_cigar.append((index, i, o))
        index += i
    return converted_cigar


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
    nuc_set = set(nucleotides)
    nucleotide_counts = {x: nucleotides.count(x) for x in nuc_set}

    for x in nucleotide_counts:
        if nucleotide_counts[x] / len(nucleotides) > cutoff:
            consensus = x
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

    consensus_sequence = []
    for position in range(min(all_nucleotides), max(all_nucleotides) + 1):
        if position in all_nucleotides:
            n = majority(all_nucleotides[position], cutoff=cutoff)
            consensus_sequence.append(n)
        else:
            consensus_sequence.append("N")
    consensus = "".join(consensus_sequence)
    return consensus.replace('.', '')


def coverage(reads):
    """Returns a tuple consisting of the positions of the first and last base
    covered by the list of reads.

    """
    if not reads:
        return (0, 0)

    first, last = reads[0].get_covered_range()
    for read in reads[1:]:
        x0, x1 = read.get_covered_range()
        if x0 < first:
            first = x0
        if x1 > last:
            last = x1
    return (first, last)


def tuple_intersection(tuple1, tuple2):
    x0, x1 = tuple1
    y0, y1 = tuple2
    return x0 <= y0 <= x1 or y0 <= x0 <= y1


def overlaps(read1, read2):
    """Returns True if the two reads cover at least one base in common and
    False otherwise.

    """
    return tuple_intersection(read1.get_covered_range(), 
                              read2.get_covered_range())


def in_bounds(read, bounds):
    """Returns True if the covered range of the read is within the bounds.
    Helper function for expressed_loci.

    """
    return tuple_intersection(read.get_covered_range(), bounds)


def expressed_loci(reads):
    """Returns a generator object which yields lists of overlapping reads.
    
    """
    locus = []
    bounds = (0, 0)

    for read in reads:
        if locus and not in_bounds(read, bounds):
            yield locus
            locus = []
            bounds = (0, 0)
        locus.append(read)

        b0, b1 = bounds
        r0, r1 = read.get_covered_range()
        p1 = read.pnext

        bounds = (min([x for x in [b0, r0] if x != 0]), max(b1, r1, p1))

    yield locus
