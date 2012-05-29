class UnmappedReadError(ValueError):
    """The exception raised when attempting an illegal operation on an unmapped
    read. A consensus sequence cannot be derived from an unmapped read, for
    example.

    """
    pass


class Read():
    """A sam-format sequence read."""
    def __init__(self, qname, flag, rname, pos, mapq, cigar, rnext, pnext,
                 tlen, seq, qual, tags=[]):
        self.qname = str(qname)
        self.flag = int(flag)
        self.rname = str(rname)
        self.pos = int(pos)
        self.mapq = int(mapq)
        self.cigar = str(cigar)
        self.rnext = str(rnext)
        self.pnext = int(pnext)
        self.tlen = int(tlen)
        self.seq = str(seq)
        self.qual = str(qual)
        self.tags = [str(x) for x in tags]


class Alignment():
    """A sam-format sequence alignment"""
    def __init__(self, head, reads):
        self.head = head
        self.reads = reads


def parse_sam_read(string):
    """Takes a string in SAMfile format and returns a Read object."""
    fields = string.strip().split()
    return Read(fields[0], fields[1], fields[2], fields[3], fields[4],
                fields[5], fields[6], fields[7], fields[8], fields[9],
                fields[10], tags=fields[11:])


def read_sam(samfile):
    """Creates an Alignment object from a correctly formatted SAM file"""
    headlines=[]
    reads=[]
    with open(samfile) as f:
        for line in f:
            if line.startswith("@"):
                headlines.append(line)
            elif line:
                reads.append(parse_sam_read(line))
    return Alignment(head="".join(headlines), reads=reads)


def majority(nucleotides, cutoff=0.5):
    """Given a collection of strings, returns the majority rule consensus among
    them. If there is no majority above the cutoff fraction, returns "N".

    """
    for i in nucleotides:
        if len([x for x in nucleotides if x == i]) >= cutoff:
            return i
    return "N"


def consensus(reads, cutoff=0.5):
    """Returns the consensus sequence of a collection of reads"""
    all_nucleotides = {}
    for read in reads:
        if read.pos == 0:
            raise UnmappedReadError
        else:
            index = read.pos
            for nuc in read.seq:
                all_nucleotides.setdefault(index, [])
                all_nucleotides[index].append(nuc)
                index += 1
    consensus_sequence = ""
    for position in range(min(all_nucleotides), max(all_nucleotides)+1):
        try:
            consensus_sequence += \
                majority(all_nucleotides[position], cutoff=cutoff)
        except KeyError:
            consensus_sequence += 'N'
    return consensus_sequence
