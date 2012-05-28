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
