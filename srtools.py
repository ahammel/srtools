class Read():
    """A sam-format sequence read."""
    def __init__(self, qname, flag, rname, pos, mapq, cigar, rnext, pnext,
                 tlen, seq, qual, tags=[]):
        self.qname = qname
        self.flag = flag
        self.rname = rname   
        self.pos = pos   
        self.mapq = mapq   
        self.cigar = cigar   
        self.rnext = rnext   
        self.pnext = pnext  
        self.tlen = tlen   
        self.seq = seq   
        self.qual = qual   
        self.tags=[]
