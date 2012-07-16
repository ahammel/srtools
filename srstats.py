import srtools
import sys
from collections import Counter

def summary_statistics(reads):
    """Returns a dictionary of summary statistics of the reads. The keys are:

        "rnames":           a Counter of the rnames of the reads
        "flags":            a Counter of the bitflags of the reads
        "cigars":           a Counter of the string representations of the 
                            cigars of the reads
        "gc":               the average GC content of the sequences.
        "read_count":       the number of sam reads
        "hash":             a Counter of the qnames of the reads
        
        "rnames_mapped":    a summary of the number of reads with each rname 

    """ 
    summary = {"rnames":Counter(),
               "flags": Counter(),
               "cigars": Counter(),
               "gc": 0,
               "read_count": 0,
               "hashes": Counter()}

    for read in reads:
        summary["rnames"][read.rname] += 1
        summary["flags"][read.flag] += 1
        summary["cigars"][srtools.print_cigar(read.cigar)] += 1
        summary["gc"] += srtools.gc_content(read.seq)
        summary["read_count"] += 1
        summary["hashes"][read.qname] += 1

    summary["gc"] /= summary["read_count"]
    return summary


#Text colouring functions for pretty-printing.

def cyan(string):
    """Returns the input string wrapped in terminal codes for cyan text."""
    return "".join(["\033[96m", string, "\033[0m"])


def green(string):
    """Returns the input string wrapped in terminal codes for green text."""
    return "".join(["\033[92m", string, "\033[0m"])


def print_summary_statistics(input_file, output_file=sys.stdout):
    """Pretty-prints the summary statistics of a sam file to the specified
    output file, or to stdout if not output_file is selected.

    """
    alignment = srtools.read_sam(input_file)
    stats = summary_statistics(alignment.reads)

    f = open(output_file, "w")

    #Header
    head_string = "".join([cyan("Summary of Sam File "),
                           cyan(input_file) ])

    print("\n" + head_string + "\n", file=f)

    #Chromosome mapping
    print(green("Total Reads Mapped to Chromosomes") + "...", file=f)

    for rname in sorted(stats["rnames"]):
        freq = stats["rnames"][rname]
        rel_freq = freq / stats["read_count"]
        p_string = "{:8s}{:9d}{:>8.1%}".format(rname, freq, rel_freq)
        print(p_string, file=f)

    #Bitflags
    print("\n" + green("Total Counts of Bit Flags") + "...", file=f)
    flags = sorted(stats["flags"], key=lambda x: str(x))
    for flag in flags:
        freq = stats["flags"][flag] 
        rel_freq = freq / stats["read_count"]
        p_string = "{:<8}{:9d}{:>8.1%}".format(flag, freq, rel_freq)
        print(p_string, file=f)
   
    #Cigars
    print("\n" + green("Total Counts of Cigar Strings") + "...", file=f)
    cigars = sorted(stats["cigars"], key=lambda x: stats["cigars"][x])
    cigars.reverse()
    for cigar in cigars[:20]:
        freq = stats["cigars"][cigar]
        rel_freq = freq / stats["read_count"]
        p_string = "{:10s}{:6d}{:>8.1%}".format(cigar, freq, rel_freq)
        print(p_string, file=f)

    #GCc content
    print("\n" + green("Average % GC")+ "...", file=f)
    print(stats["gc"], file=f)

    #Total reads (pair = one read)
    print("\n" + green("Total Reads (Mate Pairs and Orphans)"), file=f)
    print(len(stats["hashes"]), file=f)

    #Total reads (pair = two reads)
    print("\n" + green("Total Reads"), file=f)
    print(stats["read_count"], file=f)
    print("", file=f)

    f.close()
