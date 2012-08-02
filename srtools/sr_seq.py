import random
import itertools
import sr_sam

class NullSequenceError(ValueError):
    """The exception raised when attempting to illegally manipulate a null
    sequence (i.e., an empty sequence or one composed entirely of N's).
    The GC content of a null sequence is undefined, for example.

    """
    pass


COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def read_fasta(fasta_file):
    """Returns a dictionary a sequence names and values from a fasta-format
    file.

    """
    seq_dict = {}
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].strip()
            elif line.strip():
                seq_dict.setdefault(name, "")
                seq_dict[name] += line.strip()
    return seq_dict


def reverse_complement(sequence):
    rc = []
    seq = list(sequence)
    while seq:
        rc.append(COMPLEMENT[seq.pop()])
    return "".join(rc)


def gc_content(sequence):
    """Returns the fraction of the sequence which consists of GC base pairs.

    """
    base_counts = {x: sequence.count(x) for x in sequence if x in "ACGT"}
    base_counts.setdefault("G", 0)
    base_counts.setdefault("C", 0)
    total = sum(base_counts.values())

    if total == 0:
        raise NullSequenceError

    gc_count = base_counts["G"] + base_counts["C"]
    return gc_count / total


def block_sequence(seq, start, n):
    """Splits a sequence into blocks of size n, prepended by the first 'start'
    items. Helper function for reading_frames.

    """
    blocks = []
    if start != 0:
        blocks.append(seq[:start])
    blocks.extend([seq[i:i + n] for i in range(start, len(seq), n)])
    return  blocks


def reading_frames(sequence):
    """Returns the six possible reading frames of the sequence."""

    frames = []

    for i in range(3):
        frames.append(block_sequence(sequence, i, 3))
        frames.append(block_sequence(reverse_complement(sequence), i, 3))

    return frames


def open_reading_frames(sequence):
    """Returns a list of the ORFs of the sequence in all six translation
    frames.

    """
    stop_codons = ["TAG", "TAA", "TGA"]
    orfs = []
    for frame in reading_frames(sequence):
        starts = [i for i, x in enumerate(frame) if x == "ATG"]
        stops = [i for i, x in enumerate(frame) if x in stop_codons]
        product = itertools.product(starts, stops)
        fr_list = ["".join(frame[a:b]) for a, b in product if a < b]
        orfs.extend(fr_list)
    return orfs


def random_sequence(length):
    """Returns a random nucleotide sequence of the specified length."""
    return "".join([random.choice("ACGT") for i in range(length)])


def randomize_sequence(seq):
    """Randomizes a sequence of nucleotides, preserving N's"""
    nucs = []
    for nucleotide in seq:
        if nucleotide == "N":
            nucs.append("N")
        else:
            nucs.append(random.choice("ACGT"))
    return "".join(nucs)
