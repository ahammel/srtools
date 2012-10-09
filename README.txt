=======
SRtools
=======

Note: this README is a work in progress. SRtools contains many features not documented here. Please have a look at the code to find these goodies.

SRtools provides parsers and classes for easily manipulating short read genetic data, including sam files, and a number of related data types. A typical application would be to filter aligned data for interesting reads. For example, the following script will print all of the reads in an alignment which have an insertion or deletion relative to the reference frame::

    from srtools import SamAlignment 

    alignment = SamAlignment("some_data.sam")

    def has_indel(read):
        operators = [o for (i, o) in read.cigar]
        return "I" in operators or "D" in operators

    for read in alignment:
        if has_indel(read):
            print(read)

Note that instances of the SamAlignment class are iterable objects which return members of the srtools Read datatype. The same effect could be accomplished with the filter_reads method::

    for read in alignment.filter_reads(has_indel):
        print(read)

A slightly more complicated example prints all the reads in the samfile which are part of an overlapping set of reads which overlaps with a gene::

    from srtools import SamAlignment, expressed_loci
    from srtools.gff import read_gff, in_features

    alignment = SamAlignment("some_data.sam")
    features = read_gff("TAIR9_genes.gff")

    chromosomes = [("1", "Chr1"), 
                   ("2", "Chr2"), 
                   ("3", "Chr3"),
                   ("4", "Chr4"),
                   ("5", "Chr5"),
                   ("Mt", "ChrM"),
                   ("Pt", "ChrC")]

    for chr_name, chr_number in chromosomes:
        alignment.rewind()

        genes = features.filter_features(lambda x: x.sequence == chr_name
                                         and x.f_type == "gene")

        reads = alignment.filter_reads(lambda x: x.rname == chr_number)
        loci = expressed_loci(reads)

        for locus in loci:
            if in_features(locus, genes):
                for read in locus:
                    print(read)

Since SAM files can run to hundreds of gigabytes, srtools does not attempt to keep them in memory. Alignments are generator objects and the ``rewind`` method restarts the generator.

Installation
===========

The latest development can be found at <https://github.com/ahammel/srtools>, while the latest stable version (with no testing scripts) lives at <https://github.com/ahammel/srtools-sacred> (you probably want the latter if you want to use srtools but not contribute). To install just do this::

    git clone https://github.com/ahammel/srtools-sacred
    cd srtools-sacred
    sudo python3 setup.py install

And you're done.

How to Contribute
=================

Think SRtools is missing  a feature? By all means, clone the testing module and code it up! I'm desperately seeking contributors. The only groundrules for contributing are:

1. All pulls to the sacred repo are at Alex's discression.
2. Nothing gets into sacred without passing a unit-test.
3. *Nothing gets into sacred without passing a unit-test.*

All unit tests should be written so as to run with either the nose or py-test modules. If you don't have experience with either of these (and can't be bothered to learn), please flag your untested code with an issue in your branch and somebody will get around to writing a test for it. 

See the LICENSE file for licensing information. Long story short: this software
is BSD-licensed.
