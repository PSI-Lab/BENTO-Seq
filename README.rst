This directory contains the bootstrap tool to estimate PSI distributions based on aligned RNA-Seq data (BAM files).
The RNAseq data must be provided as a binary, sorted, and indexed BAM file.Before using the package, symlink it using::

    python setup.py develop

or install it globally using::

    python setup.py install

The basic usage of the script is::

    bento-seq event_definitions bamfile output_file [OPTIONS]

Two sample BAM files (``TopHat2_chr21.bam`` and ``STAR_chr21.bam``) and a sample event specfication file
(``events_chr21.tab``) are provided in the ``examples/`` folder to test the bootstrap tool.

We are providing a number of genomewide sets of alternative splicing events for human (``hg19``, ``hg39``) and mouse
(``mm9``, ``mm10``) which are downloaded automatically when first used. To use them, simply use the name of the genome
assembly in place of ``event_definitions``, e.g.::

    bento-seq hg19 examples/STAR_chr21.bam examples_chr21_STAR.results

If you already have a a set of alternative splicing events that you would like to use, provide the path to the event
definitions instead::

    bento-seq examples/events_chr21.tab examples/STAR_chr21.bam examples/events_chr21_STAR.results

BENTO-Seq supports bam-files created with STAR or TopHat and BENTO-Seq will automatically detect which kind of file you
provided. E.g. to run the same task with a TopHat bam-file::

    bento-seq examples/events_chr21.tab examples/TopHat2_chr21.bam examples/events_chr21_TopHat.results

There are a number of other command line options. Run::

    bento-seq -h

to get an explanation on their usage.

