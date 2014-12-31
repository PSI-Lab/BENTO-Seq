This directory contains the bootstrap tool to estimate PSI
distributions based on aligned RNA-Seq data (BAM files). The RNAseq
data must be provided as a binary, sorted, and indexed BAM file. 

Before using the package, symlink it using

    python setup.py develop

or install it globally using
 
    python setup.py install

The basic usage of the script is

    bento-seq event_definitions bamfile output_file [OPTIONS]

There are a number of command line options. Run

    bento-seq -h

to get an explanation on their usage.

Two sample BAM files (`TopHat2_chr21.bam` and `STAR_chr21.bam`) and a
sample event specfication file (`events_chr21.tab`) are provided in
the `examples/` folder to test the bootstrap tool. Run

    bento-seq examples/events_chr21.tab examples/STAR_chr21.bam examples/events_chr21_STAR.results -1
    
or

    bento-seq examples/events_chr21.tab examples/TopHat2_chr21.bam examples/events_chr21_TopHat.results -1

to use them.
