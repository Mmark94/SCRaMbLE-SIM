#!/bin/bash

# Test SCRaMbLE-SIM and simulate long sequencing reads.

# Test SCRaMbLE-SIM on one synthetic chromosome.
python ../simulate_reads.py --reference IXR_BACnewseq.fa -ID synIXR --synthetic True -coverage 20 -reads_len 8 -sigma 3 -circular True

# Test SCRaMbLE-SIM on two synthetic chromosomes.
#python ../simulate_reads.py --reference synII_III.fasta -ID synII_III --synthetic True -coverage 20 -reads_len 8 -sigma 3

# Test SCRaMbLE-SIM on WT chr09. This chromosome does not have loxPsym sites.
# Note: --synthetic False will evaluate to True! removing --synthetic will evaluate to False
#python ../simulate_reads.py --reference BY4741_chr09.fa -ID chr09 -coverage 20 -reads_len 8 -sigma 3 --segment_size 1000 --min_segment_size 100

# Test SCRaMbLE-SIM on two WT chromosomes. These chromosomes do not have loxPsym sites.
#python ../simulate_reads.py --reference BY4741_chr02_03.fasta -ID chr02_03 -coverage 20 -reads_len 8 -sigma 3 --segment_size 1000 --min_segment_size 100
