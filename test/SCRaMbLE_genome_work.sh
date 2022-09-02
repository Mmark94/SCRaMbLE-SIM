#!/bin/bash

# Test SCRaMbLE-SIM and simulate SCRaMbLE of a synthetic chromosome.

# Test SCRaMbLE-SIM on one synthetic chromosome.
python ../SCRaMbLE_DNA_chromosomes.py --filename IXR_BACnewseq.fa -ID synIXR --synthetic True --Number_SCRaMbLE_events 15

# Test SCRaMbLE-SIM on two synthetic chromosomes.
#python ../SCRaMbLE_DNA_chromosomes.py --filename synII_III.fasta -ID synII_III --synthetic True --Number_SCRaMbLE_events 15

# Test SCRaMbLE-SIM on WT chr09. This chromosome does not have loxPsym sites.
# Note: --synthetic False will evaluate to True! removing --synthetic will evaluate to False
#python ../SCRaMbLE_DNA_chromosomes.py --filename BY4741_chr09.fa -ID chr09 --Number_SCRaMbLE_events 15 --segment_size 1000 --min_segment_size 100

# Test SCRaMbLE-SIM on two WT chromosomes. These chromosomes do not have loxPsym sites.
#python ../SCRaMbLE_DNA_chromosomes.py --filename BY4741_chr02_03.fasta -ID chr02_03 --Number_SCRaMbLE_events 15 --segment_size 1000 --min_segment_size 100
