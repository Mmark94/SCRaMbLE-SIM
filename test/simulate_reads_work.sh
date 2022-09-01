#!/bin/bash

# Test SCRaMbLE-SIM and simulate long sequencing reads.

filename="IXR_BACnewseq.fa"
ID="synIXR"
#path_SCRaMbLE-SIM=""
python ../simulate_reads.py --reference "$filename" -ID "$ID" -coverage 20 -reads_len 8 -sigma 3 -circular True
