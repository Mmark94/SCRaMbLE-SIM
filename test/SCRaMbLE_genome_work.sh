#!/bin/bash

# Test SCRaMbLE-SIM and simulate SCRaMbLE of a synthetic chromosome.

filename="IXR_BACnewseq.fa"
ID="synIXR"
#path_SCRaMbLE-SIM=""
python ../SCRaMbLE_DNA_chromosomes.py --filename "$filename" -ID "$ID" --Number_SCRaMbLE_events 15
