#!/bin/bash
# Draw plots for tree $1 and put to tree $2

        cd /mnt/pool/rhic/4/verbalune/VorticityHic
        ./main -i /mnt/pool/rhic/4/verbalune/VorticityHic/forest_11.5/50488570_${PBS_ARRAYID}.root -o plots_11.5/plot50488570_${PBS_ARRAYID}.root

