#!/bin/bash
# Convert file $1 to tree $2

	cd /mnt/pool/rhic/4/verbalune/VorticityHic
	./convert -i /mnt/pool/rhic/4/parfenovpeter/ASCII/UrQMD/11.5gev/urqmd_ascii_50540950_${PBS_ARRAYID}.f14 -o forest_11.5/50540950_${PBS_ARRAYID}.root
