#!/bin/bash
# Draw plots for tree $1 and put to tree $2

	i=$PBS_ARRAYID
	#i=2

	input=/mnt/pool/rhic/4/verbalune/VorticityHic/forest_11.5/50488570_$i.root
	output=/mnt/pool/rhic/4/verbalune/VorticityHic/plots_11.5/plot50488570_$i.root
        cd /mnt/pool/rhic/4/verbalune/VorticityHic
        ./main -i ${input} -o ${output}
	echo ${output}
