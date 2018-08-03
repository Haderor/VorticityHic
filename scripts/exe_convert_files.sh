
#!/bin/bash
# Convert file $1 to tree $2
	i=${PBS_ARRAYID}
	#i=41
	cd /mnt/pool/rhic/4/verbalune/VorticityHic
	./convert -i /mnt/pool/rhic/4/parfenovpeter/ASCII/UrQMD/11.5gev/urqmd_ascii_50540950_${i}.f14 -o forest_11.5/50540950_${i}.root
