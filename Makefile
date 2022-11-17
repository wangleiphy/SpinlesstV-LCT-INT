main.sif: main.def
	sudo singularity build main.sif main.def 

run: input/params.in input/mylattices.xml
	singularity exec main.sif mpirun -np 4 /project/SpinlesstV-LCT-INT-PBC/bin/main -a 10 -T 120 $(pwd)/input/sif_params.in

launch: main
	sudo singularity exec --writable main sh

main: main.sif
	sudo singularity build --sandbox main main.sif

clear:
	sudo rm -rf main