module load intel-compilers/19.0.3
module load gcc/6.4.0
module load python/3.6.4

export OMP_SCHEDULE="dynamic,16"

Ns=(10000 15849 25119 39811 63096 100000)
thetas=(0.1 0.2 0.3 0.4 0.5)
types=(0 1)
ncrits=(64 128 256)
runs=(1 2 3)

for N in "${Ns[@]}"
do
    for theta in "${thetas[@]}"
    do
	for typ in "${types[@]}"
	do
	    for ncrit in "${ncrits[@]}"
	    do
		for run in "${runs[@]}"
		do
		    echo "${N} ${theta} ${typ} ${ncrit} ${run}"
		    ./main --nparticles ${N} --ncrit=${ncrit} --theta=$theta --label="run${run}" --type=$typ
		done
	    done
	done
    done
done
