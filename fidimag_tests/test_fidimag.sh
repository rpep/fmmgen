for L in {5..60}
do
	for order in 5 #3 7
	do
		for theta in 0.5 0.7 0.9
		do
			python3 test_fmm.py $L $order $theta 128
		done
	done
done

