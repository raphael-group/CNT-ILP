#!/bin/bash
c=1
for s in {1..3}
do
	for n in {5,10,15,20,30,40}
	do
		for k in {4,6,8}
		do
			for u in {1,2,3}
			do
				for d in 0.2 0.4
				do
					../utils/generateC.py -c $c -s $s -n $n -k $k -u $u -d $d -o ./
				done
			done
		done
	done
done
