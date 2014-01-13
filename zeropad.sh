#!/bin/bash
for i in tmp/*.txt
do {
	num=`expr match "$i" '.*\+_\+\([0-9]\+\).*'`
	paddednum=`printf "%03d" $((10#$num))`
	mv $i ${i/$num/$paddednum}
}
done
