#!/bin/bash
pad=3
formatting="%0"$pad"d"
for i in tmp/*.txt
do {
	[[ "$i" =~ ([0-9]+)\. ]]
	num=${BASH_REMATCH[1]}
	paddednum=`printf $formatting $((10#$num))`
	echo ${i/$num/$paddednum}
	#mv $i ${i/$num/$paddednum}
}
done
