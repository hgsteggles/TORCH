#
# Plot all data in directory with plot.py
#
clear
i=0
np=$(nproc)	
for file in tmp/data2D_*.txt
	do python plot.py "-i" "$file" &
	i=`expr $i + 1`
	if [ $i -eq $np ]
		then
		wait
		i=0
	fi
done
