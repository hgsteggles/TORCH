#set terminal wxt size 1600,900
#set term postscript eps enhanced color size 8,4.5 solid "Helvetic" 14
#enhanced color solid size 8cm,6cm
#set term png
#set output 'plot.eps'

set xlabel 'time / t_{rec}'
set size ratio 1
set xrange [*:*]
set yrange [*:*]
set key

set style line 1 lt 1 lw 1 pt 1 ps 0.5 lc rgb "black"
set style line 2 lt 1 lw 1 pt 1 ps 0.5 lc rgb "blue"
set style line 3 lt 1 lw 1 pt 1 ps 0.5 lc rgb "red"
set style line 4 lt 1 lw 1 pt 1 ps 0.5 lc rgb "purple"

set multiplot layout 1,2 rowsfirst title "Comparing Integration Schemes: 100 Timesteps"
set ylabel 'IF / pc'
plot "IF.dat" u 1:3 title "Analytical" w l ls 1 ,\
	"IF.dat" u 1:2 title "Implicit" w l ls 2
set ylabel 'IF/Rs'
set yrange [0:0.2]
plot 1 title "Analytical" w l ls 1 ,\
	"IF.dat" u 1:4 title "Implicit" w l ls 2
	
unset multiplot
set nologscale x

set term wxt
