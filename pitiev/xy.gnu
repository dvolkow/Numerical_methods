set terminal postscript eps enhanced
set output "xy.eps"
set size ratio -1
unset key
set yrange [-15:30]
set xrange [-15:30]
set mxtics 5
set mytics 5
set xlabel "X" 
set ylabel "Y"
plot 'res.txt' using 1:2 with l ls 1 lt rgb 'blue'
reset

