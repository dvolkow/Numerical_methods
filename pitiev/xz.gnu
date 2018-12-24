set terminal postscript eps enhanced
set output "xz.eps"
set size ratio -1
unset key
set yrange [0:45]
set xrange [-15:30]
set mxtics 5
set mytics 5
set xlabel "X" 
set ylabel "Z"
plot 'res.txt' using 1:3 with l ls 1 lt rgb 'blue'
reset

