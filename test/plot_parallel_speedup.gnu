set term eps size 3.5,3.5
set output "parallel_speedup.eps"
#set title "Parallel speedup"
set xlabel "Nr. of threads"
set ylabel "Speedup"
set key left top
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.15
set key font ",11"
plot    'parallel_speedup_10000_100_12.dat' u 1:4 w lp lt -1 pt 6 dt 2 t "Naive, N=1e4",\
        'parallel_speedup_100000_100_12.dat' u 1:4 w lp lt -1 pt 8 dt 2 t "Naive, N=1e5",\
        'parallel_speedup_10000_100_12.dat' u 1:5 w lp lt -1 pt 6 t "Improved, N=1e4",\
        'parallel_speedup_100000_100_12.dat' u 1:5 w lp lt -1 pt 8 t "Improved, N=1e5",\
        'parallel_speedup_100000_100_12.dat' u 1:1 w l lt -1 dt 3 t "Linear"