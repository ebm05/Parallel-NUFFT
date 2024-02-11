set term eps size 3.5,3.5
set output "parallel_time.eps"
#set title "Parallel speedup"
set xlabel "Nr. of threads"
set ylabel "Time [s]"
set key right top
set key font ",11"
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.15
set logscale y
plot    'parallel_speedup_10000_100_12.dat' u 1:2 w lp lt -1 pt 6 dt 2 t "Naive, N=1e4",\
        'parallel_speedup_100000_100_12.dat' u 1:2 w lp lt -1 pt 8 dt 2 t "Naive, N=1e5",\
        'parallel_speedup_10000_100_12.dat' u 1:3 w lp lt -1 pt 6 t "Improved, N=1e4",\
        'parallel_speedup_100000_100_12.dat' u 1:3 w lp lt -1 pt 8 t "Improved, N=1e5"