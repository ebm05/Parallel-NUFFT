set term eps size 3.5,3.5
set output "serial_speedup.eps"
set xlabel "N"
set ylabel "Speedup"
set key left top
set key font ",11"
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.15
set logscale x
plot    'FGG_1000_6_16.dat'  u 1:6 w lp lt -1 pt 4 dt 1 t "FGG vs naive gridding, Msp=6",\
        'FGG_1000_12_16.dat' u 1:6 w lp lt -1 pt 4 dt 2 t "FGG vs naive griddning, Msp=12",\
        'FGG_1000_6_16.dat'  u 1:5 w lp lt -1 pt 6 dt 1 t "FGG vs naive FGG, Msp=6",\
        'FGG_1000_12_16.dat' u 1:5 w lp lt -1 pt 6 dt 2 t "FGG vs naive FGG, Msp=12"
