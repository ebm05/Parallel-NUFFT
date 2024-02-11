set term eps size 3.5,3.5
set output "serial_time.eps"
set xlabel "N"
set ylabel "Time [s]"
set logscale y
set logscale x
set key right bottom
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.15
set format y "10^{%L}"
set key font ",11"
plot    'FGG_1000_6_16.dat'  u 1:2 w lp lt -1 pt 4 dt 1 t "naive gridding, Msp=6",\
        'FGG_1000_12_16.dat' u 1:2 w lp lt -1 pt 4 dt 2 t "naive griddning, Msp=12",\
        'FGG_1000_6_16.dat'  u 1:3 w lp lt -1 pt 6 dt 1 t "naive FGG, Msp=6",\
        'FGG_1000_12_16.dat' u 1:3 w lp lt -1 pt 6 dt 2 t "naive FGG, Msp=12",\
        'FGG_1000_6_16.dat'  u 1:4 w lp lt -1 pt 8 dt 1 t "FGG, Msp=6",\
        'FGG_1000_12_16.dat' u 1:4 w lp lt -1 pt 8 dt 2 t "FGG, Msp=12"
