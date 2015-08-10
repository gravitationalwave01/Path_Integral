set terminal png
set output "dynamics.png"

set title "Dynamics of a single harmonic oscillator while the potential is varied"
set xlabel "time"
set xrange [0:.15]
#set ylabel "en"
#set yrange[-1:5]

plot "positions.dat" u 1:2 w l title "1", "positions.dat" u 1:3 w l title "2", "positions.dat" u 1:4 w l title "3", "positions.dat" u 1:5 w l title "4"


