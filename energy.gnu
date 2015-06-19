set terminal png
set output "energy.png"

set title "Energy of an ensemble of simple harmonic oscillators"
set xlabel "time"
set ylabel "energy"
set xrange [0:20]
#set yrange[-1:5]

plot "energy.dat" u 1:2 w l title "kinetic", "energy.dat" u 1:3 w l title "potential", "energy.dat" u 1:($2+$3) w l title "total" 
#plot "energy.dat" u 1:($2+$3) w l title "total" 


