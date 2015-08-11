set terminal png
set output "energy.png"

set title "Energy of an ensemble of simple harmonic oscillators"
set xlabel "time"
set ylabel "energy"
#set xrange [20:40]
#set yrange[-1:5]

#plot "energy2.dat" u 1:2 w l title "kinetic", "energy2.dat" u 1:3 w l title "ext pot", "energy2.dat" u 1:4 w l title "bead pot"
plot "energy2.dat" u 1:($2+$3+$4) w l title "total" 


