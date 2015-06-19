set terminal png
set output "fpot.png"

#set title "Energy of an ensemble of simple harmonic oscillators"
#set xlabel "time"
#set ylabel "energy"
#set xrange [0:20]
#set yrange[-1:5]

plot "fpot.dat" u 1:2 w l notitle


