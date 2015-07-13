set terminal png
set output "plot.png"

#set title "Dynamics of a single harmonic oscillator while the potential is varied"
#set xlabel "time"
#set xrange [0:.1]
#set ylabel "en"
#set yrange[-1:6]

f(x) = a*x
fit f(x) "toplot.dat" u 1:2 via a

plot "toplot.dat" u 1:2 notitle, f(x)


