# <<<<<<< how to plot binned data in gnuplot >>>>>

set terminal png
set output "bin.png"

set title "Initial displacements of all HOs in an ensemble"
set xlabel "Displacement from equilibrium"
set ylabel "occurances"
#set xrange[-10:10]

binwidth = .02
bin(x,width)=width*floor(x/width)

#plot "dist.dat" u (bin($2,binwidth)):(1.0) smooth freq with boxes title "velocities"#, "dist.dat" u (bin($3,binwidth)):(1.0) smooth freq with boxes title "positions"
plot "work.dat" u (bin($1,binwidth)):(1.0) smooth freq with boxes notitle

