# <<<<<<< how to plot binned data in gnuplot >>>>>

set terminal png
set output "bin.png"

#set title "Initial displacements of all HOs in an ensemble"
#set xlabel "Displacement from equilibrium"
#set ylabel "occurances"
#set xrange[6:8]

binwidth = .01
bin(x,width)=width*floor(x/width)
ibeta = 1.5
f(x) = 3333 * (1/(2*pi)/ibeta)**(1/2) * exp(-1/ibeta*x**2/2)

plot "xdist1.dat" u (bin($1,binwidth)):(1.0) smooth freq with boxes title "pos1", "xdist2.dat" u (bin($1,binwidth)):(1.0) smooth freq with boxes title "pos2",\
     "xdist3.dat" u (bin($1,binwidth)):(1.0) smooth freq with boxes title "pos3", "xdist4.dat" u (bin($1,binwidth)):(1.0) smooth freq with boxes title "pos4"

#plot "rwork_t0.dat" u (bin(-$1,binwidth)):(1.0) smooth freq with boxes notitle, "fwork_t0.dat" u (bin($1,binwidth)):(1.0) smooth freq with boxes notitle
#plot "work.dat" u (bin(-$1,binwidth)):(1.0) smooth freq with boxes notitle
#plot "test.dat" u (bin(-$1,binwidth)):(1.0) smooth freq with boxes notitle, f(x) w l

