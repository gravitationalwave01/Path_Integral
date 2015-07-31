set terminal png
set output "plot.png"

set title "W^a vs simulation time at varying coupling strength to an Andersen thermostat"
set xlabel "Simulation time"
set xrange [0:110]
set ylabel "W^a"
set yrange[0.85:2.2]

#set arrow from 0,1.0397 to 110,1.0397 nohead 
#set arrow from 0,1.5 to 110,1.5 nohead 


plot "<(sed -n '11,15p' tmp.dat)" u 1:3 w l lc rgb "yellow" title "1.0", \
     "<(sed -n '16,20p' tmp.dat)" u 1:3 w l lc rgb "gray" title "0.5", \
     "<(sed -n '21,25p' tmp.dat)" u 1:3 w l lc rgb "blue" title "0.1", \
     "<(sed -n '26,30p' tmp.dat)" u 1:3 w l lc rgb "red" title "0.05", \
     "<(sed -n '31,35p' tmp.dat)" u 1:3 w l lc rgb "green" title "0.01", \
     "<(sed -n '36,40p' tmp.dat)" u 1:3 w l lc rgb "black" title "0.005", \
     "<(sed -n '41,45p' tmp.dat)" u 1:3 w l lc rgb "purple" title "0.001"




