
set terminal png
set output "cini.png"
#set lmargin 5
#set rmargin 5
#set grid

### start multiplot 1x2 layout

set multiplot layout 2,1 rowsfirst
set size 1, 0.5
# grid
set grid




# --- Solition

set autoscale
set yrange [-4:4]
set title 'Soliton'
set xlabel 'x'
set ylabel 'u'
plot "5000.txt" u 1:2 title "u(t_i)" w l lc rgb "#0000FF",\
     "5001.txt" u 1:2 title "u(t)"   w l lc rgb "#04B404",\
     "5002.txt" u 1:2 title "u(t_f)" w l lc rgb "#FF000F",\

#--- Energy

set yrange [-4:17]
set title ' '
set xlabel 't'
set ylabel ' '
plot "energy.txt" u 1:2 title "E" w l lc rgb "#0000FF",\
     "energy.txt" u 1:3 title "P"  w l lc rgb "#FF0000"
#     "energy.txt" u 1:4  w l lc rgb "#FF0000"



unset multiplot

pause mouse
reset

