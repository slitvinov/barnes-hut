set term svg
set output "box.svg"
set size sq
unset key
unset border
unset xtics
unset ytics
set xrange [-1:1]
set yrange [-1:1]

plot \
     "<./info -t 0.2 -p 0.2 0.5 -o twig < data/points" w l lw 2, \
     "<./info -t 0.2 -p 0.2 0.5 -o leaf < data/points" w l lw 2, \
     "<./info -t 0.2 -p 0.2 0.5 -o force < data/points" w l lw 1.5, \
     "data/points" w p pt 7 lt 1, \
     "<echo 0.2 0.5" w p pt 7 ps 1
