#!/opt/local/bin/gnuplot
set terminal postscript eps enhanced size 3in,3in
set output 'file.eps'
set autoscale
set key font ",20"
set xtics font "Verdana,12" 
set ytics font "Verdana,12" 
set xlabel font "Verdana,20" 
set ylabel font "Verdana,20" 
#set title "Stress_all"
set xlabel "Локальный размер сетки, мм"
set ylabel "Напряжение по Мизесу, МПа" offset 2.3,0,0
set xrange [2.5:0]
set yrange [0:250]
plot "snow.txt"  using 1:2 with linespoints title "a" lw 3, "snow.txt" using 1:3 with linespoints title "b" lw 3, "snow.txt" using 1:4 with linespoints title "c" lw 3, "snow.txt" using 1:5 with linespoints title "d" lw 3
set output "1.png"
set xrange [2.5:0]
set yrange [9.51:9.52]
set ylabel "Displacement, 10e-5 m"
plot "1.txt"  using 1:2 with linespoints title "1" lw 3
set output "2.png"
set xrange [2.5:0]
set yrange [9.737:9.747]
set ylabel "Displacement, 10e-5 m"
plot "1.txt"  using 1:3 with linespoints title "1" lw 3
set output "3.png"
set xrange [2.5:0]
set yrange [9.661:9.671]
set ylabel "Displacement, 10e-5 m"
plot "1.txt"  using 1:4 with linespoints title "1" lw 3
set output "4.png"
set xrange [2.5:0]
set yrange [9.866:9.876]
set ylabel "Displacement, 10e-5 m"
plot "1.txt"  using 1:5 with linespoints title "1" lw 3
set output "5.png"
set xrange [6:0]
set yrange [9.42:9.48]
set ylabel "Displacement, 10e-5 m"
plot "2.txt"  using 1:2 with linespoints title "1" lw 3
set output "6.png"
set xrange [6:0]
set yrange [74.9:75.4]
set ylabel "Mizes stress, MPa"
plot "2.txt"  using 1:3 with linespoints title "1" lw 3
