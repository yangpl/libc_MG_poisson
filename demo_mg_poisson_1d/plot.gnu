#set terminal epslatex size 3.5,2.62 color colortext
#set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 500, 350 
set terminal postscript eps enhanced color solid font 'Helvetica' fontscale 1.5 
set size ratio 0.75 #1,1
#set border 3
#set multiplot layout 2,2 
#set output 'misfit.eps'

set title "Multigrid convergence"
set output 'fcost.eps'
set key top right
set logscale y 10
set xlabel '# MG cycles' # font "Helvetica,16"
set ylabel 'Residual' # font "Helvetica,16"
set grid
show grid
plot 'iterate.txt' title 'V(1,1), MG' with linespoints linestyle 1

set title "Multigrid convergence"
set output 'convergence_ratio.eps'
set key top right
set logscale y 10
unset logscale y
set xlabel '# iteration' # font "Helvetica,16"
set ylabel 'Data misfit' # font "Helvetica,16"
set grid
show grid
plot 'iterate.txt' using 1:3  title 'Gauss-Seidel' with linespoints linestyle 1


#pause -1 "Hit any key to continue"
