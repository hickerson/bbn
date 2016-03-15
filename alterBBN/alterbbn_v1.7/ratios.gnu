# set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 500, 350 
# set output 'fillbetween.1.png'
set style data lines
set title "Fill area between two curves" 
set logscale y 10
set xrange [ -0.04 : 0.04 ] noreverse nowriteback
set yrange [ 1e-16 : 1 ] noreverse nowriteback
#plot 'fierz-study.dat' u 1:2:3 w filledcu, '' u 1:2 lt -1 notitle, '' u 1:3 lt -1 notitle
plot 'fierz-study.dat' u 1:2:3 "%lf %lf %lf" w filledcu, \
      '' u 1:2 lt -1 notitle, '' u 1:3 lt -1 notitle 
