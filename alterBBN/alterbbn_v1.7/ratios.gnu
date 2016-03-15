# set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 500, 350 
# set output 'fillbetween.1.png'
set style data lines
set title "Fill area between two curves" 
set logscale y 10
set xrange [ -0.04 : 0.04 ] noreverse nowriteback
set yrange [ 1e-16 : 1 ] noreverse nowriteback
plot 'fierz-study.dat' u 1:2:3 "%lf %lf %lf" w filledcu, \
      '' u 1:2 lt -1 notitle, '' u 1:3 lt -1 notitle, \
 'fierz-study.dat' u 1:4:5 "%lf %lf %lf" w filledcu, \
      '' u 1:4 lt -1 notitle, '' u 1:5 lt -1 notitle, \
 'fierz-study.dat' u 1:6:7 "%lf %lf %lf" w filledcu, \
      '' u 1:8 lt -1 notitle, '' u 1:7 lt -1 notitle, \
 'fierz-study.dat' u 1:8:9 "%lf %lf %lf" w filledcu, \
      '' u 1:8 lt -1 notitle, '' u 1:9 lt -1 notitle, \
 'fierz-study.dat' u 1:10:11 "%lf %lf %lf" w filledcu, \
      '' u 1:10 lt -1 notitle, '' u 1:11 lt -1 notitle, \
 'fierz-study.dat' u 1:12:13 "%lf %lf %lf" w filledcu, \
      '' u 1:12 lt -1 notitle, '' u 1:13 lt -1 notitle
