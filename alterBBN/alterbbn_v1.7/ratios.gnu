# BBN ratios plotter
mpl_top    = 0.4 #inch  outer top margin, title goes here
mpl_bot    = 0.7 #inch  outer bottom margin, x label goes here
mpl_left   = 0.9 #inch  outer left margin, y label goes here
mpl_right  = 0.1 #inch  outer right margin, y2 label goes here
mpl_height = 1.5 #inch  height of individual plots
mpl_width  = 2.0 #inch  width of individual plots
mpl_dx     = 0.1 #inch  inter-plot horizontal spacing
mpl_dy     = 0.1 #inch  inter-plot vertical spacing
mpl_ny     = 3   #number of rows
mpl_nx     = 1   #number of columns

# calculate full dimensions
xsize = mpl_left+mpl_right+(mpl_width*mpl_nx)+(mpl_nx-1)*mpl_dx
ysize = mpl_top+mpl_bot+(mpl_ny*mpl_height)+(mpl_ny-1)*mpl_dy

# placement functions
#   rows are numbered from bottom to top
bot(n) = (mpl_bot+(n-1)*mpl_height+(n-1)*mpl_dy)/ysize
top(n)  = 1-((mpl_top+(mpl_ny-n)*(mpl_height+mpl_dy))/ysize)
#   columns are numbered from left to right
left(n) = (mpl_left+(n-1)*mpl_width+(n-1)*mpl_dx)/xsize
right(n)  = 1-((mpl_right+(mpl_nx-n)*(mpl_width+mpl_dx))/xsize)

#set terminal postscript eps enhanced color dl 2.0 size xsize,ysize "Helvetica" 28
set encoding iso_8859_1
set tics scale 1.5

#set output 'nxm_plot.eps'

set offsets
set autoscale fix
set size 1,1
set nokey

# define x-axis settings for all subplots
set xrange [-4:4]
set xlabel ''
set format x ''
set xtics pi
set mxtics 4

# start plotting
set multiplot

#-----------------------------------------------
# subplot  1-3
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for third row (top)
set tmargin at screen top(1)
set bmargin at screen bot(1)

set title 'left'

set ylabel "amplitude"
set yrange [-1.5:1.5]
set format y "%-2.1f"
set ytics mirror 1
set mytics 2

#set logscale y 10
set xrange [ -0.1 : 0.1 ] noreverse nowriteback
set yrange [ 1e-1 : 0.25 ] noreverse nowriteback
plot 'fierz-study.dat' u 1:2:3 "%lf %lf %lf" w filledcu, \
      '' u 1:2 lt -1 notitle, '' u 1:3 lt -1 notitle #; #, \
# 'fierz-study.dat' u 1:4:5 "%lf %lf %lf" w filledcu, \
#      '' u 1:5 lt -1 notitle, '' u 1:4 lt -1 notitle; #, \
# 'fierz-study.dat' u 1:6:7 "%lf %lf %lf" w filledcu, \
#      '' u 1:8 lt -1 notitle, '' u 1:7 lt -1 notitle, \
# 'fierz-study.dat' u 1:8:9 "%lf %lf %lf" w filledcu, \
#      '' u 1:8 lt -1 notitle, '' u 1:9 lt -1 notitle, \
# 'fierz-study.dat' u 1:10:11 "%lf %lf %lf" w filledcu, \
#      '' u 1:10 lt -1 notitle, '' u 1:11 lt -1 notitle, \
# 'fierz-study.dat' u 1:12:13 "%lf %lf %lf" w filledcu, \
#      '' u 1:12 lt -1 notitle, '' u 1:13 lt -1 notitle \
#;
unset multiplot
