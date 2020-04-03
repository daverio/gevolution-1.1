reset
set term postscript eps enhanced clip color font "arial,7" size 7,3
set output '_PLOTS.eps'
set format xy '%g'
set autoscale xy
data_path = './output/'
set multiplot layout 2,3
set key outside top left font "sans,7" spacing '1.26'
model1 = 'fR'
model2 = 'lcdm'
pk_num = '000'
style = 'smooth bezier'
zero_plot_settings = "w l lt 3 dt ' -' lw 1 lc 'black' notitle"
#####################
set logscale x
var = 'T00'
file1 = data_path.model1.'_pk'.pk_num.'_'.var.'.dat'
file2 = data_path.model2.'_pk'.pk_num.'_'.var.'.dat'
pl1 = '<(join '.file1.' '.file2.')'
plot 0 @zero_plot_settings, pl1 u 1:(($2-$6)/$6) @style w l lc 'blue' lw 3 title '{/Symbol d}P_{/Symbol r}/P_{/Symbol r}'
unset logscale
#####################
set logscale xy
var = 'phi'
use = 'using 1:($2*$1**2)'
file1 = data_path.model1.'_pk'.pk_num.'_'.var.'.dat'
file2 = data_path.model2.'_pk'.pk_num.'_'.var.'.dat'
plot file1 @use @style w l lc 'blue' lw 3 title '{/Symbol F}_{f(R)}', file2 @use @style w l lc 'red' lw 3 title '{/Symbol F}_{{/Symbol L}CDM}'
unset logscale
#####################
set logscale x
var = 'B'
use = 'using 1:($2*$1**2)'
file1 = data_path.model1.'_pk'.pk_num.'_'.var.'.dat'
file2 = data_path.model2.'_pk'.pk_num.'_'.var.'.dat'
pl1 = '<(join '.file1.' '.file2.')'
plot 0 @zero_plot_settings, pl1 u 1:(($2-$6)/$6) @style w l lc 'blue' lw 3 title '{/Symbol d}P_B/P_B'
unset logscale
#####################
set logscale xy
var = 'deltaR'
file1 = data_path.model1.'_pk'.pk_num.'_'.var.'.dat'
file2 = data_path.model2.'_pk'.pk_num.'_'.var.'.dat'
file3 = data_path.model1.'_pk'.pk_num.'_zeta.dat'
plot\
 file1 @style w l lc 'blue' lw 3 title '{/Symbol d}R_{f(R)}'\
, file2 @style w l lc 'red' lw 3 title '{/Symbol d}R_{{/Symbol L}CDM}'\
, file3 @style w l lc 'green' lw 3 title '{/Symbol z}_{f(R)}'
unset logscale
#####################
set logscale xy
var = 'chi'
file1 = data_path.model1.'_pk'.pk_num.'_'.var.'.dat'
file2 = data_path.model2.'_pk'.pk_num.'_'.var.'.dat'
file3 = data_path.model1.'_pk'.pk_num.'_xi.dat'
plot\
 file1 @style w l lc 'blue' lw 3 title '{/Symbol c}_{f(R)}'\
, file2 @style w l lc 'red' lw 3 title '{/Symbol c}_{{/Symbol L}CDM}'\
, file3 @style w l lc 'green' lw 3 lt 3 dt ' .' title '{/Symbol x}_{f(R)}'
unset logscale
#####################
set logscale x
set xrange[:1.]
plot 0 @zero_plot_settings \
, '<(join '.data_path.model1.'_background.dat '.data_path.model2.'_background.dat)' u 3:(($4-$10)/$10) w l lc 'red' lw 3 title '{/Symbol d}H/H'\
, '<(join '.data_path.model1.'_background.dat '.data_path.model2.'_background.dat)' u 3:(($5-$11)/$11) w l lc 'blue' lw 3 lt 3 dt '-' title '{/Symbol d}R/R'
unset logscale
#####################
unset multiplot
