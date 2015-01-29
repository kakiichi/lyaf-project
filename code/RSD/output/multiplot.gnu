set ter postscript enhanced color
set output "RSDs.eps"

set pm3d map
#set palette rgbformulae 33,13,10
set palette rgbformulae 22,13,-31
set contour base
set cntrparam level incremental 2.0, 2.5, 30.0
set xtics font "Times-Roman,20"
set ytics font "Times-Roman,20"
set cbtics font "Times-Roman,20"
unset key
#unset surface
set xlabel "{/Symbol s} [h^{-1}cMpc]" font "Times-Roman,22" offset 0.0,-0.5
set ylabel "{/Symbol p} [h^{-1}cMpc]" font "Times-Roman,22" offset -2.0,0.0
set cblabel "{/Symbol x}_s({/Symbol s},{/Symbol p})" font "Times-Roman,22" offset 2.0,0.0

set xrange [-4:4]
set yrange [-4:4]
#set cbrange [0:30]

do for [i=1:20] {
   set linetype i  lc rgb 'black' lw 1.5
}
set style increment userstyle
#set termoption solid

set multiplot
#####
set origin 0.0,0.0
set size square 0.5,0.5
#set rmargin at screen 0.10

set cbrange [0:100]
set title "inflow model" font "Times-Roman,15"
splot 'RSD_r03p32_rc0p0_s200_v400.dat' u 1:2:3 

#####
set origin 0.0,0.5
set size square 0.5,0.5
#set rmargin at screen 0.10

set cbrange [0:100]
set title "inflow+outflow model" font "Times-Roman,15"
splot 'RSD_r03p32_rc0p0_s200_v400_vout400_rout1p0.dat' u 1:2:3 

#####
set origin 0.5,0.0
set size square 0.5,0.5
unset rmargin
#set lmargin at screen 0.10

set cbrange [0:30]
set title "inflow+local ionization model " font "Times-Roman,15"
splot 'RSD_r03p32_rc0p5_s200_v400.dat' u 1:2:3 

#####
set origin 0.5,0.5
set size square 0.5,0.5
#set lmargin at screen 0.10

set cbrange [0:20]
set title "inflow+outflow+local ionization model" font "Times-Roman,15"
splot 'RSD_r03p32_rc0p5_s200_v400_vout400_rout1p0.dat' u 1:2:3 


