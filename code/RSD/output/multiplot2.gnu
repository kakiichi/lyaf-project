set ter postscript enhanced color
set output "RSDs_1-4-5.eps"

set pm3d map
#set palette rgbformulae 33,13,10
set palette rgbformulae 22,13,-31
set contour base
#set cntrparam level incremental 2.0, 2.5, 30.0
set cntrparam level discrete 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75
set xtics font "Times-Roman,20"
set ytics font "Times-Roman,20"
set cbtics font "Times-Roman,20"
unset key
#unset surface
set xlabel "{/Symbol s} [h^{-1}cMpc]" font "Times-Roman,22" offset 0.0,-0.5
set ylabel "{/Symbol p} [h^{-1}cMpc]" font "Times-Roman,22" offset -2.0,0.0
set cblabel "log_{10} {/Symbol x}_s({/Symbol s},{/Symbol p})" font "Times-Roman,22" offset 2.0,0.0

set xrange [-4:4]
set yrange [-4:4]
#set cbrange [0:30]

do for [i=1:20] {
   set linetype i  lc rgb 'black' lw 1.5
}
set style increment userstyle
#set termoption solid

set size 1.4,0.45
set multiplot layout 1,3 
#set bmargin 5
#set multiplot
#####
set origin 0.0,0.0
set size square 0.5,0.5
#set rmargin at screen 0.10

#set cbrange [-1:2]
set cbrange [-0.6:1.6]
set title "v_{in}=143km/s (fid.)" font "Times-Roman,15"
splot 'RSD_01.out' u 3:4:(log10($5)) 

#####
set origin 0.45,0.0
set size square 0.5,0.5
#set rmargin at screen 0.10

#set cbrange [-1:2]
set cbrange [-0.6:1.6]
set title "v_{in}=100km/s" font "Times-Roman,15"
splot 'RSD_04.out' u 3:4:(log10($5))

#####
set origin 0.90,0.0
set size square 0.5,0.5
#unset rmargin
#set lmargin at screen 0.10

set cbrange [-0.6:1.6]
set title "v_{in}=200km/s" font "Times-Roman,15"
splot 'RSD_05.out' u 3:4:(log10($5))
