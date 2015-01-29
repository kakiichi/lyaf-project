set ter postscript enhanced color
set output "RSD.eps"

set size square 0.8,0.8
set xtics font "Times-Roman,20"
set ytics font "Times-Roman,20"
unset key
#unset surface
set xlabel "{/Symbol s} [h^{-1}cMpc]" font "Times-Roman,25" offset 0.0,-0.5
set ylabel "{/Symbol p} [h^{-1}cMpc]" font "Times-Roman,25" offset -2.0,0.0
set cblabel "{/Symbol x}_s({/Symbol s},{/Symbol p})" font "Times-Roman,25" offset 2.0,0.0

set xrange [-4:4]
set yrange [-4:4]
#set cbrange [0:30]

do for [i=1:20] {
   set linetype i  lc rgb 'black' lw 1.5
}
set style increment userstyle
#set termoption solid

splot 'RSD.dat' u 1:2:3 
