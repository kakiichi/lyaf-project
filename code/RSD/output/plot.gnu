set ter postscript enhanced color
set output "RSD_02.eps"

set size square 0.8,0.8
set pm3d map
#set palette rgbformulae 33,13,10
set palette rgbformulae 22,13,-31
#set palette defined ( 0 '#000090',\
#                      1 '#000fff',\
#                      2 '#0090ff',\
#                      3 '#0fffee',\
#                      4 '#90ff70',\
#                      5 '#ffee00',\
#                      6 '#ff7000',\
#                      7 '#ee0000',\
#                      8 '#7f0000')
set contour base
set cntrparam level incremental 0, 2.0, 30.0
set xtics font "Times-Roman,20"
set ytics font "Times-Roman,20"
set cbtics font "Times-Roman,20"
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

splot 'RSD_02.out' u 3:4:5
