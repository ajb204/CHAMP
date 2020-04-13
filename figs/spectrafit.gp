reset
set term post eps enh solid color 20
set output 'figs/test.outInt3_a.eps'
set title 'Fitting spectrum: Int3_a'
set label "MZ per point: 4.86" at graph 0.02,graph 0.95
set label "chi2/dof:     33.88" at graph 0.02,graph 0.9
set label "ave Error:    5.82%" at graph 0.02,graph 0.85
set xlabel 'm/z'
set xrange [*:*]
set yrange [-20:110]
plot 'out/test.outInt3_a' u 1:2 ti 'raw' w li lt 1,\
'' u 1:3 ti 'fitted' w li lt 2,\
'' u 1:(($2-$3)-10) ti 'difference' w li lt 3
unset label
reset
set term post eps enh solid color 20
unset key
set output 'figs/test.inpInt3_a.eps'
set title 'Fitted distribution: Int3_a'
set pm3d
unset key
set xlabel '[Sub]y'
set ylabel '[sHSP]x'
set palette defined (0'white',1'blue')
set pm3d map
set size square
splot 'out/test.inpInt3_a' u 2:1:3
unset label
reset
set term post eps enh solid color 20
set output 'figs/test.out.comp.Int3_a.eps'
set title 'Fitted mass spectrum: Int3_a'
set xlabel 'm/z'
set xrange [1400.000000:21500.000000]
set view map
set cblabel 'Oligomer'
set label "Mode:      sim " at graph 0.02,graph 0.95
set label "Zfudge:    0.00 " at graph 0.02,graph 0.9
set label "adduction: 1.00e-03 " at graph 0.02,graph 0.85
set label "MRes:      1000.00 " at graph 0.02,graph 0.8
set label "Resfudge:  1.00e-12 " at graph 0.02,graph 0.75
set label "Zwidth:    1.30e+00 " at graph 0.02,graph 0.7
splot \
'out/test.outInt3_a' u 1:3:(1) ti 'full' w li lt 2,\
'out/test.out.comp.Int3_a' i 0 u 1:($3-0*5):(1) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 1 u 1:($3-1*5):(2) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 2 u 1:($3-2*5):(3) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 3 u 1:($3-3*5):(4) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 4 u 1:($3-4*5):(5) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 5 u 1:($3-5*5):(6) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 6 u 1:($3-6*5):(7) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 7 u 1:($3-7*5):(8) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 8 u 1:($3-8*5):(9) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 9 u 1:($3-9*5):(10) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 10 u 1:($3-10*5):(11) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 11 u 1:($3-11*5):(12) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 12 u 1:($3-12*5):(13) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 13 u 1:($3-13*5):(14) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 14 u 1:($3-14*5):(15) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 15 u 1:($3-15*5):(16) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 16 u 1:($3-16*5):(17) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 17 u 1:($3-17*5):(18) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 18 u 1:($3-18*5):(19) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 19 u 1:($3-19*5):(20) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 20 u 1:($3-20*5):(21) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 21 u 1:($3-21*5):(22) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 22 u 1:($3-22*5):(23) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 23 u 1:($3-23*5):(24) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 24 u 1:($3-24*5):(25) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 25 u 1:($3-25*5):(26) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 26 u 1:($3-26*5):(27) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 27 u 1:($3-27*5):(28) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 28 u 1:($3-28*5):(29) noti w li lc palette,\
 'out/test.out.comp.Int3_a' i 29 u 1:($3-29*5):(30) noti w li lc palette
 