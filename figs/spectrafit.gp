reset
set term post eps enh solid color 20
set output 'figs/test.outInt3_b.eps'
set title 'Fitting spectrum: Int3_b'
set label "MZ per point: 4.86" at graph 0.02,graph 0.95
set label "chi2/dof:     3.09" at graph 0.02,graph 0.9
set label "ave Error:    1.76%" at graph 0.02,graph 0.85
set xlabel 'm/z'
set xrange [*:*]
set yrange [-20:110]
plot 'out/test.outInt3_b' u 1:2 ti 'raw' w li lt 1,\
'' u 1:3 ti 'fitted' w li lt 2,\
'' u 1:(($2-$3)-10) ti 'difference' w li lt 3
unset label
reset
set term post eps enh solid color 20
unset key
set output 'figs/test.inpInt3_b.eps'
set title 'Fitted distribution: Int3_b'
set pm3d
unset key
set xlabel '[Sub]y'
set ylabel '[sHSP]x'
set palette defined (0'white',1'blue')
set pm3d map
set size square
splot 'out/test.inpInt3_b' u 2:1:3
unset label
reset
set term post eps enh solid color 20
set output 'figs/test.out.comp.Int3_b.eps'
set title 'Fitted mass spectrum: Int3_b'
set xlabel 'm/z'
set xrange [1400.000000:21500.000000]
set view map
set cblabel 'Oligomer'
set label "Mode:      fitty " at graph 0.02,graph 0.95
set label "Zfudge:    -1.85 (1)" at graph 0.02,graph 0.9
set label "adduction: 1.77e-02 (1)" at graph 0.02,graph 0.85
set label "MRes:      1731.98 (1)" at graph 0.02,graph 0.8
set label "Resfudge:  1.97e-05 (1)" at graph 0.02,graph 0.75
set label "Zwidth:    5.83e+00 (1)" at graph 0.02,graph 0.7
set label "shspMass:  2.02e+04 (0)" at graph 0.02,graph 0.65
splot \
'out/test.outInt3_b' u 1:3:(1) ti 'full' w li lt 2,\
'out/test.out.comp.Int3_b' i 0 u 1:($3-0*5):(1) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 1 u 1:($3-1*5):(2) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 2 u 1:($3-2*5):(3) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 3 u 1:($3-3*5):(4) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 4 u 1:($3-4*5):(5) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 5 u 1:($3-5*5):(6) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 6 u 1:($3-6*5):(7) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 7 u 1:($3-7*5):(8) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 8 u 1:($3-8*5):(9) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 9 u 1:($3-9*5):(10) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 10 u 1:($3-10*5):(11) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 11 u 1:($3-11*5):(12) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 12 u 1:($3-12*5):(13) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 13 u 1:($3-13*5):(14) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 14 u 1:($3-14*5):(15) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 15 u 1:($3-15*5):(16) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 16 u 1:($3-16*5):(17) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 17 u 1:($3-17*5):(18) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 18 u 1:($3-18*5):(19) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 19 u 1:($3-19*5):(20) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 20 u 1:($3-20*5):(21) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 21 u 1:($3-21*5):(22) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 22 u 1:($3-22*5):(23) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 23 u 1:($3-23*5):(24) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 24 u 1:($3-24*5):(25) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 25 u 1:($3-25*5):(26) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 26 u 1:($3-26*5):(27) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 27 u 1:($3-27*5):(28) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 28 u 1:($3-28*5):(29) noti w li lc palette,\
 'out/test.out.comp.Int3_b' i 29 u 1:($3-29*5):(30) noti w li lc palette
 