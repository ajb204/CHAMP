reset
set term post eps enh solid color 20
set output 'figs/test.outInt2_b.eps'
set title 'Fitting spectrum: Int2_b'
set label "MZ per point: 3.65" at graph 0.02,graph 0.95
set label "chi2/dof:     19.49" at graph 0.02,graph 0.9
set label "ave Error:    4.42%" at graph 0.02,graph 0.85
set xlabel 'm/z'
set xrange [*:*]
set yrange [-20:110]
plot 'out/test.outInt2_b' u 1:2 ti 'raw' w li lt 1,\
'' u 1:3 ti 'fitted' w li lt 2,\
'' u 1:(($2-$3)-10) ti 'difference' w li lt 3
unset label
reset
set term post eps enh solid color 20
unset key
set output 'figs/test.inpInt2_b.eps'
set title 'Fitted distribution: Int2_b'
set label "Distribution: alphaB" at graph 0.02,graph 0.950000
set label "rat  77.00" at graph 0.02,graph 0.900000
set label "mult 75.64" at graph 0.02,graph 0.900000
unset cbtics
set xlabel 'number of sHSPs'
plot 'out/test.inpInt2_b' u 1:2:(1) ti 'distibution' w boxes
unset label
reset
set term post eps enh solid color 20
set output 'figs/test.out.comp.Int2_b.eps'
set title 'Fitted mass spectrum: Int2_b'
set xlabel 'm/z'
set xrange [1400.000000:21500.000000]
set view map
set cblabel 'Oligomer'
set label "Mode:      fitty " at graph 0.02,graph 0.95
set label "Zfudge:    9.70 (1)" at graph 0.02,graph 0.9
set label "adduction: 4.68e-02 (1)" at graph 0.02,graph 0.85
set label "MRes:      18.87 (1)" at graph 0.02,graph 0.8
set label "Resfudge:  1.00e-12 (1)" at graph 0.02,graph 0.75
set label "Zwidth:    4.63e+00 (1)" at graph 0.02,graph 0.7
set label "shspMass:  2.02e+04 (0)" at graph 0.02,graph 0.65
splot \
'out/test.outInt2_b' u 1:3:(1) ti 'full' w li lt 2,\
'out/test.out.indiv.Int2_b' i 0 u 1:($3-0*5):(1) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 1 u 1:($3-1*5):(2) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 2 u 1:($3-2*5):(3) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 3 u 1:($3-3*5):(4) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 4 u 1:($3-4*5):(5) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 5 u 1:($3-5*5):(6) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 6 u 1:($3-6*5):(7) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 7 u 1:($3-7*5):(8) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 8 u 1:($3-8*5):(9) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 9 u 1:($3-9*5):(10) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 10 u 1:($3-10*5):(11) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 11 u 1:($3-11*5):(12) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 12 u 1:($3-12*5):(13) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 13 u 1:($3-13*5):(14) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 14 u 1:($3-14*5):(15) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 15 u 1:($3-15*5):(16) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 16 u 1:($3-16*5):(17) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 17 u 1:($3-17*5):(18) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 18 u 1:($3-18*5):(19) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 19 u 1:($3-19*5):(20) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 20 u 1:($3-20*5):(21) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 21 u 1:($3-21*5):(22) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 22 u 1:($3-22*5):(23) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 23 u 1:($3-23*5):(24) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 24 u 1:($3-24*5):(25) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 25 u 1:($3-25*5):(26) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 26 u 1:($3-26*5):(27) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 27 u 1:($3-27*5):(28) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 28 u 1:($3-28*5):(29) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 29 u 1:($3-29*5):(30) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 30 u 1:($3-30*5):(31) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 31 u 1:($3-31*5):(32) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 32 u 1:($3-32*5):(33) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 33 u 1:($3-33*5):(34) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 34 u 1:($3-34*5):(35) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 35 u 1:($3-35*5):(36) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 36 u 1:($3-36*5):(37) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 37 u 1:($3-37*5):(38) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 38 u 1:($3-38*5):(39) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 39 u 1:($3-39*5):(40) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 40 u 1:($3-40*5):(41) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 41 u 1:($3-41*5):(42) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 42 u 1:($3-42*5):(43) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 43 u 1:($3-43*5):(44) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 44 u 1:($3-44*5):(45) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 45 u 1:($3-45*5):(46) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 46 u 1:($3-46*5):(47) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 47 u 1:($3-47*5):(48) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 48 u 1:($3-48*5):(49) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 49 u 1:($3-49*5):(50) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 50 u 1:($3-50*5):(51) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 51 u 1:($3-51*5):(52) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 52 u 1:($3-52*5):(53) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 53 u 1:($3-53*5):(54) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 54 u 1:($3-54*5):(55) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 55 u 1:($3-55*5):(56) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 56 u 1:($3-56*5):(57) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 57 u 1:($3-57*5):(58) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 58 u 1:($3-58*5):(59) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 59 u 1:($3-59*5):(60) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 60 u 1:($3-60*5):(61) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 61 u 1:($3-61*5):(62) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 62 u 1:($3-62*5):(63) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 63 u 1:($3-63*5):(64) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 64 u 1:($3-64*5):(65) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 65 u 1:($3-65*5):(66) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 66 u 1:($3-66*5):(67) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 67 u 1:($3-67*5):(68) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 68 u 1:($3-68*5):(69) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 69 u 1:($3-69*5):(70) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 70 u 1:($3-70*5):(71) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 71 u 1:($3-71*5):(72) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 72 u 1:($3-72*5):(73) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 73 u 1:($3-73*5):(74) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 74 u 1:($3-74*5):(75) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 75 u 1:($3-75*5):(76) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 76 u 1:($3-76*5):(77) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 77 u 1:($3-77*5):(78) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 78 u 1:($3-78*5):(79) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 79 u 1:($3-79*5):(80) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 80 u 1:($3-80*5):(81) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 81 u 1:($3-81*5):(82) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 82 u 1:($3-82*5):(83) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 83 u 1:($3-83*5):(84) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 84 u 1:($3-84*5):(85) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 85 u 1:($3-85*5):(86) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 86 u 1:($3-86*5):(87) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 87 u 1:($3-87*5):(88) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 88 u 1:($3-88*5):(89) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 89 u 1:($3-89*5):(90) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 90 u 1:($3-90*5):(91) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 91 u 1:($3-91*5):(92) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 92 u 1:($3-92*5):(93) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 93 u 1:($3-93*5):(94) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 94 u 1:($3-94*5):(95) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 95 u 1:($3-95*5):(96) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 96 u 1:($3-96*5):(97) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 97 u 1:($3-97*5):(98) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 98 u 1:($3-98*5):(99) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 99 u 1:($3-99*5):(100) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 100 u 1:($3-100*5):(101) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 101 u 1:($3-101*5):(102) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 102 u 1:($3-102*5):(103) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 103 u 1:($3-103*5):(104) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 104 u 1:($3-104*5):(105) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 105 u 1:($3-105*5):(106) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 106 u 1:($3-106*5):(107) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 107 u 1:($3-107*5):(108) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 108 u 1:($3-108*5):(109) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 109 u 1:($3-109*5):(110) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 110 u 1:($3-110*5):(111) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 111 u 1:($3-111*5):(112) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 112 u 1:($3-112*5):(113) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 113 u 1:($3-113*5):(114) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 114 u 1:($3-114*5):(115) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 115 u 1:($3-115*5):(116) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 116 u 1:($3-116*5):(117) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 117 u 1:($3-117*5):(118) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 118 u 1:($3-118*5):(119) noti w li lc palette,\
 'out/test.out.indiv.Int2_b' i 119 u 1:($3-119*5):(120) noti w li lc palette
 