


#./raw/adjust_go.com # extract out only the region of interest in the raw input files

#g++ prog/disty.c -lgslcblas -lgsl -Wno-write-strings  -O3 -o disty.out 

#./bin/champ_Darwin input.txt
./bin/champ_Darwin ./raw/Idt1.inp
gnuplot figs/spectrafit.gp
chmod +x figs/make_summary.com
./figs/make_summary.com
rm latex.log
rm summary.tex