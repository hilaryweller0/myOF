#!/bin/bash -e

# plots of convergence of the semi-implicit Monge-Ampere method for a range of 
# resolutions and refinements. Convergence is defined as:
# (max(equidistribution) - min(equidistribution))/max(equidistribution)

outFile=plots/MAconvergence.eps
pens=("1,black," "1,black,5_5:" "1,black,2_4_10_4:"
      "1,blue," "1,blue,5_5:" "1,blue,2_4_10_4:"
      "1,red," "1,red,5_5:" "1,red,2_4_10_4:"
      "1,cyan," "1,cyan,5_5:" "1,cyan,2_4_10_4:"
      "1,magenta," "1,magenta,5_5:" "1,magenta,2_4_10_4:")
legends=("3 X2" "3 X4" "3 X8"
         "4 X2" "4 X4" "4 X8"
         "5 X2" "5 X4" "5 X8"
         "6 X2" "6 X4" "6 X8"
         "7 X2" "7 X4" "7 X8")
col=(2)
colx=1
xlabel='Iterations'
ylabel='(max-min)/max of equidistribution'
xmin=1
xmax=40
dx=10
ymin=1e-5
ymax=1
dy=10
legPos=x12/9
nSkip=1
projection=X15c/10cl
gv=1

inputFiles=()
baseCase=MongeAmpere

for base in 3 4 5 6 7; do
    for res in 2 4 8; do
        awk '{print $1, ($5-$4)/($5+1e-12)}' $baseCase/$base/$res/MAconvergence.dat \
            > $baseCase/$base/$res/MAconvergence2.dat
        inputFiles=(${inputFiles[*]} $baseCase/$base/$res/MAconvergence2.dat)
    done
done

. gmtPlot


