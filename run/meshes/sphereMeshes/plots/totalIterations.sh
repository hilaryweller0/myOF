#!/bin/bash -e

# calculates the total number of inner and outer iterations to reach a solution
# of the Monge-Ampere solution for each resolution

baseCase=MongeAmpere

# find the number of cells and the total number of iterations for each run and
# store the result in plots/X?.dat
rm -f plots/X?.dat
for base in 3 4 5 6 7; do
    for res in 2 4 8; do
        nIts=`grep GAMG $baseCase/$base/$res/log \
               | awk 'BEGIN {sum=0} {sum += $15} END {print sum}'`
        nCells=`grep nFaces $baseCase/$base/$res/0/polyMesh/boundary | head -1 \
                | awk '{print $2}' | awk -F';' '{print $1}'`
        echo $nCells $nIts >> plots/X$res.dat
    done
done

outFile=plots/totalIterations.eps
inputFiles=(plots/X2.dat plots/X4.dat plots/X8.dat)


pens=("1,black," "1,black,5_5:" "1,black,2_4_10_4:")
symbols=("x10p" "x10p" "x10p")
spens=("1,black," "1,black," "1,black,")
legends=("X2" "X4" "X8")
col=(2)
colx=1
xlabel='Number of cells'
ylabel='Total number of iterations'
xmin=0
xmax=45000
dx=5000
ymin=20
ymax=225
dy=20
legPos=x12/2
nSkip=0
projection=X15c/10c
gv=0

. gmtPlot


