#!/bin/bash -e

# inputs common to each graph

baseCase=5

inputFiles=($baseCase/2/errorDiags.dat \
            $baseCase/4/errorDiags.dat)

legends=("x2" "x4")

pens=("2pt,blue" "2pt,red")
#xlabel="time (days)"
xmin=0
xmax=5
dx=1
ddx=1
dy=1
ddy=3
xscale=/86400
yscale=*1
projection=X9c/6cl
legPos=x3.2/2.1
gv=1

# for error plots
ymax=0.1
ymin=1e-5

# l2(h)
#ylabel="l@-2@-(h)"
outFile=$baseCase/l2h.eps
col=2
. gmtPlot

# linf(h)
#ylabel="l@-@~\\245@~@-(h)"
outFile=$baseCase/lih.eps
col=3
. gmtPlot

# l2(u)
#ylabel="l@-2@-(u)"
outFile=$baseCase/l2u.eps
col=15
. gmtPlot

# linf(u)
#ylabel="l@-@~\\245@~@-(u)"
outFile=$baseCase/liu.eps
col=16
. gmtPlot

# for conservation plots
xlabel="time (days)"
inputFiles=(${inputFiles[*]} ${inputFiles[*]})
pens=("2pt,blue" "2pt,red" "2pt,blue,10_5:0" "2pt,red,10_5:0")
yscale=("*(1)" "*(1)" "*(-1)" "*(-1)")

# energy
#ylabel="normalised energy change"
outFile=$baseCase/energy.eps
col=4
ymax=2e-6
ymin=1e-10
. gmtPlot

# enstrophy
#ylabel="normalised enstrophy change"
outFile=$baseCase/enstrophy.eps
col=5
ymax=2e-4
ymin=1e-9
. gmtPlot

