#!/bin/bash -e

# inputs common to each graph

#scheme=asymH_CLUST50PV_midhH
#scheme=Dubos_CLUST50PV_midh
scheme=$1

inputFiles=(HRbucky/5_Cd/save/dt1800_${scheme}/errorDiags.dat \
            cubeCd/24x24_eq/save/dt1800_${scheme}/errorDiags.dat \
            diamondCubeCd/17x17_eq/save/dt1800_${scheme}/errorDiags.dat)

legends=("Hexagonal" \
         "Cubed sphere" \
         "Diamondised cube")

pens=("2pt,purple" "2pt,blue" "2pt,red")
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
gv=0

# for error plots
ymax=0.1
ymin=1e-5

# l2(h)
#ylabel="l@-2@-(h)"
outFile=comparePlots/l2h_${scheme}Cd.eps
col=2
. gmtPlot

# linf(h)
#ylabel="l@-@~\\245@~@-(h)"
outFile=comparePlots/lih_${scheme}Cd.eps
col=3
. gmtPlot

# l2(u)
#ylabel="l@-2@-(u)"
outFile=comparePlots/l2u_${scheme}Cd.eps
col=15
. gmtPlot

# linf(u)
#ylabel="l@-@~\\245@~@-(u)"
outFile=comparePlots/liu_${scheme}Cd.eps
col=16
. gmtPlot

# for conservation plots
xlabel="time (days)"
inputFiles=(${inputFiles[*]} ${inputFiles[*]})
pens=("2pt,purple" "2pt,blue" "2pt,red"
      "2pt,purple,10_5:0" 
      "2pt,blue,10_5:0" "2pt,red,10_5:0")
yscale=("*(1)" "*(1)" "*(1)"
        "*(-1)" "*(-1)" "*(-1)")
gv=1

# energy
#ylabel="normalised energy change"
outFile=comparePlots/energy_${scheme}Cd.eps
col=4
ymax=2e-6
ymin=1e-10
. gmtPlot

# enstrophy
#ylabel="normalised enstrophy change"
outFile=comparePlots/enstrophy_${scheme}Cd.eps
col=5
ymax=2e-4
ymin=1e-9
. gmtPlot

