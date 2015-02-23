#!/bin/bash -e

# inputs common to each graph
legends=("perp error" "1st order")

pens=("2pt,black" 
      "1pt,black,16_16:0" "1pt,black,16_16:0")
symbols=("x8p" "x0p")
ymax=7
ymin=5e-2
ylabel="maximum error (m/s)"

xlabel="dofs"
#xlabel="number of primal mesh cells"
xmin=500
xmax=1.5e6
#xmin=500
#xmax=1e5
dx=1
ddx=1
dxg=1
dy=1
ddy=1
dyg=1
projection=X9cl/6cl
gv=0

# all data files
gridTypes=(HRbucky HRbucky HRbucky HRbucky HRbucky HRbucky HRbucky)

gridDirs=(HRbucky/2
          HRbucky/3
          HRbucky/4
          HRbucky/5
          HRbucky/6
          HRbucky/7
          HRbucky/8)

inputFiles=(comparePlots/gridType_HRbucky.dat
            comparePlots/1stOrder.dat)

rm -f comparePlots/gridType_*.dat
id=0
while [ $id -ne ${#gridDirs[*]} ]; do
    dir=${gridDirs[$id]}
    ndofs=`grep startFace $dir/constant/polyMesh/boundary | tail -1 | \
           awk '{print $2}' | awk -F';' '{print $1}'`
    perpError=`grep perp $dir/testPerpLog | head -1 | awk '{print $5}'`
    perpErrorR=`grep perp $dir/testPerpLog | tail -1 | awk '{print $6}'`
    echo $ndofs $perpError $perpErrorR >> comparePlots/gridType_${gridTypes[$id]}.dat
    let id=$id+1
done

# order of accuarcy curves for ndofs
y2=`echo | awk '{print 0.5*'$ymax'}'`
x1=`echo | awk '{print '$xmin'*1.2}'`
x2=`echo | awk '{print '$xmax'*0.8}'`
y11=`echo | awk '{print '$y2'*('$x1'/'$x2')^0.5}'`
y12=`echo | awk '{print '$y2'*('$x1'/'$x2')}'`
echo -e $x1 $y2\\n$x2 $y11 > comparePlots/1stOrder.dat

col=(2 2)
outFile=HRbucky/plots/perpError.eps
legPos=x0.2/2.3
. gmtPlot

col=(3 2)
outFile=HRbucky/plots/perpErrorR.eps
legPos=x4.2/6
. gmtPlot

