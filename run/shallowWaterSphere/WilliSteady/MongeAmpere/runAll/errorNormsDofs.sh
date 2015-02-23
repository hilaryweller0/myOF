#!/bin/bash
# -e

# inputs common to each graph
legends=("x1" "x2" "x4" "1st/2nd order")

pens=("black" "blue" "red"
      "black,16_16:0" "black,16_16:0")
symbols=("+8p" "x8p" "c8p" "x0p" "x0p")
ymax=1
ymin=1e-5

xlabel="degrees of freedom"
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
legPos=x0.2/2.3
gv=1

# all data files
gridTypes=(x1 x1 x1 x1
           x2 x2 x2 x2 x2
           x4 x4 x4 x4 x4)

dir1=../../../../../hilary-2.2.2/run/shallowWaterSphere/WilliSteady/HRbucky

gridDirs=($dir1/4/save/dt1800_diagonalH_CLUSTvort_CLUSTh
          $dir1/5/save/dt900_diagonalH_CLUSTvort_CLUSTh
          $dir1/6/save/dt450_diagonalH_CLUSTvort_CLUSTh
          $dir1/7/save/dt225_diagonalH_CLUSTvort_CLUSTh
          3/2 4/2 5/2 6/2 7/2
          3/4 4/4 5/4 6/4 7/4)

time=432000

inputFiles=(comparePlots/gridType_x1.dat
            comparePlots/gridType_x2.dat
            comparePlots/gridType_x4.dat
            comparePlots/1stOrder.dat comparePlots/2ndOrder.dat)

rm -f comparePlots/gridType_*.dat
id=0
while [ $id -ne ${#gridDirs[*]} ]; do
    dir=${gridDirs[$id]}
    ndofs=`grep startFace $dir/constant/polyMesh/boundary | tail -1 | \
           awk '{print $2}' | awk -F';' '{print $1}'`
    l2h=`tail -1 $dir/errorDiags.dat | grep $time | awk '{print $2}'`
    lih=`tail -1 $dir/errorDiags.dat | grep $time | awk '{print $3}'`
    l2u=`tail -1 $dir/errorDiags.dat | grep $time | awk '{print $15}'`
    liu=`tail -1 $dir/errorDiags.dat | grep $time | awk '{print $16}'`
    echo $ndofs $l2h $lih $l2u $liu >> comparePlots/gridType_${gridTypes[$id]}.dat
    let id=$id+1
done

# order of accuarcy curves for ndofs
y2=`echo | awk '{print 0.5*'$ymax'}'`
x1=`echo | awk '{print '$xmin'*1.2}'`
x2=`echo | awk '{print '$xmax'*0.8}'`
y11=`echo | awk '{print '$y2'*('$x1'/'$x2')^0.5}'`
y12=`echo | awk '{print '$y2'*('$x1'/'$x2')}'`
echo -e $x1 $y2\\n$x2 $y11 > comparePlots/1stOrder.dat
echo -e $x1 $y2\\n$x2 $y12 > comparePlots/2ndOrder.dat

col=(2 2 2 2 2)
outFile=comparePlots/l2h_dofs.eps
. gmtPlot

col=(3 3 3 2 2)
outFile=comparePlots/lih_dofs.eps
. gmtPlot

col=(4 4 4 2 2)
outFile=comparePlots/l2u_dofs.eps
. gmtPlot

col=(5 5 5 2 2)
outFile=comparePlots/liu_dofs.eps
. gmtPlot

#rm -f comparePlots/gridType_*.dat

