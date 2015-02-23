#!/bin/bash -e

# inputs common to each graph
legends=("Orthogonal HR grid" "Centroidal hexagonal" \
         "Centroidal cubed sphere" \
         "Centroidal diamonds" "1st order")

pens=("2pt,black" "2pt,purple" "2pt,blue" "2pt,red"
      "1pt,black,16_16:0" "1pt,black,16_16:0")
symbols=("x8p" "c8p" "s8p" "d8p" "x0p")
ymax=1
ymin=1e-5
ylabel="error"

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
yscale=/38
projection=X9cl/6cl
gv=1

# all data files
gridTypes=(HRbucky HRbucky HRbucky HRbucky HRbucky HRbucky
      HRbucky_nonOrthog HRbucky_nonOrthog HRbucky_nonOrthog
       HRbucky_nonOrthog HRbucky_nonOrthog HRbucky_nonOrthog
      cubeC cubeC cubeC cubeC cubeC cubeC cubeC cubeC cubeC
      diamondCubeC diamondCubeC diamondCubeC diamondCubeC diamondCubeC
       diamondCubeC diamondCubeC diamondCubeC)

gridDirs=(HRbucky/3
          HRbucky/4
          HRbucky/5
          HRbucky/6
          HRbucky/7
          HRbucky/8
          HRbucky/3_nonOrthog
          HRbucky/4_nonOrthog
          HRbucky/5_nonOrthog
          HRbucky/6_nonOrthog
          HRbucky/7_nonOrthog
          HRbucky/8_nonOrthog
          cubeC/6x6_eq
          cubeC/12x12_eq
          cubeC/17x17_eq
          cubeC/24x24_eq
          cubeC/32x32_eq
          cubeC/48x48_eq
          cubeC/72x72_eq
          cubeC/144x144_eq
          cubeC/288x288_eq
          diamondCubeC/6x6_eq
          diamondCubeC/12x12_eq
          diamondCubeC/17x17_eq
          diamondCubeC/24x24_eq
          diamondCubeC/32x32_eq
          diamondCubeC/48x48_eq
          diamondCubeC/72x72_eq
          diamondCubeC/144x144_eq)

inputFiles=(comparePlots/gridType_HRbucky.dat
            comparePlots/gridType_HRbucky_nonOrthog.dat
            comparePlots/gridType_cubeC.dat
            comparePlots/gridType_diamondCubeC.dat
            comparePlots/1stOrder.dat)

rm -f comparePlots/gridType_*.dat
id=0
while [ $id -ne ${#gridDirs[*]} ]; do
    dir=${gridDirs[$id]}/save/*Dubos*
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

col=(2 2 2 2 2)
outFile=comparePlots/perpError.eps
legPos=x0.2/2.3
. gmtPlot

col=(3 3 3 3 2)
outFile=comparePlots/perpErrorR.eps
legPos=x4.2/6
. gmtPlot

