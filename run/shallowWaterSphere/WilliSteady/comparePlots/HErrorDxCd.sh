#!/bin/bash -e

# inputs common to each graph
legends=("Hexagonal" \
         "Cubed sphere" \
         "Diamondised cube" "1st order")

pens=("2pt,purple" "2pt,blue" "2pt,red"
      "1pt,black,16_16:0" "1pt,black,16_16:0")
symbols=("c8p" "s8p" "d8p" "x0p")
ymax=1
ymin=1e-7
#ylabel="error"

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
legPos=x0.1/1.8
gv=0

# all data files
gridTypes=( HRbucky_Cd HRbucky_Cd HRbucky_Cd
       HRbucky_Cd HRbucky_Cd HRbucky_Cd
      cubeCd cubeCd cubeCd cubeCd cubeCd cubeCd cubeCd cubeCd cubeCd
      diamondCubeCd diamondCubeCd diamondCubeCd diamondCubeCd diamondCubeCd
       diamondCubeCd diamondCubeCd diamondCubeCd)

gridDirs=(HRbucky/3_Cd
          HRbucky/4_Cd
          HRbucky/5_Cd
          HRbucky/6_Cd
          HRbucky/7_Cd
          HRbucky/8_Cd
          cubeCd/6x6_eq
          cubeCd/12x12_eq
          cubeCd/17x17_eq
          cubeCd/24x24_eq
          cubeCd/32x32_eq
          cubeCd/48x48_eq
          cubeCd/72x72_eq
          cubeCd/144x144_eq
          cubeCd/288x288_eq
          diamondCubeCd/6x6_eq
          diamondCubeCd/12x12_eq
          diamondCubeCd/17x17_eq
          diamondCubeCd/24x24_eq
          diamondCubeCd/32x32_eq
          diamondCubeCd/48x48_eq
          diamondCubeCd/72x72_eq
          diamondCubeCd/144x144_eq)

Htypes=(asymmetric Dubos)

inputFiles=(comparePlots/gridType_HRbucky_Cd.dat
            comparePlots/gridType_cubeCd.dat
            comparePlots/gridType_diamondCubeCd.dat
            comparePlots/1stOrder.dat)

for Htype in ${Htypes[*]}; do
    rm -f comparePlots/gridType_*.dat
    id=0
    while [ $id -ne ${#gridDirs[*]} ]; do
        dir=${gridDirs[$id]}/save/*${Htype}*
        ndofs=`grep startFace $dir/constant/polyMesh/boundary | tail -1 | \
               awk '{print $2}' | awk -F';' '{print $1}'`
        HError=`grep 'Maximum H error' $dir/testPerpLog | head -1 | awk '{print $5}'`
        echo $ndofs $HError >> comparePlots/gridType_${gridTypes[$id]}.dat
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
    outFile=comparePlots/HError${Htype}Cd.eps
    . gmtPlot
done

