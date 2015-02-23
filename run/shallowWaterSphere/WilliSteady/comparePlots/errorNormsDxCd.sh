#!/bin/bash -e

#scheme=asymH_CLUST50PV_midhH
#scheme=Dubos_CLUST50PV_midh
scheme=$1

# inputs common to each graph
legends=("Hexagonal" \
         "Cubed sphere" \
         "Diamondised cube" "1st/2nd order")

pens=("2pt,purple" "2pt,blue" "2pt,red"
      "1pt,black,16_16:0" "1pt,black,16_16:0")
symbols=("x8p" "c8p" "s8p" "d8p" "x0p" "x0p")
ymax=0.2
ymin=1e-5

xlabel=""
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
gridTypes=( HRbucky_Cd HRbucky_Cd HRbucky_Cd
       HRbucky_Cd HRbucky_Cd HRbucky_Cd
      cubeCd cubeCd cubeCd cubeCd cubeCd cubeCd cubeCd cubeCd cubeCd
      diamondCubeCd diamondCubeCd diamondCubeCd diamondCubeCd diamondCubeCd
       diamondCubeCd diamondCubeCd diamondCubeCd)

gridDirs=(HRbucky/3_Cd/save/dt3600_${scheme}
          HRbucky/4_Cd/save/dt3600_${scheme}
          HRbucky/5_Cd/save/dt1800_${scheme}
          HRbucky/6_Cd/save/dt900_${scheme}
          HRbucky/7_Cd/save/dt450_${scheme}
          HRbucky/8_Cd/save/dt225_${scheme}
          cubeCd/6x6_eq/save/dt3600_${scheme}
          cubeCd/12x12_eq/save/dt3600_${scheme}
          cubeCd/17x17_eq/save/dt2400_${scheme}
          cubeCd/24x24_eq/save/dt1800_${scheme}
          cubeCd/32x32_eq/save/dt1200_${scheme}
          cubeCd/48x48_eq/save/dt900_${scheme}
          cubeCd/72x72_eq/save/dt600_${scheme}
          cubeCd/144x144_eq/save/dt300_${scheme}
          cubeCd/288x288_eq/save/dt150_${scheme}
          diamondCubeCd/6x6_eq/save/dt3600_${scheme}
          diamondCubeCd/12x12_eq/save/dt2400_${scheme}
          diamondCubeCd/17x17_eq/save/dt1800_${scheme}
          diamondCubeCd/24x24_eq/save/dt1200_${scheme}
          diamondCubeCd/32x32_eq/save/dt900_${scheme}
          diamondCubeCd/48x48_eq/save/dt600_${scheme}
          diamondCubeCd/72x72_eq/save/dt450_${scheme}
          diamondCubeCd/144x144_eq/save/dt225_${scheme})

time=432000

inputFiles=(comparePlots/gridType_HRbucky_Cd.dat
            comparePlots/gridType_cubeCd.dat
            comparePlots/gridType_diamondCubeCd.dat
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
outFile=comparePlots/l2h_dofs_${scheme}Cd.eps
. gmtPlot

col=(3 3 3 2 2)
outFile=comparePlots/lih_dofs_${scheme}Cd.eps
. gmtPlot

col=(4 4 4 2 2)
outFile=comparePlots/l2u_dofs_${scheme}Cd.eps
. gmtPlot

col=(5 5 5 2 2)
outFile=comparePlots/liu_dofs_${scheme}Cd.eps
. gmtPlot

