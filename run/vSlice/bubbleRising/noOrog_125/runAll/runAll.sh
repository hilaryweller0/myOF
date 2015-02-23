# running
exnerFoamH > log 2>&1 &
exnerFoam > log 2>&1 &

# post processing
time=1000
#sumFields $time thetaDiff $time theta 0 theta_init -scale0 1 -scale1 -1
#plotPatchData -time $time theta
plotPatchData -time $time theta_BF02
gv $time/theta_BF02.eps &
writeuvw -time $time Uf
plotPatchData -time $time w_BF02
gv $time/w_BF02.eps &

makecpt -Cjet -T300.001/302.001/0.1 > constant/gmtDicts/white_jet.cpt
plotPatchData -time $time thetaNNS
gv $time/thetaNNS.eps &

gmtPlot constant/gmtPlots/plotw.gmt



