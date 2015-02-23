# !/bin/bash -e

# Plot the initial residuals as a function of iteration number for all of the 
# MongeAmpere meshes and plot boostLaplacian

#inputs same for both graphs
cases=(Fig4.4_V0 Fig4.5_V0 Fig4.4_V1 Fig4.5_V1)
colx=1
legends=("FD, ring mesh" "FD, bell mesh" "geom, ring mesh" "geom, bell mesh")
pens=("1p,black" "1p,blue" "1p,red" "1p,green")
xlabel=''
xlabel='Iteration number'
xmin=1
xmax=1000
dx=200
nSkip=1
gv=0

#plots for residual
col=11
ylabel='Initial residual'
ymin=1e-8
ymax=1
dy=10
dyg=0
legPos=x3/5
projection=X7.5c/5cl

inputFiles=()
for case in ${cases[*]}; do
    inputFiles=(${inputFiles[*]} $case/MAconvergence.dat)
done
# plot for residuals
outFile=plots/MAresids.eps
. gmtPlot
pstitle $outFile
#gv $outFile &

# plots for boostLaplacian
pens=("3p,black" "3p,blue" "1p,red" "1p,green")
col=10
ylabel='1+@~a@~'
ymin=1
ymax=5
dy=1
dyg=0
legPos=x3/2.3
projection=X7.5c/5c
inputFiles=()
for case in ${cases[*]}; do
    inputFiles=(${inputFiles[*]} $case/MAconvergence.dat)
done
# plot for boostLaplacian
outFile=plots/boostLaplacian.eps
. gmtPlot
pstitle $outFile
gv $outFile &

