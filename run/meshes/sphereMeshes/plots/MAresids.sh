# !/bin/bash -e

# Plot the initial residuals as a function of iteration number for all of the 
# MongeAmpere meshes

cases=(MongeAmpereV0/4 MongeAmpereV0/5 MongeAmpereV0/6 MongeAmpereV0/7
       MongeAmpereV1/4 MongeAmpereV1/5 MongeAmpereV1/6 MongeAmpereV1/7)

refs=(2 4 8 16)

# inputs the same for all graphs
colx=1
legends=()
for ref in ${refs[*]}; do
    legends=(${legends[*]} X$ref)
done
pens=("1p,black" "1p,blue" "1p,red" "1p,green")
xlabel=''
#xlabel='Iteration number'
xmin=1
xmax=1000
dx=200
nSkip=1
gv=0

#plots for residual
col=11
#ylabel='Initial residual'
ymin=1e-8
ymax=1
dy=10
dyg=0
legPos=x4/5
projection=X7.5c/5cl

for case in ${cases[*]}; do
    inputFiles=()
    for ref in ${refs[*]}; do
        inputFiles=(${inputFiles[*]} $case/$ref/MAconvergence.dat)
    done
    
    # plot for residuals
    outFile=$case/MAresiduals.eps
    . gmtPlot
    pstitle $outFile
    #gv $outFile &
done

# plots for boostLaplacian
col=10
#ylabel='1+@~a@~'
ymin=1
ymax=100
dy=2
dyg=0
legPos=x5/2
projection=X7.5c/5cl
for case in ${cases[*]}; do
    inputFiles=()
    for ref in ${refs[*]}; do
        inputFiles=(${inputFiles[*]} $case/$ref/MAconvergence.dat)
    done
    
    # plot for boostLaplacian
    outFile=$case/boostLaplacian.eps
    . gmtPlot
    pstitle $outFile
    gv $outFile &
done

