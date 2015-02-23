filesToPack=`find src/ExnerTheta src/finiteVolume src/Hops src/sampling \
                  applications/solvers/ExnerFoam/ExnerFoam \
                  applications/solvers/ExnerFoam/ExnerFoamH \
                  applications/utilities/preProcessing/setExnerBalanced \
                  applications/utilities/preProcessing/setExnerBalancedH \
                  applications/utilities/preProcessing/createSpongeLayer \
                  applications/utilities/preProcessing/mapFields \
                  applications/utilities/mesh/add2dMountain \
                  applications/utilities/postProcessing/globalSum \
                  applications/utilities/postProcessing/plotPatchData_gmt4.5 \
                  applications/utilities/postProcessing/plotPatchData_gmt5.1 \
                  applications/utilities/postProcessing/sumFields \
                  applications/utilities/postProcessing/thermoVars \
                  -name "*.H" -o -name "*.C" -o -name files -o -name options | \
                  sed -e '\@Make[.A-Za-z]*/[^/]*/@d' -e '\@/lnInclude@d'`

tar czpf ExnerFoam.tar.gz README $filesToPack
scp ExnerFoam.tar.gz sws02hs@oak.rdg.ac.uk:public_html/AtmosFOAM/ExnerFoam.tar.gz
