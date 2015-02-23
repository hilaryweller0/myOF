#!/bin/bash -e

case=save/SIg_H
#gedit -s $case/system/controlDict $case/system/fvSchemes $case/system/fvSolution
exnerFoamNonOrthogIMEXEX -case $case > $case/log 2>&1 &
#tail -f $case/log

case=save/SI_H
#gedit -s $case/system/controlDict $case/system/fvSchemes $case/system/fvSolution
exnerFoamNonOrthogIMEXEX -case $case > $case/log 2>&1 &
#tail -f $case/log

case=save/SIg_gradp
#gedit -s $case/system/controlDict $case/system/fvSchemes $case/system/fvSolution
exnerFoamIMEXEX -case $case > $case/log 2>&1 &
#tail -f $case/log

case=save/SI_gradp
#gedit -s $case/system/controlDict $case/system/fvSchemes $case/system/fvSolution
exnerFoamIMEXEX -case $case > $case/log 2>&1 &
#tail -f $case/log

