#!/bin/bash -e

case=save/dt20_cubicUpCPCFit_H_CN
exnerFoamH -case $case > $case/log 2>&1 &
tail -f $case/log

case=save/dt20_cubicUpCPCFit_dpdx_CN
exnerFoam -case $case > $case/log 2>&1 &
tail -f $case/log

