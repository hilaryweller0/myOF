#!/bin/bash -e

./comparePlots/errorNorms.sh DubosH_CLUSTPV_CLUSTh
./comparePlots/errorNorms.sh asymmetricH_CLUSTPV_CLUSTh

./comparePlots/errorNormsCd.sh DubosH_CLUSTPV_CLUSTh
./comparePlots/errorNormsCd.sh asymmetricH_CLUSTPV_CLUSTh

./comparePlots/errorNormsDx.sh asymmetricH_CLUSTPV_CLUSTh
./comparePlots/errorNormsDx.sh DubosH_CLUSTPV_CLUSTh

./comparePlots/errorNormsDxCd.sh asymmetricH_CLUSTPV_CLUSTh
./comparePlots/errorNormsDxCd.sh DubosH_CLUSTPV_CLUSTh

./comparePlots/perpErrorDx.sh
./comparePlots/perpErrorDxCd.sh

./comparePlots/HErrorDx.sh
./comparePlots/HErrorDxCd.sh

