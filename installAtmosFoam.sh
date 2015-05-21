#!/bin/bash -e

# Instructions for installing AtmosFOAM on ubuntu

# This installation won't run as a script. You will need
# to copy line by line onto the command line
# You will need admin previligdes for sudo commands
# The file also includes documentation

# installing openfoam
# instructions from:
# http://www.openfoam.org/download/ubuntu.php
VERS=$(lsb_release -cs)
sudo sh -c "echo deb http://www.openfoam.org/download/ubuntu $VERS main > /etc/apt/sources.list.d/openfoam.list" 
sudo apt-get update
sudo apt-get install openfoam231


#** To use OpenFOAM please add
#**
#**    . /opt/openfoam231/etc/bashrc
#**
#** To your ~/.bashrc
. /opt/openfoam231/etc/bashrc

# pre-requisites for installing gmt
sudo apt-get install cmake libnetcdf-dev

# installing gmt (5.1.?)
# instructions at
# http://gmt.soest.hawaii.edu/projects/gmt/wiki/Download
# to get a sufficiently recent version, it needs to be compiled from sources
mkdir -p $HOME/gmt
cd $HOME/gmt
wget ftp://ftp.geologi.uio.no/pub/gmt/gmt-5.1.1-src.tar.bz2 \
     ftp://ftp.geologi.uio.no/pub/gmt/gshhg-gmt-2.3.4.tar.gz
     ftp://ftp.geologi.uio.no/pub/gmt/README.GMT
bzip2 -dc gmt-5.1.1-src.tar.bz2 | tar xvf -
cd gmt-5.1.1-src
mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make -j1
sudo make -j1 install


#** To use gmt please add the following lines to .bashrc
#** 
#**    export GMTU=$HOME/gmt/usr
#**    export GMTC=$GMTU/colours/special


# installing git
sudo apt-get install git

# install gv (for viewing postscript files)
sudo apt-get install gv

# install pdfcrop
sudo apt-get install texlive-extra-utils

# installing AtmosFOAM
mkdir -p $HOME/OpenFOAM
cd $HOME/OpenFOAM
git clone https://github.com/hertzsprung/AtmosFOAM AtmosFOAM-$WM_PROJECT_VERSION
cd AtmosFOAM-$WM_PROJECT_VERSION

#** add the following line to your $HOME/.bashrc
#** 
#** export ATMOSFOAM_SRC=$HOME/$WM_PROJECT/AtmosFOAM-$WM_PROJECT_VERSION/src
#** 

export ATMOSFOAM_SRC=$HOME/$WM_PROJECT/AtmosFOAM-$WM_PROJECT_VERSION/src

# and then build AtmosFOAM:
./Allwmake



# To create and run an OpenFOAM tutorial case
mkdir -p $FOAM_RUN/tutorials/incompressible/icoFoam
cp -r $FOAM_TUTORIALS/incompressible/icoFoam/cavity \
      $FOAM_RUN/tutorials/incompressible/icoFoam
cd $FOAM_RUN/tutorials/incompressible/icoFoam/cavity
blockMesh
icoFoam

# Instructions on how to use gmtFoam (for plotting results)
# from the OpenFOAM case where you want to use it (assuming this is 
# $FOAM_RUN/tutorials/incompressible/icoFoam/cavity)
cd $FOAM_RUN/tutorials/incompressible/icoFoam/cavity
cp -r $HOME/OpenFOAM/AtmosFOAM-$WM_PROJECT_VERSION/applications/utilities/postProcessing/gmtFoam/gmtDicts constant

# and to finally run gmtFoam to create a mesh plot and a results plot
gmtFoam -time 0 mesh
gv 0/mesh.pdf &
gmtFoam -time 0.5 pU
gv 0.5/pU.pdf &

# Have a look at the gmt dictionaries in the constant/gmtDicts directory
# too see you to plot things differently


# not pre-requesites but other things I installed to help

sudo apt-get install gedit-plugins openssh-server #\
#    gv psutils gedit-plugins gedit-latex-plugin

# fix for gedit (which you may not need)
# sudo chown -R $USER:$USER ~/.cache



