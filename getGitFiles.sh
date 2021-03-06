#!/bin/bash -e
#

# Find all .C, .H, .sh, Makefile, files and options files to put in git repository

exts=".C .H .sh Dict"
fileNames="Makefile files options runAll gmtDicts fvSchemes fvSolution README"
excludes="lnInclude linux ~"

files=''

for ext in $exts; do
    newFiles=`find . -name "*$ext"`
    files="$files $newFiles"
done

for name in $fileNames; do
    newFiles=`find . -name $name`
    files="$files $newFiles"
done

for ex in $excludes; do
    filesTmp=`for file in $files; do echo $file; done | grep -v $ex`
    files=$filesTmp
done

git add $files


# Other git things needed (in this order)
# git init
# ./getGitFiles.sh
# git commit -m 'All my developments in 2.3.0'
# git remote add myOF https://github.com/hilaryweller0/myOF
# git pull https://github.com/hilaryweller0/myOF
# git push myOF master


# how to clone and rename:
git clone https://github.com/hilaryweller0/myOF Joshua-2.3.1

# how to pull AtmosFOAM
git pull https://github.com/hertzsprung/AtmosFOAM

