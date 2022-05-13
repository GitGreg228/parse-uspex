#!/bin/sh

for dir_name in `ls -d */`
do
    echo $dir_name
    cd $dir_name
    sbatch script.sh
    cd ..
done
