#!/bin/sh

for dir_name in `ls -d */`
do
    cd $dir_name
    if [ -d phonopy/ ]; then
        cd phonopy/
    if [ -f total_dos.pdf ]; then
        if [ -f thermal_properties.yaml ]; then
            echo deleting disps in $dir_name
            rm -rf disp-*/
        fi
        else 
            echo work is not finished in $dir_name
        fi
        cd ..
    fi
    cd ..
done
    