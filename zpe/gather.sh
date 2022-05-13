#!/bin/sh

for dir_name in `ls -d */`
do
    cd $dir_name
    stat=$dir_name": "
    if [ -d phonopy/ ]; then 
        cd phonopy/
        phonstat=""
        if [ ! -f total_dos.pdf ]; then
            for disp in `ls -d disp-*/`
            do
                cd $disp
                if [ -f log ]; then
                    err=`cat err`
                    log=`tail -n 1 log | awk '{print $2}'`
                    if [ "$log" == "F=" ]; then
                        problems=`grep 'very serious problems' log`
                        if [ "$problems" == "" ]; then
                            phonstat=$phonstat""
                            else
                            phonstat=$phonstat"problems in $disp "
                            #sed -i "s/NPAR = 4/NPAR = 8/g" INCAR
                            #sbatch script.sh
                        fi
                    else
                        if [ "$err" != "" ]; then
                            phonstat=$phonstat" error in $disp"
                        else
                            phonstat=$phonstat"$disp still working "
                        fi
                    fi
                else
                    phonstat=$phonstat"$disp not yet started "
                fi
                cd ..
            done
        fi
        cd ..
    else
        stat=$stat"no phonopy dir created yet"
    fi
    if [ "$phonstat" == "" ]; then
        phonstat="all done"
        cd phonopy
        if [ ! -f FORCE_SETS ]; then
            disp_s=`ls -d disp-*/ | head -n 1`
            disp_f=`ls -d disp-*/ | tail -n 1`
            disp_s="${disp_s/disp-/}"
            disp_f="${disp_f/disp-/}"
            disp_s="${disp_s////}"
            disp_f="${disp_f////}"
            disp_range=disp-\{${disp_s}..${disp_f}\}
            cat>phon.sh<<!
#!/bin/sh
phonopy -f $disp_range/vasprun.xml
!
            chmod 777 phon.sh
            ./phon.sh
        fi
        if [ ! -f mesh.conf ]; then
                atom_name=`sed '6q;d' POSCAR`
                cat>mesh.conf <<!
ATOM_NAME = $atom_name
DIM = 2 2 2
MP = 31 31 31
!
        fi
        if [ ! -f total_dos.pdf ]; then
                phonopy -s -p mesh.conf --tolerance=2e-1
        fi
        if [ ! -f thermal_properties.yaml ]; then
                phonopy -t -p mesh.conf --tolerance=2e-1 --tmax 2000 > thermal.log
        fi
        cd ..
    fi
    printf "$stat $phonstat \n"
    cd ..
done

