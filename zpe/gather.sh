#!/bin/sh

for dir_name in `ls -d */`
do
    cd $dir_name
    stat=$dir_name": "
    if [ -d phonopy/ ]; then
        cd phonopy/
        echo $PWD
        phonstat=""
        if [ ! -f total_dos.pdf ]; then
            for disp in `ls -d disp-*/`
            do
                cd $disp
                if [ -f log ]; then
                    err=`cat err`
                    err=$err`grep 'BAD TERMINATION' log`
                    err=$err`grep 'I REFUSE TO CONTINUE WITH THIS SICK JOB' log`
                    log=`tail -n 1 log | awk '{print $2}'`
                    if [ "$log" == "F=" ]; then
                        problems=`grep 'very serious problems' log`
                        if [ "$problems" == "" ]; then
                            phonstat=$phonstat""
                        else
                            phonstat=$phonstat"problems in $disp "
                            #sed -i "s/NPAR = 4/NPAR = 8/g" INCAR
                            #sed -i "s/NPAR = 8/NPAR = 16/g" INCAR
                            #sed -i "s/NPAR = 16/NPAR = 8/g" INCAR
                            #sed -i "s/SBATCH -n 8/SBATCH -n 16/g" script.sh
                            #sed -i "s/SBATCH -n 16/SBATCH -n 32/g" script.sh
                            #sbatch script.sh
                        fi
                    else
                        if [ "$err" != "" ]; then
                            phonstat=$phonstat" error in $disp"
                            #sed -i "s/NPAR = 4/NPAR = 8/g" INCAR
                            #sed -i "s/NPAR = 8/NPAR = 16/g" INCAR
                            #sed -i "s/NPAR = 16/NPAR = 8/g" INCAR
                            #sed -i "s/SBATCH -n 8/SBATCH -n 16/g" script.sh
                            #sed -i "s/SBATCH -n 16/SBATCH -n 32/g" script.sh
                            #sbatch script.sh
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
            echo "creating FORCE_SETS"
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
        dim=`grep dim phonopy_disp.yaml | awk -F\" '{print $2}'`
        tol=`grep symmetry_tolerance: phonopy_disp.yaml | tail -n 1 | awk -F\" '{print $2}'`
        if [ ! -f mesh.conf ]; then
                echo "creating mesh.conf"
                atom_name=`sed '6q;d' POSCAR`
                cat>mesh.conf <<!
ATOM_NAME = $atom_name
DIM = $dim
MP = 31 31 31
!
        fi
        if [ ! -f total_dos.pdf ]; then
                echo "creating total_dos.pdf"
                phonopy -s -p mesh.conf --tolerance=$tol
        fi
        if [ ! -f thermal_properties.yaml ]; then
                echo "creating thermal_properties.yaml"
                phonopy -t -p mesh.conf --tolerance=$tol --tmax 2000 > thermal.log
        fi
        if [ ! -f pdos.conf ]; then
                echo "creating pdos.conf"
                atom_name=`sed '6q;d' POSCAR`
                cat>pdos.conf <<!
ATOM_NAME = $atom_name
DIM = $dim
PDOS = AUTO
!
        fi
        if [ ! -f partial_dos.pdf ]; then
                echo "creating partial_dos.pdf"
                phonopy -s -p pdos.conf --tolerance=$tol
        fi
	if [ ! -f band.conf ]; then
                echo "creating band.conf"
                atom_name=`sed '6q;d' POSCAR`
                cat>band.conf <<!
ATOM_NAME = $atom_name
DIM = $dim
BAND = AUTO
!
        fi
        if [ ! -f phonon_band.dat ]; then
                 echo "creating phonon_band.dat"
                 phonopy -p band.conf --tolerance=$tol
        fi
        if [ ! -f band.pdf ]; then
                 phonopy -p -s band.conf --tolerance=$tol
        fi
        if [ ! -f band-pdos.conf ]; then
                echo "creating band-pdos.conf"
                atom_name=`sed '6q;d' POSCAR`
                cat>band-pdos.conf <<!
ATOM_NAME = $atom_name
DIM = $dim
BAND = AUTO
PDOS = AUTO
!
        fi
        if [ ! -f band_dos.pdf ]; then
                 phonopy -p -s band-pdos.conf --tolerance=$tol
        fi
        if [ ! -f phonon_band.pdf ]; then
                echo "creaing phonon_band.pdf"
                sumo-phonon-bandplot --dim $dim --height 6 --width 9 --dos total_dos.dat --symprec=$tol
        fi
        rm mesh.yaml
        cd ..
    fi
    printf "$stat $phonstat \n"
    cd ..
done

