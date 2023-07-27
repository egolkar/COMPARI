#! /bin/bash
### embedded SGE parameters ###
#$ -S /bin/bash
#$ -cwd 
#$ -o /cbica/home/parkerwi/comp_space/COMPARI/sge/${JOB_ID}.${JOB_NAME}.stdout
#$ -e /cbica/home/parkerwi/comp_space/COMPARI/sge/${JOB_ID}.${JOB_NAME}.stderr
#$ -l short

module load python/anaconda/DiCIPHR

scriptsdir=$(dirname $0)
ldelta=0.01293
bdelta=0.03666
shorder=12
flip=""
memory=48g
usage() {
    cat << EOF 
Usage: COMPARI-pipeline.sh -i dwi.nii -m mask.nii -o outputDir [options]

Options:
    -b  bvals   Text file containing bvals 
    -r  bvecs   Text file containing bvecs 
    -l  int     SH order (Default: 12)
    -d  float   Little delta time (seconds)  (Default: 0.01293)
    -D  float   Big Delta time  (seconds)  (Default: 0.03666)
    -f  x|y|z   Flip bvecs along provided axis. (Default: none)
EOF
echo $* 1>&2
exit 
}
while getopts ":ho:i:m:b:r:l:d:D:f:M:" OPT
do
    case $OPT in
        h) # help
            usage
            ;;
        o) 
            o=$OPTARG
            ;;
        i) 
            dwi=$OPTARG
            ;;
        m) 
            mask=$OPTARG
            ;;
        b) 
            bvals=$OPTARG
            ;;
        r) 
            bvecs=$OPTARG
            ;;
        l) 
            shorder=$OPTARG
            ;;
        d) 
            ldelta=$OPTARG
            ;;
        D) 
            bdelta=$OPTARG
            ;;
        f) 
            flip=$OPTARG
            ;;
        M)  memory=$OPTARG 
            ;;
        *) # getopts issues an error message
            echo "UNHANDLED OPTION" 1>&2 
            usage
            ;;
    esac
done

if [ -z "$dwi" ] || [ -z "$o" ] || [ -z "$mask" ]; then 
    usage "Provide all required arguments (-i, -o, -m)"
fi 
if [ -z "$bvals" ]; then 
    bvals="${dwi%.nii.gz}.bval"
fi 
if [ -z "$bvecs" ]; then 
    bvecs="${dwi%.nii.gz}.bvec"
fi 
flip="--${flip}flip -1"
s=$(basename $o)
#output dir 
mkdir -pv $o/logs

# get slices of mask that have data 
get_mask_slices() {
    python << EOF 
from numpy import any
from nibabel import load 
dat = load("$mask").get_fdata()
slices = list(map(str,filter(lambda i: any(dat[:,i,:]), range(dat.shape[1]))))
print(*slices)
EOF
}
mask_slices=$(get_mask_slices) 

sleep 1 && qsub -l short -l h_vmem=8g -N "COMPARI_pre-${s}" \
    -o $o/logs/\$JOB_NAME.\$JOB_ID.stdout \
    -e $o/logs/\$JOB_NAME.\$JOB_ID.stderr \
    $scriptsdir/wrapper.sh ./COMPARI_preprocessing.py -i $dwi \
            -m $mask -b $bvals -r $bvecs \
            -d $ldelta -s $bdelta -o $o

pickle="$o/pickle.pkl"
hold=""
for i in $mask_slices; 
do
    d=$o/DWIs/DWI${i}.nii.gz
    m=$o/Masks/Mask${i}.nii.gz
    sleep 1 && qsub -hold_jid "COMPARI_pre-${s}" -l h_vmem=$memory -N "COMPARI_fit-${s}_${i}" \
        -o $o/logs/\$JOB_NAME.\$JOB_ID.stdout \
        -e $o/logs/\$JOB_NAME.\$JOB_ID.stderr \
        $scriptsdir/wrapper.sh ./COMPARI_Fractions_FODs.py \
                -i $d -o $o/results/$i -m $m -b $bvals -r $bvecs -p $pickle \
                $flip -d $ldelta -s $bdelta -l $shorder
    hold="COMPARI_fit-${s}_${i},$hold"
done 

sleep 1 && qsub -hold_jid "$hold" -N "COMPARI_post_${s}" -l short -l h_vmem=8g \
    -o $o/logs/\$JOB_NAME.\$JOB_ID.stdout \
    -e $o/logs/\$JOB_NAME.\$JOB_ID.stderr \
        $scriptsdir/wrapper.sh ./merg_COMPARI.py $s $(dirname $o) $mask $shorder
