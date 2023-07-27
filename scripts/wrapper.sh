#! /bin/bash
### embedded SGE parameters ###
#$ -S /bin/bash
#$ -cwd 
#$ -l hostname=!compute-fed1.bicic.local 
#$ -l hostname=!compute-fed2.bicic.local 
#$ -l hostname=!compute-fed3.bicic.local 
#$ -l hostname=!compute-fed4.bicic.local 
#$ -l hostname=!compute-fed5.bicic.local 

echo "Executing on: $(hostname)" | tee -a /dev/stderr
echo "Executing in: $(pwd)" | tee -a /dev/stderr
echo "Executing at: $(date)" | tee -a /dev/stderr
echo "Executing   : $0" | tee -a /dev/stderr 
echo "Arguments   : $*" | tee -a /dev/stderr

module load python/anaconda/DiCIPHR

echo python $* 

python $*
