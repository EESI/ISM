#!/bin/bash
### tell SGE to use bash for this script
#$ -S /bin/bash
### execute the job from the current working directory, i.e. the directory in which the qsub command is given
#$ -cwd
### join both stdout and stderr into the same file
#$ -j y
### set email address for sending job status
#$ -M XXX@XXX.EDU
### project - basically, your research group name with "Grp" replaced by "Prj"
#$ -P rosenclassPrj
### select parallel environment, and number of job slots
#$ -pe shm 1
### request 15 min of wall clock time "h_rt" = "hard real time" (format is HH:MM:SS, or integer seconds)
#$ -l h_rt=48:00:00
### a hard limit 16 GB of memory per slot - if the job grows beyond this, the job is killed
#$ -l h_vmem=16G
### want node with at least 15 GB of free memory per slot
#$ -l m_mem_free=14G
### select the queue all.q, using hostgroup @intelhosts
#$ -q all.q

. /etc/profile.d/modules.sh
### These four modules must ALWAYS be loaded
module load shared
module load proteus
module load sge/univa
module load gcc
module load python/anaconda3

# before running the pipeline, run ```for f in $DIR/*.gz ; do gunzip -c "$f" > "${f%.*}" ; done``` to uncompress the files

# reference=${absolute path to the reference data (e.g., covid-19-reference.fasta)}
# input=${absolute path to the fasta data (e.g., gisaid.fasta)}
# output=${absolute path to the output data (e.g., gisaid_filtered.fasta)}



start=`date +%s`
python sequence_filtering.py $input $reference $output
end=`date +%s`
runtime=$((end-start))
echo "Filter time is $runtime"