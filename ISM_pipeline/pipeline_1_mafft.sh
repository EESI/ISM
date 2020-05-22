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
#$ -l h_rt=180:00:00
### a hard limit 16 GB of memory per slot - if the job grows beyond this, the job is killed
#$ -l h_vmem=62G
### want node with at least 15 GB of free memory per slot
#$ -l m_mem_free=61G
#$ -l vendor=intel
### select the queue all.q, using hostgroup @intelhosts
#$ -q long.q

. /etc/profile.d/modules.sh
### These four modules must ALWAYS be loaded
module load shared
module load proteus
module load sge/univa
module load gcc
module load mafft/7.453

# LOG=${absolute path to the log file (e.g., PROGRESS_2020.output)}
# input=${absolute path to the filtered data (e.g., gisaid_filtered.fasta)}
# output=${absolute path to the alignment output (e.g., mafft_2020.output)}

start=`date +%s`
mafft --thread 1 --progress $LOG --auto $input > $output
end=`date +%s`
runtime=$((end-start))
echo "Alignment time is $runtime"