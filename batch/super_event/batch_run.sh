#!/bin/bash

format='+%Y/%m/%d-%H:%M:%S'

date $format

job_num=$(($SLURM_ARRAY_TASK_ID))

filelist=$lists_dir/$(ls $lists_dir | sed "${job_num}q;d")

cd $output_dir
mkdir -p $job_num
cd $job_num

module load /cvmfs/vae.gsi.de/centos7/modules/linux-centos7-x86_64/gcc-8.1.0-gcc-4.8.5-oyp4lmr

echo "loading " $ownroot
source $ownroot

echo "executing $build_dir/super_event -i $filelist -t hades_analysis_tree -n -1 -o output.root --cuts-macro Hades/AuAu1.23.C --pdg-code $pdg_code --group-events 20"

$build_dir/super_event -i $filelist -t hades_analysis_tree -n -1 -o output.root --cuts-macro Hades/AuAu1.23.C --pdg-code $pdg_code --group-events 5

date $format
echo JOB FINISHED!