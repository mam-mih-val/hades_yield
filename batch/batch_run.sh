#!/bin/bash

format='+%Y/%m/%d-%H:%M:%S'

date $format

job_num=$(($SLURM_ARRAY_TASK_ID))

filelist=$lists_dir/$(ls $lists_dir | sed "${job_num}q;d")

cd $output_dir
mkdir -p $job_num
cd $job_num

while read line; do
    echo $line >> list.txt
done < $filelist
echo >> list.txt

module load /cvmfs/vae.gsi.de/centos7/modules/linux-centos7-x86_64/gcc-8.1.0-gcc-4.8.5-oyp4lmr

echo "loading " $ownroot
source $ownroot

echo "executing $build_dir/yield -i list.txt -t hades_analysis_tree -n -1 -o yield.root --cuts-macro Hades/AuAu1.23.C --pdg-code $pdg_code --acceptance-protons=/lustre/nyx/hades/user/mmamaev/hades_yield/src/acceptance_2212_pT_theta.root"

$build_dir/yield -i list.txt -t hades_analysis_tree -n -1 -o yield.root --cuts-macro Hades/AuAu1.23.C --pdg-code $pdg_code --acceptance-protons=/lustre/nyx/hades/user/mmamaev/hades_yield/src/acceptance_2212_pT_theta.root

date $format
echo JOB FINISHED!