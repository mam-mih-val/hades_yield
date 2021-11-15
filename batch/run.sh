#!/bin/bash

file_list=$1
output_dir=$2
pdg_code=$3

ownroot=/lustre/nyx/hades/user/mmamaev/install/root-6.18.04-centos7-cxx17/bin/thisroot.sh

current_dir=$(pwd)
partition=main
time=8:00:00
build_dir=/lustre/nyx/hades/user/mmamaev/hades_yield/build/src

lists_dir=${output_dir}/lists/
log_dir=${output_dir}/log/

mkdir -p $output_dir
mkdir -p $log_dir
mkdir -p $lists_dir

split -l 1 -d -a 3 --additional-suffix=.list "$file_list" $lists_dir

n_runs=$(ls $lists_dir/*.list | wc -l)

job_range=1-$n_runs

echo file list=$file_list
echo output_dir=$output_dir
echo log_dir=$log_dir
echo lists_dir=$lists_dir
echo n_runs=$n_runs
echo job_range=$job_range

sbatch -J Yield -p $partition -t $time -a $job_range -e ${log_dir}/%A_%a.e -o ${log_dir}/%A_%a.o --export=output_dir=$output_dir,file_list=$file_list,ownroot=$ownroot,lists_dir=$lists_dir,build_dir=$build_dir,pdg_code=$pdg_code -- /lustre/nyx/hades/user/mmamaev/hades_yield/batch/batch_run.sh
