#!/bin/bash

#MOAB -l nodes=1:ppn=4
#MOAB -l walltime=3:00:00:00
#MOAB -l pmem=8gb
#MOAB -m a
#MOAB -N Mouse.Brain.RPF22.tRPF22
#MOAB -o R1.Mouse.Brain.RPF22.tRPF22.${MOAB_JOBID}.out
#MOAB -e R1.Mouse.Brain.RPF22.tRPF22.${MOAB_JOBID}.out


export pathScripts=/beegfs/home/hd/hd_hd/hd_fd141/scripts

#### Mapping
## RP
# sh RiboSeqPipelineMapping.sh : Species : Tissue : Library id : Data type : Replicate number : Number of CPUs
sh ${pathScripts}/RiboSeqPipelineMapping.sh Mouse Brain RPF22 RP R1 4

## TR
# sh RiboSeqPipelineMapping.sh : Species : Tissue : Library id : Data type : Replicate number : Number of CPUs
sh ${pathScripts}/RiboSeqPipelineMapping.sh Mouse Brain tRPF22 TR R1 4


#### Processing
## RP
# sh RiboSeqPipelineProcessing.sh : Species : Tissue : Library id : Data type : Replicate number : Number of CPUs
sh ${pathScripts}/RiboSeqPipelineProcessing.sh Mouse Brain RPF22 RP R1 4

## TR
# sh RiboSeqPipelineProcessing.sh : Species : Tissue : Library id : Data type : Replicate number : Number of CPUs
sh ${pathScripts}/RiboSeqPipelineProcessing.sh Mouse Brain tRPF22 TR R1 4
