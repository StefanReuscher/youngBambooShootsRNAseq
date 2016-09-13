#! /bin/bash

# checks if a $1 is set, if not quit
if [ -z "$1" ]; then
	echo "You must enter a tissue_replicate as first argument !"
	exit
fi

# create a copy of the jobscript template
cp ~/job_scripts/PH_refbase_transcriptome/jobscriptTemplate.sh ~/job_scripts/PH_refbase_transcriptome/PH_"$1".sh

# replace TISSUE with $1
sed -i -e "s/TISSUE/$1/g" ~/job_scripts/PH_refbase_transcriptome/PH_"$1".sh

# submit to queue
qsub -l short -l s_vmem=16G -l mem_req=16G ./job_scripts/PH_refbase_transcriptome/PH_"$1".sh
