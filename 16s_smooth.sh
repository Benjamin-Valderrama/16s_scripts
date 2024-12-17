#!/bin/bash

# TO-DO:
#
# 1. PARSE ARGUMENTS OF DADA2 IN THE DADA2.R SCRIPT

# Default values
current_wd="`pwd`"
study_folder=""
accession_number=""
run_all=false
run_download=false
run_dada2=false
run_picrust2=false
run_modules=false

# Function to display script usage
function display_usage() {
    echo "Usage: $0 -s|--study_folder STUDY_FOLDER [-r|--run_all] [-n|--accession_number ACCESSION_NUMBER]"
    echo "	[--run_download] [--run_bowtie] [--run_woltka] [--run_gbms] [-h|--help]"
    echo "Optional arguments:"
    echo "  -h, --help               Display this help message."
    echo ""
    echo "Required arguemnts:"
    echo "  -s, --study_folder       Specify the name for the study included in this meta-analysis. (Required)"
    echo "  -n, --accession_number   Specify the accession number of the raw data at the ENA (Only used if --run_all or --run_download are provided)."
    echo "  -a, --arguments          Path to file with arguemnts used on each software along the workflow"
    echo ""
    echo "Workflow arguments:"
    echo "  -r, --run_all            Run all steps of the workflow."
    echo "  --run_download           Run data download [uses fastq-dl]."
    echo "  --run_dada2              Run reads quality check and alignment [uses dada2 in R]."
    echo "  --run_picrust2           Run taxonomic and functional profilling [uses PICRUSt2]."
    echo "  --run_modules            Run module abundance and coverage calculation [uses OmixerRpm in R]."
    exit 1
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        -s|--study_folder)
            study_folder="$2"
            shift
            shift
            ;;
        -r|--run_all)
            run_all=true
            shift
            ;;
        -n|--accession_number)
            accession_number="$2"
            shift
            shift
            ;;
#	-a|--arguments)
#	    arguments="$2"
#	    shift
#	    shift
#	    ;;
        --run_download)
            run_download=true
            shift
            ;;
        --run_dada2)
            run_dada2=true
            shift
            ;;
        --run_picrust2)
            run_picrust2=true
            shift
            ;;
        --run_modules)
            run_modules=true
            shift
            ;;
        -h|--help)
            display_usage
            ;;
        *)
            echo "Unknown option: $1"
            display_usage
            ;;
    esac
done

# Check if required flag is provided
if [ -z "$study_folder" ]; then
    echo "Study folder is required."
    display_usage
fi

# -1. setting up the study folder
if [ ! -d ${study_folder} ]; then
	echo "PROGRESS -- Creating study folder : ${study_folder}"

	mkdir ${study_folder}
	mkdir ${study_folder}/00.rawdata
	mkdir ${study_folder}/nohups
else
	echo "PROGRESS -- The folder '${study_folder}' already exists. Moving to the next step ..."
fi

# 0. fastq-dl: download data from ENA.
if [ "$run_all" = true ] || [ "$run_download" = true ]; then
    if [ -z "$accession_number" ]; then
        echo "Accession number is required for downloading raw data."
        exit 1
    fi

    # DOWNLOAD DATA
    echo "PROGRESS -- Download raw data from ENA. Project accession number : ${accession_number}."
    bash /home/bvalderrama/scripts/smooth_analysis/16s_scripts/fastqdl.sh ${current_wd}/${study_folder} ${accession_number} > ${study_folder}/nohups/download.out
fi


# 1. dada2: preparing count table of taxas.
if [ "$run_all" = true ] || [ "$run_dada2" = true ]; then
    mkdir ${study_folder}/01.dada2_o
    mkdir ${study_folder}/02.picrust2
    mkdir ${study_folder}/02.picrust2/input

    # run dada2
    echo "PROGRESS -- Performing taxonomic profiling with DADA2."
    #source activate rbase41
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate rbase44
    Rscript /home/bvalderrama/scripts/smooth_analysis/16s_scripts/dada2.R ${current_wd}/${study_folder} > ${current_wd}/${study_folder}/nohups/dada2.out
    micromamba deactivate
    #conda deactivate

    # save parameter used for dada2
#    cp ${arguments} ${current_wd}/${study_folder}/nohups/workflow.args
fi


# 2. PICRUSt2: inference of the genomic content based on 16s.
if [ "$run_all" = true ] || [ "$run_picrust2" = true ]; then

    # RUN PICRUSt2 USING FILES PRODUCED IN THE PREVIOUS STEP
    echo "PROGRESS -- Performing functional inference with PICRUSt2"
#    eval "$(micromamba shell hook --shell bash)"
#    micromamba activate picrust2
    bash /home/bvalderrama/scripts/smooth_analysis/16s_scripts/PICRUSt2.sh ${current_wd}/${study_folder}
fi


# 3. OmixerRpm: calculate modules
if [ "$run_all" = true ] || [ "$run_modules" = true ]; then

    # CALCULATE THE MODULES USING THE FUNCTIONAL ANNOTATION GENERATED ABOVE
    echo "PROGRESS -- Calculating modules using the KO-based functional profiling."
    mkdir ${study_folder}/03.modules

    bash /home/bvalderrama/scripts/smooth_analysis/16s_scripts/run_modules.sh ${current_wd}/${study_folder}/02.picrust2/output/KO_metagenome_out ${current_wd}/${study_folder}/03.modules -m GBMs,GMMs > ${current_wd}/${study_folder}/nohups/omixer.out
fi

echo "PROGRESS -- WGS primary analysis finished."
