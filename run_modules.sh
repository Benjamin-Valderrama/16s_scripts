#!/bin/bash

# Empty arrays to store elements and positional arguments
modules=()
positional_args=()

# Parse command line options and positional arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -m)
            shift
            IFS=',' read -ra modules <<< "$1"
            ;;
        *)
            positional_args+=("$1")
            ;;
    esac
    shift
done


# Extract the input and output positional arguments
input="${positional_args[0]}"
output="${positional_args[1]}"


#source /home/miniconda/miniconda3/bin/activate rbase41
eval "$(micromamba shell hook --shell bash)"
micromamba activate rbase44

# Iterate through the elements using a for loop and echo each element
for module in "${modules[@]}"; do

    # Calculate stratified modules
    if [ -e "$input/pred_metagenome_contrib.tsv.gz" ]; then
        echo "RUN STRATIFIED ANALYSIS"
        Rscript /home/bvalderrama/scripts/16s_scripts/stratified_omixer.R $input $output $module
    fi

    # Calculate unstratified modules
    if [ -e "$input/pred_metagenome_unstrat.tsv.gz" ]; then
        echo "RUN UNSTRATIFIED ANALYSIS"
        Rscript /home/bvalderrama/scripts/16s_scripts/unstratified_omixer.R $input $output $module
    fi

done
