#!/bin/bash

# Print the chain swap percentage for the given run
# e.g., chain-swap-percent ./bamm -c divcontrol.txt --numberOfGenerations 1000
# Can use all BAMM options except --outName and any option
#     that changes the name of an output file

set -e    # Exit on any error
set -u    # Report unbound variables

if [ $# = 0 ]
then
    echo Usage: chain-swap-percent [bamm_path] [bamm_options]
fi

outName=chain-swap-percent-$RANDOM

runInfoFileName=run_info.txt
mcmcFileName=mcmc_out.txt
eventDataFileName=event_data.txt
chainSwapFileName=chain_swap.txt

fileNameOptions="--runInfoFilename $runInfoFileName \
    --mcmcOutfile $mcmcFileName \
    --eventDataOutfile $eventDataFileName \
    --chainSwapFileName $chainSwapFileName"

echo 'Running BAMM to calculate chain swap percentage...'
$@ --outName $outName $fileNameOptions > /dev/null

realChainSwapFileName=${outName}_$chainSwapFileName

Rscript -e "cat(mean(read.csv(\"$realChainSwapFileName\")\$swapAccepted)); \
    cat(\"\\\n\")"

rm $outName*
