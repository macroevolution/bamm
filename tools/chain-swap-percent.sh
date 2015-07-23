#!/bin/bash

# Print the chain swap percentage, accepted, and total for the given run
# and for each of the specified number of chains, delta T, and swap period

set -e    # Exit on any error
set -u    # Report unbound variables

if [ $# = 0 ]
then
    echo Prints the percentage of chain swaps accepted, the number accepted,
    echo and the number proposed \(space-delimited\) after a 20% burn-in.
    echo "Usage: chain-swap-percent.sh --numberOfChains <value1> [<value2> ...]"
    echo "                             --deltaT <value1> [<value2> ...]"
    echo "                             --swapPeriod <value1> [<value2> ...]"
    echo "                             --run <bamm_path> <bamm_args>"
    exit
fi

# Last option type for desired values
option=none

# Lists of desired values (empty at first)
numberOfChains=""
deltaT=""
swapPeriod=""
run=""

for arg in $@
do
    if [ $arg = "--numberOfChains" ]
    then
        option=numberOfChains
    elif [ $arg == "--deltaT" ]
    then
        option=deltaT
    elif [ $arg == "--swapPeriod" ]
    then
        option=swapPeriod
    elif [ $arg == "--run" ]
    then
        option=run
    else
        if [ $option = numberOfChains ]
        then
            numberOfChains="$numberOfChains $arg"
        elif [ $option = deltaT ]
        then
            deltaT="$deltaT $arg"
        elif [ $option = swapPeriod ]
        then
            swapPeriod="$swapPeriod $arg"
        elif [ $option = run ]
        then
            run="$run $arg"
        fi
    fi
done

# Basic error-checking
if [ -z "$numberOfChains" ]
then
    echo "ERROR: No --numberOfChains specified"
    exit
fi

if [ -z "$deltaT" ]
then
    echo "ERROR: No --deltaT specified"
    exit
fi

if [ -z "$swapPeriod" ]
then
    echo "ERROR: No --swapPeriod specified"
    exit
fi

if [ -z "$run" ]
then
    echo "ERROR: No --run specified"
    exit
fi

# Print formatted header
printf "%13s%13s%13s%13s%13s%13s\n" \
    nChains deltaT swapPeriod swapPercent swapAccepted swapProposed

# Run for every combination of number of chains, delta T, and swap period
for nChains in $numberOfChains
do
    for dT in $deltaT
    do
        for swapP in $swapPeriod
        do
            outName=chain-swap-percent-$RANDOM

            runInfoFileName=run_info.txt
            mcmcFileName=mcmc_out.txt
            eventDataFileName=event_data.txt
            chainSwapFileName=chain_swap.txt

            fileNameOptions="--runInfoFilename $runInfoFileName \
                --mcmcOutfile $mcmcFileName \
                --eventDataOutfile $eventDataFileName \
                --chainSwapFileName $chainSwapFileName"

            nChainsOption="--numberOfChains $nChains"
            dTOption="--deltaT $dT"
            swapPOption="--swapPeriod $swapP"

            $run --outName $outName $fileNameOptions \
                $nChainsOption $dTOption $swapPOption > /dev/null

            realSwapFileName=${outName}_$chainSwapFileName

            percent=$(Rscript -e "data = read.csv(\"$realSwapFileName\"); \
                                  data = tail(data, 0.8 * nrow(data)); \
                                  data = data[which(data$rank_1 == 1),]; \
                                  cat(mean(data\$swapAccepted))")

            accepted=$(Rscript -e "data = read.csv(\"$realSwapFileName\"); \
                                   data = tail(data, 0.8 * nrow(data)); \
                                   data = data[which(data$rank_1 == 1),]; \
                                   cat(sum(data\$swapAccepted))")

            proposed=$(Rscript -e "data = read.csv(\"$realSwapFileName\"); \
                                   data = tail(data, 0.8 * nrow(data)); \
                                   data = data[which(data$rank_1 == 1),]; \
                                   cat(nrow(data))")

            # Print current options
            printf "%13s%13s%13s%13s%13s%13s\n" \
                $nChains $dT $swapP $percent $accepted $proposed

            rm ${outName}*
        done
    done
done
