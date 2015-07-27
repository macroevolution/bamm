# Function to calculate percent of chain swaps that are successful,
# given different values of deltaT, swapPeriod and number of chains.
# If multiple values are provided, all combinations of parameter values
# will be run. The results are returned in a table.

# Navigate to appropriate directory in R.
# The call to terminal will operate in the directory set by R.

# bammPath = the path to the bamm binary;
#     leave as 'bamm' if the program is installed in the bin folder
# controlfile = name of the controlfile
# nChains = number of chains, can provide a vector of values
# deltaT = can provide a vector of values for deltaT
# swapPeriod = can provide a vector of values for swapPeriod
# nGenerations = number of generations used for evaluation
# burnin = proportion of generations to discard, default is 20%
# deleteTempFiles = boolean, should results files be automatically removed;
#     default is TRUE

chainSwapPercent <- function(bammPath = 'bamm', controlfile, nChains = 4,
    deltaT, swapPeriod = 1000, nGenerations = 10000, burnin = 0.2,
    deleteTempFiles = TRUE)
{
    if (Sys.which(bammPath) == '') {
        stop('BAMM executable not found. Please check bammPath.')
    }
    
    if (burnin >= 1) {
        stop("Burn-in should be between 0 and 1.")
    }

    if (length(nGenerations) > 1) {
        stop(paste('You may only provide multiple values for nChains,',
             'deltaT and swapPeriod.'))
    }

    if (!is.logical(deleteTempFiles)) {
        stop('deleteTempFiles must be either TRUE or FALSE.')
    }

    if (!controlfile %in% list.files()) {
        stop('Control file not found in working directory.')
    }
    
    res <- data.frame()
    pb <- txtProgressBar(min = 0,
        max = length(nChains) * length(deltaT) * length(swapPeriod), style = 3)
    counter <- 0 #to add increment to progress bar
    
    for (i in 1:length(nChains)) {
        for (j in 1:length(deltaT)) {
            for (k in 1:length(swapPeriod)) {
                id <- paste('TMP', round(runif(1)*100000), sep='_')
                call <- paste(bammPath, '-c', controlfile,
                    '--numberOfGenerations', as.integer(nGenerations),
                    '--numberOfChains', nChains[i],
                    '--swapPeriod', swapPeriod[k],
                    '--deltaT', deltaT[j],
                    '--outName', id,
                    '--simulatePriorShifts', '0', sep=' ')
                if (.Platform$OS.type != 'windows') {
                    system(command = call, intern = FALSE, wait = TRUE,
                        ignore.stdout=TRUE, ignore.stderr=TRUE)
                }
                if (.Platform$OS.type == 'windows') {
                    system(command = call, intern = FALSE, wait = TRUE,
                        ignore.stdout=TRUE, ignore.stderr=TRUE,
                        invisible = TRUE, show.output.on.console = FALSE)
                }
                data <- read.csv(
                    paste(id, 'chain_swap.txt', sep='_'))
                data <- data[which(data$rank_1 == 1), 'swapAccepted']
                data <- tail(data, (1 - burnin) * length(data))
                resVec <- c(nChains[i], deltaT[j], swapPeriod[k],
                    mean(data), sum(data), length(data))
                res <- rbind(res, resVec)
                if (deleteTempFiles) {
                	file.remove(list.files(pattern = id))
                }
                counter <- counter + 1
                setTxtProgressBar(pb, counter)
            }
        }
    }

    cat('\n')
    colnames(res) <-
        c('nChains', 'deltaT', 'swapP', 'percent', 'accepted', 'proposed')

    return(res)
}
