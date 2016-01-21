# DO NOT USE THIS CODE unless you understand its problems.
# See the README.md file for details.

# Copyright 2011-2016 Martin and Melanie von Gagern
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(random.pivot)                   # https://github.com/gagern/random_pivot
library(vegan)                          # designdist, commsimulator
library(plyr)                           # ddply, rename

# Does a complete computation using given data and counts.
# If you want to try the effect of different bin counts, or different choices
# for C-Scores zero or one, re-using the same null model will save time.
# In that case, this function can be used as a guideline.
RePairs <- function(filename,           # will be used to store randomized data
                    species,            # a factor describing species obs.
                    sites,              # a factor with the corresponding sites
                    num.rand,           # the number of randomized matrices
                    num.bins)           # the number of bins to use
{
    rrcs.to.file(filename, species, sites, numrand)
    C.R <- rrcs.from.file(filename)
    return(calc.sig.pairs(C.R, nBins=num.bins, zeroBin=T, oneBin=T))
}

# Compute significant pairs using a classification of C-Scores into bins.
#
# Input:
#  - pairs: a data.frame, describing original data and randomizations.
#           Each row corresponds to one of the (#species choose 2)
#           possible species pairs, with the columns as follows:
#    * pairs$i: the index (level number) of the first species of the pair
#    * pairs$j: the index of the second species of the pair
#    * pairs$orig.cscore: the original Checkerboard Score of this pair
#    * pairs$cscore:      a matrix, where each column corresponds to
#                         one randomized null model
#  - nBins: the number of bins to use for classification of the C-Scores
#  - zeroBin: whether to use a dedicated bin for cscore == 0
#  - oneBin: whether to use a dedicated bin for cscore == 1
#  - pairCL: significance level used to compute the upper and lower confidence
#            limit for individual pairs, will be split 50:50 for two-sided test
#  - binCL: significance level used to compute the upper confidence limit
#           for each bin, which in turn is used to choose pairs from it
#
# Output: a list containing the following named elements:
#  - data.seg: data for the segregated pairs, like input pairs but also with
#    * data.seg$bin: the bin used for this pair in the actual observation
#    * data.seg$mean: average null model C-Score
#    * data.seg$sd: null model C-Score standard deviation
#    * data.seg$uppercl: upper pairCL/2 quantile over all null models
#    * data.seg$lowercl: lower pairCL/2 quantile over all null models
#    * data.seg$ord: index of pair in bin if ordered by distance from mean
#  - data.agg: data for the aggregated pairs, like data.seg
#  - nBins: number of bins, just for reference
#  - ocsb: number of pairs in each bin for the original C-Scores
#  - sig.bins: indices of bins with counts significantly above expected
#  - bin.means: average number of pairs in each bin
#  - bin.uppercl: upper cl for the number of pairs in each bin
calc.sig.pairs <- function(pairs, nBins,
                           zeroBin = TRUE, oneBin = TRUE,
                           pairCL = 0.05, binCL = pairCL/2) {

    # helper function to assign given C-Scores to bins
    scale <- nBins - zeroBin - oneBin
    inBins <- function(score) {
        bin <- ceiling(score * scale); # Pairs 1.1 uses floor here
        if (zeroBin) {
            bin <- bin + 1
        }
        else {
            bin[bin == 0] <- 1
        }
        if (oneBin) {
            bin[score >= 1] <- nBins
        }
        return(bin)
    }

    # classify pairs by observed C-Scores into bins
    orig.bins <- inBins(pairs$orig.cscores)

    # bins: one row per bin, collecting various data
    # bins$bin: index identifying this bin
    bins <- data.frame(bin = 1:nBins)
    # bins$obs: number of actually observed pairs in each bin
    bins$obs <- tabulate(orig.bins, nBins)
    # bins$runs: matrix of pair counts for randomized null models
    bins$runs <- apply(inBins(pairs$cscores), 2, tabulate, nBins)
    # bins$exp: expected number of pairs, mean over all null models
    bins$exp <- apply(bins$runs, 1, mean)
    # bins$upperCL: upper binCL quantile of observed pair count spectrum
    bins$upperCL <- apply(bins$runs, 1, atRelPos, 1 - binCL)
    # bins$surplus: number of actually observed pairs in excess of upperCL
    bins$surplus <- pmax(0, bins$obs - bins$upperCL)

    # data: Like pairs, but with some additional variables
    reservedNames <- c("bin", "mean", "sd", "uppercl", "lowercl",
                       "cscore", "cscores", "ord")
    data <- cbind(pairs[setdiff(names(pairs), reservedNames)],
                  bin = orig.bins,
                  mean = apply(pairs$cscores, 1, mean),
                  sd = apply(pairs$cscores, 1, sd),
                  uppercl = apply(pairs$cscores, 1, atRelPos, 1 - pairCL/2),
                  lowercl = apply(pairs$cscores, 1, atRelPos, pairCL/2)
                 )
    data <- transform(data, ord = (orig.cscores - mean)^2)

    # surplus: those pairs in each bin which are above the CL of that
    # bin, and which are most significant in that bin according to ord.
    surplus <- ddply(data, .(bin), function(x) {
        # The 2010 paper suggests ordering by distance from CL, and other
        # places can be interpreted as suggesting an ordering directly by
        # observed cscores. We order by distance from mean here.
        idxs <- order(-x$ord, x$i, x$j)
        idxs <- idxs[seq.int(length.out = bins$surplus[x$bin[1]])]
        x[idxs, ]
    })
    surplus <- rename(surplus, c("orig.cscores" = "cscore"))
    surplus <- surplus[setdiff(names(surplus), "cscores")]
    agg <- subset(surplus, cscore < lowercl)
    seg <- subset(surplus, cscore > uppercl)
    return(list(data.seg = seg,
                data.agg = agg,
                nBins = nBins,
                ocsb = bins$obs,
                sig.bins = which(bins$surplus > 0),
                bin.means = bins$exp,
                bin.uppercl = bins$upperCL))
}

# RRCS: Repeated Randomizes C-Score computation, writing result to file.
# Generate a number of null models and incrementally write them to a file.
# Writing them to file can help to keep the memory load low during this
# CPU-intensive task, so that the machine can still be used for e.g. office
# work.  The data can be loaded into memory later on for a shorter time.
# The file is written in a binary format and gzip-compressed to reduce size.
# The data contains no floating point numbers, so the content is exact,
# and not subject to any rounding decisions or float computation precision.
rrcs.to.file <- function(fileName, rows, cols, repeats) {
    rows <- as.factor(rows)		# coerce to factor unless it already is
    cols <- as.factor(cols)		# coerce to factor unless it already is
    f <- gzfile(paste0(fileName, ".gz"), "wb")
    writeBin("rrcs2file-v2", f)         # magic header to identify these files
    writeBinFactor(rows, f)             # write rows and cols in order, so the
    writeBinFactor(cols, f)             #  original matrix can be restored.
    m <- pivot(rows, cols, TRUE, FALSE) # build a logical matrix
    x <- tcrossprod(m)                  # count co-occurrences of species pairs
    x <- as.dist(x)                     # discard one triangle and diagonal
    x <- as.integer(x)                  # linearize the remaining counts
    writeBin(x, f, size=2)              #  and write them to the file
    writeBin(as.integer(repeats), f, size=2) # write number of repetitions
    start <- proc.time()                # time used for progress feedback
    rep <- 0
    while (rep < repeats) {
        rep <- rep + 1
        cat("Repetition ", rep, " of ", repeats, " .", sep="")
        m <- random.pivot(rows, cols)   # initial randomization
        cat(".")
        m <- random.swap(m)             # swaps for better distribution
        cat(". ")
        x <- tcrossprod(m)              # count and save co-occurrences as above
        x <- as.dist(x)
        x <- as.integer(x)
        writeBin(x, f, size=2)
        now <- proc.time()
        cat(fmttime(now[3]-start[3]), " so far, estimate another ",
            fmttime((now[3]-start[3])/rep*(repeats-rep)), " to completion\n",
            sep="")
    }
    close(f)                            # done, close the output file
}

rrcs.from.file <- function(fileName) {
    f <- gzfile(paste0(fileName, ".gz"), "rb")
    format <- readBin(f, character(0))
    facsize <- 4                        # we use 32-bit numbers for factors
    if (format == "rrcs2file-v1") {     # we had a different version at first
        facsize <- 2                    # which used 16-bit numbers instead
        format <- "rrcs2file-v2"        # but was the same otherwise
    }
    stopifnot(format == "rrcs2file-v2") # anything else is not our format
    rows <- readBinFactor(f, facsize)
    cols <- readBinFactor(f, facsize)
    n <- nlevels(rows)
    baseFactor <- levels(rows)
    baseFactor <- factor(baseFactor, baseFactor)
    ones <- as.integer(table(rows))     # number of observations per species
    m <- matrix(1:n, nrow=n, ncol=n)    # compute pair members in a way which is
    pairs <- data.frame(i=as.integer(as.dist(t(m))), #  compatible with as.dist
                        j=as.integer(as.dist(m)))
    pairs$fi <- baseFactor[pairs$i]     # turn indices into factor
    pairs$fj <- baseFactor[pairs$j]
    pairs$ki <- ones[pairs$i]           # turn indices into observation counts
    pairs$kj <- ones[pairs$j]
    npairs <- nrow(pairs)               # should be n choose 2
    denom <- (pairs$ki*pairs$kj)        # denominator used for normalization
    k <- readBin(f, integer(0), npairs, size=2) # actual co-occurrence counts
    pairs$orig.cscores <- (pairs$ki - k)*(pairs$kj - k)/denom
    repeats <- readBin(f, integer(0), size=2)           # number of null models
    k <- readBin(f, integer(0), npairs*repeats, size=2) # simulated counts
    pairs$cscores <- matrix((pairs$ki - k)*(pairs$kj - k)/denom,
                            nrow=npairs, ncol=repeats)
    stopifnot(length(readBin(f, integer(0), size=1)) == 0) # expect end of file
    close(f)
    return(pairs)
}

# Write a factor to a binary file
writeBinFactor <- function(fac, con, size=4) {
    writeBin(as.integer(nlevels(fac)), con, size=size)
    writeBin(levels(fac), con)
    writeBin(as.integer(length(fac)), con, size=size)
    writeBin(as.integer(fac), con, size=size)
}

# Reverse the write process above
readBinFactor <- function(con, size=4) {
    nl <- readBin(con, integer(0), size=size)
    lvls <- readBin(con, character(0), nl)
    raw <- factor(lvls, lvls)
    len <- readBin(con, integer(0), size=size)
    iv <- readBin(con, integer(0), len, size=size)
    return (raw[iv])
}

# Turn a pair of factors with equal length into a single matrix.
# In other words, convert from long to wide format.

# values may be a vector of the same length (i.e. a third column of
# the long format) or a single value to be used for all matching entries.
# fillwith will be used to fill any matrix cells not defined by the data.
pivot <- function(rows, cols, values, fillwith=NA) {
    rf <- as.factor(rows)
    cf <- as.factor(cols)
    rl <- levels(rf)
    cl <- levels(cf)
    m <- matrix(fillwith, length(rl), length(cl))
    m[matrix(c(as.numeric(rf), as.numeric(cf)), ncol=2)] <- values
    rownames(m) <- rl
    colnames(m) <- cl
    return(m)
}

# Perform random swaps in the logical matrix m.
# Uses functions from the vegan library.
random.swap <- function(m, repeatScale=1) {
    nr <- nrow(m)
    nc <- ncol(m)
    # By default, the number of trials is chosen such that the expected
    # number of successful swaps equals the size of the matrix. The size
    # of the matrix equals nr*nc, whereas the number of trials per
    # successful swap is estimated by the number of potential swap
    # locations, i.e. choose(nr,2)*choose(nc,2), divided by the number
    # of actually possible swaps in the original matrix. The latter is
    # computed using designdist. As this estimate is based on the
    # original input matrix, it might be biased if that matrix has an
    # extremely low or high number of checkerboards.
    initialCheckerboards <- sum(designdist(m, "(A-J)*(B-J)", "binary"))
    trials <- nr*nr*(nr - 1)*nc*nc*(nc - 1)/4/initialCheckerboards
    trials <- round(repeatScale*trials)
    m <- commsimulator(m, method="tswap", thin=trials)
    return (m != 0)                     # Turn result back into a logical matrix
}

# Compute quantiles of a randomized distribution.
#
# Input:
#  - x:     a Vector of numerical data
#  - level: a relative position between 0 and 1; may also be a vector
#
# Output:
#  - those entries of x which are located at the positions identified by level,
#    if the vector has been sorted first. 0 is the smallest element, 1 the
#    largest, with everything in between interpolated linearly.
#
# Pairs 1.1 uses a reverse order here, with some asymmetry due to
# uncompensated 1-based indexing.
atRelPos <- function(x, level) {
    p <- round((length(x) - 1)*level + 1)
    return (sort(x, partial = p)[p])
}

# Format a time, given in seconds
fmttime <- function(t) {
    t <- round(t)
    return(sprintf("%d:%02d", t %/% 60, t %% 60))
}

# Try the above code on some sample data (not provided), to see whether
# everything works and there are no calls to undefined functions or similar.
demo.run <- function() {
    demo.nRep <- 5000
    demo.nBins <- 22
    demo.data <- read.table("demo.txt", header=TRUE, sep="\t")
    demo.species <- as.factor(demo.data[[1]])
    demo.sites <- as.factor(demo.data[[2]])
    rrcs.to.file("demo", demo.species, demo.sites, demo.nRep)
    demo.rrcs <- rrcs.from.file("demo")
    demo.sig <- calc.sig.pairs(demo.rrcs, demo.nBins, zeroBin=T, oneBin=T)
    invisible(demo.sig)                 # invisible in case this result is BIG
}
