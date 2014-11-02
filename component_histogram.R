#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

data <- read.table(args[1], row.names=1)

pdf(args[2])

hist(data, main="Sizes of weakly connected components", xlab="size (bp)")

hist(data, main='Sizes of WCCs - up to 100', xlab="size (bp)", breaks="FD")

dev.off()