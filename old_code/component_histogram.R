#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

data <- read.table(args[1], row.names=1)

pdf(args[2])

hist(data[,1], main="Sizes of weakly connected components", xlab="size (bp)")

hist(data[,1], main='Sizes of WCCs - 50bp to 200bp', xlab="size (bp)", breaks="FD", xlim=c(50,200))

dev.off()