#!/usr/bin/env Rscript

args<-commandArgs(TRUE)

if (length(args) == 0) {
    stop("Error, require params: normal_cells.matrix");
}

matrix.file = args[1]

pdf(paste0(matrix.file, '.dispersion_estimation.pdf'))

library(edgeR)
library(fitdistrplus)
library(infercnv)

expr.matrix = read.table(matrix.file)


## estimate dropout params
mean_vs_p0_table = infercnv:::.get_mean_vs_p0_from_matrix(expr.matrix)
logistic_params = infercnv:::.get_logistic_params(mean_vs_p0_table)

iterations=1
dispersion_params = c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)

resultset=matrix(0, ncol=3, nrow=iterations*length(dispersion_params))
colnames(resultset) = c('target', 'before_Zinf', 'after_Zinf')


row = 0


for (common.dispersion in dispersion_params) {
    message(sprintf("Exploring common.dispersion set at: %g", common.dispersion)) 
    for (iter in 1:iterations) {
        message(sprintf("\titer: %d", iter))
        
        row = row + 1

        ## simulate w/o zero-inflation
        sim_counts = infercnv:::.get_simulated_cell_matrix(mean_vs_p0_table$m, NULL, 100, common_dispersion=common.dispersion)

        ## estimate common disp from these data:
        design <- matrix(1, ncol(sim_counts), 1)


        disps <- edgeR::estimateDisp(sim_counts, design = design)
                                        #print(sprintf("estimated disp before dropouts: %g", disps$common.dispersion))

        resultset[row,1] <- common.dispersion
        resultset[row,2] <- disps$common.dispersion
        

        ## include zero-inflation
        sim_counts = infercnv:::.get_simulated_cell_matrix(mean_vs_p0_table$m, mean_vs_p0_table, 100,
                                                           common_dispersion=common.dispersion)
        
        
        disps <- edgeR::estimateDisp(sim_counts, design = design)
        resultset[row,3] <- disps$common.dispersion   

    }

    
}


resultset = as.data.frame(resultset)
print(resultset)
write.table(resultset, file=paste0(matrix.file, ".dispersion_estimation.dat"), quote=F, sep="\t")

## examples:
## 10x:  0.221 + 1.05 * (true_dispersion)  # colon single sample
##       0.223 + 1.05 * (true_dipersion)   # multiple colon samples

## smrtSeq: 0.95 + 1.56 * (true_dispersion)   # oligodendro
##          1.073 + 1.628 * (true_dispersion) # melanoma


res.lm = lm(resultset[,3] ~ resultset[,1])

print(res.lm)

coeff  = res.lm$coefficients
intercept = coeff[1]
slope = coeff[2]

plot(resultset[,1], resultset[,3], main=sprintf("y=%g + %g * x", intercept, slope), col='green')
points(resultset[,1], resultset[,2])


