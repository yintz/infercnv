#!/usr/bin/env Rscript

library(HoneyBADGER)
library(infercnv)

data(gexp) ## tumor cells, dim:  [6082,75]
data(ref) ## reference, length: 6082



raw.data = cbind(gexp, data.frame('GTEX'=ref))

cell.annots = data.frame(cell=colnames(gexp), type='tumor')
cell.annots = rbind(cell.annots, data.frame(cell='GTEX', type='normal'))

write.table(raw.data, file="hb.example.matrix", quote=F, sep="\t")
write.table(cell.annots, file='hb.example.cell_annots', quote=F, sep="\t", col.names=F, row.names=F)

