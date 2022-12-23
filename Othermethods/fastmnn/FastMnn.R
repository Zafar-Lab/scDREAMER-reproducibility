library(batchelor)
library(Seurat)
library(sceasy)
sobj = readRDS("/home/krushna/Documents/Data_integration/Othermethods/FastMnn/Human_Mouse_fastmnn.rds")
batch = 'batch'

sobj = RenameAssays(object = sobj, originalexp = 'RNA')


expr <- GetAssayData(sobj, slot = "data")
sce <- fastMNN(expr, batch = sobj@meta.data[[batch]])
corrected_data <- assay(sce, "reconstructed")

sobj <- SetAssayData(sobj, slot = "data", new.data = as.matrix(corrected_data))
sobj@reductions['X_emb'] <- CreateDimReducObject(reducedDim(sce, "corrected"), key = 'fastmnn_')

sceasy::convertFormat(sobj, from="seurat", to="anndata",
                      outFile='/home/krushna/Documents/Data_integration/Othermethods/FastMnn/Human_Mouse.h5ad')

