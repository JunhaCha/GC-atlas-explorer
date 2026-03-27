setwd("/Users/junhacha/Documents/Playground/seurat-shiny-explorer")
shiny::runApp()


#clean metadata for all seurat objects
library(Seurat)

all <- readRDS('../seurat_merged_TME_malignant_final_umap_app_slim.rds')
malig <- readRDS('../seurat_cancercells_final_app_slim.rds')

cd8t <- readRDS('../seurat_CD8T_final2_app_slim.rds')
cd8t$final_celltype <- as.character(cd8t@active.ident)
cd8t$final_group <- paste(cd8t$rev_pathological_subtype, cd8t$rev_condition, sep = '_')

cd4t <- readRDS('../seurat_CD4T_final2_app_slim.rds')
cd4t$final_celltype <- as.character(cd4t@active.ident)
cd4t$final_group <- paste(cd4t$rev_pathological_subtype, cd4t$rev_condition, sep = '_')

stromal <- readRDS('../seurat_Stromal_final_app_slim.rds')
stromal$final_celltype <- as.character(stromal@active.ident)
stromal$final_group <- paste(stromal$rev_pathological_subtype, stromal$rev_condition, sep = '_')

epi <- readRDS('../seurat_epithelial_normal_final_final_app_slim.rds')
epi$final_celltype <- as.character(epi@active.ident)
epi$final_group <- paste(epi$rev_pathological_subtype, epi$rev_condition, sep = '_')

myeloid <- readRDS('../seurat_Mye_final2_app_slim.rds')
myeloid$final_celltype <- as.character(myeloid@active.ident)
myeloid$final_group <- paste(myeloid$rev_pathological_subtype, myeloid$rev_condition, sep = '_')


#filter to only relevant metadata
all@meta.data <- all@meta.data[,c('nCount_RNA', 'nFeature_RNA', 'percent.mito', 'nCount_SCT', 'nFeature_SCT', 'Phase','sampleID','patientID',
                                  'age','sex','rev_condition', 'rev_pathological_subtype', 'rev_molecular_subtype', 'rev_stage', 'final_celltype', 'final_group')]

malig@meta.data <- malig@meta.data[,c('nCount_RNA', 'nFeature_RNA', 'percent.mito', 'nCount_SCT', 'nFeature_SCT', 'Phase','sampleID','patientID',
                                      'age','sex','rev_condition', 'rev_pathological_subtype', 'rev_molecular_subtype', 'rev_stage')]

cd8t@meta.data <- cd8t@meta.data[,c('nCount_RNA', 'nFeature_RNA', 'percent.mito', 'nCount_SCT', 'nFeature_SCT', 'Phase','sampleID','patientID',
                                    'age','sex','rev_condition', 'rev_pathological_subtype', 'rev_molecular_subtype', 'rev_stage', 'final_celltype','final_group')]

cd4t@meta.data <- cd4t@meta.data[,c('nCount_RNA', 'nFeature_RNA', 'percent.mito', 'nCount_SCT', 'nFeature_SCT', 'Phase','sampleID','patientID',
                                    'age','sex','rev_condition', 'rev_pathological_subtype', 'rev_molecular_subtype', 'rev_stage', 'final_celltype', 'final_group')]

stromal@meta.data <- stromal@meta.data[,c('nCount_RNA', 'nFeature_RNA', 'percent.mito', 'nCount_SCT', 'nFeature_SCT', 'Phase','sampleID','patientID',
                                          'age','sex','rev_condition', 'rev_pathological_subtype', 'rev_molecular_subtype', 'rev_stage', 'final_celltype', 'final_group')]

epi@meta.data <- epi@meta.data[,c('nCount_RNA', 'nFeature_RNA', 'percent.mito', 'nCount_SCT', 'nFeature_SCT', 'Phase','sampleID','patientID',
                                  'age','sex','rev_condition', 'rev_pathological_subtype', 'rev_molecular_subtype', 'rev_stage', 'final_celltype', 'final_group')]

myeloid@meta.data <- myeloid@meta.data[,c('nCount_RNA', 'nFeature_RNA', 'percent.mito', 'nCount_SCT', 'nFeature_SCT', 'Phase','sampleID','patientID',
                                          'age','sex','rev_condition', 'rev_pathological_subtype', 'rev_molecular_subtype', 'rev_stage', 'final_celltype', 'final_group')]


saveRDS(all,'../seurat_merged_TME_malignant_final_umap_app_slim.rds')
saveRDS(malig,'../seurat_cancercells_final_app_slim.rds')
saveRDS(cd8t,'../seurat_CD8T_final2_app_slim.rds')
saveRDS(cd4t,'../seurat_CD4T_final2_app_slim.rds')
saveRDS(stromal,'../seurat_Stromal_final_app_slim.rds')
saveRDS(epi,'../seurat_epithelial_normal_final_final_app_slim.rds')
saveRDS(myeloid,'../seurat_Mye_final2_app_slim.rds')