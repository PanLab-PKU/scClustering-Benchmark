#################################
#                               #
#       Clustering Survey       #
#                               #
#       R 4.3.3 (2024-02-29)    #
#       env -surveyR43          #
#                               #
#       by Ziying Huang         #
#       03/10/2024              #
#                               #
#################################

# definition -------------------------------------------------------------------
rm(list=ls())

all_sets <- c(
    "E-MTAB-3321","Klein",
    "GSE189120",
    "GSE171993",
    "Zeisel","PBMC", #previous
    "fig11981034","GSE137537",
    "GSE109816","E-MTAB-6701/Placenta","E-MTAB-6701/Decidua",
    "GSE125970/Rectum","GSE125970/Colon","GSE125970/Ileum", #other tissue
    "Muris/FACS/Pancreas","Muris/FACS/Fat","Muris/drop/muscle","Muris/drop/lung","Sapiens/Tongue","Sapiens/Prostate", #downsampling
    "NHP",
    "NM/14","NM/19","NM/20","NM/29","NM/39","NM/41","NM/57","NM/58","NM/60","NM/63",
    "NG/A13","NG/A30","NG/E1","NG/E2","NG/E3",
    "Cao"
)

sets <- c("E-MTAB-6701/Placenta")

set.seed(9999)

dir_use <- "/home/ziyinghuang/project/survey/dataset/" 
dir_out <- "/home/ziyinghuang/project/survey/output/"

library(Seurat)
################################# Seurat ---------------------------------------
func <- "Seurat"
# function ---------------------------------------------------------------------


runSeurat <- function(sets,
                      min.cell = 1,
                      min.feature = 1,
                      nnmethod = "annoy",
                      annoy.met = "euclidean",
                      dim = 1:10,
                      met = "matrix",
                      algor = 1
){
    
    memory <- gc()
    errors <- list()
    
    for (set in sets) {
        
        time <- data.frame(Dataset = character(), 
                           user.self = numeric(), 
                           sys.self = numeric(), 
                           elapsed = numeric(), 
                           user.child = numeric(), 
                           sys.child = numeric(), 
                           stringsAsFactors = FALSE)
        
        wd <- paste0(dir_use,set)
        data <- Read10X(paste0(wd,"/"))
        out <- paste0(dir_out,set)
        
        if (!dir.exists(out)) {dir.create(out, recursive = T)}
        
        labels <- read.csv(paste0(wd,"/labels.csv"))
        num_labels <- length(levels(factor(labels[,2])))
        
        
        r_min <- 0.02
        r_max <- 2
        r <- 0.02
        clusters <- 1
        
        while (length(clusters) != num_labels) {
            
            if (length(clusters) > num_labels) {
                r_max <- r
            } else if (length(clusters) < num_labels) {
                r_min <- r 
            }
            
            r <- (r_min + r_max)/2
            
            
            pbmc <- CreateSeuratObject(counts = data, project = "test",
                                       min.cells = min.cell, min.features = min.feature)
            pbmc <- NormalizeData(pbmc)
            pbmc <- FindVariableFeatures(pbmc)
            all.genes <- rownames(pbmc)
            pbmc <- ScaleData(pbmc, features = all.genes)
            pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
            pbmc <- FindNeighbors(pbmc,
                                  nn.method = nnmethod,
                                  annoy.metric = annoy.met,
                                  dims = dim)
            pbmc <- FindClusters(pbmc,
                                 resolution = r,
                                 method = met,
                                 algorithm = algor)
            
            clusters <- unique(pbmc$seurat_clusters)
            
            print(paste0(set," resolution: ",r))
            print(paste0(set," origin labels: ", num_labels))
            print(paste0(set," labels: ",length(unique(clusters))))
            
            res_chose <- r
        }
        
        
        tryCatch({
            mem_before <- gc()
            rownames(mem_before) <- paste0(rownames(mem_before),"_pre")
            
            elapsed_time <- system.time({
                wd <- paste0(dir_use,set)
                data <- Read10X(paste0(wd,"/"))
                out <- paste0(dir_out,set)
                
                if (!dir.exists(out)) {dir.create(out, recursive = T)}
                
                labels <- read.csv(paste0(wd,"/labels.csv"))
                num_labels <- length(levels(factor(labels[,2])))
                
                pbmc <- CreateSeuratObject(counts = data, project = "test",
                                           min.cells = min.cell, min.features = min.feature)
                pbmc <- NormalizeData(pbmc)
                pbmc <- FindVariableFeatures(pbmc)
                all.genes <- rownames(pbmc)
                pbmc <- ScaleData(pbmc, features = all.genes)
                pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
                pbmc <- FindNeighbors(pbmc,
                                      nn.method = nnmethod,
                                      annoy.metric = annoy.met,
                                      dims = dim)
                pbmc <- FindClusters(pbmc,
                                     resolution = res_chose,
                                     method = met,
                                     algorithm = algor)
            })
            
            clusters <- pbmc$seurat_clusters
            embed <- pbmc@reductions$pca
            
            write.csv(clusters,paste0(out,"/cluster_",func,"_res_",res_chose,"_",Sys.Date(),".csv"))
            saveRDS(embed,paste0(out,"/embedding_",func,"_res_",res_chose,"_",Sys.Date(),".rds"))
            
            mem_after <- gc()
            rownames(mem_after) <- paste0(rownames(mem_after),"_post")
            mem <- rbind(mem_before,mem_after)
            rownames(mem) <- paste0(rownames(mem),"_",set)
            memory <- rbind(memory,mem)
            
            elapsed_time <- as.data.frame(elapsed_time)
            time <- data.frame(Dataset = set, 
                               user.self = elapsed_time[1,], 
                               sys.self = elapsed_time[2,], 
                               elapsed = elapsed_time[3,], 
                               user.child = elapsed_time[4,], 
                               sys.child = elapsed_time[5,])
            
            write.csv(time,paste0(out,"/time_",func,"_res_",res_chose,"_",Sys.Date(),".csv"))
            
            print(paste0(set," finished."))
            print(paste0(set," origin labels: ", num_labels))
            print(paste0(set," labels: ",length(unique(clusters))))
            print(paste0(set," resolution: ",res_chose))
            
        }, error = function(e) {
            errors[[length(errors) + 1]] <- list(set = set, error = conditionMessage(e))
        })
        
    }
    warnings()
    
    write.csv(memory,paste0(dir_out,"/memory_",func,"_",Sys.Date(),".csv"))
    print(errors)
    

}

runSeurat(sets)

















