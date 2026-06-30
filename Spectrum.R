#################################
#                               #
#       Clustering Survey       #
#                               #
#       R 4.3.3 (2024-02-29)    #
#       env -surveyR43          #
#                               #
#       by Ziying Huang         #
#       08/10/2024              #
#                               #
#################################

# definition -------------------------------------------------------------------
rm(list=ls())

# all_sets <- c(
#     "E-MTAB-3321","Klein",
#     "GSE189120",
#     "GSE171993",
#     "Zeisel","PBMC", #previous
#     "fig11981034","GSE137537",
#     "GSE109816","E-MTAB-6701/Placenta","E-MTAB-6701/Decidua",
#     "GSE125970/Rectum","GSE125970/Colon","GSE125970/Ileum", #other tissue
#     "Muris/FACS/Pancreas","Muris/FACS/Fat","Muris/drop/muscle","Muris/drop/lung","Sapiens/Tongue","Sapiens/Prostate", #downsampling
#     "NHP",
#     "NM/14","NM/19","NM/20","NM/29","NM/39","NM/41","NM/57","NM/58","NM/60","NM/63",
#     "NG/A13","NG/A30","NG/E1","NG/E2","NG/E3",
#     "Cao"
# )


# all_sets <- c(
#     "GSE189120",
#     "GSE137537",
#     "E-MTAB-6701/Placenta","E-MTAB-6701/Decidua",
#     "NM/29",
#     "NG/A13",
#     "NG/A30"
# )

all_sets <- c("E-MTAB-6701/Placenta")

set.seed(9999)

dir_use <- "/home/ziyinghuang/project/survey/dataset/" 
dir_out <- "/home/ziyinghuang/project/survey/output/"

library(Seurat)
################################# Spectrum -------------------------------------
library(Spectrum)
func <- "Spectrum"
# function ---------------------------------------------------------------------
runSpectrum <- function(sets,
                        large = FALSE) {
    
    time <- data.frame(Dataset = character(), 
                       user.self = numeric(), 
                       sys.self = numeric(), 
                       elapsed = numeric(), 
                       user.child = numeric(), 
                       sys.child = numeric(), 
                       stringsAsFactors = FALSE)
    
    memory <- gc()
    errors <- list()
    
    
    for (set in sets) {
        

        tryCatch({
            mem_before <- gc()
            rownames(mem_before) <- paste0(rownames(mem_before),"_pre")
            
            elapsed_time <- system.time({
                wd <- paste0(dir_use,set)
                data <- Read10X(paste0(wd,"/"))
                out <- paste0(dir_out,set)
                
                clusters <- c()
                
                if (!dir.exists(out)) {dir.create(out, recursive = T)}
                
                labels <- read.csv(paste0(wd,"/labels.csv"))
                num_labels <- length(levels(factor(labels[,2])))
                
                pbmc <- CreateSeuratObject(counts = data, project = "test",
                                           min.cells = 3, min.features = 200)
                pbmc <- NormalizeData(pbmc)
                
                data_norm <- as.matrix(GetAssayData(pbmc, layer = "data"))
                
                print(paste0(set, " dim: ",dim(data)))
                print(paste0(set, " normalized dim: ",dim(data_norm)))
                
                large <- TRUE
                
                if (!large) {
                    res <- Spectrum(data_norm,
                                    method = 3, fixk = num_labels)
                    clusters <- res$assignments

                } else {
                    cells <- dim(pbmc)[2]
                    res <- Spectrum(data_norm,
                                    method = 3, fixk = num_labels,
                                    FASP=TRUE,FASPk=round(cells*1/10))
                    clusters <- res$allsample_assignments
                }
            })
            
            embed <- res$eigensystem
            
            write.csv(clusters,paste0(out,"/cluster_",func,"_",Sys.Date(),".csv"))
            saveRDS(embed,paste0(out,"/embedding_",func,"_",Sys.Date(),".rds"))
            
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
            
            write.csv(time,paste0(out,"/time_",func,"_",Sys.Date(),".csv"))
            
            print(paste0(set," finished."))
            print(paste0(set," origin labels: ", num_labels))
            print(paste0(set," labels: ",length(unique(clusters))))
            
        }, error = function(e) {
            errors[[length(errors) + 1]] <- list(set = set, error = conditionMessage(e))
        })

    }
    warnings()
    
    write.csv(memory,paste0(dir_out,"/memory_",func,"_",Sys.Date(),".csv"))
    print(errors)
    return(clusters)
}

library(R.utils)
timeout_limit <- 24*60*60

for (sets in all_sets) {
    withTimeout({
        clusters <- runSpectrum(sets)
    }, timeout = timeout_limit)
    
    if (is.null(clusters)) {
        withTimeout({
            runSpectrum(sets, large = TRUE)
        }, timeout = timeout_limit)
    }
}































