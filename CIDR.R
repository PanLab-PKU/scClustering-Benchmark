#################################
#                               #
#       Clustering Survey       #
#                               #
#       R 4.3.3 (2024-02-29)    #
#       env -surveyR43          #
#                               #
#       by Ziying Huang         #
#       24/09/2024              #
#                               #
#################################

# definition -------------------------------------------------------------------
rm(list=ls())

# sets <- c("Cao","sim-DENOVO","sim-PBMC")

# sets <- c(
#     # "E-MTAB-3321","Klein",
#     # "GSE189120","GSE171993","Zeisel","PBMC", #previous
#     # "fig11981034","GSE137537","GSE109816",
#     # "E-MTAB-6701/Placenta","E-MTAB-6701/Decidua","GSE125970/Rectum","GSE125970/Colon","GSE125970/Ileum", #other tissue
#     # "Muris/FACS/Pancreas","Muris/FACS/Fat","Muris/drop/muscle","Muris/drop/lung","Sapiens/Tongue","Sapiens/Prostate", #downsampling
#     # "NHP","NM/14","NM/19","NM/20","NM/29","NM/39","NM/41","NM/57","NM/58","NM/60","NM/63",
#     # "NG/A13","NG/A30","NG/E1","NG/E2","NG/E3"
# )

sets <- c("E-MTAB-6701/Placenta")

set.seed(9999)

dir_use <- "/home/ziyinghuang/project/survey/dataset/" 
dir_out <- "/home/ziyinghuang/project/survey/output/"

library(Seurat)
################################# CIDR -----------------------------------------
library(cidr)
func <- "CIDR"
# function ---------------------------------------------------------------------
runCIDR <- function(sets) {
    
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
                
                if (!dir.exists(out)) {dir.create(out, recursive = T)}
                
                labels <- read.csv(paste0(wd,"/labels.csv"))
                num_labels <- length(levels(factor(labels[,2])))
                
                cidr <- scDataConstructor(as.matrix(data))
                cidr <- determineDropoutCandidates(cidr)
                cidr <- wThreshold(cidr)
                cidr <- scDissim(cidr)
                cidr <- scPCA(cidr, plotPC = FALSE)
                cidr <- nPC(cidr)
                cidr <- scCluster(cidr, nCluster = num_labels)
                
                clusters <- cidr@clusters
                embed <- cidr@PC
            })
            
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

runCIDR(sets)





