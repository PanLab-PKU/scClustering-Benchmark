#################################
#                               #
#       Clustering Survey       #
#                               #
#       R 4.3.3 (2024-02-29)    #
#       env -surveyR43          #
#                               #
#       by Ziying Huang         #
#       21/10/2024              #
#                               #
#################################

# definition -------------------------------------------------------------------
rm(list=ls())

# sets <- c("E-MTAB-3321","Klein",
#           "GSE189120","GSE171993","Zeisel","PBMC", #previous
#           "fig11981034","GSE137537","GSE109816","E-MTAB-6701/Placenta","E-MTAB-6701/Decidua","GSE125970/Rectum","GSE125970/Colon","GSE125970/Ileum", #other tissue
#           "Muris/FACS/Pancreas","Muris/FACS/Fat","Muris/drop/muscle","Muris/drop/lung","Sapiens/Tongue","Sapiens/Prostate", #downsampling
#           "NHP","NM/14","NM/19","NM/20","NM/29","NM/39","NM/41","NM/57","NM/58","NM/60","NM/63",
#           "NG/A13","NG/A30","NG/E1","NG/E2","NG/E3","Cao")

# sets <- c("Klein","GSE189120","Zeisel",
#           "GSE109816","E-MTAB-6701/Placenta","E-MTAB-6701/Decidua",
#           "Muris/FACS/Pancreas","Muris/FACS/Fat","Muris/drop/muscle","Muris/drop/lung", #downsampling
#           "NHP",
#           "NM/39","NM/63",
#           "NG/A13","NG/E1","NG/E3")

sets <- c("E-MTAB-6701/Placenta")

set.seed(9999)

dir_use <- "/home/ziyinghuang/project/survey/dataset/" 
dir_out <- "/home/ziyinghuang/project/survey/output/"

library(Seurat)
################################# Monocle3 -------------------------------------
library(monocle3)
func <- "Monocle3"
# function ---------------------------------------------------------------------

runMonocle3 <- function(sets,
                        pp.method = "PCA",
                        scale = FALSE,
                        norm = "log",
                        rd.method = "UMAP",
                        umap.met = "cosine",
                        umap.nn = "annoy",
                        cls.method = "leiden"){
    
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
        
        wd <- paste0(dir_use,set)
        data <- Read10X(paste0(wd,"/"))
        out <- paste0(dir_out,set)
        
        if (!dir.exists(out)) {dir.create(out, recursive = T)}
        
        labels <- read.csv(paste0(wd,"/labels.csv"))
        num_labels <- length(levels(factor(labels[,2])))

        r_min <- 1e-40
        r_max <- 1e-20
        r <- 1e-40
        clusters <- 1
        
        gene_annotation <- data.frame(gene_short_name = rownames(data))
        rownames(gene_annotation) <- rownames(data)
        
        cds <- new_cell_data_set(as.matrix(data),
                                 gene_metadata = gene_annotation)
        cds <- preprocess_cds(cds,
                              method = pp.method,
                              scaling = scale,
                              norm_method = norm)
        cds <- reduce_dimension(cds,
                                reduction_method = rd.method,
                                preprocess_method = pp.method,
                                umap.metric = umap.met,
                                umap.nn_method = umap.nn)
        
        r <- r_min
        cds <- cluster_cells(cds,
                             reduction_method = rd.method,
                             cluster_method = cls.method,
                             resolution = r)
        
        clusters <- unique(cds@clusters$UMAP$clusters)
        print(paste0(set," min resolution: ",r))
        print(paste0(set," origin labels: ", num_labels))
        print(paste0(set," labels min: ",length(unique(clusters))))
        min <- length(unique(clusters))
        
        if (min == num_labels) {res_chose = r}
        
        r <- r_max
        cds <- cluster_cells(cds,
                             reduction_method = rd.method,
                             cluster_method = cls.method,
                             resolution = r)
        clusters <- unique(cds@clusters$UMAP$clusters)
        print(paste0(set," max resolution: ",r))
        print(paste0(set," origin labels: ", num_labels))
        print(paste0(set," labels max: ",length(unique(clusters))))
        max <- length(unique(clusters))
        
        if (max == num_labels) {res_chose = r}
        
        if (min < num_labels & max > num_labels) {
            
            r_min <- 1e-20
            r_max <- 1
            r <- 1e-20
            clusters <- 1
            
            while (length(clusters) != num_labels) {

                if (length(clusters) > num_labels) {
                    r_max <- r
                } else if (length(clusters) < num_labels) {
                    r_min <- r 
                }
                
                r <- (r_min + r_max) / 2
                
                cds <- cluster_cells(cds,
                                     reduction_method = rd.method,
                                     cluster_method = cls.method,
                                     resolution = r)
                
                clusters <- unique(cds@clusters$UMAP$clusters)
                print(paste0(set," resolution: ",r))
                print(paste0(set," origin labels: ", num_labels))
                print(paste0(set," labels: ",length(unique(clusters))))
                 
                res_chose = r
            }
            
            tryCatch({
                mem_before <- gc()
                rownames(mem_before) <- paste0(rownames(mem_before),"_pre")
                
                elapsed_time <- system.time({
                    gene_annotation <- data.frame(gene_short_name = rownames(data))
                    rownames(gene_annotation) <- rownames(data)
                    
                    cds <- new_cell_data_set(as.matrix(data),
                                             gene_metadata = gene_annotation)
                    cds <- preprocess_cds(cds,
                                          method = pp.method,
                                          scaling = scale,
                                          norm_method = norm)
                    cds <- reduce_dimension(cds,
                                            reduction_method = rd.method,
                                            preprocess_method = pp.method,
                                            umap.metric = umap.met,
                                            umap.nn_method = umap.nn)
                    cds <- cluster_cells(cds,
                                         reduction_method = rd.method,
                                         cluster_method = cls.method,
                                         resolution = res_chose)
                })
                
                clusters <- cds@clusters$UMAP$clusters
                embed <- cds@int_metadata$reduce_dim_metadata
                
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
        } else {
            print(paste0(set, " resolution out of bound."))
        }

    }
    warnings()
    
    write.csv(memory,paste0(dir_out,"/memory_",func,"_",Sys.Date(),".csv"))
    print(errors)
}

runMonocle3(sets)

















