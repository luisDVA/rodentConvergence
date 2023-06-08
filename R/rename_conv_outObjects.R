# give names to convergent regime model objects
library(fs)
library(glue)
library(purrr)
library(l1ou)

# PC1cr 
# relevant RDS files
PC1convRDSs <- dir_ls("2021_conv/out/",regexp = "PC1_conv_regs_cr")
PC1_conv_regs_cr <- flatten(map(PC1convRDSs,readRDS))
PC1_conv_regs_cr
names(PC1_conv_regs_cr)

# 12 to 50
remtree_inds <- c(34,78,40,88,61,69,77,47,100,71,55,42,17,7,16,82,86,87,99,
                  65,36,74,80,68,49,70,48,62,15,10,52,33,31,94,75,12,5,96,29)
# 51 to 100
trees_ready <- c(9,11,18,19,25,38,43,60,67,72,89,
                 34,78,40,88,61,69,77,47,100,71,55,42,
                 17,7,16,82,86,87,99,65,36,74,80,68,
                 49,70,48,62,15,10,52,33,31,94,75,12,5,96,29)
pending_trees <- which(1:100 %in% trees_ready)
pending_trees

names(PC1_conv_regs_cr)[12:50] <- glue("2021_conv/out/emodPC1_cr{tree}_bic_pca.rds",tree=remtree_inds)
names(PC1_conv_regs_cr)[51:100] <- glue("2021_conv/out/emodPC1_cr{tree}_bic_pca.rds",tree=pending_trees)
names(PC1_conv_regs_cr)

PC1_conv_regs_cr
saveRDS(object = PC1_conv_regs_cr,file = "2021_conv/out/processed_l1ouModels/PC1_convRegs_crAll.rds")

# PC2cr 
# relevant RDS files
PC2convRDSs <- dir_ls("2021_conv/out/",regexp = "PC2_conv_regs_cr")
PC2_conv_regs_cr <- flatten(map(PC2convRDSs,readRDS))
names(PC2_conv_regs_cr)
# repeat as for PC1
names(PC2_conv_regs_cr)[12:50] <- glue("2021_conv/out/emodPC2_cr{tree}_bic_pca.rds",tree=remtree_inds)
names(PC2_conv_regs_cr)[51:100] <- glue("2021_conv/out/emodPC2_cr{tree}_bic_pca.rds",tree=pending_trees)
names(PC2_conv_regs_cr)
PC2_conv_regs_cr
saveRDS(object = PC2_conv_regs_cr,file = "2021_conv/out/processed_l1ouModels/PC2_convRegs_crAll.rds")

# PC3cr 
# relevant RDS files
PC3convRDSs <- dir_ls("2021_conv/out/",regexp = "PC3_conv_regs_cr")
# object with all trees
PC3_conv_regs_cr <- readRDS(PC3convRDSs[4])
#PC3_conv_regs_cr <- flatten(map(PC3convRDSs,readRDS))
names(PC3_conv_regs_cr)
# repeat as for PC1
names(PC3_conv_regs_cr)[12:50] <- glue("2021_conv/out/emodPC3_cr{tree}_bic_pca.rds",tree=remtree_inds)
names(PC3_conv_regs_cr)[51:100] <- glue("2021_conv/out/emodPC3_cr{tree}_bic_pca.rds",tree=pending_trees)
names(PC3_conv_regs_cr)
PC3_conv_regs_cr
saveRDS(object = PC3_conv_regs_cr,file = "2021_conv/out/processed_l1ouModels/PC3_convRegs_crAll.rds")

# Multivariate cr
MultconvRDSs <- dir_ls("2021_conv/out/",regexp = "mult_conv")
MultconvRDSs
conv_mult_cr <- readRDS(MultconvRDSs[1])
conv_mult_cr
names(conv_mult_cr) <- glue("2021_conv/out/emod_mult_cr{tree}_bic_pca.rds",tree=1:100)
names(conv_mult_cr)
saveRDS(object = conv_mult_cr,file = "2021_conv/out/processed_l1ouModels/mult_convRegs_crAll.rds")

# Multivariate ext
MultconvRDSs <- dir_ls("2021_conv/out/",regexp = "mult_conv")
MultconvRDSs
conv_mult_ext <- readRDS(MultconvRDSs[2])
conv_mult_ext
names(conv_mult_ext) <- glue("2021_conv/out/emod_mult_ext{tree}_bic_pca.rds",tree=1:100)
names(conv_mult_ext)
saveRDS(object = conv_mult_ext,file = "2021_conv/out/processed_l1ouModels/mult_convRegs_extAll.rds")

# EXTERNAL
# PC1ext
# relevant RDS files
PC1convRDSsext <- dir_ls("2021_conv/out/",regexp = "PC1_conv_regs_ext")
PC1convRDSsext
PC1_conv_regs_ext <- flatten(map(PC1convRDSsext,readRDS))
PC1_conv_regs_ext
PC1_conv_regs_ext[[60]]$Y

names(PC1_conv_regs_ext)<- glue("2021_conv/out/emodPC1_ext{tree}_bic_pca.rds",tree=1:100)
names(PC1_conv_regs_ext)
saveRDS(object = PC1_conv_regs_ext,file = "2021_conv/out/processed_l1ouModels/PC1_convRegs_extAll.rds")

# PC2ext
# relevant RDS files
## PC2 ext renamed during model fitting

## PC3ext
# relevant RDS files
PC3convRDSsext <- dir_ls("2021_conv/out/",regexp = "PC3_conv_regs_ext")
PC3convRDSsext
PC3_conv_regs_ext <- flatten(map(PC3convRDSsext,readRDS))
PC3_conv_regs_ext
names(PC3_conv_regs_ext)
names(PC3_conv_regs_ext)<- glue("2021_conv/out/emodPC3_ext{tree}_bic_pca.rds",tree=1:100)
names(PC3_conv_regs_ext)
saveRDS(object = PC3_conv_regs_ext,file = "2021_conv/out/processed_l1ouModels/PC3_convRegs_extAll.rds")

## ALL PLCS 
ALLPCSconvRDS <- dir_ls("2021_conv/out/",regexp = "AllPCS_conv")
ALLPCSconvRDS
ALLPCs_conv_regs <- readRDS(ALLPCSconvRDS)
ALLPCs_conv_regs
names(ALLPCs_conv_regs)

# phylo
rodtrees <- read.nexus("2021_conv/phylo/rodent_trees265.nex")
# morphology
allPCS <- read.csv("2021_conv/data/allPCs_fm.csv")
# PCS as matrix
allPCSmat <- as.data.frame(allPCS[2:6])
row.names(allPCSmat) <- allPCS$sp
allPCSmat <- as.matrix(allPCSmat)

batchadjust <- function(treeph) {
   adjust_data(tree=treeph,Y=allPCSmat,)
}

adjusted <- map(rodtrees,batchadjust)
adjtrees <- adjusted |> map('tree')
names(adjtrees) <- c(1:100)

modtrees <- map(ALLPCs_conv_regs,'tree')

adjtrees

find_matchingTree <- function(modeltree){
  whichtree <- map(adjtrees,~identical(.,modeltree))
  names(whichtree[whichtree==T])
}

forModnames <- map_chr(modtrees,find_matchingTree)
forModnames

names(ALLPCs_conv_regs) <- glue("2021_conv/out/emodALLPCs_{tree}_bic_pca.rds",tree=forModnames)
names(ALLPCs_conv_regs)
ALLPCs_conv_regs
saveRDS(object = ALLPCs_conv_regs, file = "2021_conv/out/processed_l1ouModels/ALLPCs_convRegs_All.rds")
