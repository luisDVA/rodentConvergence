# PC1-5 all measurements
library(analogsea)
library(magrittr)
library(future)
library(parallelly)
library(furrr)
library(l1ou)


# Path to private SSH key that matches key on DigitalOcean
ssh_private_key_file <- "/home/luisd/.ssh/id_rsa.pub"

#analogsea::sizes()
droplets()

# Set up remote machines --------------------------------------------------
# Create droplet with Docker pre-installed
# "s-4vcpu-8gb", has 4 CPUs and 8 GB of RAM
#droplet1 <- docklet_create(size = "s-4vcpu-8gb")

# Pull the docker image with the environment for this project
# NB: Wait for a minute before running this so that Docker is ready to
# run on the remote machines
#droplet(droplet1$id) %>% 
#  docklet_pull("luisdva/l1outest")


# Create snapshot
# droplet(droplet1$id) %>% 
#   droplet_power_off() %>% 
#   droplet_snapshot(name = "l1ou_ready") %>% 
#   droplet_power_on()

# Create a new droplet based on this snapshot. This new computer will already
# have rocker/tidyverse on it
images(private = TRUE)
# You can see a list of available snapshots and get the name/id with
# images(private = TRUE)
droplet1 <- docklet_create(image = "106746335",size = "s-4vcpu-8gb")
droplet2 <- droplet_create(image = "106746335",size = "s-4vcpu-8gb")
droplet3 <- docklet_create(image = "106746335",size = "s-4vcpu-8gb")
droplet4 <- docklet_create(image = "106746335",size = "s-4vcpu-8gb")

droplet5 <- docklet_create(image = "106746335",size = "s-4vcpu-8gb")
droplet6 <- docklet_create(image = "106746335",size = "s-4vcpu-8gb")
droplet7 <- docklet_create(image = "106746335",size = "s-4vcpu-8gb")
droplet8 <- docklet_create(image = "106746335",size = "s-4vcpu-8gb")
droplet9 <- docklet_create(image = "106746335",size = "s-4vcpu-8gb")
droplet10 <- docklet_create(image = "106746335",size = "s-4vcpu-8gb")

# Get IP addresses
ip1 <- droplet(droplet1$id)$networks$v4[[1]]$ip_address
ip2 <- droplet(droplet2$id)$networks$v4[[1]]$ip_address
ip3 <- droplet(droplet3$id)$networks$v4[[1]]$ip_address
ip4 <- droplet(droplet4$id)$networks$v4[[1]]$ip_address
ip5 <- droplet(droplet5$id)$networks$v4[[1]]$ip_address
ip6 <- droplet(droplet6$id)$networks$v4[[1]]$ip_address
ip7 <- droplet(droplet7$id)$networks$v4[[1]]$ip_address
ip8 <- droplet(droplet8$id)$networks$v4[[1]]$ip_address
ip9 <- droplet(droplet9$id)$networks$v4[[1]]$ip_address
ip10 <- droplet(droplet10$id)$networks$v4[[1]]$ip_address

ips <- c(ip1, ip2)
ips <- c(ip1, ip2,ip3,ip4)#,ip4,ip5,ip6,ip7,ip8,ip9)
ips <- c(ip1, ip2,ip3,ip4,ip5)
ips <- c(ip1, ip2,ip3,ip4,ip5,ip6,ip7,ip8,ip9,ip10)

# Make remote cluster -----------------------------------------------------
# Command to run on each remote machine
# The script loads the custom  image
# --net=host allows it to communicate back to this computer
rscript <- c("sudo", "docker", "run", "--net=host", 
             "luisdva/l1outest", "Rscript")

workers <- c(ips)
workers

# Connect and create a cluster
customcluster <- parallelly::makeClusterPSOCK(
  ips,
  # User name; DO droplets use root by default
  user = "root",
  
  # Use private SSH key registered with DO
  rshopts = c(
    "-o", "ServerAliveInterval=1100", #1400
    "-o", "ServerAliveCountMax=3",
    "-o", "StrictHostKeyChecking=no",
    "-o", "IdentitiesOnly=yes",
    "-i", ssh_private_key_file
  ),
  
  rscript = rscript,
  
  dryrun = FALSE
)

customcluster


# Use the cluster of computers as the backend for future fns
future::plan(cluster, workers = customcluster)

# morphology
allPCS <- read.csv("2021_conv/data/allPCs_fm.csv")

# phylo
rodtrees <- read.nexus("2021_conv/phylo/rodent_trees265.nex")


# PCS as matrix
allPCSmat <- as.data.frame(allPCS[2:6])
row.names(allPCSmat) <- allPCS$sp
allPCSmat <- as.matrix(allPCSmat)


allPCSmat
# shift config batch, PC1-5
est_shifts_phyl <- function(treeph) {
  adjdata <- adjust_data(tree=treeph,Y=allPCSmat)
  eModel <- estimate_shift_configuration(adjdata$tree, adjdata$Y, max.nShifts = 22,
                                         criterion = "BIC",nCores = 4)
  eModel
}

#est_shifts_phyl(rodtreessub[[]])

treesblck <- rodtrees[c(9,11,18,19,25,38,43,60,67)] #,72,89)]
treesblck <- rodtrees[c(72,89)]
#treesblck <- rodtrees[c(9,11,18,19,25,38,43,60,67)] #,72,89)]
#treesblcksub <- rodtreessub[c(9,11,18)]# ,18)]#,19,25,38,43,60,67,72,89)]

allPCSshifts_first5 <- future_map(rodtrees[1:5],est_shifts_phyl)
allPCSshifts_next20 <- future_map(rodtrees[6:25],est_shifts_phyl)
allPCSshifts_n25to30 <- future_map(rodtrees[26:30],est_shifts_phyl)
allPCSshifts_n30to35 <- future_map(rodtrees[31:35],est_shifts_phyl)
allPCSshifts_n35to40 <- future_map(rodtrees[36:40],est_shifts_phyl)
allPCSshifts_n40to45 <- future_map(rodtrees[41:45],est_shifts_phyl)
#allPCSshifts_n45to50 <- future_map(rodtrees[45:50],est_shifts_phyl)
allPCSshifts_n45to50 <- map(rodtrees[46:50],est_shifts_phyl)

allPCSshifts_n50to55 <- future_map(rodtrees[51:55],est_shifts_phyl)
allPCSshifts_n55to60 <- future_map(rodtrees[56:60],est_shifts_phyl)
allPCSshifts_n61to65 <- future_map(rodtrees[61:65],est_shifts_phyl)
allPCSshifts_n66to70 <- future_map(rodtrees[66:70],est_shifts_phyl)
allPCSshifts_n71to75 <- future_map(rodtrees[71:75],est_shifts_phyl)
allPCSshifts_n76to80 <- future_map(rodtrees[76:80],est_shifts_phyl)
allPCSshifts_n81to85 <- future_map(rodtrees[81:85],est_shifts_phyl)
allPCSshifts_n85to90 <- future_map(rodtrees[86:90],est_shifts_phyl)
allPCSshifts_n91to95 <- future_map(rodtrees[91:95],est_shifts_phyl)
allPCSshifts_n96to100 <- future_map(rodtrees[96:100],est_shifts_phyl)

allpcs95 <- mget(ls(pattern="allPCSs"))
allpcs95 <- flatten(allpcs95)
saveRDS (object = allpcs95,file = "2021_conv/out/emodAllPCS_first95trees_bic_pca.rds")
allpcs95 <- readRDS(file = "2021_conv/out/emodAllPCS_first95trees_bic_pca.rds")




all100PCStrees <- c(allpcs95,allPCSshifts_n25to30)
all100PCStrees[6:10]
all100PCStrees[96:100]
identical(all100PCStrees[6:10],all100PCStrees[96:100])
identical(all100PCStrees[6:10],all100PCStrees[16:20])

all95pending <- all100PCStrees[c(1:5,11:100)]

allPCS100trees <- c(all95pending,allPCSshifts_n45to50)
allPCS100trees %>% map_dbl('score') %>% unique() %>% length()


saveRDS(object = allPCS100trees,file = "2021_conv/out/emodAllPCS_100trees_bic_pca.rds")
allPCS100trees <- readRDS("2021_conv/out/emodAllPCS_100trees_bic_pca.rds")
write
library(purrr)

adjusteddats <-rodtrees %>% map(~adjust_data(tree=.x,Y=allPCSmat))

treesadj <- adjusteddats %>% map('tree')

purrr::map(mami,"score")

map2(map(mami,"tree"),treesadj, ~identical(.x,.y))

identical(treesadj[[25]],mami[[25]]$tree)

#saveRDS(object = allPCSshifts_first5,file = "2021_conv/out/emodAllPCS_first5trees_bic_pca.rds")
#saveRDS(object = allPCSshifts_next20,file = "2021_conv/out/emodAllPCS_6to25trees_bic_pca.rds")

allPCSshifts_first5 <- readRDS("2021_conv/out/emodAllPCS_first5trees_bic_pca.rds")
allPCSshifts_next20 <- readRDS("2021_conv/out/emodAllPCS_6to25trees_bic_pca.rds")

allPCS_conv_regs_all100 <- 
  future_map(allPCS100trees,estimate_convergent_regimes,nCores=4)

allPCS_conv_regs_all100
saveRDS(object = allPCS_conv_regs_all100,file = "2021_conv/out/AllPCS_conv_regs_all100trees.rds")

stopCluster(customcluster)
rm(customcluster)
droplets()
droplet_delete(droplet1)
droplet_delete(droplet2)
droplet_delete(droplet3)
droplet_delete(droplet4)
droplet_delete(droplet5)
droplet_delete(droplet6)
droplet_delete(droplet7)
droplet_delete(droplet8)
droplet_delete(droplet9)
droplet_delete(droplet10)

# lagartijas
data(lizard.tree, lizard.traits)
# here lizard.traits already has row names:

lizard <- adjust_data(lizard.tree, lizard.traits[,1])

lista_arbs <- list(lizard.tree,lizard.tree,lizard.tree)
class(lista_arbs) <- "multiPhylo"
lista_arbs

est_shifts_phyl <- function(treeph) {
  adjdata <- adjust_data(tree=treeph,Y=lizard.traits[,1])
  eModel <- estimate_shift_configuration(adjdata$tree, adjdata$Y, 
                                         criterion = "BIC",nCores = 4)
  eModel
}

future_map(lista_arbs,est_shifts_phyl)


eModel <- estimate_shift_configuration(lizard$tree, lizard$Y)
eModel












# fit ou model
eModel <- estimate_shift_configuration(cranial$tree, cranial$Y, 
                                       criterion = "BIC",nCores = 3)

debug(l1ou:::select_best_solution)

x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE, FALSE, FALSE, TRUE))
lapply(x, FUN = quantile, probs = 1:3/4)
future_lapply(x, FUN = quantile, probs = 1:3/4)


droplets()


objs <- 
c("allPCSshifts_first5", "allPCSshifts_next20", "allPCSshifts_n25to30", "allPCSshifts_n30to35", "allPCSshifts_n35to40", "allPCSshifts_n40to45", "allPCSshifts_n45to50", "allPCSshifts_n50to55", "allPCSshifts_n55to60", "allPCSshifts_n61to65", "allPCSshifts_n66to70", "allPCSshifts_n71to75", "allPCSshifts_n76to80", "allPCSshifts_n81to85", "allPCSshifts_n85to90", "allPCSshifts_n91to95", "allPCSshifts_n96to100")

sort(objs)
mget()

modtrees <- map(mami,"tree")

map(modtrees,map(treesadj,~identical(.x,.y)))

findcorresponding <- function(adj_tree){
which(map_lgl(map(mami,"tree"),~identical(.x,adj_tree)))
}

xuales <- map(treesadj,findcorresponding)
names(xuales) <- 1:100
xuales
chu <- xuales %>% map(paste0,collapse=",") %>% tibble::enframe() %>% 
  dplyr::mutate(name=as.numeric(name)) %>% 
  tidyr::unnest_auto(value)

jurge <- tibble::enframe(map(treesadj,findcorresponding)) %>% 

