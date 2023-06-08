# PC1 cranial sensitivity
library(analogsea)
library(magrittr)
library(future)
library(parallelly)
library(furrr)
library(l1ou)
library(fs)
library(purrr)
library(stringr)

# Path to private SSH key that matches key on DigitalOcean
ssh_private_key_file <- "/home/luisd/.ssh/id_rsa.pub"

#analogsea::sizes()

# Set up remote machines --------------------------------------------------
# Create droplet with Docker pre-installed
# "s-4vcpu-8gb", has 4 CPUs and 8 GB of RAM
#droplet1 <- docklet_create(size = "s-4vcpu-8gb")

# Pull the docker image with the environment for this project
# NB: Wait for a minute before running this so that Docker is ready to
# run on the remote machines
# droplet(droplet1$id) %>% 
#   docklet_pull("luisdva/l1outest")


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
droplet1 <- docklet_create(image = "105106896",size = "s-4vcpu-8gb")
droplet2 <- droplet_create(image = "105106896",size = "s-4vcpu-8gb")

droplet3 <- docklet_create(image = "105106896",size = "s-4vcpu-8gb")
droplet4 <- docklet_create(image = "105106896",size = "s-4vcpu-8gb")
droplet5 <- docklet_create(image = "105106896",size = "s-4vcpu-8gb")
droplet6 <- docklet_create(image = "105106896",size = "s-4vcpu-8gb")
droplet7 <- docklet_create(image = "105106896",size = "s-4vcpu-8gb")
droplet8 <- docklet_create(image = "105106896",size = "s-4vcpu-8gb")
droplet9 <- docklet_create(image = "105106896",size = "s-4vcpu-8gb")

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

ips <- c(ip1, ip2, ip3)#,ip4,ip5,ip6,ip7,ip8,ip9)
ips <- c(ip1, ip2,ip3,ip4,ip5,ip6,ip7,ip8,ip9)

# Make remote cluster -----------------------------------------------------
# Command to run on each remote machine
# The script loads the custom  image
# --net=host allows it to communicate back to this computer
rscript <- c("sudo", "docker", "run", "--net=host", 
             "luisdva/l1outest", "Rscript")

workers <- c(ips)

# Connect and create a cluster
customcluster <- parallelly::makeClusterPSOCK(
  ips,
  # User name; DO droplets use root by default
  user = "root",
  
  # Use private SSH key registered with DO
  rshopts = c(
    "-o", "ServerAliveInterval=1400",
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
extPCS <- read.csv("2021_conv/data/extPCs_fm.csv")
crPCS <- read.csv("2021_conv/data/crPCs_fm.csv")
# phylo
rodtrees <- read.nexus("2021_conv/phylo/rodent_trees265.nex")


# PCS as matrix
crPCSmat <- as.data.frame(crPCS[2:4])
row.names(crPCSmat) <- crPCS$sp
crPCSmat <- as.matrix(crPCSmat)

# shift config batch, PC1
est_shifts_phyl <- function(treeph) {
  adjdata <- adjust_data(tree=treeph,Y=crPCSmat[,1])
  eModel <- estimate_shift_configuration(adjdata$tree, adjdata$Y, max.nShifts = 55,
                                         criterion = "BIC",nCores = 4)
  eModel
}

#est_shifts_phyl(rodtreessub[[]])

treesblck <- rodtrees[c(9,11,18)]#,19,25,38,43,60,67)] #,72,89)]
#treesblck <- rodtrees[c(9,11,18,19,25,38,43,60,67)] #,72,89)]
#treesblck <- rodtrees[c(72,89)]
#treesblck <- rodtrees[c(9,11,18,19,25,38,43,60,67)] #,72,89)]
#treesblcksub <- rodtreessub[c(9,11,18)]# ,18)]#,19,25,38,43,60,67,72,89)]

crPC1shiftsrem_3_10 <- future_map(treesblck,est_shifts_phyl)
crPC1shiftsrem_3_20 <- future_map(treesblck,est_shifts_phyl)
crPC1shiftsrem_3_25 <- future_map(treesblck,est_shifts_phyl)
crPC1shiftsrem_3_30 <- future_map(treesblck,est_shifts_phyl)
crPC1shiftsrem_3_35 <- future_map(treesblck,est_shifts_phyl)
crPC1shiftsrem_3_40 <- future_map(treesblck,est_shifts_phyl)
crPC1shiftsrem_3_45 <- future_map(treesblck,est_shifts_phyl)
crPC1shiftsrem_3_50 <- future_map(treesblck,est_shifts_phyl)
crPC1shiftsrem_3_55 <- future_map(treesblck,est_shifts_phyl)
crPC1shiftsrem_3_60 <- future_map(treesblck,est_shifts_phyl)
crPC1shiftsrem_3_70 <- future_map(treesblck,est_shifts_phyl)
crPC1shiftsrem_3_80 <- future_map(treesblck,est_shifts_phyl)
crPC1shiftsrem_3_90 <- future_map(treesblck,est_shifts_phyl)
crPC1shiftsrem_3_100 <- future_map(treesblck,est_shifts_phyl)
crPC1shiftsrem_3_110 <- future_map(treesblck,est_shifts_phyl)

crPC1shiftsrem_3_10
crPC1shiftsrem_3_20
crPC1shiftsrem_3_25
crPC1shiftsrem_3_30
crPC1shiftsrem_3_35
crPC1shiftsrem_3_40
crPC1shiftsrem_3_45
crPC1shiftsrem_3_50
crPC1shiftsrem_3_55
crPC1shiftsrem_3_60
crPC1shiftsrem_3_70
crPC1shiftsrem_3_80
crPC1shiftsrem_3_90
crPC1shiftsrem_3_100
crPC1shiftsrem_3_110

saveRDS(object = crPC1shiftsrem_3_20,file = "2021_conv/out/sens_emodPC1_cr2trees20_bic_pca.rds")
saveRDS(object = crPC1shiftsrem_3_25,file = "2021_conv/out/sens_emodPC1_cr2trees25_bic_pca.rds")
saveRDS(object = crPC1shiftsrem_3_30,file = "2021_conv/out/sens_emodPC1_cr2trees30_bic_pca.rds")
saveRDS(object = crPC1shiftsrem_3_35,file = "2021_conv/out/sens_emodPC1_cr2trees35_bic_pca.rds")
saveRDS(object = crPC1shiftsrem_3_40,file = "2021_conv/out/sens_emodPC1_cr2trees40_bic_pca.rds")
saveRDS(object = crPC1shiftsrem_3_45,file = "2021_conv/out/sens_emodPC1_cr2trees45_bic_pca.rds")
saveRDS(object = crPC1shiftsrem_3_50,file = "2021_conv/out/sens_emodPC1_cr2trees50_bic_pca.rds")
saveRDS(object = crPC1shiftsrem_3_55,file = "2021_conv/out/sens_emodPC1_cr2trees55_bic_pca.rds")
saveRDS(object = crPC1shiftsrem_3_60,file = "2021_conv/out/sens_emodPC1_cr2trees60_bic_pca.rds")
saveRDS(object = crPC1shiftsrem_3_70,file = "2021_conv/out/sens_emodPC1_cr2trees70_bic_pca.rds")
saveRDS(object = crPC1shiftsrem_3_80,file = "2021_conv/out/sens_emodPC1_cr2trees80_bic_pca.rds")
saveRDS(object = crPC1shiftsrem_3_90,file = "2021_conv/out/sens_emodPC1_cr2trees90_bic_pca.rds")
saveRDS(object = crPC1shiftsrem_3_100,file = "2021_conv/out/sens_emodPC1_cr2trees100_bic_pca.rds")
saveRDS(object = crPC1shiftsrem_3_110,file = "2021_conv/out/sens_emodPC1_cr2trees110_bic_pca.rds")

stopCluster(customcluster)
customcluster
droplet_delete(droplet1)
droplet_delete(droplet2)
droplet_delete(droplet3)
droplet_delete(droplet4)
droplet_delete(droplet5)
droplet_delete(droplet6)
droplet_delete(droplet7)
droplet_delete(droplet8)
droplet_delete(droplet9)
droplets()



library(dplyr)

# relevant RDS files
resRDSsSensPC1 <- dir_ls("2021_conv/out/",regexp = "sens_emodPC1")
resRDSsSensPC1
# PC1
PC1_l1oumodsSens <- map(resRDSsSensPC1,readRDS)

seq(20,55,5)
names(PC1_l1oumodsSens) <- c(seq(20,55,5),seq(60,110,10))
names(PC1_l1oumodsSens)
PC1_l1oumodsSens[[1]][[1]]$score

BICs <- 
map_depth(PC1_l1oumodsSens,2,'score') %>% 
  tibble::enframe() %>% 
  tidyr::unnest_longer(value) %>% 
  dplyr::select(nmax=name,BIC=value) %>% 
  mutate(tree=rep(1:3,14))

BICs

nshiftsdf <- 
  map_depth(PC1_l1oumodsSens,2,'nShifts') %>% 
  tibble::enframe() %>% 
  tidyr::unnest_longer(value) %>% 
  dplyr::select(nmax=name,nShifts=value) %>% 
  mutate(tree=rep(1:3,14))
rm(nshfts)
PC1sens <- left_join(nshiftsdf,BICs)
PC1sens

library(ggplot2)

PC1sens %>% filter(tree==1) %>% 
   ggplot(aes(nmax,nShifts))+geom_point()

ggplot(PC1sens,aes(x=as.numeric(nmax),nShifts))+
  geom_point()+
  geom_line()+
  facet_wrap(~tree,scales = "free")

ggplot(PC1sens,aes(x=as.numeric(nmax),BIC))+
  geom_point()+
  geom_line()+
  facet_wrap(~tree,scales = "free")


