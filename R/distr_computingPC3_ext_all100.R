# PC3 external, all 100 trees
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
#images(private = TRUE)
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
droplet10 <- docklet_create(image = "105106896",size = "s-4vcpu-8gb")

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
ip10 <- droplet(droplet9$id)$networks$v4[[1]]$ip_address

ips <- c(ip1, ip2,ip3,ip4,ip5,ip6,ip7,ip8,ip9,ip10)

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
# phylo
rodtrees <- read.nexus("2021_conv/phylo/rodent_trees265.nex")


# PCS as matrix
extPCSmat <- as.data.frame(extPCS[2:4])
row.names(extPCSmat) <- extPCS$sp
extPCSmat <- as.matrix(extPCSmat)

# shift config batch, PC1
est_shifts_phyl <- function(treeph) {
  adjdata <- adjust_data(tree=treeph,Y=extPCSmat[,3])
  eModel <- estimate_shift_configuration(adjdata$tree, adjdata$Y, max.nShifts = 26,
                                         criterion = "BIC",nCores = 4)
  eModel
}

treesblck <- rodtrees

extPC3shifts_first50 <- future_map(treesblck[1:50],est_shifts_phyl)
extPC3shifts_first50
saveRDS(object = extPC3shifts_first50,file = "2021_conv/out/emodPC3_ext_first50trees_bic_pca.rds")

extPC3shifts_second50 <- future_map(treesblck[51:100],est_shifts_phyl)
extPC3shifts_second50
saveRDS(object = extPC3shifts_second50,file = "2021_conv/out/emodPC3_ext_second50trees_bic_pca.rds")
#extPC3shifts_second50 <- readRDS("2021_conv/out/emodPC3_ext_second50trees_bic_pca.rds")

PC3ext_conv_regs_first50 <- future_map(extPC3shifts_first50,estimate_convergent_regimes,nCores=4)
PC3ext_conv_regs_first50
PC3ext_conv_regs_second50 <- future_map(extPC3shifts_second50,estimate_convergent_regimes,nCores=4)
PC3ext_conv_regs_second50

saveRDS(object = PC3ext_conv_regs_first50,file = "2021_conv/out/PC3_conv_regs_ext-first50trees.rds")
saveRDS(object = PC3ext_conv_regs_second50,file = "2021_conv/out/PC3_conv_regs_ext-second50trees.rds")

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
droplet_delete(droplet10)
droplets()
