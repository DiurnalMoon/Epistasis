Geno <- readRDS("/Volumes/Nothing/Pig/Chip/Data/Trimmed_1cM_F2Geno-1.rds")
# Trimmed_1cM_F2Geno-1.rds 
# Contains an nx(2p) matrix with Chip data 
# trimmed down to ~ 1 SNP per 1 cM interval
# With both additive & dominant effect scaled according to 
# our beloved generalized 2 allelic model
# Column names of the matrix are:
# chromosome_location-A (for additive effect) or 
# chromosome_location-D (for dominant effect)

add <- Geno[, grep('-A', colnames(Geno))]
dom <- Geno[, grep('-D', colnames(Geno))]

n.main <- 50 # Number of main effect loci
n.size <- 150 # Number of epistatic effect loci
n.orlp <- round(n.main * 0.5) # Number of loci with both main & AA effect
p.neg <- c(1/3, 0.1, 0.9, 1/2, 1/2) # Percentage of negative effect size
v.D <- 0.15 # Variance component for dominant effect
v.AA <- 0.2 # Variance component for add-by-add effect
var_cop <- c(0.25, v.D, v.AA, 0.1, 0.1, 0.2)
powerrate=1 # rate of power law
exp=5 # rate of effect size distribution

dataDir <- "/Volumes/Nothing/Working/Epistasis/Testdata/"
# Directory for saving the simulated data.
scenario <- paste(n.size, n.main, v.D, v.AA, powerrate, sep='-' )
# Set the name describing parameters in each scenarios
# Simulated data and all results are better to be saved with the name 
# "Something_scenario_Sseed.rds"
# Example: EffectIndex_150-50-0.4-0.1-1_S1.rds
# contains all information of simulated 'true' genetic architecture
# with a network of 150 nodes following powarlaw with rate 1; 
# 50 main effect QTLs; with
# variance components of D & AA being 0.4 & 0.1.

# Function Simu_pheno
# saves simulated data to files under the dataDir
# returens phenotypic values, a vector of length = # individual
# Two packages needed: igraph MASS
pheno <- Simu_pheno(add=add, dom=dom, n.main=n.main, n.size=n.size, 
                   n.orlp=n.orlp, p.neg=p.neg,
                   dataDir=dataDir, scenario=scenario,
                   var_cop=var_cop, seed=1,
                   powerrate=1, exp=1)
