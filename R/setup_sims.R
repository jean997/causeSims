#'@title Set up simulation directory
#'@description A helper function that prepares a working directory for running simulations.
#'@details This function prepares a working directory to replicate the simulations presented in Morrison et al 2019
#'(https://www.biorxiv.org/content/10.1101/682237v1). After running this function you will be able to execute simulation
#'experiments using DSC using the commands
#'
#'  nohup dsc --replicate 100 --host config.yml -c 4 pwr.dsc > pwr_dsc.out &
#'
#'  nohup dsc --replicate 100 --host config.yml -c 4 false_positives.dsc > fp_dsc.out &
#'
#'  You will need to modify the config.yml file to work with your compute cluster.
#'@export
setup_sims <- function(){
  #Download DSC files
  cat("Downloading DSC files\n")
  download.file(url = "https://raw.githubusercontent.com/jean997/causeSims/master/dsc_files/power.dsc" , destfile="power.dsc")
  download.file(url = "https://raw.githubusercontent.com/jean997/causeSims/master/dsc_files/false_positives.dsc" , destfile="false_positives.dsc")

  #Download data
  if(!dir.exists("data/")) system("mkdir data")
  cat("Downloading LD Data\n")
  download.file(url="https://zenodo.org/record/3235780/files/chr19_snpdata_hm3only.RDS?download=1", destfile = "data/chr19_snpdata_hm3only.RDS")
  download.file(url="https://zenodo.org/record/3235780/files/evd_list_chr19_hm3.RDS?download=1", destfile="data/evd_list_chr19_hm3.RDS")

  cat("Downloading config file. You will need to modify this for your setup\n")
  download.file(url="https://raw.githubusercontent.com/jean997/causeSims/master/dsc_files/config.yml", destfile="config.yml")
}
