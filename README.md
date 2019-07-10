causeSims: A package to help with simulating summary statistic data and testing MR methods
======



## Installation 

```{r}
devtools::install_github("jean997/causeSims")
```

You will need to install these R packages manually following the instructions in the links provided. 

+ CAUSE (http://github.com/jean997/cause) 
+ gsmr (http://cnsgenomics.com/software/gsmr/)
+ MR-PRESSO (https://github.com/rondolab/MR-PRESSO)

`causeSims` has other dependencies that should install automatically.
If you have trouble installing `RcppZiggurat` you may need to install GSL. On Linux use

```
sudo apt-get install libgsl-dev
```
Please make sure you are using `mixsqp-0.1-97` which is currently the version on CRAN, rather than the development version. 


## Generate a data set 

The main data simulation function, `sum_stats` generates a data frame containing GWAS summary statistics for two traits. Data are simulated using LD information that must be provided by the user. We use LD estimates from the 1000 Genomes CEU population on chromosome 19. These can be downloaded [here](https://zenodo.org/record/3235780) or with the following R commands

```
if(!dir.exists("data/")) system("mkdir data")
cat("Downloading LD Data\n")
download.file(url="https://zenodo.org/record/3235780/files/chr19_snpdata_hm3only.RDS?download=1", destfile = "data/chr19_snpdata_hm3only.RDS")
download.file(url="https://zenodo.org/record/3235780/files/evd_list_chr19_hm3.RDS?download=1", destfile="data/evd_list_chr19_hm3.RDS")
```

We recomend putting the two files in a directory called `data` (the R commands above do this for you). This will make it easy to run this example and also to run more simulations using the DSC file in this repository. 

LD data is passed to `sum_stats` as two objects. The first, `snps` is a data frame with variant information. The data frame must have at least these three columns with the following names:

+ `AF` allele frequency
+ `SNP` snp name
+ `region_id` Region ID corresponding to `evd_list`

The second object, `evd_list` is a list of Eigen decompositions for varaint correlation matrices. The elements of `evd_list` can be created using the `eigen` function and are lists with elements `values` and `vectors`. The correlation matrix specified by `evd_list[[i]]` is the correlation matrix for variants with `region_id==i` in the variant information data frame. Thus
we should have `sum(snps$region_id == i) == nrow(evd_list[[i]]$vectors)` for all `i`. The `sum_stats` function assumes that variants in the variant information data frame appear in order.

On four cores it takes about 20 seconds to make one data set:
```{r}
library(causeSims)
snps <- readRDS("data/chr19_snpdata_hm3only.RDS") 
evd_list <- readRDS("data/evd_list_chr19_hm3.RDS")

set.seed(1)
dat <- sum_stats(snps, evd_list,
                  n_copies = 30,
                  n1 = 50000, n2=50000, h1=0.3, h2=0.3,
                  neffect1 = 1000, neffect2 = 1000,
                  gamma = 0, eta = 0.3, q = 0.3,
                  cores = 4, ld_prune_pval_thresh = 1e-3, r2_thresh = 0.1)


dim(dat)
#[1] 584700     18

head(dat)
#     beta_hat_1  beta_hat_2        seb1        seb2 ld_b1 ld_b2 b1 b2 beta_hat_1_nold beta_hat_2_nold          snp region_id
#1 -0.0008899914 0.004018979 0.006471832 0.006471832     0     0  0  0    -0.001196115    0.0099583122  rs8100066.1         1
#2  0.0017467308 0.002086462 0.006350853 0.006350853     0     0  0  0    -0.006669769   -0.0002480333  rs8102615.1         1
#3 -0.0005518635 0.003663230 0.006471832 0.006471832     0     0  0  0    -0.007202967   -0.0117717884  rs8105536.1         1
#4 -0.0151432368 0.009930815 0.013253119 0.013253119     0     0  0  0     0.004230415   -0.0007057890  rs2312724.1         1
#5 -0.0014797713 0.003971240 0.006457704 0.006457704     0     0  0  0    -0.003104569    0.0022862145  rs1020382.1         1
#6  0.0012875652 0.004577477 0.006659461 0.006659461     0     0  0  0     0.003908199   -0.0005627174 rs12459906.1         1
#         AF ldscore_hm3 rep   p_value p_value_nold ld_prune
#1 0.3939394    7.503727   1 0.8906216    0.8533713    FALSE
#2 0.4545455    7.029503   1 0.7832864    0.2936187    FALSE
#3 0.3939394    7.505135   1 0.9320455    0.2657206    FALSE
#4 0.9393939    1.887018   1 0.2531977    0.7495737    FALSE
#5 0.3989899    7.366604   1 0.8187537    0.6306912    FALSE
#6 0.3434343    7.516243   1 0.8466898    0.5572949    FALSE
```

## Run some MR methods

Now you can analyze the data set with many methods that have wrappers built into the package. The data we just produced contains
summary statistics generated both with and without LD. To run the no LD version, use the option `no_ld = TRUE` in all of the wrappers below. All of the analysis wrappers in `causeSims` are very simple so looking at the code may be the fastest way to clear up any questions about what they are doing. 

CAUSE: This takes about 4 minutes for parameter estimation and 45 seconds for model fitting (your times might vary).
We do these steps separately to facilitate experimenting with parameters.

```{r}
params <- cause_params_sims(dat, null_wt = 10, no_ld=FALSE)
cause_res <- cause_sims(dat, params, no_ld = FALSE)
summary(cause_res)
```


LCV: LCV is not currently distributed as a package so we coppied the code into this package. Code was originally downloaded from https://github.com/lukejoconnor/LCV.
This takes about 40 seconds. By default, if the data simulated with no LD are used, LCV is run with `crosstrait.intercept = 0` and `ldsc.intercept = 0`. If the data simulated with LD, we use `crosstrait.intercept = 1` and `ldsc.intercept = 1`. This behavior can be modified by changing the `intercepts` argument. 

```{r}

lcv_res <- lcv_sims(dat,no_ld = FALSE, intercepts = TRUE, sig.thresh = 30)

lcv_res
```

GSMR: GSMR requires access to LD data. We use to true correlation struction that the data was simulated from which is more accurate than might be available in general. GSMR also requires a $p$-value threshold used for determining which variants to include. 

```{r}
gsmr_res <- gsmr_sims(dat, evd_list, p_val_thresh  = 5e-8, no_ld = FALSE)
```


The other methods included all have the same function arguments and return format (a list with a $p$-value and an estimate). They all require a $p$-value threshold. 

```{r}
ivw_res <- ivw(dat, p_val_thresh=5e-8, no_ld = FALSE)

egger_res <- egger(dat, p_val_thresh=5e-8, no_ld = FALSE)

mrpresso_res <- mrpresso(dat, p_val_thresh=5e-8, no_ld = FALSE)

wm_res <- weighted_median(dat, p_val_thresh=5e-8, no_ld = FALSE)

twas_res <- twas(dat, p_val_thresh=5e-8, no_ld = FALSE)

pwm_res <- pweighted_median(dat, p_val_thresh=5e-8, no_ld = FALSE)

```

## Use DSC to run many simulations

We have included two dsc files in the `dsc_files` subdirectory of this repository. These can be used to run lots of simulations experimenting with parameters. You can find more instructions on doing this [here](https://jean997.github.io/cause/simulations.html). 

