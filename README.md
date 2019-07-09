causeSims: A package to help with simulating summary statistic data and testing MR methods
======



## Installation 

```{r}
devtools::install_github("jean997/causeSims")
```

You will need 

+ CAUSE (http://github.com/jean997/cause) 
+ gsmr (http://cnsgenomics.com/software/gsmr/)
+ MR-PRESSO (https://github.com/rondolab/MR-PRESSO)

and some other packages available through CRAN that should install automattically.
If you have trouble installing `RcppZiggurat` you may need to install GSL. On Linux use

```
sudo apt-get install libgsl-dev
```
Please use `mixsqp-0.1-97` which is currently the version on CRAN. 


## Generate a data set 

We are using some LD data that can be downloaded [here.](https://zenodo.org/record/3235780)
We recomend putting these files in a directory called `data`. This will make it easy to run this example and also to run more simulations using the DSC file in this repository.
The data consist of a data frame with snp information. This data frame has at least these three columns

+ AF allele frequency
+ SNP snp name
+ region_id region_id corresponding to `evd_list`

`evd_list` is a list of eigen decompositions for correlation matrices (can be created with `eigen`). 
These are lists with elements `values` and `vectors`. `evd_list[[i]]` should give the correlation matrix for SNPs with `region_id==i`. Thus
we should have `sum(snps$region_id == i) == nrow(evd_list[[i]]$vectors)`.

On four cores it takes about 20 seconds to make one data set. 

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
#[1] 584700     19
```

## Run some MR methods

Now you can analyze the data set with many methods that have wrappers built into the package. The data we just produced contains
summary statistics generated both with and without LD. To run the no LD version, use the option `no_ld = TRUE` in all of the wrappers below. All of the analysis wrappers in `causeSims` are very simple so looking at the code may be the fastest way to clear up any questions about what they are doing. We have written them only to facilitate easily analyzing simulated data and switching between the data generated with and without LD. 

CAUSE. This takes about 4 minutes for parameter estimation and 45 seconds for model fitting (for one data set on Jean's laptop). 
We do these steps separately to facilitate experimenting with parameters.

```{r}
params <- cause_params_sims(dat, null_wt = 10, no_ld=FALSE)
cause_res <- cause_sims(dat, params, no_ld = FALSE)
summary(cause_res)
```


LCV: LCV is not currently distributed as a package so we coppied the code into this package. Code was originally downloaded from https://github.com/lukejoconnor/LCV.
This takes about 40 seconds. 

```{r}

lcv_res <- lcv_sims(dat,no_ld = FALSE, intercepts = TRUE, sig.thresh = 30)

lcv_res
```

GSMR: GSMR needs access to LD data

```{r}
gsmr_res <- gsmr_sims(dat, evd_list, p_val_thresh  = 5e-8, no_ld = FALSE)
```


Many other methods that all have the same function arguments and return format

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

