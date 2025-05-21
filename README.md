# DegreeAFT Quick‑Start Guide
A lightweight **R** package for **network‑regularized accelerated failure‑time (AFT) survival modelling**.  
DegreeAFT implements inverse‑degree so that biological network topology guides feature selection, improving interpretability and robustness in high‑dimensional omics studies. It also provides an option to use Laplacian penalty.

---

## Installation
```r
# install.packages("devtools")
devtools::install_github("yliu38/DegreeAFT")
```

---

## Worked Example (Kidney Renal Clear‑Cell Carcinoma)
### Load libaries
```r
library(DegreeAFT)
library(tidyverse); library(Matrix); library(survival)
library(STRINGdb);  library(igraph); library(WGCNA)
library(survAUC);   library(future.apply)
```

### 1  Load the demo data
```r
load(file="kirc.RData")

exp_kirc_agg <- kirc$expr        # gene expression (p × n)
sur_kirc     <- kirc$survival    # data.frame(sample, OS.time, OS)
marker_kirc  <- kirc$marker      # vector of marker genes
```

### 2  Build a co‑expression network (WGCNA)
```r
wgcna_net <- build_wgcna_net(exp_kirc_agg[marker_kirc, ])
```
`wgcna_net` includes:
* `genes`    – gene symbols present in the network  
* `A`        – adjacency matrix  
* `DegInv`   – inverse‑degree diagonal matrix  
* `Lap`      – graph Laplacian  

### 3  Hyper‑parameter tuning via 5‑fold CV
Prepare the tuning grid and align samples:
```r
grid <- expand.grid(
  lam1 = c(0.005, 0.01, 0.02),
  lam2 = c(0.05, 0.10, 0.20, 0.40)
)

samples <- intersect(colnames(exp_kirc_agg), sur_kirc$sample)
sur_df  <- sur_kirc[match(samples, sur_kirc$sample), ]
time    <- sur_df$OS.time
status  <- sur_df$OS
X_full  <- t(exp_kirc_agg[wgcna_net$genes, samples])  # n × p
```

#### 3.1  Inverse‑degree penalty
```r
wgcna_net$DegInv <- as.matrix(wgcna_net$DegInv)   # ensure base matrix
cv_deg <- cv_select(
  time, status, X_full,
  wgcna_net$DegInv, prox_deg,
  grid, K = 5
)
```

#### 3.2  Laplacian penalty (reference)
```r
cv_lap <- cv_select(
  time, status, X_full,
  wgcna_net$Lap, prox_lap,
  grid, K = 5
)

bind_rows(
  InverseDegree = cv_deg$cv_metrics,
  Laplacian     = cv_lap$cv_metrics,
  .id = "Method")
```

### 4  Fit the final model & inspect selected genes
```r
best <- cv_deg$best_lambda

beta <- prox_deg(
  time, status, X_full,
  wgcna_net$DegInv,
  best$lam1, best$lam2
)

sel <- which(beta != 0)

tibble(
  gene     = wgcna_net$genes[sel],
  coef     = beta[sel],
  isolated = isol_status(wgcna_net$genes[sel], wgcna_net$A)
) %>% 
  arrange(desc(abs(coef))) %>% 
  print(n = 20)
```
The table lists all genes with non‑zero coefficients, sorted by effect size, and flags those that are network‑isolated.

---

## Citation
If you use **DegreeAFT** in published work, please cite:

```
Liu Y et al. Network‑regularized accelerated failure‑time models improve biomarker discovery in cancer. bioRxiv, 2025.
```

---
