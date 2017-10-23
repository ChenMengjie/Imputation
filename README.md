
### Installation

**Imputation** relies on the following R packages: **Rcpp**, **RcppArmadillo**, **quadprog**, **glmnet**. All packagess are hosted on CRAN. 
  ```R
  install.packages ("Rcpp")
  install.packages ("RcppArmadillo")
  install.packages ("quadprog")
  install.packages ("glmnet")
  ```

**scPoissonGamma** can be installed from github directly as follows:

  ```R
  install.packages ("devtools")
  library(devtools)
  install_github("ChenMengjie/Imputation")
  ```
  
  
### Examples

```R
library(Imputation)
data(GSE75748_sc_time_course) # gene expression is a gene * cell count matrix
gene.expression <- gene.expression[1:5000, 1:200] # test run on a subset of genes
system.time(imp.res <- try(Imputation3_cpp(gene.expression, percentage.cutoff = 0.1, num = 3000, percentage.samples = 0.8, minbool = FALSE))) 
str(imp.res)
```
**imputed** and **imputed.by.gene** save the matrices after imputing zero values using sample relationship and both sample/gene relationship, respectively. **predicted** and **predicted.by.gene** save the matrices after predicting all entries (including nonzero entries) using sample relationship and both sample/gene relationship, respectively. 

### Author

**Mengjie Chen** (UChicago)
