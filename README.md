
### Installation

**Imputation** relies on the following R packages: **Rcpp**, **RcppArmadillo**, **quadprog**, **glmnet**. All packagess are hosted on CRAN. 
  ```R
  install.packages ("Rcpp")
  install.packages ("RcppArmadillo")
  install.packages ("quadprog")
  install.packages ("glmnet")
  ```

**Imputation** can be installed from github directly as follows:

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
system.time(imp.res <- try(Imputation1_cpp(gene.expression, percentage.cutoff = 0.1, num = 3000, percentage.samples = 0.8, minbool = FALSE))) # function using lasso for variable screening 
str(imp.res)
system.time(imp.res2 <- try(Imputation1_cpp_with_elasticnet(gene.expression, 0.1, 5000, FALSE, alpha = 0.5))) # function using elastic net for variable screening

```

**imputed** saves the matrix after imputing zero values. **predicted** save the matrix after predicting all entries (including nonzero entries). 

### Author

**Mengjie Chen** (UChicago)
