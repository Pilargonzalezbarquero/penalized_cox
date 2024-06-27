## penalized_cox: Comparing Lasso and Adaptive Lasso in High-Dimensional Data
This repository contains all the code needed to implement the formulas and the simulations used in the following paper: 

* **Title**: "Comparing Lasso and Adaptive Lasso in High-Dimensional Data: A Genetic Survival Analysis in Triple-Negative Breast Cancer"
* **Authors :** Pilar Gonzalez-Barquero, Rosa E. Lillo, Alvaro Mendez-Civieta
* **Institutions**: uc3m-Santander Big Data Institute, Universidad Carlos III de Madrid, Department of Statistics, Universidad Carlos III de Madrid and Department of Biostatistics, Columbia University.
* **Keywords**: Survival analysis, Cox regression, penalization, lasso, adaptive lasso, weight calculation.




### Abstract
This study aims to evaluate the performance of Cox regression with lasso penalty and adaptive lasso penalty in high-dimensional settings. Variable selection methods are necessary in this context to reduce dimensionality and make the problem feasible. Several weight calculation procedures for adaptive lasso are proposed to determine if they offer an improvement over lasso, as adaptive lasso addresses its inherent bias. These proposed weights are based on principal component analysis, ridge regression, univariate Cox regressions and random survival forest (RSF). The proposals are evaluated in simulated datasets. 

A real application of these methodologies in the context of genomic data is also carried out. The study consists of determining the variables, clinical and genetic, that influence the survival of patients with triple-negative breast cancer (TNBC), which is a type breast cancer with low survival rates due to its aggressive nature.

### License
This project is licensed under the terms of the MIT license. See the [LICENSE](LICENSE) file for license rights and limitations.
