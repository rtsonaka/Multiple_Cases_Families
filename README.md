# Multiple_Cases_Families
Marginal genetic effects estimation in family and twin studies using random-effects models


* RcodeFams.R contains code for the analysis of multiple-cases family data as descibed in Tsonaka etal. (2015). 
* MainFuns.R and SupportFuns.R contain all the supporting functions.
* Dataset.txt contains data from a simulated dataset of 48 multiple-cases families selected such that they contain at least one case.

This code has been used to run the simulation study descibed in Tsonaka et.al. (2015). It runs only for families with equal size and it can handle any kinship matrix (twins, sibships, etc). 
Results can be derived for: 
  1. the marginalized mixed-effects models (marginalized = TRUE),  
  2. mixed-effects logistic regression (marginalized = FALSE),
  3. without ascertainment, i.e. randomly sampled families (ascert = FALSE),
  4. with ascertainment (i.e. at least 1 or at least 2 cases, ascert = TRUE).

Reference: Tsonaka, R., van der Woude, D. and Houwing-Duistermaat, J. (2015). 
            Marginal genetic effects estimation in family and twin studies using random-effects models. 
            Biometrics, 71, 1130 - 1138.

