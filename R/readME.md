``` r
library(doParallel)
```

    ## Loading required package: foreach

    ## Loading required package: iterators

    ## Loading required package: parallel

``` r
source("C:/Users/dasil/Dropbox/code/mlVAR/mlVAR_helpers.R")
```

Using a longitudinal dataset from Intensive Longitudinal Methods: An
Introduction to Diary and Experience Sampling Research (Bolger &
Laurenceau), we’ll randomly create some missingness to illustrate some
of the functions in this repo.

``` r
simData <- bmlm::BLch9

vars <- colnames(simData[, 3:5])

simDATA_split <- split.data.frame(simData, simData$id)

nTime <- unique(sapply(simDATA_split, nrow))

nPerson <- length(simDATA_split)

iList <- list()

for (i in 1:nPerson){
  
  jList <- list()
  
  for (j in 1:length(vars)) {
    
    perMiss <- sample(seq(.05,.5,.05),1)
    
    nMiss<- perMiss * nTime
    
    missInds <- sort(sample(1:nTime, nMiss, replace = FALSE ))
    
    tempVec <- simDATA_split[[i]][, vars[j] ] 
    
    tempVec[missInds]<- NA
    
    jList[[j]] <-  tempVec
    
  }
  
  tempDat <- do.call("cbind", jList)
  
  colnames(tempDat) <- vars
  
  tempDat <- data.frame(tempDat, ID = simDATA_split[[i]]$id, day = 1:nTime )
  
  iList[[i]]<- tempDat
  
}

simData_cmbnd <- do.call("rbind", iList)
```

Over 50% of the rows in the dataset contain at least 1 missing value and
about a quarter of the data is missing in total.

``` r
rowMiss <- apply(simData_cmbnd, 1, anyNA)

propMiss <- round(sum(rowMiss == TRUE)/ length(rowMiss), 2)

propMiss
```

    ## [1] 0.59

``` r
apply(simData_cmbnd[, 1:3], 2, function(x) sum(is.na(x))/ length(x))
```

    ##   fwkstrs    fwkdis   freldis 
    ## 0.2571429 0.2623810 0.2742857

Impute the missing data using Amelia. This is about the most simple
imputation one can do with Amelia. Along with speed, amelia’s strength
is it’s ability to incorporate temporal information into the imputation
model along with a host of other prior information. Further, Amelia has
been shown to possess both solid imputation accuracy and model
prediction accuracy relative to other imputation methods (kim et al,
2019).

``` r
a_out <- Amelia::amelia(simData_cmbnd, cs = "ID", ts = "day" , intercs = TRUE , m = 10, p2s = FALSE)

imps <- a_out$imputations
```

We’ll now use the function “par\_mlVAR00” this calls the main function
“mlVAR00” to fit network models and store the resulting relevant
information need to create estimates from the imputed datasets (Rubin,
1987).

``` r
doParallel::registerDoParallel(7)

bootImp <- foreach(i=1:length(imps)) %dopar% par_mlVAR00(dat = imps[[i]], scale = TRUE, variables = vars, ID = "ID", rfStructure = c("correlated", "correlated"), timeArgs = list(NULL, NULL, FALSE))
```

Combine the estimates according to rubins rule (rubin, 1987).

``` r
combined_ests <- rubin_combine(bootImp, m =length(imps))

list(temporal = combined_ests$temporal$t_value, contemporaneous = combined_ests$contemporaneous_network$t_value, betweenSubjects = combined_ests$`between-subjects_network`$t_value)
```

    ## $temporal
    ##            fwkstrs     fwkdis    freldis
    ## fwkstrs -1.1202206  0.3815585 -0.7380703
    ## fwkdis   0.2133221 -1.3327593  0.1440999
    ## freldis -0.9339473 -0.5856300 -1.8123485
    ## 
    ## $contemporaneous
    ##          fwkstrs   fwkdis  freldis
    ## fwkstrs      NaN 3.446373 2.532986
    ## fwkdis  3.284127      NaN 5.489678
    ## freldis 2.432979 5.292629      NaN
    ## 
    ## $betweenSubjects
    ##           fwkstrs    fwkdis   freldis
    ## fwkstrs       NaN 0.1001936 0.9718311
    ## fwkdis  0.3530003       NaN 2.2148289
    ## freldis 1.0072917 2.3718595       NaN

Compare the test-statistics from the imputed data to the complete
dataset

``` r
completeData <- mlVAR00(dat = simData, scale = TRUE, variables = vars, ID = "id", temporal = "correlated", contemporaneous = "correlated")
```

    ## Check output for LMER convergence! Convergence issues may arise with a fully specified random effect stucture for a relatively larger number of variables

    ## 
    ##  estimating temporal and between-subject networks 
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |======================                                           |  33%
      |                                                                       
      |===========================================                      |  67%
      |                                                                       
      |=================================================================| 100%
    ## 
    ##  estimating contemporaneous network 
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |======================                                           |  33%
      |                                                                       
      |===========================================                      |  67%
      |                                                                       
      |=================================================================| 100%

``` r
list(temporal = completeData$results$temporal$`T-value`, contemporaneous = completeData$results$contemporaneous$`T-value`, betweenSubjects = completeData$results$`between-subjects`$`T-value`)
```

    ## $temporal
    ##            fwkstrs     fwkdis   freldis
    ## fwkstrs -1.6211098  1.9629190 -1.169135
    ## fwkdis  -0.8208542 -1.0999510 -0.190265
    ## freldis -0.7635176 -0.2474368 -1.940117
    ## 
    ## $contemporaneous
    ##          fwkstrs   fwkdis  freldis
    ## fwkstrs      NaN 4.622962 4.980959
    ## fwkdis  4.163088      NaN 4.911895
    ## freldis 4.774537 4.777671      NaN
    ## 
    ## $betweenSubjects
    ##           fwkstrs    fwkdis  freldis
    ## fwkstrs       NaN 0.6789197 0.900753
    ## fwkdis  0.8738897       NaN 2.437498
    ## freldis 0.4763302 2.2686996      NaN

Finally, temporal information can be included in the model through the
“time” arguments. We can now again fit another model - this time
including time and time^2 as additional fixed effects along with their
respective random slopes.

``` r
completeData_time <- mlVAR00(simData, variables = vars, ID = "id", temporal = "correlated", contemporaneous = "correlated" , scale = TRUE, timeVar = "time", timePoly = 2, timeRandom = TRUE )
```

    ## Check output for LMER convergence! Convergence issues may arise with a fully specified random effect stucture for a relatively larger number of variables

    ## 
    ##  estimating temporal and between-subject networks 
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |======================                                           |  33%
      |                                                                       
      |===========================================                      |  67%
      |                                                                       
      |=================================================================| 100%
    ## 
    ##  estimating contemporaneous network 
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |======================                                           |  33%
      |                                                                       
      |===========================================                      |  67%
      |                                                                       
      |=================================================================| 100%

``` r
list(temporal = completeData_time$results$temporal$`T-value`, contemporaneous = completeData_time$results$contemporaneous$`T-value`, betweenSubjects = completeData_time$results$`between-subjects`$`T-value`)
```

    ## $temporal
    ##            fwkstrs     fwkdis   freldis
    ## fwkstrs -2.5707598  2.1582631 -1.125736
    ## fwkdis  -0.8497286 -1.5653791 -0.260191
    ## freldis -0.7069282 -0.1609755 -2.747096
    ## 
    ## $contemporaneous
    ##          fwkstrs   fwkdis  freldis
    ## fwkstrs      NaN 4.579648 5.032454
    ## fwkdis  4.181363      NaN 4.826621
    ## freldis 4.832927 4.696426      NaN
    ## 
    ## $betweenSubjects
    ##          fwkstrs    fwkdis   freldis
    ## fwkstrs      NaN 0.6800027 0.8652408
    ## fwkdis  0.764811       NaN 2.5366225
    ## freldis 1.229973 2.2030348       NaN
