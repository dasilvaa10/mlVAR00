
##Contents

Work in progress:

A set of tools building off the mlVAR package (https://github.com/SachaEpskamp/mlVAR), a package that fits multilevel vector autoregressive models - models that are extremely useful for analyzing intensively collected longitudinal data. The functions here contain a bit of added functionality for dealing with missing data via multiple imputation in addition to providing more flexibility when dealing with the effect of "time".

##Demo

#Create some missing data

Using a longitudinal dataset from Intensive Longitudinal Methods: An Introduction to Diary and Experience Sampling Research (Bolger &Laurenceau), we’ll randomly create some missingness to illustrate some of the functions in this repo.

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

Over 50% of the rows in the dataset contain at least 1 missing value and about a quarter of the data is missing in total.

``` r
rowMiss <- apply(simData_cmbnd, 1, anyNA)

propMiss <- round(sum(rowMiss == TRUE)/ length(rowMiss), 2)

propMiss
```

    ## [1] 0.56

``` r
apply(simData_cmbnd[, 1:3], 2, function(x) sum(is.na(x))/ length(x))
```

    ##   fwkstrs    fwkdis   freldis 
    ## 0.2400000 0.2466667 0.2352381

#Impute with Amelia

Impute the missing data using Amelia. This is one of the most simple imputations one can do with Amelia. Along with speed, amelia’s strength is it’s ability to incorporate temporal information into the imputation model along with a host of other prior information. Further, Amelia has been shown to possess both solid imputation accuracy and model prediction accuracy relative to other imputation methods (kim et al, 2019).

``` r
a_out <- Amelia::amelia(simData_cmbnd, cs = "ID", ts = "day" , intercs = TRUE , m = 10, p2s = FALSE, lags = vars, leads = vars)

imps <- a_out$imputations
```

#Analyze each imputed dataset

We’ll now use the function “par_mlVAR00” this calls the main function“mlVAR00” to fit network models and store the resulting relevant information need to create estimates from the imputed datasets (Rubin,1987).

``` r
doParallel::registerDoParallel(7)

bootImp <- foreach(i=1:length(imps)) %dopar% par_mlVAR00(dat = imps[[i]], scale = TRUE, variables = vars, ID = "ID", rfStructure = c("correlated", "correlated"), timeArgs = list(NULL, NULL, FALSE))
```

#Combine the estimates

``` r
combined_ests <- rubin_combine(bootImp, m =length(imps))

list(temporal = combined_ests$temporal$t_value, contemporaneous = combined_ests$contemporaneous_network$t_value, betweenSubjects = combined_ests$`between-subjects_network`$t_value)
```

    ## $temporal
    ##             fwkstrs     fwkdis    freldis
    ## fwkstrs -1.7516728  0.8370902 -0.8235370
    ## fwkdis  -0.1369372 -1.7100384 -0.1525515
    ## freldis -0.5589783 -0.3972492 -1.3976761
    ## 
    ## $contemporaneous
    ##          fwkstrs   fwkdis  freldis
    ## fwkstrs      NaN 3.926828  4.027920
    ## fwkdis   3.759270    NaN   5.197376
    ## freldis  3.983772 5.065068      NaN

    ## $betweenSubjects
    ##            fwkstrs   fwkdis   freldis
    ## fwkstrs        NaN 1.045810 0.5115741
    ## fwkdis  1.10200309      NaN 2.3210565
    ## freldis 0.05466549 2.538301       NaN

#Compare the test-statistics from the imputed data to the complete dataset
In this 3 variable model we are looking at 21 different associates in total over 3 different networks (9 temporal, 6 contemporaneous, 6 between-subject). We can see that the conclusions drawn would be the exact same in the imputed data relative to the complete data in 20/21 pairings, the only expection being a weak temporal relationship where "fwkdis" predicts "fwkstrs" dissapearing.

``` r
completeData <- mlVAR00(dat = simData, scale = TRUE, variables = vars, ID = "id", temporal = "correlated", contemporaneous = "correlated")
```

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


#Including temporal information
Finally, temporal information can be included in the model through the “time” arguments. We can now again fit another model - this time including time and time^2 as additional fixed effects along with their respective random slopes.

``` r
completeData_time <- mlVAR00(simData, variables = vars, ID = "id", temporal = "correlated", contemporaneous = "correlated" , scale = TRUE, timeVar = "time", timePoly = 2, timeRandom = TRUE )
```

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
