
# Contents

Work in progress:

A set of tools building off the mlVAR package (https://github.com/SachaEpskamp/mlVAR), a package that fits multilevel vector autoregressive models. These models are extremely useful for analyzing continuous intensively collected longitudinal data. Extending gaussian graphical models to muiltilevel data, they allow you to extract temporal (lagged relationships), contemporaneous (relationships on a given measurement period), and between-subjects (trait level relationships) networks from the data. The functions here contain a bit of added functionality for dealing with missing data via multiple imputation in addition to providing more flexibility when dealing with the effect of "time".

#Demo

## Create some missing data with the 'ml_missing' function

Using a longitudinal dataset from Intensive Longitudinal Methods: An Introduction to Diary and Experience Sampling Research (Bolger & Laurenceau), we’ll randomly create some missingness to illustrate some of the functions in this repo. The variables "fwkstrs", "fwkdis", "freldis" relate to work stress, work dissatisfaction, and relationship dissatisfaction respectively.

``` r
dat <- bmlm::BLch9

datMiss <- ml_missing(dat = dat, variables = c("fwkstrs", "fwkdis", "freldis"), 
                    temporalVar = "time", ID = "id", ranSeq = seq(.05, .7, .05))

```

We can see that around 35% percent of the data are missing

``` r
apply(datMiss[, 1:3], 2, function(x) sum(is.na(x)) / nrow(datMiss[, 1:3]))
```

    ##   fwkstrs    fwkdis   freldis 
    ## 0.3428571 0.3490476 0.3600000

In total, 73% of the rows contain at least one missing observation

``` r
sum(apply(datMiss, 1, anyNA)) / nrow(datMiss)
```

    ## [1] 0.73

When we call "format_mlVAR" to create the lagged relationships, the total number of rows containing missing data increases.

``` r
datMiss_formatted <- format_mlVAR(datMiss, scale = FALSE, timeVar = "time")

sum(apply(datMiss_formatted, 1, anyNA)) / nrow(datMiss)
```

    ## [1] 0.8661905

## Impute with Amelia

One of the strengths of mixed models are their robustness to missing data. However, in my experience, the "shift" created by looking at lagged relationships results in too much missing data to ingnore. Thus, we impute the missing data using Amelia. Below is one of the most simple imputations one can do with Amelia. Along with speed, Amelia’s strength is it’s ability to incorporate temporal information into the imputation model along with a host of other prior information. Further, Amelia has been shown to possess both solid imputation accuracy and model prediction accuracy relative to other imputation methods (kim et al, 2019).

``` r
a_out <- Amelia::amelia(datMiss, cs = "ID", ts = "time" , intercs = TRUE , m = 10, p2s = FALSE)

imps <- a_out$imputations
```

## Analyze each imputed dataset

We’ll now use the function “par_mlVAR00” this calls the main function “mlVAR00” to fit network models and store the resulting relevant information needed to create estimates from the imputed datasets (Rubin,1987).

``` r
doParallel::registerDoParallel(7)

bootImp <- foreach(i=1:length(imps)) %dopar% par_mlVAR00(dat = imps[[i]], scale = TRUE, 
														variables = c("fwkstrs", "fwkdis", "freldis"), ID = "ID", 
                                                        rfStructure = c("correlated", "correlated"), 
                                                        timeArgs = list(NULL, NULL, FALSE))
```

## Combine the estimates using mlVAR_rubin

``` r
combined_ests <- mlVAR_rubin(bootImp, m = length(imps))

list(temporal = combined_ests$temporal$t_value, 
	between = combined_ests$`between-subjects`$t_value,  
	contemporaneous = combined_ests$contemporaneous$t_value)
```

    ## $temporal
    ##            fwkstrs     fwkdis   freldis
    ## fwkstrs -1.0329728  0.8769986 -0.514470
    ## fwkdis  -0.1525965 -1.4344470 -0.126235
    ## freldis -0.3164594 -0.6359068 -1.860104
    ## 
    ## $between
    ##           fwkstrs    fwkdis  freldis
    ## fwkstrs       NaN 0.0579413 1.108926
    ## fwkdis  0.2333846       NaN 2.234443
    ## freldis 1.1261756 2.1729906      NaN
    ## 
    ## $contemporaneous
    ##          fwkstrs   fwkdis  freldis
    ## fwkstrs      NaN 4.299106 4.943368
    ## fwkdis  4.168510      NaN 4.471184
    ## freldis 4.892881 4.700831      NaN

## Compare the test-statistics from the imputed data to the complete dataset
In this 3 variable model we are looking at 21 different associations in total (from a series of univariate multilevel models) over 3 different networks (9 temporal, 6 contemporaneous, 6 between-subject). We can see that the conclusions drawn would be the exact same in the imputed data relative to the complete data in 20/21 pairings, the only expection being a weak temporal relationship where "fwkdis" predicts "fwkstrs" disappearing.

``` r
completeData <- mlVAR00(dat = dat, scale = TRUE, variables = c("fwkstrs", "fwkdis", "freldis"), ID = "id", 
                        temporal = "correlated", contemporaneous = "correlated")
```

    ## Check output for LMER convergence! Convergence issues may arise 
       with a fully specified random effect stucture for a relatively larger number of variables

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
list(temporal = completeData$results$temporal$`T-value`, 
 	betweenSubjects = completeData$results$`between-subjects`$`T-value`,
 	contemporaneous = completeData$results$contemporaneous$`T-value`)
```

    ## $temporal
    ##            fwkstrs     fwkdis   freldis
    ## fwkstrs -1.6211098  1.9629190 -1.169135
    ## fwkdis  -0.8208542 -1.0999510 -0.190265
    ## freldis -0.7635176 -0.2474368 -1.940117
    ## 
    ## $betweenSubjects
    ##           fwkstrs    fwkdis  freldis
    ## fwkstrs       NaN 0.6789197 0.900753
    ## fwkdis  0.8738897       NaN 2.437498
    ## freldis 0.4763302 2.2686996      NaN
    ## 
    ## $contemporaneous
    ##          fwkstrs   fwkdis  freldis
    ## fwkstrs      NaN 4.622962 4.980959
    ## fwkdis  4.163088      NaN 4.911895
    ## freldis 4.774537 4.777671      NaN


## Visualizing the networks
The function "prep4graph" takes in either a list of combined imputed data or the results of a normal "mlVAR00" analysis and prepares it to be input into the "qgraph" package for visualization. Here, we'll graph the imputed data (top) and the complete data (bottom).

``` r
toGraph <- prep4graph(combined_ests, type  = "imputed")

toGraph_complete <- prep4graph(completeData, type = "complete")
```


![Imputed](https://raw.githubusercontent.com/dasilvaa10/mlVAR00/master/img/imp_graph.png)

![Complete](https://raw.githubusercontent.com/dasilvaa10/mlVAR00/master/img/complete_graph.png)


## Including temporal information
Finally, temporal information can be included in the model through the “time” arguments. We can now  fit another model - this time including time and time^2 as additional fixed effects along with their respective random slopes.

``` r
completeData_time <- mlVAR00(dat = dat, scale = TRUE, variables = c("fwkstrs", "fwkdis", "freldis"), ID = "id", 
                            temporal = "correlated", contemporaneous = "correlated", 
                            timeVar = "time", timePoly = 2, timeRandom = TRUE)
```

``` r
list(temporal = completeData_time$results$temporal$`T-value`, 
	betweenSubjects = completeData_time$results$`between-subjects`$`T-value`, 
	contemporaneous = completeData_time$results$contemporaneous$`T-value`)
```

    ## $temporal
    ##            fwkstrs     fwkdis   freldis
    ## fwkstrs -2.5707598  2.1582631 -1.125736
    ## fwkdis  -0.8497286 -1.5653791 -0.260191
    ## freldis -0.7069282 -0.1609755 -2.747096
    ## 
    ## $betweenSubjects
    ##          fwkstrs    fwkdis   freldis
    ## fwkstrs      NaN 0.6800027 0.8652408
    ## fwkdis  0.764811       NaN 2.5366225
    ## freldis 1.229973 2.2030348       NaN
    ## 
    ## $contemporaneous
    ##          fwkstrs   fwkdis  freldis
    ## fwkstrs      NaN 4.579648 5.032454
    ## fwkdis  4.181363      NaN 4.826621
    ## freldis 4.832927 4.696426      NaN

