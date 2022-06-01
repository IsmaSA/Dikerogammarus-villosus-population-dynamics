
# Script for Tracking a voracious killer: Dikerogammarus villosus invasion dynamics across Europe

##load packages
suppressMessages({
  library(dplyr, quiet = TRUE, warn.conflicts = FALSE)
  library(reshape, quiet = TRUE, warn.conflicts = FALSE)
  library(metafor, quiet = TRUE, warn.conflicts = FALSE)
  library(ggplot2)
  library(tidyr)  
  library(stringr)
  library(readxl)
  library(vegan)
  library(sjPlot)
  library(MuMIn); eval(metafor:::.MuMIn)
})

###############   Load data    #####################
##read in species data
df <- read_excel("D.villosus.xlsx")

# In case of missing data in time series; we filled the missing years with NAs
new.row <- head(Time_series1[NA,], 5)    #Create new row
new.row["year"] <- c(2007,2013,2014,2015,2016)
new.row["year"] <- 2006
# assign the year without data

#We extracted the slope of each time using the Mann-Kendall trend test
TS = ts(df$Abundance,
        frequency=1,
        start= 1994)

res <-  MannKendall(TS)
print(res)
summary(res)
#we repeated this process for each time series. 

##Meta-regression of Dikerogammarus villosus
df1<- df %>% mutate(const=1)
res <- rma.mv(S_Dv, Var_Dv,  random = ~ site_id | const, data= df1)
res
forest(res, showweights = T, order="obs", slab= df1$site_id)

#We repeated this for each community metrics (i.e. Abundance, Richness, Shannon diversity and Evenness Pielou)

##Cheked bias
#Funnel plot
funnel(res, shade = c("white","gray55", "gray75"), refline = 0, legend=F) #We repeated this for each meta-regression model

#Egger test
resid <- rstandard(res)
eggers <- regtest(x = resid$resid, sei =sqrt(df$Var_Dv), model = "lm") 
eggers  


## Generalise Lineal Model 

#Before run this model, first we checked the collinearity using the corvif function
# Functions from Zuur book on mixed models
# at https://highstat.com/index.php/mixed-effects-models-and-extensions-in-ecology-with-r
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}

corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  #correlation part
  cat("Correlations of the variables\n\n")
  tmp_cor <- cor(dataz,use="complete.obs")
  print(tmp_cor)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}

# Run a model selection
#Using as example Raw-Evennes trend, adding all covariates
res1 <- glmulti(Eveness ~    Middle_point+ 
                  Years.Total+ 
                  Biogeo+
                  N+
                  Elevation+
                  Precipitation + 
                  slope_preci+
                  Temperature + 
                  slope_temp+
                  Mean_precipita + 
                  Mean_temp + 
                  DistanceKm +
                  slope + 
                  Tmin + 
                  Slope_Tmin+
                  Slope_Tmax+
                  Tmax+
                  S_Dv,
                data=df,
                level=1, method="ML",crit="aicc", confsetsize=128)
print(res1)
top <- weightable(res1)
top <- top[top$aicc <= min(top$aicc) + 1,]
top #We selected the model with the lowest AIC

#This process was repeated for each model 

#Each moderator was obtained following the methods described in the manuscript
#Bio_region: Refer to Biogeographical region
#Years_Tota: The total number of each time series
#Elevation: Elevation of the stream/river
#East: Coordinates
#Precipitation: Mean daily precipitation
#Temperature: Mean daily temperature
#DistanceKm: Distance to the next Dams
#Middle_point: Middle year of each time series
#Ecossytem: Type of ecosystems (i.e. stream or river)

colnames(df)

model_Dv <- glm(S_Dv ~ Bio_region +  Years_Total + Elevation  + East + Precipitation + Temperature + 
                  DistanceKm + Middle_point+ Ecosystem, data=df)

summary(model_Dv)










