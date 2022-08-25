
# Script for Tracking a killer shrimp: Dikerogammarus villosus invasion dynamics across Europe

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

##Meta-regression 

df1<- df[!duplicated(df$site_id), ] 
df1<- df1 %>% mutate(const=1)

res <- rma.mv(S_Dv, Var_Dv,  random = ~ E+N | const,  struct="SPGAU", data= df1)

res
forest(res, showweights = T, order="obs", slab= df1$site_id)


res1 <- rma.mv(S_Abun, Var_Abun,  random = ~ E+N | const, struct="SPGAU",
               mods=~ S_Dv +Middle_point , data= df1, control=list(maxiter=1000))
res1
forest(res1, showweights = T, order="obs", slab= df1$site_id)

res2 <- rma.mv(S_Rich,Var_Rich, random = ~ E+N | const ,
               mods=~ S_Dv +Middle_point, struct="SPGAU",
               data= df1, control=list(maxiter=1000))
res2
forest(res2, showweights = T, order="obs", slab= df1$site_id)


res3 <- rma.mv(S_Diver, Var_Diver, random = ~ E+N | const,
               mods=~ S_Dv +Middle_point, struct="SPGAU",
               data= df1, control=list(maxiter=1000))
res3
forest(res3, showweights = T, order="obs", slab= df1$site_id)

res4 <- rma.mv(S_Turn, Var_Turn, random = ~ E+N | const, struct="SPGAU",mods=~ S_Dv +
                 Middle_point, data= df1, control=list(maxiter=1000))
res4
forest(res4, showweights = T, order="obs", slab= df1$site_id)

res5 <- rma.mv(S_Eve, Var_Eve, random = ~ E+N | const, struct="SPGAU",mods=~ S_Dv +
                 Middle_point, data= df1, control=list(maxiter=1000))
res5
forest(res5, showweights = T, order="obs", slab= df1$site_id)



#Example of Egger test (Repeated for each model)
res
resid = rstandard(res)
eggers <- regtest(x = resid$resid, sei =sqrt(df1$Var_Dv), model = "lm")
eggers 

#Funnel plot
funnel(res, shade = c("white","gray55", "gray75"), refline = 0, legend=F) #We repeated this for each meta-regression model


#Heterogeneity I2
W <- diag(1/df1$Var_Dv)
X <- model.matrix(res)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2 <- 100 * res$sigma2 / (sum(res$sigma2) + (res$k-res$p)/sum(diag(P)))
names(I2) <- c("Easy")
round(I2, 2)
cat("Total I2: ", round(sum(I2), 2))




## Generalised Lineal Model 

#Before run this model, first we checked the collinearity using the corvif function
# Functions from Zuur et al., 2009
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
#Ecossytem: Type of ecosystems (i.e. stream or large river)



colnames(df)
df1 <- df[!duplicated(df$site_id),]
moderators <- c("S_Dv", "E","N","Elevation","Slope_preci","Slope_temp","Mean_precipita",
                "Mean_temp","DistanceKm","slope","Avg_Tmin","Slope_Tmin","Slope_Tmax","Avg_Tmax","Middle_point")

corvif(mutate_if(df1[, moderators], is.factor, as.numeric))
#VIF =5

#Rate of change of D. villosus trend
m1 <- glm(S_Dv ~  
            Elevation+
            Slope_preci+
            Slope_temp+
            Mean_precipita +
            Mean_temp+
            DistanceKm+
            slope+
            Avg_Tmax+
            Avg_Tmin+
            Slope_Tmax+
            Middle_point+
            Biogeo+
            Ecosystem,
          data = df1, na.action="na.fail")

summary(m1)


p<- plot_model(m1, vline.color = "black", sort.est = TRUE,show.values = TRUE, value.offset = .3)


p <- plot_model(m1, vline.color = "red", transform = NULL)
p<-plot_model(m1, type = "eff", terms = "Middle_point")
p + theme_sjplot()
p

moderators <- c("Proportion", "E","N","Elevation","Precipitation","Temperature",
                "DistanceKm","slope","Middle_point","Tmax","Tmin","year","Avg_Tmax","Avg_Tmin","Mean_precipita",
                "Mean_temp")


corvif(mutate_if(df[, moderators], is.factor, as.numeric))



## Relative abundance model (i.e. Dominance).
m3 <- glm(Proportion ~  
            Elevation+
            Precipitation  +
            Temperature    +
            Mean_temp+
            DistanceKm+
            slope+
            Tmin           +
            Avg_Tmax       +
            Tmax           +
            year           +
            Biogeo+
            Ecosystem,
          data = df,family = "quasibinomial",  na.action="na.fail")

summary(m3)

p<- plot_model(m3, vline.color = "black", sort.est = TRUE, show.values = TRUE, value.offset = .3,
               transform = NULL)
p<-plot_model(m3, type = "eff", terms = "year")
p + theme_sjplot()








