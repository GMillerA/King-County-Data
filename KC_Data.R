#########King County Spatial EDA##########
###Load Libraries####
library(geoR)
library(spdep)
library(tidyverse)
library(boot)
library(readr)
library(Hmisc)
library(glmnet) # fit glms with elastic net
library(gstat) 
###Load Data###
kc_data <- read_csv("~/King County House Data/kc_house_data.csv")

######EDA#####
head(kc_data, 5)
#Get rid of ID and Date First
df <- kc_data[,-c(1,2)]
#Compute Correlation matrix for df
res <- cor
round(res, 2)
#Plot Correlation matrix(df)
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = res, col = col, symm = TRUE)
#Price, bathrooms, grade, sqft_above, sqft_living15, sqft_living all correlated

#visualize pairwise correlations
pairs(df[,c('price', 'bathrooms', 'grade', 'sqft_above', 'sqft_living15')])

#Feature selection from those four with price as the target
# convert input variables to a matrix
input <- df[,c('bathrooms', 'grade', 'sqft_above', 'sqft_living15')]
input <- as.matrix(input)
# get a vector with our output variable
output <- df$price

cv_fit <- cv.glmnet(input, output, family = "gaussian")

# get coefficents for the best model
coef(cv_fit, s = "lambda.min")

# get a (non-sparse) matrix of the coefficents for the best model
coef_matrix <- coef(cv_fit, s = "lambda.min") %>% 
  as.matrix()


# get the variables with a coefficent that's not 0 
variables <- row.names(coef_matrix)[coef_matrix != 0] %>% # get variables w/ non-zero intercepts
  setdiff("(Intercept)") #remove the intercept (if it's included)

# this variable has just our selected features. 
# print the first few
head(variables)

# turn our list of formulas into a variable
variables_selected <- paste(variables, collapse="+")
formula <- paste("price ~ ",variables_selected,sep = "") %>%
  as.formula()

# fit a glm model
model <- glm(formula, # formula
             data = df, # dataset
             family = ("gaussian")) # fit model

# output plots in a 2 x 2 grid 
par(mfrow = c(2,2)) 

# diagnostic plots
plot(model)

# take a closer look at our model
summary(model)

#Plot
ggplot(df, aes(x = bathrooms + grade + sqft_above +sqft_living15, y = price)) + # draw a plot
  geom_point()  + # add points
  geom_smooth(method = "glm", # plot a regression...
              method.args = list(family = "gaussian")) 

#Model performs somewhat well but there are some major outliers with very high prices
#Now to check for spatial effects


####Clean and Convert to Geodata
df.geo <- df %>%
  select(lat, long, price)
#Test for spatial autocorrelation
coordinates(df.geo) <- ~long+lat

nb <- knn2nb(knearneigh(coordinates(df.geo), k = 3))
corrI <- sp.correlogram(nb, df.geo$price, order = 5, method = 'I', zero.policy = TRUE)
plot(corrI)
corrC <- sp.correlogram(nb, df.geo$price, order = 5, method = 'C', zero.policy = TRUE)
plot(corrC)
correl <- sp.correlogram(nb, df.geo$price, order = 5, method = 'corr', zero.policy = TRUE)
plot(correl)
#Significant spatial autocorrelation at multiple spatial scales. Probably reflects different
#neighborhoods

#Local Moran's I and Autocorrelation Plots

#Jitter data to avoid duplicates

df.gd <- as.geodata(df.geo)
df.gd1 <- jitterDupCoords(df.gd, max = .01)
par(mfrow = c(1,1))
plot.gd(df.gd1)
plot(df.gd1)

df.vg<-variog(df.gd1, max.dist = .6)  # emperical semivariogram
plot(df.vg,pch=16, col="darkblue",cex=1.2)

### examine spatial trend
df.lm <- lm(price ~ long+lat, data=df.geo)
summary (df.lm)
##remove trend
# predict trend surface
trend.gd <- as.geodata(cbind(df.gd1$coords,predict(df.lm, newdata=data.frame(df.gd1$coords))))
plot.gd(trend.gd)

#make copy of data and subtract trend
mydat2.gd<-df.gd1   #make copy
mydat2.gd$data <- mydat2.gd$data - trend.gd$data  # subtract trend
plot.gd(mydat2.gd)

###examine residual data in mydat2.gd ###
qqnorm(mydat2.gd$data)

mydat2.vg<-variog(mydat2.gd, max.dist=.2)  # emperical semivariogram
plot(mydat2.vg,pch=16, col="darkblue",cex=1.2) 

#ini.cov.pars = partial sill, range
mydat2.fit <- variofit(mydat2.vg, ini.cov.pars=c(1.0e+11, .1), nugget=1.0e+6, cov.model="gaussian", weights="cressie")
plot(mydat2.vg)
lines(mydat2.fit)

### regression - kriging  ###
min(mydat2.gd$coords[,1])
#[1] 47.1559
max(mydat2.gd$coords[,1])
#[1] 47.78269
min(mydat2.gd$coords[,2])
#[1] -122.519
max(mydat2.gd$coords[,2])
#[1] -121.315
locations <- expand.grid(-122.52:-121.32, 47.16:47.78)  # make grid for predictions
krige.mod <- krige.control(type.krige="ok", obj.model=mydat2.fit)
mydat2.ok <- krige.conv(mydat2.gd, loc=locations, krige=krige.mod)    # krige residuals
image(mydat2.ok, loc=locations, col=rev(heat.colors(12)))
points(mydat2.gd, pch=19, add.to.plot=T)

summary(mydat2.ok$predict)
length(mydat2.ok$predict)

trend <- predict(df.lm, newdata=data.frame(east=locations[,1],north=locations[,2]))  # predict trend surface 
summary(trend)

mydat.rk<-mydat2.ok
mydat.rk$predict<-mydat.rk$predict+trend
summary(mydat.rk$predict)
image(mydat.rk, loc=locations, col=rev(heat.colors(12))) 
points(mydat.gd, pch=19, add.to.plot=T)
