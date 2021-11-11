##############
## Script file for running HGAMS and calculating relative transmissibility
##############

####Authored by B. Swallow

if(!require(tidyverse)){install.packages("tidyverse");require(tidyverse)}
if(!require(mgcv)){install.packages("mgcv");require(mgcv)}
if(!require(ggplot2)){install.packages("ggplot2");require(ggplot2)}
if(!require(broom)){install.packages("broom");require(broom)}
if(!require(rgdal)){install.packages("rgdal");require(rgdal)}
if(!require(emmeans)){install.packages("emmeans");require(emmeans)}




#Read in COG data and local boundaries for LTLAs
data <- read.delim("data/data.tsv")
bounds <- read.csv("data/Local_Authority_Districts_(December_2020)_UK_BUC.csv", stringsAsFactors=TRUE)
data$day<-as.numeric(as.Date(data$WeekEndDate)-as.Date(data$WeekEndDate[1])+1)

#Select only Alpha, Delta and B.1.177 variants; sum together Delta variants 
# by summing counts over each LTLA and date
d2 <- data[data$Lineage=="B.1.617.2"|data$Lineage=="B.1.1.7"|data$Lineage=="B.1.177"|data$Lineage=="AY.4.2",]
d2A <- data[data$Lineage=="B.1.177"|data$Lineage=="B.1.1.7",]
d2d <- data[data$Lineage=="AY.4"|data$Lineage=="B.1.617.2",]
d2d = d2d%>%group_by(WeekEndDate,LTLA)%>%summarise(Lineage="Delta",Count=sum(Count),day=day)
d2 = rbind(d2A,d2d)
d2=d2%>%mutate(LTLA = as.factor(LTLA))
d2=d2%>%mutate(Lineage = as.factor(Lineage))

#Match the LTLA ID to the corresponding lat/long from the boundary data
d2$lat <- bounds$LAT[match(d2$LTLA,bounds$LAD20CD)]
d2$lon <- bounds$LONG[match(d2$LTLA,bounds$LAD20CD)]
d2<-d2%>%group_by(day, Lineage)
sums = d2%>%group_by(day,lat,lon)%>%summarise(sum=sum(Count))




### HGAMs model fitting #####


#### NB: Can also use parallel versions using below and replacing "gam" function with "bam"
# specifying "cluster=cl" as an extra agrugment in the gam functions below
# Further speedups can be obtained for some models using 'discrete=T' argument with 'method=fREML'
# See help(bam) for more details
#require(parallel)  
#nc <- 2   ## cluster size, set for example portability
#if (detectCores()>1) { ## no point otherwise
#  cl <- makeCluster(nc) 
#  ## could also use makeForkCluster, but read warnings first!
#} else cl <- NULL

#
### Model S
cov_modS <- gam(Count ~ t2(lat, lon, day, Lineage, bs=c("tp", "tp", "tp", "re"),
                           k=c(10, 10, 10, 6), m=2, full=TRUE),
                data=d2, method="REML", family="nb")

#### Model I
cov_modI3 <- gam(Count ~ Lineage + te(lat, lon, day, by=Lineage,
                                      bs=c("tp", "tp", "tp"), k=c(10, 10, 10), m=2),
                 data=d2, method="REML", family="nb")


gam.check(cov_modS)
gam.check(cov_modI)

k.check(cov_modS)
k.check(cov_modI)

summary(cov_modS)
summary(cov_modI)

### Model G
cov_modG <- gam(Count ~ t2(lat, lon, day, bs=c("tp", "tp", "tp"), k=c(10, 10, 10)),
                data=d2, method="REML", family="nb")

### Model GS
cov_modGS <- bam(Count ~ t2(lat, lon, day, bs=c("tp", "tp", "tp"),
                            k=c(10, 10), m=2) +
                   t2(lat, lon, day, Lineage, bs=c("tp", "tp", "tp", "re"),
                      k=c(10, 10, 10, 6), m=2, full=TRUE),
                 data=d2, method="REML", family="nb")
### Model GI
cov_modGI <- gam(Count ~ Lineage+
                   t2(lat, lon, day, bs=c("tp", "tp", "tp"), k=c(10, 10, 10), m=2) +
                   te(lat, lon, day, by=Lineage, bs= c("tp", "tp", "tp"),
                      k=c(10, 10, 10), m=1),
                 data=d2, method="REML", family=poisson(link = "log"))


gam.check(cov_modG)
gam.check(cov_modGS)
gam.check(cov_modGI)

k.check(cov_modG)
k.check(cov_modGS)
k.check(cov_modGI)

AIC(cov_modS,cov_modI,cov_modG,cov_modGS,cov_modGI)



### Calculate (code adapted from Davies et al., 2021, Science)

M.from.delta_r = function (delta_r, g=5.5) { 
  delta_R = exp(delta_r*g)
  return( delta_R ) 
}

Rt.from.r = function(r, g=4.7, sigma=2.9) {
  k <- (sigma / g)^2
  Rt <- (1 + k * r * g)^(1 / k)
  return(Rt) }

M.from.delta_r_df = function (df, g1=5.5, g2=3.6, 
                              coln=c("M1","M1.LCL","M1.UCL",
                                     "M2","M2.LCL","M2.UCL")) { 
  df_num = df[,which(unlist(lapply(df, is.numeric))), drop=F]
  df_nonnum = df[,which(!unlist(lapply(df, is.numeric))), drop=F]
  df_out1 = apply(df_num, 2, function (delta_r) M.from.delta_r(delta_r, g1))
  if (class(df_out1)[1]=="numeric") df_out1=as.data.frame(t(df_out1), check.names=F)
  df_out2 = apply(df_num, 2, function (delta_r) M.from.delta_r(delta_r, g2))
  if (class(df_out2)[1]=="numeric") df_out2=as.data.frame(t(df_out2), check.names=F)
  df_out = data.frame(df_out1, df_out2, check.names=F)
  if (!is.null(coln)) colnames(df_out) = coln
  return( data.frame(df_nonnum, df_num, df_out, check.names=F) )
}

# Specify specific lat/longs of all LTLAs
ull <- na.omit(unique(d2[,c("lat","lon")]))

output <- NULL

# For each LTLA calculate relative multiplicative growth rate factor and 
for(i in 1:nrow(ull)){

  #Get trends for alpha v B.1.177
  mfit_emtrends <- emtrends(cov_modGI, ~ Lineage*lat*lon, "day", mode="latent",at=list(lat=unlist(ull[i,1]),lon=unlist(ull[i,2]),day=seq(as.numeric(as.Date("2020-11-07",format="%Y-%m-%d")-as.Date(data$WeekEndDate[1])),as.numeric(as.Date("2021-01-30",format="%Y-%m-%d")-as.Date(data$WeekEndDate[1])), by=2)))
  #Get trends for delta v alpha
  mfit_emtrends2 <- emtrends(cov_modGI, ~ Lineage*lat*lon, "day", mode="latent",at=list(lat=unlist(ull[i,1]),lon=unlist(ull[i,2]),day=seq(as.numeric(as.Date("2021-04-24",format="%Y-%m-%d")-as.Date(data$WeekEndDate[1])),as.numeric(as.Date("2021-07-12",format="%Y-%m-%d")-as.Date(data$WeekEndDate[1])), by=2)))
  
  #Calculate difference between trends for each comparison
  mfit_contrasts = data.frame(as.data.frame(contrast(mfit_emtrends, method="trt.vs.ctrl", ref=1, reverse=TRUE, adjust="sidak")),
                              as.data.frame(confint(contrast(mfit_emtrends, method="trt.vs.ctrl", ref=1, reverse=TRUE, adjust="sidak")))[,c("lower.CL","upper.CL")])
  mfit_contrasts2 = data.frame(as.data.frame(contrast(mfit_emtrends2, method="trt.vs.ctrl", ref=3, reverse=TRUE, adjust="sidak")),
                               as.data.frame(confint(contrast(mfit_emtrends2, method="trt.vs.ctrl", ref=3, reverse=TRUE, adjust="sidak")))[,c("lower.CL","upper.CL")])
  mfit_contrasts3 = rbind(mfit_contrasts[1,],mfit_contrasts2[1,])
  colnames(mfit_contrasts3)[which(colnames(mfit_contrasts3) %in% c("estimate","lower.CL","upper.CL"))] = 
    c("delta_r","delta_r.lower.CL","delta_r.upper.CL")
  
  mfit_contrasts = data.frame(mfit_contrasts3[,c("contrast","SE","t.ratio","p.value")],
                              M.from.delta_r_df(mfit_contrasts3[,c("delta_r","delta_r.lower.CL","delta_r.upper.CL")]))
  
  
  output <- rbind(output,mfit_contrasts) 

}

# Output estimate of mulp. rate assuming 5.5 day generation time
out2comp <- matrix(output$M1,nrow=2,byrow=T)


latlon <- expand.grid(unlist(ull[,1]),unlist(ull[,2]))
latlon <- latlon[order(latlon[,2]),]

llind <- sapply(1:311,function(i)which(latlon[,1]==unlist(ull[i,1])&latlon[,2]==unlist(ull[i,2])))


#### PLOTS #####

# Read in England shape file; downloaded from https://borders.ukdataservice.ac.uk/easy_download_data.html?data=England_ct_2011
my_spdf <- readOGR( 
  dsn= "england_ct_2011.shp" , 
  verbose=FALSE
)

summary(my_spdf)

#Transform projection into lat/long to match LTLA data
my_spdf <- spTransform(my_spdf, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

spdf_fortified <- tidy(my_spdf)


# Plot shparefile
names(latlon) <- c("lat","lon")
plt = ggplot() +
  geom_polygon(data = spdf_fortified, aes( x = long, y = lat, group = group), fill="#69b3a2", color="white")

#Add points for LTLAs where Alpha vs B.1.177 is greater than 1.5 i.e. 50% increase
plt2 = plt +
  geom_point(data = latlon[llind,][which(out2comp[1,]>1.5),], aes(x = lon, y = lat, size = out2comp[1,which(out2comp[1,]>1.5)],col= out2comp[1,which(out2comp[1,]>1.5)]), 
             shape = 20) +
  scale_fill_gradient(name="Multiplicative Rt",low = "yellow", high = "red", na.value = NA,aesthetics = "colour")  +
  labs(title="Alpha vs B.1.177")+ guides(size="none") +
  theme_void() 
ggsave('alphavssumcontr.pdf',plt2)

#Add points for LTLAs where Delta vs Alpha is greater than 1.5 i.e. 50% increase
plt3 = plt +
  geom_point(data = latlon[llind,][which(out2comp[2,]>1.5),], aes(x = lon, y = lat, size = out2comp[2,which(out2comp[2,]>1.5)],col= out2comp[2,which(out2comp[2,]>1.5)]), 
             shape = 20) +
  scale_fill_gradient(name="Multiplicative Rt",low = "yellow", high = "red", na.value = NA,aesthetics = "colour")  +
  labs(title="Delta vs Alpha")+ guides(size="none") +
  theme_void() 

ggsave('deltavsalphacontr.pdf',plt3)


