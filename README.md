# GBTM and k-means analysis of longitudinal data

  __Paul Schneider__  
  2018-04-30  
<br><br>


Step-by-Step tutorial for cluster analysis of longitudinal data using the R-packages crimCV and traj  
Simply copy-paste the code chunks or look at the code only version: 
_link_   


# Setup packages and functions

The chunk below installs and loads the required packages and also load some slightly *modified* functions from both packages.

```r
url = "https://raw.githubusercontent.com/bitowaqr/traj/master/Setup_n_functions.R"
source(url)
```

# Prepare the data

```r
# prepare data

# Data needs the following format:
  # rows = cases
  # columns = observations at time points
  # no time variables needed

# Enter file path or select manually
  test.data = read.csv("https://raw.githubusercontent.com/bitowaqr/traj/master/test.data.csv")
  
# Otherwise select your own data:
  # test.data = read.csv(file.choose())
  
# reducing the size of the demo data set
  # test.data = test.data[1:100,] # looking at a subset 

# Missing values need to be coded as negative values for crimCV
  test.data[is.na(test.data)] = -0.00001

# data frame to matrix
  df = as.matrix(test.data[,-1])
  rownames(df) = test.data$ID

# resulting data set
  kable(head(df),format="markdown")
```



|        t.1|        t.2|         t.3|        t.4|        t.5|       t.6|
|----------:|----------:|-----------:|----------:|----------:|---------:|
|  707.77330|  707.77330| 10651.47980|  247.44559|  247.44559|  224.3058|
| 3336.66827| 4892.62165|  8621.54143| 5698.54852| 5846.25382| 6487.2372|
|   -0.00001|   -0.00001|    -0.00001|  621.12652| 3312.48985| 9529.3954|
| 1957.65656|  763.94261|   949.09004| 1656.94409| 1561.56394|  375.1921|
|   -0.00001|   -0.00001|    -0.00001|   -0.00001|   -0.00001| 1822.2672|
|  257.52852|  257.52852|   257.52852|  257.52852|  257.52852|  257.5285|


# GBTM CLUSTERING: crimCV

__GBTM modeling using crimCV:__
  
  1. Set a grid to evaluate which model specifications best fit the data
  2. Choose a model and retrieve data
  3. Evaluate group mean estimates and averages
  

## set the grid of models to evaluate

```r
# Set: 
  n.cluster = c(1,2,3) # 1.which k's to evaluate, 
  p.poly = c(1,2,3) # 2. which p's to evaluate, 
  rcv = T # 3. do you want to run cross validation? T/F
```

## model evaluation

```r
# retrieve AIC, BIC and CV error for each of the models
points.for.AIC.plot = points.for.BIC.plot = points.for.cv.plot = data.frame(x=NA,value=NA,cluster=NA)
for(c in 1:length(cv.eval.list)){
  for(p in 1:length(cv.eval.list[[c]])){
    

    tryCatch({
        cv = ifelse(!is.null(cv.eval.list[[c]][[p]]$cv),cv.eval.list[[c]][[p]]$cv,NA)
        points.for.cv.plot = rbind(points.for.cv.plot,
                                   data.frame(x=p,value=cv ,cluster=c))
    }, error =function(e){})
    
    
    tryCatch({
        aic = ifelse(!is.null(cv.eval.list[[c]][[p]]$AIC),cv.eval.list[[c]][[p]]$AIC,NA)
        points.for.AIC.plot = rbind(points.for.AIC.plot,
                                    data.frame(x=p,value=aic,cluster=c))
    }, error =function(e){})
    
        tryCatch({
        bic = ifelse(!is.null(cv.eval.list[[c]][[p]]$BIC),cv.eval.list[[c]][[p]]$BIC,NA)
        points.for.BIC.plot = rbind(points.for.BIC.plot,
                                    data.frame(x=p,value=bic,cluster=c))
        }, error =function(e){})
    
  }}

points.for.AIC.plot = points.for.AIC.plot[-1,]  
points.for.BIC.plot = points.for.BIC.plot[-1,] 
points.for.cv.plot = points.for.cv.plot[-1,]

AIC.gbtm.plot = ggplot(points.for.AIC.plot) + 
  geom_line(aes(x=x,y=value,col=as.factor(cluster))) +
  #geom_point(aes(x=x,y=value,col=as.factor(cluster))) +
  geom_text(aes(x=x,y=value,col=as.factor(cluster),label=cluster),size=5) +
  xlab("Polynomial") +
  ylab("AIC")

BIC.gbtm.plot = ggplot(points.for.BIC.plot) + 
  geom_line(aes(x=x,y=value,col=as.factor(cluster))) +
  # geom_point(aes(x=x,y=value,col=as.factor(cluster))) +
  geom_text(aes(x=x,y=value,col=as.factor(cluster),label=cluster),size=5) +
  xlab("Polynomial") +
  ylab("BIC")

cv.error.gbtm.plot = ggplot(points.for.cv.plot) + 
  geom_line(aes(x=x,y=value,col=as.factor(cluster))) +
  geom_text(aes(x=x,y=value,col=as.factor(cluster),label=cluster),size=5) +
  #geom_point(aes(x=x,y=value,col=as.factor(cluster))) +
  xlab("Polynomial")  +
  ylab("LOOCV Absolute error") 

plot.legend = get_legend(AIC.gbtm.plot + theme(legend.position = "bottom"))

model.eval.plot = 
  plot_grid(
    plot_grid(cv.error.gbtm.plot + theme(legend.position = "none"),
                            BIC.gbtm.plot + theme(legend.position = "none"),
                            AIC.gbtm.plot + theme(legend.position = "none"),
                            ncol=3),
    plot.legend,nrow = 2,rel_heights = c(10,1))

# plot model evaluation
model.eval.plot
```

![](https://github.com/bitowaqr/traj/raw/master/figure-html/unnamed-chunk-5-1.png)<!-- -->

## Set the parameters for your model of choice


```r
# Select k and p 
k.set = 3
p.set = 3

# plot details
y.axis.label = "Cost per patient month"
x.axis.label = "Months before death"
plot.title = "Cost per patient month over the last 6 months before death"
```

## Retrive data from your model of choice


```r
# retrieve the final model

  # select model
  ind.k = which(grepl(k.set,names(cv.eval.list)))
  ind.p = which(grepl(p.set,names(cv.eval.list[[ind.k]])))
  # retrieve participants membership
  gbtm.members = data.frame(ID =  rownames(df),
                                 cluster = apply(summary(cv.eval.list[[ind.k]][[ind.p]]),
                                                 1,
                                                 function(x)which(x ==max(x))))

  members.per.cluster = data.frame(table(gbtm.members$cluster))
```

## Plot *estimated* trajectories


```r
# estimated trajectories

modelled.list = plot(cv.eval.list[[ind.k]][[ind.p]],size=1,plot=F)
modelled.list$time = modelled.list$time 

  model.plot.modelled = 
    ggplot(modelled.list) +
    geom_line(aes(x=time,y=value,col=cluster)) +
    scale_y_continuous(name=y.axis.label) +
    scale_x_continuous(name=x.axis.label) +
    ggtitle(plot.title) +
    scale_color_manual(lab=paste("Group ",members.per.cluster$Var1," (n=",members.per.cluster$Freq,")",sep=""),
                       values=c(2,3,4),
                       name="Estimated group trajectories") +
    theme_minimal()
  
  model.plot.modelled
```

![](https://github.com/bitowaqr/traj/raw/master/figure-html/unnamed-chunk-8-1.png)<!-- -->

## Retrieve group function terms


```r
# retrieve model function terms with intercept and * for p < .05
  long.test.dat = melt(df)
  names(long.test.dat)  = c("ID","time","value")
  long.test.dat$time = as.numeric(gsub("t.","",long.test.dat$time))
  long.test.dat = merge(long.test.dat,gbtm.members,"ID")
  long.test.dat$cluster = as.factor(long.test.dat$cluster)

  polynomial.model.results = summary(lm(value ~ -1+poly(time,p.set):cluster+cluster, long.test.dat))
  model.spec = round(polynomial.model.results$coefficients[,1],2)
  sig.model.specification = ifelse(polynomial.model.results$coefficients[,4]<0.05,"*"," ")
  model.spec = paste(model.spec,sig.model.specification,sep="")
  model.spec = formatC(model.spec)
  model.spec = matrix(data=model.spec, ncol=p.set+1)
  
  colnames(model.spec) = c("Intercept",paste("Polynomial",1:p.set))
  rownames(model.spec) = c(paste("Group",1:k.set))
  
  kable(model.spec,format="markdown")
```



|        |Intercept |Polynomial 1 |Polynomial 2 |Polynomial 3 |
|:-------|:---------|:------------|:------------|:------------|
|Group 1 |2760.21*  |133725.18*   |13492.78*    |-5078.66     |
|Group 2 |4251.83*  |96114.93*    |-47743.14*   |2492.84      |
|Group 3 |685.78*   |36158.03*    |-10915.68    |-1298.57     |

## Plot *average* group trajectories


```r
# Retrieve observed group trajectories
  long.test.dat$value[long.test.dat$value<0] = NA
  long.test.dat.means = aggregate(value ~ cluster + time, long.test.dat, function(x) mean( x , na.rm = T))
  
  model.plot.from.data = ggplot(long.test.dat.means) +
    geom_line(aes(x=time,y=value,col=cluster)) +
    scale_y_continuous(name=y.axis.label) +
    scale_x_continuous(name=x.axis.label) +
    ggtitle(plot.title) +
    scale_color_manual(lab=paste("Group ",members.per.cluster$Var1," (n=",members.per.cluster$Freq,")",sep=""),
                       values=c(2,3,4),
                       name="Observed group trajectories") +
    theme_minimal()
  
  model.plot.from.data 
```

![](https://github.com/bitowaqr/traj/raw/master/figure-html/unnamed-chunk-10-1.png)<!-- -->



## Setup for combined plot


```r
# give names to clusters? 
cluster.names = paste(
  c("First cluster", 
    "Second",
    "Third"),
    " (n=",format(as.numeric(members.per.cluster$Freq),digits=1),")",sep="")


# set a y limits to have all sub plots on the same scale?
set.y.limit = c(0,15000)

# Group average overview and individual trajectories

pop.average.traj = aggregate(value ~time, long.test.dat, function(x) mean(x,na.rm=T))
model.plot.modelled.plus.pop.average = 
    ggplot() +
  geom_line(data=modelled.list,aes(x=time,y=value,col=cluster,linetype="Estimated")) +
  geom_line(data=pop.average.traj,aes(x=time,y=value,col="Total",linetype="Average")) +  
    scale_y_continuous(name=y.axis.label) +
    scale_x_continuous(name=x.axis.label) +
    ggtitle(plot.title) +
    scale_color_manual(lab=c(paste("Group ",members.per.cluster$Var1," (n=",members.per.cluster$Freq,")",sep=""),
                             paste("Total", " (n=",sum(members.per.cluster$Freq),")",sep="")),
                       values=c(2:(k.set+1),1),
                       name="Estimated group trajectories") +
  scale_linetype_manual(lab=c("Average","Estimated"),values = c(2,1), name="")+
  guides(color = guide_legend(order = 1),
      linetype = guide_legend(order=0)) +
    theme_minimal()


  individual.plot.list = list()
  times.ex = unique(long.test.dat$time)
  for(i in 1:length(unique(long.test.dat$cluster))){
    individual.plot.list[[i]] = 
      ggplot() +
      theme_minimal() +
          geom_line(data=long.test.dat[long.test.dat$cluster==i,],
                    aes(x=time,y=value,group=ID,linetype="Average"),col=i+1,alpha=0.3,size=0.4) +
      geom_line(data=long.test.dat.means[long.test.dat.means$cluster==i,],
                    aes(x=time,y=value,linetype="Average"),col=i+1,size=1) +
       geom_line(data=modelled.list[modelled.list$cluster==i,],
                    aes(x=time,y=value,col=cluster,linetype="Estimate"),col=i+1,size=1) +
        scale_linetype_manual(lab=c("Average","Estimated"),values = c(2,1), name="") +
      theme(legend.position = "none") +
      ylab("") +
            # ggtitle(plot.titles.for.mega[i]) +
            coord_cartesian(ylim=set.y.limit) 
  }
  # individual.plot.list[[length(individual.plot.list)+1]] = get.legend
  
  gbtm.mega.plot = 
    plot_grid(model.plot.modelled.plus.pop.average+
                coord_cartesian(ylim=set.y.limit),
              plot_grid(plotlist=individual.plot.list,ncol=round(k.set/2,0)),ncol=1,rel_heights = c(2,3))
  
  gbtm.mega.plot
```

![](https://github.com/bitowaqr/traj/raw/master/figure-html/unnamed-chunk-12-1.png)<!-- -->




# k-means clustering

## k-means of what?

```r
?step1measures # info shows the 24 measurements on which k-means is performend
```

The 24 measures are:

  1. Range
  2. Mean-over-time*
  3. Standard deviation (SD)
  4. Coefficient of variation (CV)
  5. Change
  6. Mean change per unit time
  7. Change relative to the first score
  8. Change relative to the mean over time
  9. Slope of the linear model*
  10. R^2: Proportion of variance explained by the linear model
  11. Maximum of the first differences
  12. SD of the first differences
  13. SD of the first differences per time unit
  14. Mean of the absolute first differences*
  15. Maximum of the absolute first differences
  16. Ratio of the maximum absolute difference to the mean-over-time
  17. Ratio of the maximum absolute first difference to the slope
  18. Ratio of the SD of the first differences to the slope
  19. Mean of the second differences
  20. Mean of the absolute second differences
  21. Maximum of the absolute second differences
  22. Ration of the maximum absolute second difference to the mean-over-time
  23. Ratio of the maximum absolute second difference to mean absolute first difference
  24. Ratio of the mean absolute second difference to the mean absolute first difference



## Automated function for plotting results 

from `step1measures` --> `step2factors` --> `step3clusters`


```r
# Takes a 'data frame' and creates traj cluster analysis with mean and individual plots
    traj.k.mean = function( fill.data.matrix = test.data[,-1],
                            ID = test.data$ID ,
                            nclusters = NULL # set number of clusters, or NULL to let R decide
                       )
      {
      times = dim(fill.data.matrix)[2]
      colnames(fill.data.matrix) = paste("x",1:dim(fill.data.matrix)[2],sep="")
      time.mat =  matrix(data=rep(1:dim(fill.data.matrix)[2],each=dim(fill.data.matrix)[1]),
                         nrow=dim(fill.data.matrix)[1],
                         ncol=dim(fill.data.matrix)[2])
      colnames(time.mat) = paste("t",1:dim(fill.data.matrix)[2],sep="")
      cluster.data = cbind(ID = ID, fill.data.matrix,time.mat)
      
      fill.data.matrix = cbind(ID,fill.data.matrix)
      time.mat = cbind(ID,time.mat)
      
      cat("\n Clustering needas at least 5 observation per case... Cases with fewer observations are being removed...! \n")
      s1 = step1measures(Data = fill.data.matrix, 
                         Time = time.mat,
                         ID = T )
      s2all = step2factors(s1)
      s3all = step3clusters(s2all, nclusters = nclusters)  # 5 clusters
      # clust.build.plot(s3all,y.max.lim =y.max.lim)
      k.membership.table = data.frame(table(s3all$clusters$cluster))
      k.membership = s3all$clusters
      k.k = unique(k.membership$cluster)
      
      
      cluster.plot.df = data.frame(cluster=NA,value=NA,time=NA)
      for(i in k.k){
        mean.per.cluster = colMeans(fill.data.matrix[fill.data.matrix$ID %in% k.membership$ID[k.membership$cluster==i],-1],na.rm=T)
        cluster.name = paste(i," (",length(fill.data.matrix[fill.data.matrix$ID %in% k.membership$ID[k.membership$cluster==i],1]),")",sep="")
        cluster.plot.df=rbind(cluster.plot.df,
                              data.frame(cluster=cluster.name,
                                         value=mean.per.cluster,
                                         time = 1:length(mean.per.cluster)))
      }
      cluster.plot.df = cluster.plot.df[-1,]
      cluster.plot.df$time = as.numeric(cluster.plot.df$time)
      n.per.cluster = as.numeric(table(cluster.plot.df$cluster))
      cluster.plot.df$cluster = as.factor(cluster.plot.df$cluster)
      
      mean.cluster.plot = 
        ggplot(cluster.plot.df) + 
        geom_line(aes(x=time,y=value,col=cluster)) +
        guides(color=guide_legend("Cluster (n)"))  # add guide properties by aesthetic
      
      # Individual plot
      indiv.data = merge(fill.data.matrix,k.membership,"ID")
      indiv.data = melt(indiv.data,id.vars=c("ID","cluster"))
      indiv.data$variable = as.character(indiv.data$variable)
      indiv.data$variable = gsub("x","",indiv.data$variable)
      indiv.data$variable = as.numeric(indiv.data$variable)
      indiv.data$cluster = as.factor(indiv.data$cluster)
      
      lev.clust = levels(indiv.data$cluster)
      n.clust = as.numeric(by(indiv.data$ID,indiv.data$cluster,function(x)length(unique(x))))
      names.clust = data.frame(cluster = lev.clust,n = paste(lev.clust," (",n.clust,")",sep=""))
      indiv.data = merge(indiv.data,names.clust,by="cluster")
      indiv.data$ID = as.factor(indiv.data$ID)

      individual.cluster.plot = 
        ggplot(indiv.data,aes(x=variable,y=value,col= cluster)) + 
            geom_line(aes(group=ID),alpha=0.3,size=0.6) +
            stat_summary(aes(group = cluster), fun.y = mean, geom = 'line', size=1, alpha=1) +
            guides(color=guide_legend("Cluster")) 
      
      
      cluster.analysis = list(
        s1 = s1,
        s2all = s2all,
        s3all = s3all,
        membership = k.membership,
        mean.cluster.plot=mean.cluster.plot,
        individual.cluster.plot = individual.cluster.plot)
      
      return(cluster.analysis)
    }
```

## k-24-means cluster analysis 


```r
cluster.analysis.full = traj.k.mean( fill.data.matrix = test.data[,-1], ID = test.data$ID ,nclusters = NULL)
```

```
## 
##  Clustering needas at least 5 observation per case... Cases with fewer observations are being removed...! 
## [1] "Correlation of m5 and m6 : 1"
## [1] "Correlation of m12 and m13 : 1"
## [1] "Correlation of m17 and m18 : 1"
## [1] "m6 is removed because it is perfectly correlated with m5"  
## [2] "m13 is removed because it is perfectly correlated with m12"
## [1] "Computing reduced correlation e-values..."
```

![](https://github.com/bitowaqr/traj/raw/master/figure-html/unnamed-chunk-15-1.png)<!-- -->

```r
# results
cluster.analysis.full$s3
```

```
## Number of observations:  558 
## 
## Cluster distribution:
## 
##   1   2 
## 531  27 
## 
## Measures with max.loading in factors:  m5 m10 m14 m16
## 
## If you report these results, please cite:
## Sylvestre MP, et al. (2006). Classification of patterns of delirium severity scores over time in an elderly population. 
## International Psychogeriatrics,18(4), 667-680. doi:10.1017/S1041610206003334.
```

```r
cluster.analysis.full$mean.cluster.plot
```

![](https://github.com/bitowaqr/traj/raw/master/figure-html/unnamed-chunk-15-2.png)<!-- -->

```r
cluster.analysis.full$individual.cluster.plot
```

![](https://github.com/bitowaqr/traj/raw/master/figure-html/unnamed-chunk-15-3.png)<!-- -->

## k-means clusters within GBTM clusters


```r
# k-means clusters within GBTM clusters
unique.GBTM.clusters = unique(long.test.dat$cluster)

clusters.within.clusters.plots = list()
set.y.limit = set.y.limit 

for(k in unique.GBTM.clusters){
  temp.pat = unique(long.test.dat$ID[long.test.dat$cluster==k])
  temp.fill.data = test.data[test.data$ID %in% temp.pat,-1]
  temp.k.means = traj.k.mean( fill.data.matrix = temp.fill.data, ID = temp.pat ,nclusters = NULL)
  clusters.within.clusters.plots[[k]] = temp.k.means$individual.cluster.plot + coord_cartesian(ylim=set.y.limit) 
}
```

```
## 
##  Clustering needas at least 5 observation per case... Cases with fewer observations are being removed...! 
## [1] "Correlation of m5 and m6 : 1"
## [1] "Correlation of m12 and m13 : 1"
## [1] "Correlation of m17 and m18 : 1"
## [1] "m6 is removed because it is perfectly correlated with m5"  
## [2] "m13 is removed because it is perfectly correlated with m12"
## [1] "Computing reduced correlation e-values..."
```

![](https://github.com/bitowaqr/traj/raw/master/figure-html/unnamed-chunk-16-1.png)<!-- -->

```
## 
##  Clustering needas at least 5 observation per case... Cases with fewer observations are being removed...! 
## [1] "Correlation of m5 and m6 : 1"
## [1] "Correlation of m12 and m13 : 1"
## [1] "m6 is removed because it is perfectly correlated with m5"  
## [2] "m13 is removed because it is perfectly correlated with m12"
## [1] "Computing reduced correlation e-values..."
```

![](https://github.com/bitowaqr/traj/raw/master/figure-html/unnamed-chunk-16-2.png)<!-- -->

```
## 
##  Clustering needas at least 5 observation per case... Cases with fewer observations are being removed...! 
## [1] "Correlation of m5 and m6 : 1"
## [1] "Correlation of m12 and m13 : 1"
## [1] "Correlation of m17 and m18 : 0.999"
## [1] "m6 is removed because it is perfectly correlated with m5"  
## [2] "m13 is removed because it is perfectly correlated with m12"
## [1] "Computing reduced correlation e-values..."
```

![](https://github.com/bitowaqr/traj/raw/master/figure-html/unnamed-chunk-16-3.png)<!-- -->

## Overiew plot


```r
plot_grid(plotlist = clusters.within.clusters.plots)
```

![](https://github.com/bitowaqr/traj/raw/master/figure-html/unnamed-chunk-17-1.png)<!-- -->

## References
  
  Jason D. Nielsen (2013). crimCV: Group-Based Modelling of Longitudinal Data. R package version 0.9.3.
  https://CRAN.R-project.org/package=crimCV
  
  Sylvestre MP, et al. (2006). Classification of patterns of delirium severity scores over time in an elderly population.
  International Psychogeriatrics, 18(4), 667-680. doi:10.1017/S1041610206003334.

  Leffondree, K. et al. (2004). Statistical measures were proposed for identifying longitudinal patterns of change in
  quantitative health indicators. Journal of Clinical Epidemiology, 57, 1049-1062. doi : 10.1016/j.jclinepi.2004.02.012.



## end



