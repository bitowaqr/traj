# setup
# load and install libraries
  required_packages<-c("crimCV","foreign","ggplot2","cowplot","reshape2","dplyr","traj","cowplot","knitr") 

  pft_packages <- function(package){
      for(i in 1:length(package)){
        if(eval(parse(text=paste("require(",package[i],")")))==0) {
          install.packages(package)}}
      return (eval(parse(text=paste("require(",package,")"))))}
  
  pft_packages(required_packages)

# some convenient functions


load.traj.data = function(ID.index = NULL,
                            time.points = NULL,
                            path = NULL){
    
    if(is.null(path)) {
      cat("\n \n Select your dta or csv data file \n")
      cat("Works with data organized like this (with or withour ID) \n \n")
      print(head(sample_traj_data)[1:3,1:4])
      
      path = file.choose(file.path())
    }
    
    ending = substr(path,nchar(path)-2,nchar(path))
    if(ending == "csv"){
      dat = read.csv(path)
    }
    if(ending == "dta"){
      dat = read.dta(path)
    }
    
    # Set ID column
    cat("\n ")
    if(is.null(ID.index)){
    ID.index = readline(prompt="Is the first column an ID variable? y/n  " )}
    if(ID.index %in% c("yes","y","Yes","Y","N","n","No","no")){
      if(ID.index == "yes" |ID.index == "Yes" | ID.index == "y" | ID.index == "Y"){
        ID = dat[,1]
        dat =dat[,-1]
      }
      
      if(ID.index == "no" | ID.index == "No" | ID.index == "n" | ID.index == "N"){
        cat("\n ID variable created")
        ID = 1:length(dat[,1])
      }
    } else{stop("\n Invalid input <-  Types y or n")}
    
    # Set time points
    if(is.null(time.points)){
    time.points = readline(prompt="How many observations per case?  " )}
    time.points = as.numeric(time.points)
    if(time.points!=round(time.points)) stop(" Time points must be an integer")
    if(time.points > length(dat)) stop(" Data does not have enough columns")
    mat = dat[,1:time.points]
    
    # make data set
    traj_data  = cbind(ID,mat)
    names(traj_data) = c("ID",paste("t",1:time.points,sep=""))
    cat("\n Data set with dimensions: cases=",dim(traj_data)[1],"; timepoints=",dim(traj_data)[2]-1," created!",sep="")
    
    # missing data must be replaced with negative values
    na.sum = sum(is.na(traj_data))
    if(na.sum >0){
      traj_data[is.na(traj_data)] = -0.001
      cat("\n",na.sum,"missing values are being replaces with negative values!")
    }
    
    return(traj_data)
  }  
  
  


# Plot model fit stats
     fit_eval = function(eval.stats.data = eval.stats){
        eval.stats.data$p = as.integer(eval.stats.data$p)
        types = unique(eval.stats.data$type)
        eval.data.list = list(); eval.plot.list = list()
        for(i in types){
          eval.data.list[[i]] = eval.stats.data[which(eval.stats.data$type==i),]
          eval.plot.list[[i]] = ggplot(eval.data.list[[i]]) +   
          geom_line(aes(x=p,y=value,col=as.factor(k))) +
          #geom_point(aes(x=x,y=value,col=as.factor(cluster))) +
          geom_text(aes(x=p,y=value,col=as.factor(k),label=k),size=5) +
          xlab("Polynomial") +
          ylab(paste(i)) +
            scale_color_manual(name="Number of groups",values=unique(eval.data.list[[i]]$k)) +
            scale_x_continuous(breaks = c(min(eval.data.list[[i]]$p),max(eval.data.list[[i]]$p)))
          }
        return(list(eval.data.list = eval.data.list,
                    eval.plot.list = eval.plot.list))
        }
# fit multiple crimCV GBTM models over k and p grid
fit.gbtm = function(data = traj_data,
           evaluate.k.groups = NULL,
           evaluate.p.polynomials = NULL,
           run.loocv = NULL,
           method = NULL,
           init.proc = 20,
           Risk.set = NULL,
           ...){
    
    if(is.null(evaluate.k.groups)) {
      evaluate.k.groups = readline(prompt="Which number of groups to evaluate? For multiple, type: 1,2,...,k  " )
      evaluate.k.groups <- as.integer(strsplit(evaluate.k.groups, ",")[[1]])
      }
    if(is.null(evaluate.p.polynomials)){ 
      evaluate.p.polynomials =  readline(prompt="Which degrees of polynomials to evaluate? For multiple, type: 1,2,...,p  " )
      evaluate.p.polynomials <- as.integer(strsplit(evaluate.p.polynomials, ",")[[1]])
    }
    if(is.null(run.loocv)) run.loocv = as.logical(readline(prompt="Run leave-one-out cross-validation? Type True or False  " ))
    if(is.null(method)) method = readline(prompt="Which method to use - type 'ZIP' or 'ZIPt'?  " )
    
    
    cv.eval.list = list()
    eval.stats = data.frame(k=NA,p=NA,value=NA,type=NA)
    index = 0
    for(k in evaluate.k.groups){
      for(p in evaluate.p.polynomials){
        
        index = index + 1
        cv.eval.list[[index]] = list()
        attributes(cv.eval.list[[index]])$Groups = k
        attributes(cv.eval.list[[index]])$Poly = p
        names(cv.eval.list)[index] = paste("k",k,"p",p,sep="")
        
        cat("\n Fitting GBTM for k=", k, " and p=",p,"... \n",sep="")
        
        tryCatch({
          crimCV.model.temp = crimCV(Dat = as.matrix(data[,-1]),
                                     ng = k,
                                     dpolyp = p,
                                     dpolyl = p,
                                     rcv = run.loocv,     
                                     model = method,
                                     init = init.proc,
                                     Risk = Risk.set
          )
          cv.eval.list[[index]] = crimCV.model.temp
          
          tryCatch({
            
            AIC.temp = crimCV.model.temp$AIC
            eval.stats = rbind(eval.stats,data.frame(k=k,p=p,value=AIC.temp,type="AIC"))
            
            BIC.temp = crimCV.model.temp$BIC
            eval.stats = rbind(eval.stats,data.frame(k=k,p=p,value=BIC.temp,type="BIC"))
            
            if(with(crimCV.model.temp,exists("cv"))){
              RCV.temp = crimCV.model.temp$cv  
              eval.stats = rbind(eval.stats,data.frame(k=k,p=p,value=RCV.temp,type="CV"))
            }
            
            
            
          }, error = function(e){cat("Something went wrong with eval. stats. of k=",k," p=",p," - Continue")}
          )
        }, error = function(e){cat("Something went wrong while fitting k=",k," p=",p," - Continue")}
        )
        
      }
    }
    
    eval.stats = eval.stats[-1,]
    
    tryCatch({eval.stats = fit_eval(eval.stats.data = eval.stats)}, error = function(e){cat("")})
    
    
    return(list(cv.eval.list = cv.eval.list,
                eval.stats = eval.stats))
  }
  
 # Select model for further analysis
set.model = function(models.list = fitted.gbtm$cv.eval.list,
                       set.k = NULL,
                       set.p = NULL){
    cat("\n Select the model you want to further evaluate:")
    if(is.null(set.k)) set.k = readline(prompt="How many groups?  " )
    if(is.null(set.p)) set.p = readline(prompt="What degree of polynomial?  " )
    set.available = grepl(paste("k",set.k,"p",set.p,sep=""),names(models.list))
    if(sum(set.available) == 0 ) stop(" This combination of groups and polynomials is not available.")
    select.model = models.list[set.available][[1]]
    return(select.model)
  }

# get model terms from data

    get.model.terms = function(model = select.model, 
                               data = traj_data){
      
      membership = data.frame(ID = data$ID,
                                                 group = apply(summary(model),1,function(x)which(x == max(x))))
      traj_data_long = melt(data,id.vars="ID")
      names(traj_data_long)  = c("ID","time","value")
      traj_data_long$time = as.numeric(gsub("t","",traj_data_long$time))
      traj_data_long = merge(traj_data_long,membership,"ID")
      traj_data_long$group = as.factor(traj_data_long$group)
      traj_data_long$value[traj_data_long$value<0] = NA
      
      p.poly = dim(model$beta)[1]-1
      k.group = dim(model$beta)[2]
      
      model.spec = matrix(data=NA, ncol=p.poly+1,nrow=k.group)
      colnames(model.spec) = c("Intercept",paste("Polynomial",1:p.poly))
      rownames(model.spec) = c(paste("Group",1:k.group))
      
      for(k in 1:k.group){
        temp = summary(lm(value ~ 1+poly(time,p.poly,raw=T), traj_data_long[traj_data_long$group==k,]))
        temp.coef = round(temp$coefficients[,1],2)
        significance = ifelse(temp$coefficients[,4]<0.05,"*"," ")
        significance.coef = paste(temp.coef,significance,sep="")
        temp.input = formatC(significance.coef)
        model.spec[k,] = temp.input
      }
     
      return(model.spec)
    }
                                                             
                                                             
 # retrieve probabilistic and deterministic group membership
    gbtm.members = function(model = select.model,
                            data = traj_data){
      # probabilistic membership
      prob.membership = summary(model)
      # deterministic membership
      deter.membership = data.frame(ID = data$ID,
                                    group = apply(prob.membership,1,function(x)which(x == max(x))))
      # Group membership frequency table
      deter.membership.table = as.data.frame(table(deter.membership$group))
      names(deter.membership.table) = c("group","n_members")
      member_list = list(prob.membership = prob.membership,
                         deter.membership = deter.membership,
                         deter.membership.table = deter.membership.table)
      return(member_list)
    }

# plot observed average traj per group
    
  plot.mean.per.group = function(model = select.model, 
                                 data = traj_data,
                                 y.axis.label = "Y-axis label",
                                 x.axis.label = "X-axis label",
                                 plot.title = "Plot title",
                                 plot.total = T){
    membership = data.frame(ID = data$ID,
                            group = apply(summary(model),1,function(x)which(x == max(x))))
    n_members = data.frame(table(membership$group))
    names(n_members) = c("group","n_members")
    
    traj_data_long = melt(data,id.vars="ID")
    names(traj_data_long)  = c("ID","time","value")
    traj_data_long$time = as.numeric(gsub("t","",traj_data_long$time))
    traj_data_long = merge(traj_data_long,membership,"ID")
    traj_data_long$group = as.factor(traj_data_long$group)
    
    p.poly = as.numeric(attributes(model)$p)
    k.group = as.numeric(attributes(model)$k)
    
    traj_data_long$value[traj_data_long$value<0] = NA
    long.dat.means = aggregate(value ~ group + time, traj_data_long, function(x) mean( x , na.rm = T))
    pop.average.traj =  aggregate(value ~ time, traj_data_long, function(x) mean( x , na.rm = T))
    
    model.plot.from.data = ggplot(long.dat.means) +
      geom_line(aes(x=time,y=value,col=group)) +
      scale_color_manual(lab=paste("Group ",n_members$group," n=",n_members$n_members," (",round(n_members$n_members/sum(n_members$n_members),2)*100,"%)",sep=""),
                         values=as.numeric(as.character(n_members$group))+1,
                         name="Average group trajectories") +
      theme_minimal()
    if(plot.total ==T){
      
      
      model.plot.from.data = 
        model.plot.from.data +
        geom_line(data = pop.average.traj, aes(x=as.numeric(time),y=value,group="Total",col="Total")) +
        scale_color_manual(lab=c(paste("Group ",n_members$group," n=",n_members$n_members," (",round(n_members$n_members/sum(n_members$n_members),2)*100,"%)",sep=""),
                                 paste("Total n=",sum(n_members$n_members)," (100%)")),
                           values=c(as.numeric(as.character(n_members$group))+1,1),
                           name="Average group trajectories") 
    }
    
    return(model.plot.from.data)
  }

                  
# plot group traj overview and individual groups
                             

plot.gbtm.groups = function(data = traj_data,
                             model = select.model,
                             set.y.limit = NULL,  # c(0,40)
                             individual.group.names = NULL, # give names
                             plot.total = T,
                             plot_overview = T,
                             n_row = NULL,
                             y.axis.label = "y axis label",
                             x.axis.label = "x axis label",
                             plot.title = ""){
    
    membershiplist = gbtm.members(data = data,model = model)
    membership = membershiplist$deter.membership
    membershiptable = membershiplist$deter.membership.table
    
    # modeled data
    modelled.list = plot(model,size=1,plot=F) 
    
    # long data
    traj_data_long = melt(data,id.vars="ID")
    names(traj_data_long)  = c("ID","time","value")
    traj_data_long$time = as.numeric(gsub("t","",traj_data_long$time))
    traj_data_long = merge(traj_data_long,membership,"ID")
    traj_data_long$group = as.factor(traj_data_long$group)
    
    # population average
    pop.average.traj = aggregate(value ~time, traj_data_long, function(x) mean(x,na.rm=T))
    
    #
    p.poly = as.numeric(attributes(model)$p)
    k.group = as.numeric(attributes(model)$k)
    
    traj_data_long$value[traj_data_long$value<0] = NA
    long.dat.means = aggregate(value ~ group + time, traj_data_long, function(x) mean( x , na.rm = T))
    
    ######################## HALTED ############################model.plot.modelled.plus.pop.average = 
    
    
    model.plot.modelled.plus.pop.average =  
      ggplot() +
      geom_line(data=modelled.list,aes(x=time,y=value,col=cluster,linetype="Estimated")) +
      scale_y_continuous(name=y.axis.label) +
      scale_x_continuous(name=x.axis.label) +
      ggtitle(plot.title) +
        scale_color_manual(lab=c(paste("Group ",membershiptable[,1]," n=",membershiptable[,2]," (",round(membershiptable[,2]/sum(membershiptable[,2]),2)*100,"%)",sep="")),
                           values=c(as.numeric(as.character(membershiptable[,1]))+1,1),
                           name="Estimated group trajectories") +
      scale_linetype_manual(lab=c("Estimated"),values = c(1), name="")+
      guides(color = guide_legend(order = 1),
             linetype = guide_legend(order=0)) +
      theme_minimal()
    
    
    if(plot.total == T) {
      model.plot.modelled.plus.pop.average = 
        model.plot.modelled.plus.pop.average + 
            geom_line(data=pop.average.traj,aes(x=time,y=value,col="Total",linetype="Average")) +
        scale_color_manual(lab=c(paste("Group ",membershiptable[,1]," n=",membershiptable[,2]," (",round(membershiptable[,2]/sum(membershiptable[,2]),2)*100,"%)",sep=""),
                                 paste("Total n=",sum(membershiptable[,2])," (100%)",sep="")),
                           values=c(as.numeric(as.character(membershiptable[,1]))+1,1),
                           name="Estimated group trajectories") +
        scale_linetype_manual(lab=c("Average","Estimated"),values = c(2,1), name="")}
    
    
    if(is.null(individual.group.names)) individual.group.names = paste("Group ",membershiptable[,1]," n=",membershiptable[,2]," (",round(membershiptable[,2]/sum(membershiptable[,2]),2)*100,"%)",sep="")
    
    individual.plot.list = list()
    times.ex = unique(traj_data_long$time)
    for(i in 1:length(unique(traj_data_long$group))){
      individual.plot.list[[i]] = 
        ggplot() +
        theme_minimal() +
        geom_line(data=traj_data_long[traj_data_long$group==i,],
                  aes(x=time,y=value,group=ID,linetype="Estimate"),col=i+1,alpha=0.3,size=0.4) +
        geom_line(data=long.dat.means[long.dat.means$group==i,],
                  aes(x=time,y=value,linetype="Average"),col=i+1,size=1) +
        geom_line(data=modelled.list[modelled.list$cluster==i,],
                  aes(x=time,y=value,col=cluster,linetype="Estimate"),col=i+1,size=1) +
        scale_linetype_manual(lab=c("Average","Estimated"),values = c(2,1), name="") +
        theme(legend.position = "none") +
        ylab("") +
        xlab("") +
        ggtitle(individual.group.names[i]) +
        coord_cartesian(ylim=set.y.limit) 
    }
    # individual.plot.list[[length(individual.plot.list)+1]] = get.legend
    
    
    if(is.null(n_row)) n_row = round(length(membershiptable[,1])/2,0)
    
    if(plot_overview == T){
    p1 = 
      plot_grid(model.plot.modelled.plus.pop.average+
                  coord_cartesian(ylim=set.y.limit),
                plot_grid(plotlist=individual.plot.list,
                          nrow=n_row),
                ncol=1,rel_heights = c(2,3))
    } else {
      p1 = 
                  plot_grid(plotlist=individual.plot.list,
                            nrow=n_row)
    }
    return(p1)
}

                             
                             
                             
                             
                             
                             
                             
                             

# load adjusted crimCV functions
  
  # crimCV IMPROVEMENTS:
    # update plot.zip --> makes better plots!
    # update dmZIPt --> doesnt print cross validation trail
    # update summary.dmZIP --> doesnt print stuff
  
    plot.dmZIPt = function (x, size=1,plot=T...) {
      prob <- x$prob
      Xb <- x$X %*% x$beta
      lambda <- exp(Xb)
      p <- exp(-x$tau * t(Xb))
      p <- t(p)
      p <- p/(1 + p)
      mu <- (1 - p) * lambda
      tt <- 1:nrow(mu)
      mu = melt(mu)
      names(mu) = c("time","cluster","value")
      mu$cluster = as.factor(mu$cluster)
      if(plot==T){
      ggplot(mu) +
        geom_line(aes(x=time, y=value,col=cluster),size=size)
      } else {
        return(mu)
      }
    }
      plot.dmZIP <-
        function (x, size=1,plot=T,...) {
        prob <- x$prob
        Xb <- x$X %*% x$beta
        Zg <- x$Z %*% x$gamma
        lambda <- exp(Xb)
        p <- exp(Zg)
        p <- p/(1 + p)
        mu <- (1 - p) * lambda
        tt <- 1:nrow(mu)
        mu = melt(mu)
        names(mu) = c("time","cluster","value")
        mu$cluster = as.factor(mu$cluster)
        table = as.numeric(table(mu$cluster))
        if(plot==T){
      ggplot(mu) +
        geom_line(aes(x=time, y=value,col=cluster),size=size)
      } else {
        return(mu)
      }
    }
    dmZIPt = 
      function (Dat, X, ng, rcv = TRUE, init = 20, Risk = NULL) 
        {
        initpr <- function(Dat, X, Risk, ni, no, npp, ng, npop) {
            pparam <- rep(0, (npp + 2) * ng * npop)
            pllike <- rep(0, npop)
            Datm <- Dat
            Datm[Datm < 0] <- 0
            Dmu <- apply(Datm, 2, mean)
            Dmu <- rep(Dmu, ni)
            Dmu <- matrix(Dmu, no, ni)
            Dmu <- t(Dmu)
            Dmu[Dat > -0.1] <- 0
            Datm <- Datm + Dmu
            Frtr <- .Fortran("r_dmzipt_init_param", as.double(X), 
                as.double(Datm), as.double(Risk), pparam = as.double(pparam), 
                pllike = as.double(pllike), as.integer(ni), as.integer(no), 
                as.integer(npp), as.integer(ng), as.integer(npop))
            out <- NULL
            out$pparam <- matrix(Frtr$pparam, npop, (npp + 2) * ng)
            out$pllike <- Frtr$pllike
            out
        }
        ni <- nrow(Dat)
        no <- ncol(Dat)
        npp <- ncol(X)
        if (nrow(X) != no) 
            stop("dmZIPt: Incorrect input nrows(X) must equal ncol(Dat)!")
        if (is.null(Risk)) {
            Risk <- matrix(1, ni, no)
        }
        nn <- ni * no
        npr <- (npp + 1) * ng + ng
        cat("|----Initialize----:")
        flush.console()
        npop <- init * ng
        outi <- initpr(Dat, X, Risk, ni, no, npp, ng, npop)
        if (is.null(outi)) {
            cat(" [FAILED] \n")
            return(outi)
        }
        cat(" [DONE] \n")
        flush.console()
        pparam <- matrix(outi$pparam, npop, (npp + 2) * ng)
        bparam <- NULL
        bllike <- NULL
        llike <- -1e+20
        cat("|-----Fitting------:")
        flush.console()
        for (i in 1:npop) {
            param <- pparam[i, ]
            param <- matrix(param, npp + 2, ng)
            beta <- param[1:npp, ]
            tau <- param[npp + 1, ]
            prob <- param[npp + 2, ]
            ggt <- matrix(0, ni, ng)
            Info <- matrix(0, npp * ng + 2 * ng - 1, npp * ng + 2 * 
                ng - 1)
            tFrtr <- .Fortran("r_dmzipt", as.double(X), as.double(Dat), 
                as.double(Risk), ggt = as.double(ggt), prob = as.double(prob), 
                beta = as.double(beta), tau = as.double(tau), llike = as.double(0), 
                Info = as.double(Info), as.integer(nn), as.integer(ni), 
                as.integer(no), as.integer(npp), as.integer(ng), 
                err = as.integer(0))
            if (tFrtr$err != 0) 
                next
            if (tFrtr$llike > llike) {
                Frtr <- tFrtr
                llike <- tFrtr$llike
            }
            if (i == 1) {
                beta <- matrix(tFrtr$beta, npp, ng)
                tau <- exp(tFrtr$tau)
                prob <- tFrtr$prob
                param <- rbind(beta, tau, prob)
                param <- c(param)
                bparam <- cbind(bparam, param)
                bllike <- c(bllike, tFrtr$llike)
            }
            else {
                if (!(min(abs(tFrtr$llike - bllike)) < 1e-05)) {
                    beta <- matrix(tFrtr$beta, npp, ng)
                    tau <- exp(tFrtr$tau)
                    prob <- tFrtr$prob
                    param <- rbind(beta, tau, prob)
                    param <- c(param)
                    bparam <- cbind(bparam, param)
                    bllike <- c(bllike, tFrtr$llike)
                }
            }
        }
        if (length(bllike) > 1) {
            rllike <- order(-bllike)
            bllike <- bllike[rllike]
            bparam <- bparam[, rllike]
        }
        out <- NULL
        if (Frtr$err != 0) {
            cat(" [FAILED] \n")
            return(out)
        }
        cat(" [DONE] \n")
        flush.console()
        out$beta <- matrix(Frtr$beta, npp, ng)
        out$tau <- exp(Frtr$tau)
        out$prob <- Frtr$prob
        out$gwt <- matrix(Frtr$ggt, ni, ng)
        out$Info <- matrix(Frtr$Info, npp * ng + 2 * ng - 1, npp * 
            ng + 2 * ng - 1)
        llike <- Frtr$llike
        out$llike <- llike
        out$AIC <- -2 * llike + 2 * (ng * (npp + 1) + ng - 1)
        out$BIC <- -2 * llike + (ng * (npp + 1) + ng - 1) * log(nn)
        out$bparam <- bparam
        out$bllike <- bllike
        out$ni <- ni
        out$ng <- ng
        out$X <- X
        if (rcv) {
            cat("|---CV/Jackknife---: \n")
            flush.console()
            param <- matrix(0, npr, ni)
            parm <- c(c(out$beta), out$tau, out$prob)
            iseq <- 1:ni
            imask <- rep(TRUE, ni)
            nn <- (ni - 1) * no
            cvllike <- rep(0, ni)
            cvDat <- matrix(0, ni, no)
            cvDat2 <- matrix(0, ni, no)
            Info <- matrix(0, npp * ng + 2 * ng - 1, npp * ng + 2 * 
                ng - 1)
            errcnt <- 0
            for (i in 1:ni) {
                # cat("                   : i =", i, "\n")
                flush.console()
                imask[i] <- FALSE
                indi <- iseq[imask]
                tDat <- Dat[indi, ]
                tRisk <- Risk[indi, ]
                ggt <- out$gwt[indi, ]
                beta <- out$beta
                tau <- log(out$tau)
                if (ng > 1) {
                    prob <- out$prob
                }
                else {
                    prob <- 1
                }
                cFrtr <- .Fortran("r_dmzipt", as.double(X), as.double(tDat), 
                    as.double(tRisk), ggt = as.double(ggt), prob = as.double(prob), 
                    beta = as.double(beta), tau = as.double(tau), 
                    llike = as.double(0), Info = as.double(Info), 
                    as.integer(nn), as.integer(ni - 1), as.integer(no), 
                    as.integer(npp), as.integer(ng), err = as.integer(0))
                if (cFrtr$err != 0) {
                    param[, i] <- c(out$beta, out$tau, out$prob)
                    cvDat[i, ] <- NA
                    errcnt <- errcnt + 1
                    imask[i] <- TRUE
                    next
                }
                beta <- cFrtr$beta
                tau <- exp(cFrtr$tau)
                prob <- cFrtr$prob
                param[, i] <- c(beta, tau, prob)
                beta <- matrix(beta, npp, ng)
                Xb <- X %*% beta
                Zg <- -tau * t(Xb)
                Xb <- log(Risk[i, ]) + Xb
                Zg <- t(Zg)
                lam <- exp(Xb)
                qq <- exp(Zg)
                qq <- qq/(1 + qq)
                mu <- (1 - qq) * lam
                mu2 <- mu
                llike = log(exp(Zg) + exp(-exp(Xb)))
                DD <- Dat[i, ]
                DD[DD < 0] <- 0
                llike[DD > 0.5, ] <- 0
                tllike = c(DD) * Xb - exp(Xb) - lgamma(DD + 1)
                tllike[DD < 0.5, ] <- 0
                llike = llike + tllike
                llike = llike - log(1 + exp(Zg))
                llike[Dat[i, ] < 0, ] <- 0
                ggt <- apply(llike, 2, sum)
                ggt <- log(prob) + ggt
                rscl <- max(ggt)
                ggt <- exp(ggt - rscl)
                cvllike[i] <- sum(ggt)
                ggt <- ggt/cvllike[i]
                cvllike[i] <- log(cvllike[i]) + rscl
                for (j in 1:ng) {
                    mu[, j] <- ggt[j] * mu[, j]
                    mu2[, j] <- prob[j] * mu[, j]
                }
                cvDat[i, ] <- apply(mu, 1, sum)
                cvDat2[i, ] <- apply(mu2, 1, sum)
                imask[i] <- TRUE
            }
            cat("|---CV/Jackknife---: [DONE] \n")
            flush.console()
            cv <- mean(abs(Dat - cvDat), na.rm = TRUE)
            cvunc <- mean(abs(Dat - cvDat2), na.rm = TRUE)
            cv2 <- mean((Dat - cvDat)^2, na.rm = TRUE)
            dparm <- ni * parm - (ni - 1) * param
            jVar <- var(t(dparm))/ni
            out$cv <- cv
            out$cvunc <- cvunc
            out$cv2 <- cv2
            out$lcv <- -2 * sum(cvllike)
            out$param <- param
            out$jVar <- jVar
            out$errcnt <- errcnt
        }
        class(out) <- c("dmZIPt", class(out))
        out
    }
    dmZIP = 
      function (Dat, X, Z, ng, rcv = TRUE, init = 20, Risk = NULL) 
        {
        initpr <- function(Dat, X, Z, Risk, ni, no, npp, npl, ng, 
            npop) {
            pparam <- rep(0, (npp + npl + 1) * ng * npop)
            pllike <- rep(0, npop)
            Datm <- Dat
            Datm[Datm < 0] <- 0
            Dmu <- apply(Datm, 2, mean)
            Dmu <- rep(Dmu, ni)
            Dmu <- matrix(Dmu, no, ni)
            Dmu <- t(Dmu)
            Dmu[Dat > -0.1] <- 0
            Datm <- Datm + Dmu
            Frtr <- .Fortran("r_dmzip_init_param", as.double(X), 
                as.double(Z), as.double(Datm), as.double(Risk), pparam = as.double(pparam), 
                pllike = as.double(pllike), as.integer(ni), as.integer(no), 
                as.integer(npp), as.integer(npl), as.integer(ng), 
                as.integer(npop))
            out <- NULL
            out$pparam <- matrix(Frtr$pparam, npop, (npp + npl + 
                1) * ng)
            out$pllike <- Frtr$pllike
            out
        }
        ni <- nrow(Dat)
        no <- ncol(Dat)
        npp <- ncol(X)
        npl <- ncol(Z)
        if (nrow(X) != nrow(Z)) 
            stop("dmZIP: Incorrect input nrows(X) must equal ncol(Z)!")
        if (nrow(X) != no) 
            stop("dmZIP: Incorrect input nrows(X) must equal ncol(Dat)!")
        if (is.null(Risk)) {
            Risk <- matrix(1, ni, no)
        }
        nn <- ni * no
        npr <- (npp + npl) * ng + ng
        cat("|----Initialize----:")
        flush.console()
        npop <- init * ng
        outi <- initpr(Dat, X, Z, Risk, ni, no, npp, npl, ng, npop)
        if (is.null(outi)) {
            cat(" [FAILED] \n")
            return(outi)
        }
        cat(" [DONE] \n")
        flush.console()
        pparam <- matrix(outi$pparam, npop, (npp + npl + 1) * ng)
        bparam <- NULL
        bllike <- NULL
        llike <- -1e+20
        cat("|-----Fitting------:")
        flush.console()
        for (i in 1:npop) {
            param <- pparam[i, ]
            param <- matrix(param, npp + npl + 1, ng)
            beta <- param[1:npp, ]
            gamma <- param[(npp + 1):(npp + npl), ]
            prob <- param[npp + npl + 1, ]
            ggt <- matrix(0, ni, ng)
            Info <- matrix(0, (npp + npl) * ng + ng - 1, (npp + npl) * 
                ng + ng - 1)
            tFrtr <- .Fortran("r_dmzip", as.double(X), as.double(Z), 
                as.double(Dat), as.double(Risk), ggt = as.double(ggt), 
                prob = as.double(prob), beta = as.double(beta), gamma = as.double(gamma), 
                llike = as.double(0), Info = as.double(Info), as.integer(nn), 
                as.integer(ni), as.integer(no), as.integer(npp), 
                as.integer(npl), as.integer(ng), err = as.integer(0))
            if (tFrtr$err != 0) 
                next
            if (tFrtr$llike > llike) {
                Frtr <- tFrtr
                llike <- tFrtr$llike
            }
            if (i == 1) {
                beta <- matrix(tFrtr$beta, npp, ng)
                gamma <- matrix(tFrtr$gamma, npl, ng)
                prob <- tFrtr$prob
                param <- rbind(beta, gamma, prob)
                param <- c(param)
                bparam <- cbind(bparam, param)
                bllike <- c(bllike, tFrtr$llike)
            }
            else {
                if (!(min(abs(tFrtr$llike - bllike)) < 1e-05)) {
                    beta <- matrix(tFrtr$beta, npp, ng)
                    gamma <- matrix(tFrtr$gamma, npl, ng)
                    prob <- tFrtr$prob
                    param <- rbind(beta, gamma, prob)
                    param <- c(param)
                    bparam <- cbind(bparam, param)
                    bllike <- c(bllike, tFrtr$llike)
                }
            }
        }
        if (length(bllike) > 1) {
            rllike <- order(-bllike)
            bllike <- bllike[rllike]
            bparam <- bparam[, rllike]
        }
        out <- NULL
        if (Frtr$err != 0) {
            cat(" [FAILED] \n")
            return(out)
        }
        cat(" [DONE] \n")
        flush.console()
        out$beta <- matrix(Frtr$beta, npp, ng)
        out$gamma <- matrix(Frtr$gamma, npl, ng)
        out$prob <- Frtr$prob
        out$gwt <- matrix(Frtr$ggt, ni, ng)
        out$Info <- matrix(Frtr$Info, (npp + npl) * ng + ng - 1, 
            (npp + npl) * ng + ng - 1)
        llike <- Frtr$llike
        out$llike <- llike
        out$AIC <- -2 * llike + 2 * (ng * (npp + npl) + ng - 1)
        out$BIC <- -2 * llike + (ng * (npp + npl) + ng - 1) * log(nn)
        out$bparam <- bparam
        out$bllike <- bllike
        out$ni <- ni
        out$ng <- ng
        out$X <- X
        out$Z <- Z
        if (rcv) {
            cat("|---CV/Jackknife---: \n")
            flush.console()
            param <- matrix(0, npr, ni)
            parm <- c(c(out$beta), c(out$gamma), out$prob)
            iseq <- 1:ni
            imask <- rep(TRUE, ni)
            nn <- (ni - 1) * no
            cvllike <- rep(0, ni)
            cvDat <- matrix(0, ni, no)
            Info <- matrix(0, (npp + npl) * ng + ng - 1, (npp + npl) * 
                ng + ng - 1)
            errcnt <- 0
            for (i in 1:ni) {
                # cat("                   : i =", i, "\n")
                flush.console()
                imask[i] <- FALSE
                indi <- iseq[imask]
                tDat <- Dat[indi, ]
                tRisk <- Risk[indi, ]
                ggt <- out$gwt[indi, ]
                beta <- out$beta
                gamma <- out$gamma
                if (ng > 1) {
                    prob <- out$prob
                }
                else {
                    prob <- 1
                }
                cFrtr <- .Fortran("r_dmzip", as.double(X), as.double(Z), 
                    as.double(tDat), as.double(tRisk), ggt = as.double(ggt), 
                    prob = as.double(prob), beta = as.double(beta), 
                    gamma = as.double(gamma), llike = as.double(0), 
                    Info = as.double(Info), as.integer(nn), as.integer(ni - 
                      1), as.integer(no), as.integer(npp), as.integer(npl), 
                    as.integer(ng), err = as.integer(0))
                if (cFrtr$err != 0) {
                    param[, i] <- c(out$beta, out$gamma, out$prob)
                    cvDat[i, ] <- NA
                    errcnt <- errcnt + 1
                    imask[i] <- TRUE
                    next
                }
                beta <- cFrtr$beta
                gamma <- cFrtr$gamma
                prob <- cFrtr$prob
                param[, i] <- c(beta, gamma, prob)
                beta <- matrix(beta, npp, ng)
                gamma <- matrix(gamma, npl, ng)
                Xb <- X %*% beta
                Zg <- Z %*% gamma
                Xb <- log(Risk[i, ]) + Xb
                lam <- exp(Xb)
                qq <- exp(Zg)
                qq <- qq/(1 + qq)
                mu <- (1 - qq) * lam
                llike = log(exp(Zg) + exp(-exp(Xb)))
                DD <- Dat[i, ]
                DD[DD < 0] <- 0
                llike[DD > 0.5, ] <- 0
                tllike = c(DD) * Xb - exp(Xb) - lgamma(DD + 1)
                tllike[DD < 0.5, ] <- 0
                llike = llike + tllike
                llike = llike - log(1 + exp(Zg))
                llike[Dat[i, ] < 0, ] <- 0
                ggt <- apply(llike, 2, sum)
                ggt <- log(prob) + ggt
                rscl <- max(ggt)
                ggt <- exp(ggt - rscl)
                cvllike[i] <- sum(ggt)
                ggt <- ggt/cvllike[i]
                cvllike[i] <- log(cvllike[i]) + rscl
                for (j in 1:ng) {
                    mu[, j] <- ggt[j] * mu[, j]
                }
                cvDat[i, ] <- apply(mu, 1, sum)
                imask[i] <- TRUE
            }
            cat("|---CV/Jackknife---: [DONE] \n")
            flush.console()
            cv <- mean(abs(Dat - cvDat), na.rm = TRUE)
            cv2 <- mean((Dat - cvDat)^2, na.rm = TRUE)
            dparm <- ni * parm - (ni - 1) * param
            jVar <- var(t(dparm))/ni
            out$cv <- cv
            out$cv2 <- cv2
            out$lcv <- -2 * sum(cvllike)
            out$param <- param
            out$jVar <- jVar
            out$errcnt <- errcnt
        }
        class(out) <- c("dmZIP", class(out))
        out
      }
    
    summary.dmZIP = function (object, ...) 
      {
        # cat("Summary \n")
        # cat("======= \n \n")
        # cat("AIC: ", object$AIC, "\n")
        # cat("BIC: ", object$BIC, "\n")
        if (!is.null(object$cv)) {
            #cat("CVE: ", object$cv, "\n")
        }
        gwt <- object$gwt
        gwt <- round(gwt, 3)
        gwt <- gwt/apply(gwt, 1, sum)
        colnames(gwt) <- paste("Grp", 1:object$ng, sep = "")
        rownames(gwt) <- paste("Indiv", 1:object$ni, sep = "")
        #cat("\n")
        #cat("Group Membership Prob. \n")
        #cat("====================== \n \n")
        return(gwt)
    }
    
    summary.dmZIPt = function (object, ...) 
      {
        # cat("Summary \n")
        # cat("======= \n \n")
        # cat("AIC: ", object$AIC, "\n")
        # cat("BIC: ", object$BIC, "\n")
        if (!is.null(object$cv)) {
           # cat("CVE: ", object$cv, "\n")
        }
        gwt <- object$gwt
        gwt <- round(gwt, 3)
        gwt <- gwt/apply(gwt, 1, sum)
        colnames(gwt) <- paste("Grp", 1:object$ng, sep = "")
        rownames(gwt) <- paste("Indiv", 1:object$ni, sep = "")
        # cat("\n")
        # cat("Group Membership Prob. \n")
        # cat("====================== \n \n")
        return(gwt)
    }
    
      assignInNamespace("dmZIP",dmZIP ,ns="crimCV")
      assignInNamespace("dmZIPt",dmZIPt ,ns="crimCV")

      
      
# IMPROVE treaj function
        # solves an issue with highly correlated measurements which led to an error
    
    step1measures = function (Data, Time, ID = FALSE, verbose = TRUE) 
    {
      data = Data
      time = Time
      input.data = data
      input.time = time
      if (dim(data)[1] != dim(time)[1] || dim(data)[2] != dim(time)[2]) 
        stop("data and time must be the same size.")
      sample.size = dim(data)[1]
      if (ID) {
        IDvector = data[, 1]
        data = data[, -1]
        time = time[, -1]
      }
      max.num.obs = dim(data)[2]
      clean.data = matrix(ncol = max.num.obs, nrow = sample.size)
      clean.time = matrix(ncol = max.num.obs, nrow = sample.size)
      num.obs = rep(999, sample.size)
      less.than.4.obs = NULL
      for (i_sample in 1:sample.size) {
        real.obs.pos = which(!is.na(data[i_sample, ]))
        num.obs[i_sample] = length(real.obs.pos)
        clean.data[i_sample, ] = as.vector(c(unlist(data[i_sample, 
                                                         real.obs.pos]), rep(NA, max.num.obs - num.obs[i_sample])))
        clean.time[i_sample, ] = as.vector(c(unlist(time[i_sample, 
                                                         real.obs.pos]), rep(NA, max.num.obs - num.obs[i_sample])))
        if (length(real.obs.pos) < 4) 
          less.than.4.obs = c(less.than.4.obs, i_sample)
        clean.data.pos = which(!is.na(clean.data[i_sample, ]))
        if (any(is.na(clean.time[i_sample, clean.data.pos]))) 
          stop(paste("There must be a time associated to every observation. Line: ", 
                     i_sample, sep = ""))
      }
      if (!is.null(less.than.4.obs)) {
        clean.data = clean.data[-less.than.4.obs, ]
        clean.time = clean.time[-less.than.4.obs, ]
      }
      sample.size = nrow(clean.data)
      if (ID) {
        if (!is.null(less.than.4.obs)) 
          IDvector = IDvector[-less.than.4.obs]
      } else {IDvector = seq(1:sample.size)}
      data = clean.data
      time = clean.time
      output = data.frame(matrix(ncol = 25, nrow = sample.size))
      names = c("ID", "m1", "m2", "m3", "m4", "m5", "m6", "m7", 
                "m8", "m9", "m10", "m11", "m12", "m13", "m14", "m15", 
                "m16", "m17", "m18", "m19", "m20", "m21", "m22", "m23", 
                "m24")
      colnames(output) = names
      output$ID = IDvector
      for (i in 1:sample.size) {
        output$m1[i] = max(data[i, ], na.rm = TRUE) - min(data[i, 
                                                               ], na.rm = TRUE)
      }
      for (i in 1:sample.size) {
        output$m2[i] = mean(data[i, ], na.rm = TRUE)
      }
      for (i in 1:sample.size) {
        output$m3[i] = sqrt(var(data[i, ], na.rm = TRUE))
      }
      for (i in 1:sample.size) {
        output$m4[i] = 100 * output$m3[i]/output$m2[i]
      }
      for (i in 1:sample.size) {
        output$m5[i] = pastecs::last(data[i, ], na.rm = TRUE) - pastecs::first(data[i, 
                                                                                    ], na.rm = TRUE)
      }
      for (i in 1:sample.size) {
        output$m6[i] = (pastecs::last(data[i, ], na.rm = TRUE) - pastecs::first(data[i, 
                                                                                     ], na.rm = TRUE))/(pastecs::last(time[i, ], na.rm = TRUE) - 
                                                                                                          pastecs::first(time[i, ], na.rm = TRUE) + 1)
      }
      for (i in 1:sample.size) {
        output$m7[i] = (pastecs::last(data[i, ], na.rm = TRUE) - pastecs::first(data[i, 
                                                                                     ], na.rm = TRUE))/pastecs::first(data[i, ], na.rm = TRUE)
      }
      for (i in 1:sample.size) {
        output$m8[i] = (pastecs::last(data[i, ], na.rm = TRUE) - pastecs::first(data[i, 
                                                                                     ], na.rm = TRUE))/output$m2[i]
      }
      for (i in 1:sample.size) {
        b = coefficients(lm(data[i, ] ~ time[i, ]))
        output$m9[i] = b[2]
      }
      for (i in 1:sample.size) {
        model = lm(data[i, ] ~ time[i, ])
        r = resid(model)
        RSS = r %*% r
        Y = subset(data[i, ], is.na(data[i, ]) == FALSE)
        m = length(Y)
        SYY = Y %*% Y - (sum(Y)^2)/m
        SSREG = SYY[1] - RSS[1]
        output$m10[i] = SSREG/SYY
      }
      FD = matrix(nrow = sample.size, ncol = max.num.obs - 1)
      for (i in 1:sample.size) {
        for (j in 1:(max.num.obs - 1)) {
          FD[i, j] = data[i, (j + 1)] - data[i, j]
        }
      }
      for (i in 1:sample.size) {
        output$m11[i] = max(FD[i, ], na.rm = TRUE)
      }
      for (i in 1:sample.size) {
        output$m12[i] = sqrt(var(FD[i, ], na.rm = TRUE))
      }
      FDunit = matrix(nrow = sample.size, ncol = max.num.obs - 
                        1)
      for (i in 1:sample.size) {
        for (j in 1:(max.num.obs - 1)) {
          FDunit[i, j] = FD[i, j]/(time[i, j + 1] - time[i, 
                                                         j])
        }
      }
      for (i in 1:sample.size) {
        output$m13[i] = sqrt(var(FDunit[i, ], na.rm = TRUE))
      }
      for (i in 1:sample.size) {
        output$m14[i] = mean(abs(FD[i, ]), na.rm = TRUE)
      }
      for (i in 1:sample.size) {
        output$m15[i] = max(abs(FD[i, ]), na.rm = TRUE)
      }
      for (i in 1:sample.size) {
        output$m16[i] = output$m15[i]/output$m2[i]
      }
      for (i in 1:sample.size) {
        output$m17[i] = output$m15[i]/output$m9[i]
      }
      for (i in 1:sample.size) {
        output$m18[i] = output$m12[i]/output$m9[i]
      }
      SD = matrix(nrow = sample.size, ncol = max.num.obs - 2)
      for (i in 1:sample.size) {
        for (j in 1:(max.num.obs - 2)) {
          SD[i, j] = FD[i, (j + 1)] - FD[i, j]
        }
      }
      for (i in 1:sample.size) {
        output$m19[i] = mean((SD[i, ]), na.rm = TRUE)
      }
      for (i in 1:sample.size) {
        output$m20[i] = mean(abs(SD[i, ]), na.rm = TRUE)
      }
      for (i in 1:sample.size) {
        output$m21[i] = max(abs(SD[i, ]), na.rm = TRUE)
      }
      for (i in 1:sample.size) {
        output$m22[i] = output$m21[i]/output$m2[i]
      }
      for (i in 1:sample.size) {
        output$m23[i] = output$m21[i]/output$m14[i]
      }
      for (i in 1:sample.size) {
        output$m24[i] = mean(abs(SD[i, ]), na.rm = TRUE)/output$m14[i]
      }
      temp.data = output$m2[(output$m2 != 0)]
      abs.temp.data = abs(temp.data)
      if (length(temp.data) == 0) {
        mean.0 = 0.0001
      } else {mean.0 = temp.data[which.min(abs.temp.data)]/100}
      m4.na.pos = which(is.na(output$m4) | is.infinite(output$m4))
      if (length(m4.na.pos) != 0) 
        output$m4[m4.na.pos] = 100 * output$m3[m4.na.pos]/mean.0
      m8.na.pos = which(is.na(output$m8) | is.infinite(output$m8))
      if (length(m8.na.pos) != 0) {
        if (length(m8.na.pos) > 1) 
          output$m8[m8.na.pos] = (apply(data[m8.na.pos, ], 
                                        1, pastecs::last, na.rm = TRUE) - apply(data[m8.na.pos, 
                                                                                     ], 1, pastecs::first, na.rm = TRUE))/mean.0
        else output$m8[m8.na.pos] = (pastecs::last(data[m8.na.pos, ], 
                                                   na.rm = TRUE) - pastecs::first(data[m8.na.pos, ], na.rm = TRUE))/mean.0
      }
      m16.na.pos = which(is.na(output$m16) | is.infinite(output$m16))
      if (length(m16.na.pos) != 0) 
        output$m16[m16.na.pos] = output$m15[m16.na.pos]/mean.0
      m22.na.pos = which(is.na(output$m22) | is.infinite(output$m22))
      if (length(m22.na.pos) != 0) 
        output$m22[m22.na.pos] = output$m21[m22.na.pos]/mean.0
      temp.data = data[(data[, 1] != 0), ]
      abs.temp.data = abs(temp.data)
      if (nrow(temp.data) == 0) {
        y1.0 = 0.0001} else {y1.0 = temp.data[which.min(abs.temp.data)]/100 +0.0001}
      m7.na.pos = which(is.na(output$m7) | is.infinite(output$m7))
      if (length(m7.na.pos) != 0) {
        if (length(m7.na.pos) > 1) {
          output$m7[m7.na.pos] = (apply(data[m7.na.pos, ], 
                                        1, pastecs::last, na.rm = TRUE) - apply(data[m7.na.pos, 
                                                                                     ], 1, pastecs::first, na.rm = TRUE))/y1.0} else { output$m7[m7.na.pos] = (pastecs::last(data[m7.na.pos, ], 
                                                                                                                                                                             na.rm = TRUE) - pastecs::first(data[m7.na.pos, ], na.rm = TRUE))/y1.0}
      }
      form10 = vector(length = sample.size)
      for (i_test in 1:sample.size) {
        model = lm(data[i_test, ] ~ time[i_test, ])
        r = resid(model)
        RSS = r %*% r
        Y = subset(data[i_test, ], is.na(data[i_test, ]) == 
                     FALSE)
        m = length(Y)
        form10[i_test] = Y %*% Y - (sum(Y)^2)/m
        if (form10[i_test] == 0) 
          form10[i_test] = NA
      }
      if (!is.na(min(form10))) 
        syy.0 = min(form10)
      else syy.0 = 0.0001
      m10.na.pos = which(is.na(output$m10) | is.infinite(output$m10))
      if (length(m10.na.pos) != 0) {
        for (i_na in m10.na.pos) {
          model = lm(data[i_na, ] ~ time[i_na, ])
          r = resid(model)
          RSS = r %*% r
          Y = subset(data[i_na, ], is.na(data[i_na, ]) == 
                       FALSE)
          m = length(Y)
          SSREG = SYY[1] - RSS[1]
          output$m10[i_na] = SSREG/syy.0
        }
      }
      temp.data = output$m9[(output$m9 != 0)]
      abs.temp.data = abs(temp.data)
      if (length(temp.data) == 0) {slope.0 = 0.0001
      } else { slope.0 = temp.data[which.min(abs.temp.data)]/100}
      m17.na.pos = which(is.na(output$m17) | is.infinite(output$m17))
      if (length(m17.na.pos) != 0) 
        output$m17[m17.na.pos] = output$m15[m17.na.pos]/slope.0
      m18.na.pos = which(is.na(output$m18) | is.infinite(output$m18))
      if (length(m18.na.pos) != 0) 
        output$m18[m18.na.pos] = output$m12[m18.na.pos]/slope.0
      temp.data = output$m14[(output$m14 != 0)]
      abs.temp.data = abs(temp.data)
      if (length(temp.data) == 0) {
        mean.abs.0 = 0.0001} else {mean.abs.0 = temp.data[which.min(abs.temp.data)]/100}
      m23.na.pos = which(is.na(output$m23) | is.infinite(output$m23))
      if (length(m23.na.pos) != 0) 
        output$m23[m23.na.pos] = output$m21[m23.na.pos]/mean.abs.0
      m24.na.pos = which(is.na(output$m24) | is.infinite(output$m24))
      if (length(m24.na.pos) != 0) {
        abs.na.data = abs(SD[m24.na.pos, ])
        if (length(m24.na.pos) > 1) {
          output$m24[m24.na.pos] = apply(abs.na.data, 1, mean, 
                                         na.rm = TRUE)/mean.abs.0 } else {
                                           output$m24[m24.na.pos] = mean(abs.na.data, na.rm = TRUE)/mean.abs.0}
      }
      if(sum(is.infinite(output$m7))){
        output$m7[is.infinite(output$m7)] = max(output$m7[!is.infinite(output$m7)])
        
      }
      check.correlation( output[, -1], verbose)
      trajMeasures = structure(list(measurments = output, data = cbind(IDvector, 
                                                                       clean.data), time = cbind(IDvector, clean.time)), class = "trajMeasures")
      return(trajMeasures)
    }
    
    check.correlation=function (output, verbose = TRUE, is.return = FALSE) 
    {
      cor.mat = cor(output)
      mes.names = names(output)
      is.corr = FALSE
      corr.var = NULL
      for (i_row in mes.names[-length(mes.names)]) {
        i_pos = which(mes.names == i_row)
        res.names = mes.names[(i_pos + 1):length(mes.names)]
        for (i_col in res.names) {
          if (cor.mat[i_row, i_col] > 0.999) {
            corr.var = rbind(corr.var, c(i_col, i_row))
            if (verbose) {
              print(paste("Correlation of ", i_row, " and ", 
                          i_col, " : ", round(cor.mat[i_row, i_col], 
                                              3), sep = ""))
              is.corr = TRUE
            }
          }
        }
      }
      if (!is.corr && verbose) 
        print("No correlations found. That is good.")
      if (is.return) 
        return(corr.var)
    }
    
    plotMeanTraj = function (x, clust.num = NULL, ...) 
    {
      clust = x$clusters
      data = as.matrix(x$data)
      time = as.matrix(x$time)
      data = data[, -1]
      time = time[, -1]
      unique.clusters = sort(unique(clust[, 2]))
      num.clust = length(unique.clusters)
      if (is.null(clust.num)) {
        num.plot.x = round(sqrt(num.clust))
        if (num.plot.x^2 < num.clust) 
          num.plot.y = num.plot.x + 1
        else num.plot.y = num.plot.x
        par(mfrow = c(num.plot.y, num.plot.x), oma = c(0, 0, 
                                                       4, 0))
      }
      else unique.clusters = clust.num
      min.y = min(data, na.rm = T)
      max.y = max(data, na.rm = T)
      min.x = min(time, na.rm = T)
      max.x = max(time, na.rm = T)
      xlab = "time"
      ylab = "y"
      args = list(...)
      params = names(args)
      for (i_clust in unique.clusters) {
        clust.data.pos = which(clust$cluster == i_clust)
        clust.data = data[clust.data.pos, ]
        clust.time = time[clust.data.pos, ]
        unique.clust.time = as.numeric(names(table(clust.time)))
        comp.mean = data.frame(matrix(ncol = 3, nrow = 1))
        for (i_time in unique.clust.time) {
          time.pos = which(clust.time == i_time)
          data.time = clust.data[time.pos]
          mean.data = mean(data.time)
          centiles = 0
          comp.mean = rbind(comp.mean, c(mean.data, centiles[1], 
                                         centiles[2]))
        }
        comp.mean = comp.mean[-1, ]
        args[["x"]] = unique.clust.time
        args[["y"]] = comp.mean[, 1]
        if (!("xlim" %in% params)) 
          args[["xlim"]] = c(min.x, max.x)
        if (!("ylim" %in% params)) 
          args[["ylim"]] = c(min.y, max.y)
        if (!("main" %in% params)) 
          args[["main"]] = paste("Cluster ", i_clust, sep = "")
        if (!("xlab" %in% params)) 
          args[["xlab"]] = xlab
        if (!("ylab" %in% params)) 
          args[["ylab"]] = ylab
        do.call(plot, args)
        lines(unique.clust.time, comp.mean[, 1])
      }
      if (is.null(clust.num)) {
        text = "Mean for Every Cluster"
        mtext(text, side = 3, line = 0.5, outer = TRUE, cex = 1.2)
      }
      if (is.null(clust.num)) {
        cat("nope")
      }
    }

# Load sample data
sample_traj_data = read.csv(text='"ID","t 1","t 2","t 3","t 4","t 5","t 6"
                              1,707.773301494225,707.773301494225,10651.4798046765,247.445586779939,247.445586779939,224.305830021114
                              2,3336.66826966042,4892.62165068284,8621.54142701205,5698.5485159355,5846.25381703163,6487.2371862624
                              3,-1e-05,-1e-05,-1e-05,621.126521648944,3312.48984962403,9529.39539355984
                              4,1957.65655920987,763.942609669943,949.09003555274,1656.94408804365,1561.56394381288,375.192068812877
                              5,-1e-05,-1e-05,-1e-05,-1e-05,-1e-05,1822.26722128284
                              6,257.528524447414,257.528524447414,257.528524447414,257.528524447414,257.528524447414,257.528524447414
                              7,6559.25841371479,6559.25841371479,5818.86290494286,982.907818868382,979.599993221458,4003.5407701786
                              8,889.463002397763,889.463002397763,889.463002397763,4405.18455285578,1139.25695385647,375.132753856472
                              9,-1e-05,1965.94229578523,1978.9058232261,1978.9058232261,8417.46348971495,3674.1052632261
                              10,603.052962766365,3512.32309179862,2244.11157033907,7280.41848189044,14685.5075827308,5261.27561634422
                              11,6478.09917599566,4106.45517599566,1247.25117599566,272.131175995657,272.131175995657,272.131175995657
                              12,1521.45760289976,1521.45760289976,4108.70560289976,1521.45760289976,929.28238978501,7882.82672702155
                              13,-1e-05,-1e-05,77.1784251915575,1157.67637787336,6749.55706023313,3232.76248165136
                              14,151.206115813525,151.206115813525,151.206115813525,151.206115813525,151.206115813525,151.206115813525
                              15,823.006860255652,823.006860255652,823.006860255652,823.006860255652,5065.43596025565,11153.6493602557
                              16,3139.1894285543,3147.36129931965,2604.85480488479,2604.85480488479,2857.66713818256,4695.12221620555
                              17,1347.39171202082,3934.63971202082,2161.37925126946,907.138448532558,485.784559702941,336.497999702941
                              18,533.707139488881,533.707139488881,533.707139488881,1399.3304715693,9825.71306895786,666.9876944124
                              19,4084.45766220507,3414.85975852015,1236.26712699999,1883.07912699999,1163.25352651585,10662.5418941133
                              20,-1e-05,-1e-05,-1e-05,364.672706714269,9432.32971138833,370.00225138833
                              21,2445.84158629389,582.649663394786,321.797291723357,122.3219486805,122.3219486805,122.3219486805
                              22,197.30458998288,197.30458998288,197.30458998288,197.30458998288,197.30458998288,2725.42792296058
                              23,180.448521996917,433.260855294687,877.351950185872,6263.94899354773,8480.68241046741,6570.10941622248
                              24,70.9739612676056,70.9739612676056,70.9739612676056,70.9739612676056,70.9739612676056,70.9739612676056
                              25,80.6624696853688,80.6624696853688,80.6624696853688,80.6624696853688,80.6624696853688,80.6624696853688
                              26,4406.45923014359,4561.55415010124,5855.17815010124,4561.55415010124,4644.84712222822,3749.27986824147
                              27,2708.7465277154,1103.3977277154,1103.3977277154,1103.3977277154,983.46765778533,3433.65272072239
                              28,819.08262097805,819.08262097805,760.381444507462,598.953209213345,1249.30461680421,6683.57512528997
                              29,-1e-05,-1e-05,4407.85879936575,2194.03837125194,929.976704763092,6751.28470476309
                              30,4129.45347094477,3623.82880434923,3623.82880434923,11326.4481373269,15800.2457662096,12012.7462637177')
