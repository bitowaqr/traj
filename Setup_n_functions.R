# setup
# load and install libraries
  required_packages<-c("crimCV","ggplot2","cowplot","reshape2","dplyr","traj","cowplot") 

  pft_packages <- function(package){
      for(i in 1:length(package)){
        if(eval(parse(text=paste("require(",package[i],")")))==0) {
          install.packages(package)}}
      return (eval(parse(text=paste("require(",package,")"))))}
  
  pft_packages(required_packages)

# load some functions
  
  
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
