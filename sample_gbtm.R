
# GBTM TRAJ ANALYSIS USING the crimCV package
# Adjusted by Paul Schneider 2018

#----------------------------------------#
#  1. Load functions & packages          #
#----------------------------------------#

  url = "https://raw.githubusercontent.com/bitowaqr/traj/master/Setup_n_functions.R"
  source(url)

  
#----------------------------------------#
#  2. Load your data                     #
#----------------------------------------#

  traj_data = load.traj.data(ID.index = NULL, time.points = NULL, path = NULL)
  
  
#----------------------------------------#
#  3. Fit GBTM models                    #
#----------------------------------------#
  
  # You can fit your own models (which might take a while!)
  fitted.gbtm = fit.gbtm()
  
  # Or load some models I have ran and saved for your convenience:
  load(url("https://github.com/bitowaqr/traj/blob/master/sample_gbtm_results.rdata?raw=true"))
  # this loads an object calles 'fitted.gbtm' which correpsonds to:
  # fitted.gbtm = fit.gbtm(data = traj_data,
  #                        evaluate.k.groups = c(2:8),
  #                        evaluate.p.polynomials = c(1:5),
  #                        run.loocv = T,
  #                        method = "ZIP",
  #                        init.proc = 20,
  #                        Risk.set = NULL)
  
  
#----------------------------------------#
#  4. Compare fit across models          #
#----------------------------------------#
  
    # Look at individual values of individual models
    # models for 2 groups with 3 degree polynomials are stored under:
    fitted.gbtm$cv.eval.list$k2p3$AIC # AIC
    fitted.gbtm$cv.eval.list$k2p3$BIC # BIC
    fitted.gbtm$cv.eval.list$k2p3$cv # LOOCV
    
    # Look at all AIC values of all models (showing the first 6)
    head(fitted.gbtm$eval.stats$eval.data.list$AIC)
    
    # Plot AIC values of all models
    fitted.gbtm$eval.stats$eval.plot.list$AIC # or BIC or CV
    
    # Look at plots for AIC, BIC and CV in comparison 
      plot.legend = get_legend(fitted.gbtm$eval.stats$eval.plot.list[[1]] + theme(legend.position = "bottom"))
      eval.plots.legendless = lapply(fitted.gbtm$eval.stats$eval.plot.list, function(x) x + theme(legend.position = "none") + scale_x_continuous(labels=1:5))
      # IGNORE ONE OUTLIER
      eval.plots.legendless$CV = eval.plots.legendless$CV + coord_cartesian(ylim=c(4,8))
      eval.plot.overview = plot_grid(plot_grid(plotlist =  eval.plots.legendless, nrow=1),
                                     plot.legend,nrow=2,rel_heights = c(10,1))    
    eval.plot.overview  # note: missing values in CV plot maybe due to non-convergence?

    
#----------------------------------------#
#  5. Select and evaluate 1 GBTM         #
#----------------------------------------#
  
  # select your model (eg. 6 groups and 2 poly)
    select.model = set.model(models.list = fitted.gbtm$cv.eval.list, set.k = NULL, set.p = NULL)
    
    # Retrieve membership info from model
    member_list = gbtm.members(model = select.model,
                               data = traj_data)
    
    # show beginining of probabilistic membership table
    head(member_list$prob.membership)
    # show beginining of determinsitic membership table
    head(member_list$prob.membership)
    # Group membership frequency table
    member_list$deter.membership.table
    

  # retrieve model function terms with intercept and * for p < .05
    model.spec = get.model.terms(model = select.model, 
                                 data = traj_data)
    
    # print model
    model.spec$model
    
    # expected values per group (showing the first 19 rows)
    model.spec$predicted.group.values[2:19,]
    
    
#----------------------------------------#
#  6. Create fancy GBTM-traj plots       #
#----------------------------------------#
    
  ## Create complex plot: overview and individual group average traj + individual traj
    complex.plot = plot.gbtm.groups(data = traj_data,
                                    model = select.model,
                                    set.y.limit = NULL,  # c(0,40)
                                    individual.group.names = NULL, # give names
                                    plot.total = T,
                                    plot_overview = T,
                                    n_row = NULL,
                                    y.axis.label = "y axis label",
                                    x.axis.label = "x axis label",
                                    plot.title = "")
    complex.plot
    
    # create individual plots per group only
    individual.plots = plot.gbtm.groups(data = traj_data,
                                        model = select.model,
                                        set.y.limit = NULL,  # c(0,40)
                                        individual.group.names = NULL, # give names
                                        plot.total = T,
                                        plot_overview = F,
                                        n_row = 2,
                                        y.axis.label = "y axis label",
                                        x.axis.label = "x axis label",
                                        plot.title = "")

    individual.plots
    
    
    ## end
    #  ; )
    
    
    
