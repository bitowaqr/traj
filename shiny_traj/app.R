# a shiny app test
# to run 'traj' cluster analysis
# modified by Paul Schneider, 2018

 # # # to execute this app remotely, run the follwoing three lines of code:
    # if(eval(require(RCurl))==0){install.packages("RCurl")}
    # url = "https://raw.githubusercontent.com/bitowaqr/traj/master/shiny_traj/app.R"
    # eval(parse(text = RCurl::getURL(url)), envir= .GlobalEnv)


library(shiny)
library(DT)
url = "https://raw.githubusercontent.com/bitowaqr/traj/master/Setup_n_functions.R"
source(url)
load.traj_24.data_shiny = function(ID.index = NULL,time.points = NULL,input_data = raw_data,cluster_col=NULL){
  
  dat = input_data
  # Set ID column
  cat("\n ")
  
  
  if(ID.index %in% c("yes","y","Yes","Y","N","n","No","no")){
    if(ID.index == "yes" |ID.index == "Yes" | ID.index == "y" | ID.index == "Y"){
      ID = dat[,1]
      dat =dat[,-1]
    }
    
    if(ID.index == "no" | ID.index == "No" | ID.index == "n" | ID.index == "N"){
      cat("\n ID variable created")
      ID = 1:length(dat[,1])
    }
  } else{stop("\n Invalid input <-  Type y or n")}
  
  # Set time points
  
  time.points = as.numeric(time.points)
  if(time.points!=round(time.points)) stop(" Time points must be an integer")
  if(time.points > length(dat)) stop(" Data does not have enough columns")
  mat = dat[,1:time.points]
  
  # make data set
  time_mat = matrix(ncol=time.points,data=rep(1:time.points,each=dim(mat)[1]))
  traj_data  = cbind(ID,mat,time_mat)
  names(traj_data) = c("ID",
                       paste("x",1:time.points,sep=""),
                       paste("t",1:time.points,sep=""))
  
  if(cluster_col != "NULL"){
    traj_data = data.frame(traj_data[,1],input_data[,as.numeric(cluster_col)],traj_data[,-1])
    names(traj_data) = c("ID","predefined_cluster",
                         paste("x",1:time.points,sep=""),
                         paste("t",1:time.points,sep=""))
  }
  
  cat("\n Data set with dimensions: cases=",dim(traj_data)[1],"; timepoints=",dim(mat)[2]," created!",sep="")
  return(traj_data)
}
traj.k.mean_shiny = function( processed_data = traj_data,
                        nclusters_set = 1, # set number of clusters, or NULL to let R decide
                        select_criteria = "ccc",
                        set_subgroup = input$set_subgroup_input,
                        init_value = 50,
                        set_y_label ="",
                        set_x_label = "",
                        sample_per_plot = 15,
                        n_nrow = NULL){
  
  s1 = s2all = s3all = k.membership = k.membership.table= sample.combined.individual.plot = sample.individual.plot = all.individual.plot = NULL
  
  
  
  
  if(!is.null(nclusters_set)){
  if(nclusters_set!="NULL"){
    nclusters_set = as.numeric(nclusters_set)
  }else{
    nclusters_set = NULL
  }}
  
  if(!is.null(n_nrow)){
    if(n_nrow!="NULL"){
      n_nrow = as.numeric(n_nrow)
    }else{
      n_nrow = NULL
    }}
  
  if(set_subgroup != "NULL"){
    set_subgroup <- as.integer(strsplit(set_subgroup, ",")[[1]])
    processed_data = processed_data[processed_data$predefined_cluster %in% set_subgroup,]
  } 
  
  
  if(sum("predefined_cluster" %in% names(processed_data))>0){
    remove.predefined.cluster = which(names(processed_data) == "predefined_cluster")
    processed_data = processed_data[,-remove.predefined.cluster]
  }
    
  id.index = 1
  min.data.index = 2
  max.data.index = length(processed_data[,-1])/2+1
  
  min.time.index = length(processed_data[,-1])/2+2
  max.time.index = length(processed_data[,]) 
  
  ID = processed_data[,id.index]
  data.raw = processed_data[,min.data.index:max.data.index]
  time.raw = processed_data[,min.time.index:max.time.index]
  
  minimum.check = apply(data.raw,1,function(x)(sum(!is.na(x))))
  minimum.check = sum(minimum.check<5)
  
  if(minimum.check >0){
    cat("\n Clustering needs at least 5 observation per case:",minimum.check,"cases are being removed...! \n")
  }
  cat("\n Start Step 1: \n ")
  tryCatch({
    s1 = step1measures(Data = processed_data[,c(1,min.data.index:max.data.index)], 
                       Time = processed_data[,c(1,min.time.index:max.time.index)], 
                       ID = T ) 
  }, error = function(e) cat("\n Oops...Something went wrong in step 1"))
  cat("\n Start Step 2: ")
  tryCatch({
    s2all = step2factors(s1) 
  }, error = function(e) cat("\n Oops...Something went wrong in step 2. Maybe one or more columns/time points contains less than 5 valid observations?!"))
  cat("\n Start Step 3: ")
  tryCatch({
    print(nclusters_set)
    s3all = step3clusters(trajFactors = s2all, 
                          nclusters = nclusters_set,
                          criteria = select_criteria,
                          nstart = init_value)  # 5 clusters
  }, error = function(e) cat("\n Oops...Something went wrong in step 3"))
  # clust.build.plot(s3all,y.max.lim =y.max.lim)
  k.membership.table = data.frame(table(s3all$clusters$cluster))
  k.membership = s3all$clusters
  k.k = unique(k.membership$cluster)
  
  cat("\n \n Creating Plots")
  tryCatch({
    # plot all individuals Individual plot
    indiv.data = merge(processed_data[,c(1,min.data.index:max.data.index)],k.membership,"ID")
    indiv.data = melt(indiv.data,id.vars=c("ID","cluster"))
    indiv.data$variable = as.character(indiv.data$variable)
    indiv.data$variable = gsub("x","",indiv.data$variable)
    indiv.data$variable = as.numeric(indiv.data$variable)
    indiv.data$cluster = as.factor(indiv.data$cluster)
    lev.clust = levels(indiv.data$cluster)
    n.clust = as.numeric(by(indiv.data$ID,indiv.data$cluster,function(x)length(unique(x))))
    names.clust = data.frame(cluster = lev.clust,n = paste("Cluster ",lev.clust," n=",n.clust," (",round(n.clust/sum(n.clust),2)*100,"%)",sep=""))
    indiv.data = merge(indiv.data,names.clust,by="cluster")
    indiv.data$ID = as.factor(indiv.data$ID)
    
    if(is.null(n_nrow)) n_nrow = length(n.clust)
    
    all.individual.plot = 
      ggplot(indiv.data,aes(x=variable,y=jitter(value),col= cluster)) + 
      facet_wrap(~cluster,nrow = n_nrow) + 
      geom_line(aes(group=ID),alpha=0.3,size=0.5) +
      stat_summary(aes(group = cluster), fun.y = mean, geom = 'line',linetype="dashed", size=1, alpha=1,col="black") +
      scale_color_manual(values= k.k+1,labels=names.clust$n,name="Clusters") +
      theme(legend.position = "top") +
      ylab(set_y_label) +
      xlab(set_x_label)
    
    
    pre_all_overview = 
      ggplot(indiv.data,aes(x=variable,y=jitter(value))) + 
      geom_line(aes(col=ID),alpha=0.3,size=0.5) +
      stat_summary(aes(x=variable), fun.y = mean, geom = 'line',linetype="dashed", size=1, alpha=1,col="black") +
      theme(legend.position = "none") +
      ylab(set_y_label) +
      xlab(set_x_label)
    
    cluster.plot.df = data.frame(ID=NA,cluster=NA,value=NA,time=NA)
    for(i in k.k){
      temp.sample = processed_data[,c(1,min.data.index:max.data.index)][ID %in% k.membership$ID[k.membership$cluster==i],]
      if(length(temp.sample[,1])<sample_per_plot) sample_per_plot = length(temp.sample[,1])
      draw = sample(1:length(temp.sample[,1]),sample_per_plot)
      sample.data = temp.sample[draw,]
      sample.data = melt(sample.data,id.vars="ID")
      sample.data = data.frame(ID = sample.data$ID,
                               cluster = i,
                               value = sample.data$value,
                               time = sample.data$variable
      )
      cluster.plot.df = rbind(cluster.plot.df,sample.data)
    }
    cluster.plot.df = cluster.plot.df[-1,]
    cluster.plot.df$time = gsub("x","",cluster.plot.df$time)
    cluster.plot.df$time = as.numeric(as.character(cluster.plot.df$time))
    cluster.plot.df$cluster = as.factor(cluster.plot.df$cluster)
    cluster.plot.df$ID = as.factor(cluster.plot.df$ID)
    
    sample.combined.individual.plot = 
      ggplot() + 
      geom_line(data= cluster.plot.df, aes(x=time,y=jitter(value,factor = 2),col=cluster,group =ID)) +
      scale_color_manual(values= k.k+1,labels=names.clust$n,name="Clusters") +
      theme(legend.position = "top") +
      ylab(set_y_label) +
      xlab(set_x_label)
    
    sample.individual.plot = 
      ggplot(cluster.plot.df) + 
      geom_line(aes(x=jitter(time),y=jitter(value),col =ID),alpha=0.5) +
      geom_point(aes(x=jitter(time),y=jitter(value),col =ID),alpha=0.5) +
      facet_wrap(~ cluster,nrow = n_nrow) +
      stat_summary(aes(x=(time),y=(value)),geom="line",fun.y = "mean",size=2) +
      ylab(set_y_label) +
      xlab(set_x_label) +
      theme(legend.position = "none")
    
  }, error = function(e) cat("\n Oops...Something went wrong while plotting"))
  
  cat("\n Makig individual plots and statistics")
  tryCatch({
    individual.group.plots = list()
    for(l in 1:max(k.k)){
      
      individual.group.plots[[l]] = list()
      names(individual.group.plots)[l] = paste("group_",l,sep="")
      
      temp = indiv.data[indiv.data$cluster==l,]
      temp$ID = droplevels(temp$ID)
      
      individual.group.plots[[l]]$plot = 
        ggplot(temp) + 
        geom_line(aes(x=jitter(variable),y=jitter(value),col =ID),alpha=0.5) +
        geom_point(aes(x=jitter(variable),y=jitter(value),col =ID),alpha=0.5) +
        stat_summary(aes(y = value,x=variable), fun.y = "mean", geom = "line",size=2,col="black") +
        ylab(set_y_label) +
        xlab(set_x_label) +
        #scale_color_manual(values=c(rainbow(length(ID))) )+
        theme(legend.position = "none")
      
      summary.stat = data.frame(n = length(unique(temp$ID)),
                                obs = sum(!is.na(temp$value)),
                                mean= mean(temp$value,na.rm=T),
                                sd=sd(temp$value,na.rm=T),
                                median = median(temp$value,na.rm=T),
                                IQR = IQR(temp$value,na.rm=T))
      by.summary.stat = by(data.frame(value=temp$value,
                                      ID = temp$ID),temp$variable,function(x) data.frame(n = length(unique(x$ID[!is.na(x$value)])),
                                                                                         mean= mean(x$value,na.rm=T),
                                                                                         sd=sd(x$value,na.rm=T),
                                                                                         median = median(x$value,na.rm=T),
                                                                                         IQR = IQR(x$value,na.rm=T)))
      
      by.summary.stat = matrix(unlist(by.summary.stat),nrow=5)
      rownames(by.summary.stat) = c("n","mean","sd","median","IQR")
      colnames(by.summary.stat) = paste("Timepoint_",1:max(temp$variable,na.rm=T),sep="")
      by.summary.stat = t(by.summary.stat)
      individual.group.plots[[l]]$summary.stat = list(Overall=summary.stat,
                                                      per.timepoint = by.summary.stat)
    }
    
  }, error = function(e) cat("\n Oops...Something went wrong while plotting"))
  
  cluster.analysis = list(
    s1 = s1,
    s2all = s2all,
    s3all = s3all,
    membership = k.membership,
    membership.table = k.membership.table,
    sample.combined.individual.plot = sample.combined.individual.plot,
    sample.individual.plot = sample.individual.plot,
    all.individual.plot = all.individual.plot,
    individual.group.plots = individual.group.plots,
    pre_all_overview = pre_all_overview
  )
  
  return(cluster.analysis)
}

server <- function(input, output, session){
  
  
  
  myData <- reactive({
    inFile <- input$file1
    print(inFile[,1])
    if (is.null(inFile)) return(NULL)
    
    if(grepl("\\.csv",inFile$name)){
    data <- read.csv(inFile$datapath, header = TRUE)
    } 
    if(grepl("\\.dta",inFile$name)){
      data <- foreign::read.dta(inFile$datapath)
    }
    # if((grepl("\\.dta",inFile$name) | grepl("\\.dta",inFile$name))==F){
    #   stop("Sorry, this only works with .csv or .dta files.")
    # }
    data  <- load.traj_24.data_shiny(ID.index = input$id.var,
                                     time.points = input$time.points,
                                     input_data = data,
                                     cluster_col = input$cluster_col_index)
    rownames(data) = NULL
    data
  })
  
  output$contents <- DT::renderDataTable({DT::datatable(myData(),rownames=F)})
  
  output$pre_clusters <- renderTable({as.data.frame(table(myData()[,"predefined_cluster"], dnn = list("Predefined cluster")), responseName = "Frequency")})
  
  values <- reactiveValues()
  values$list <- list()
  newEntry <- observe({
    observeEvent(input$do, {
      
      traj_analysis_results = isolate(traj.k.mean_shiny(myData()[,],
                                                        nclusters = input$nclusters_input, 
                                                        set_subgroup = input$set_subgroup_input,
                                                        n_nrow = input$n_row_input))
      isolate(values$list$plot_overview <- traj_analysis_results$all.individual.plot )
      isolate(values$list$plot_pre_overview <- traj_analysis_results$pre_all_overview + ggtitle("Pre-defined Pop. mean traj"))
      isolate(values$list$individual.group.plots <- traj_analysis_results$individual.group.plots)
      isolate(values$list$membership.table <- traj_analysis_results$membership.table)
      isolate(values$list$s_stat <- traj_analysis_results[1:3])
    })
    
    output$plot_overview <- renderPlot({values$list$plot_overview}
                                       #,height = 480, width = "auto"
                                       )
    output$plot_pre_overview <- renderPlot({values$list$plot_pre_overview}
                                           #,height = 480, width = "auto"
                                           )
    output$plot_individual_group <- renderPlot({values$list$individual.group.plots[[as.numeric(input$select_group_for_plot)]][[1]]})
    output$summary_individual_group <- renderPrint({print(values$list$individual.group.plots[[as.numeric(input$select_group_for_plot)]][[2]])})
    output$freq_table <- renderTable({values$list$membership.table})
    output$step_stat_printput <- renderPrint({print(values$list$s_stat[[as.numeric(as.character(input$select_step))]])})
  })
  
  
}

write.csv(data.frame(a = 1:10, b = letters[1:10]), 'test.csv')



ui<- shinyUI(fluidPage(
  
  fluidRow(
    navbarPage(" Traj24",
               
               tabPanel("Load your data",
                        sidebarPanel(
      textInput("id.var", "ID var in first column?",value="yes"),
      textInput("time.points", "How many time pints (i.e. columns?)",value="9"),
      textInput("cluster_col_index", "Cluster col index",value="37"),
      br(),
      fileInput('file1', 'Choose CSV File',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv'))
    ),
    mainPanel(
      DT::dataTableOutput('contents')

    )
  ),
  tabPanel("Run the Traj24 analysis",
           sidebarPanel(
             textInput("nclusters_input", "How many clusters to analyze (a number or NULL)",value="3"),
             textInput("set_subgroup_input", "Analyse subgroup (if predetermined clusters are available)",value="NULL"),
             textInput("n_row_input", "How many rows should the plot have?",value="NULL"),
             actionButton("do", "Run the traj24 analysis!"),
             br(),   br(),
             "Frequency table of predefined clusters:",
             tableOutput('pre_clusters')
           ),
           mainPanel("Pre-defined Pop. mean trajectory",
                     plotOutput("plot_pre_overview"),
                     "Identified latent group trajectories",
             plotOutput("plot_overview")
                     
                     )
  ),
  tabPanel("Plot individual grpuos",
           sidebarPanel(
             textInput("select_group_for_plot", "Which group to plot",value="3"),
             br(),   br(),
             "Traj24 cluster sizes:",
             tableOutput("freq_table")
           ),
           mainPanel(plotOutput("plot_individual_group"),
                     verbatimTextOutput("summary_individual_group"))
  ),
  tabPanel("traj24 stats",
           sidebarPanel(
             selectInput("select_step", "Select statistic",choices = c(1,2,3),selected = 3)
           ),
           mainPanel(verbatimTextOutput("step_stat_printput"))
  )
)
)))

shinyApp(ui,server)
