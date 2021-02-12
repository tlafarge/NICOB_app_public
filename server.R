
# server.R
library(metafor)
library(R2jags)
# library(rmutil)

shinyServer(function(input, output, session) {
  # outputOptions(output, "bayesRPest", suspendWhenHidden = FALSE)

  source("significator.R")
  col1="#FDE725FF"
  col2="#71CF57FF"
  col3="#440154FF"
  col4="#2E6E8EFF"
  col5="#9FDA3AFF"


  output$dyn_input1<-renderUI({
    x=as.numeric(unlist(strsplit(gsub("\\s", "", input$mean),",")))
    lab=trimws( strsplit(input$lablabels,",")[[1]])
    if(length(lab)!=0){
      sanitize = !startsWith( lab,"-")
      tempvalue =  mad(x[sanitize])
    }else
    {
      tempvalue=mad(x)
    }
    numericInput("halfCauchyScaleTau",
                 label = h4("Scale for half-Cauchy prior on between laboratory variance"),
                 # value = mad(labmeansbayes()))#median(abs(labmeansbayes()-median(labmeansbayes()))))
                 value=tempvalue)

  })

  output$dyn_input2<-renderUI({
    x=as.numeric(unlist(strsplit(gsub("\\s", "", input$mean),",")))
    se=as.numeric(unlist(strsplit(gsub("\\s", "", input$se),",")))
    lab=trimws( strsplit(input$lablabels,",")[[1]])
    if(length(lab)!=0){
      sanitize = !startsWith( lab,"-")
      tempvalue = median(se[sanitize])
    }else
    {
      tempvalue=median(se)
    }
    numericInput("halfCauchyScaleSigma",
                 label = h4("Scale for half-Cauchy prior on within laboratory variances"),
                 # value =  median(labsesbayes()))
                 #value =  median(as.numeric(unlist(strsplit(gsub("\\s", "", input$se),",")))))
                value=tempvalue)

  })
  output$dyn_input1_lap<-renderUI({
  x=as.numeric(unlist(strsplit(gsub("\\s", "", input$mean),",")))
  lab=trimws( strsplit(input$lablabels,",")[[1]])
  if(length(lab)!=0){
    sanitize = !startsWith( lab,"-")
    tempvalue =  mad(x[sanitize])
  }else
  {
    tempvalue=mad(x)
  }
  numericInput("halfCauchyScaleTau_lap",
               label = h4("Scale for half-Cauchy prior on between laboratory variance"),
               # value = mad(labmeansbayes()))#median(abs(labmeansbayes()-median(labmeansbayes()))))
               value=tempvalue)

  })

  output$dyn_input2_lap<-renderUI({
    x=as.numeric(unlist(strsplit(gsub("\\s", "", input$mean),",")))
    se=as.numeric(unlist(strsplit(gsub("\\s", "", input$se),",")))
    lab=trimws( strsplit(input$lablabels,",")[[1]])
    if(length(lab)!=0){
      sanitize = !startsWith( lab,"-")
      tempvalue = median(se[sanitize])
    }else
    {
      tempvalue=median(se)
    }
    numericInput("halfCauchyScaleSigma_lap",
                 label = h4("Scale for half-Cauchy prior on within laboratory variances"),
                 # value =  median(labsesbayes()))
                 #value =  median(as.numeric(unlist(strsplit(gsub("\\s", "", input$se),",")))))
                 value=tempvalue)

  })


#######################################
  ######## Save and Load configuration file
  #######################################
  formReset = function()
  {
    inputIDs =c("lablabels", "mean", "se", "df", "units", "poolweights" )

    for(ii in inputIDs)
      session$sendInputMessage(ii,  list(value =""))

    session$sendInputMessage("seed",  list(value = "5"))
    session$sendInputMessage("coverage",  list(value ="0.95"))
    session$sendInputMessage("niters",  list(value ="250000"))
    session$sendInputMessage("nthin",  list(value ="25"))
    session$sendInputMessage("nburnin",  list(value ="50000"))
    session$sendInputMessage("DoEbootrep",  list(value ="10000"))
    updateCheckboxInput(session, "knhacheck", value = FALSE)
    updateCheckboxInput(session, "doboot", value = FALSE)
    updateCheckboxInput(session, "DoE", value = FALSE)
    updateRadioButtons(session, "LOODoE", selected = FALSE)

  }

  fileReadData = function(datapath)
  {
  	if (is.null(datapath))
  		return(NULL)
  	else{
  		fh = file(datapath, open = 'rt')
  		line = readLines(fh, n = 1)
  		savedInputs = list()
  		labManual = list()
  		while (length(line) != 0)
  		{
  			if (length(grep("=", line)) == 1)
  			{
  				key = strsplit(line, "=")[[1]][1]
  				value = paste(strsplit(line, "=")[[1]][-1], collapse = "=")

  				tempPair = list(key = value)
  				names(tempPair) = key
  				savedInputs = c(savedInputs, tempPair)
  			} else if (length(strsplit(line, ",")[[1]]) >= 2)
  			{
  				labManual = c(labManual, strsplit(line, ","))
  			}

  			line = readLines(fh, n = 1)
  		}

  		if (length(labManual)>0 && length(unique(unlist(lapply(  labManual, length  )))) != 1)
  			showNotification("Manual entry lines of different length detected, use csv format",type = "error")
  		else{
  			if (length(labManual) >= 2)
  			{
  				if (length(labManual[[1]]) == 2)
  				{
  					savedInputs$lablabels = ""
  					savedInputs$mean = paste(unlist(lapply(labManual, `[[`, 1)), collapse = ",")
  					savedInputs$se = paste(unlist(lapply(labManual, `[[`, 2)), collapse = ",")
  					savedInputs$df = ""

  				}

  				if (length(labManual[[1]]) == 3)
  				{
  					if (all(!is.na(suppressWarnings(as.numeric(
  						unlist(lapply(labManual, `[[`, 1))
  					)))))
  					{
  						savedInputs$lablabels = ""
  						savedInputs$mean = paste(unlist(lapply(labManual, `[[`, 1)), collapse = ",")
  						savedInputs$se = paste(unlist(lapply(labManual, `[[`, 2)), collapse = ",")
  						savedInputs$df = paste(unlist(lapply(labManual, `[[`, 3)), collapse = ",")
  					} else{
  						savedInputs$lablabels = paste(unlist(lapply(labManual, `[[`, 1)), collapse = ",")
  						savedInputs$mean = paste(unlist(lapply(labManual, `[[`, 2)), collapse = ",")
  						savedInputs$se = paste(unlist(lapply(labManual, `[[`, 3)), collapse = ",")
  						savedInputs$df = ""
  					}
  				}

  				if (length(labManual[[1]]) == 4)
  				{
  					savedInputs$lablabels = paste(unlist(lapply(labManual, `[[`, 1)), collapse = ",")
  					savedInputs$mean = paste(unlist(lapply(labManual, `[[`, 2)), collapse = ",")
  					savedInputs$se = paste(unlist(lapply(labManual, `[[`, 3)), collapse = ",")
  					savedInputs$df = paste(unlist(lapply(labManual, `[[`, 4)), collapse = ",")
  				}
  			}
  		}
  		close(fh)
  		inputIDs   = names(savedInputs)

  		inputvalues   = unlist(savedInputs)


  		if (length(inputIDs) >= 1)
  		{
  			formReset()

  			for (i in 1:length(inputIDs)) {
  				if (
              inputIDs[i] != "NICOB version" &&
              inputIDs[i] != "knhacheck" &&
              inputIDs[i] != "LOODoE" &&
              inputIDs[i] != "doboot" &&
              inputIDs[i] != "DoE")
  				{
  					session$sendInputMessage(inputIDs[i],  list(value = inputvalues[[i]]))
  				}
  			}


  			for(i in which(inputIDs=="knhacheck"))
  				updateCheckboxInput(session, "knhacheck", value = as.logical(inputvalues[[i]]))
  			for(i in which(inputIDs=="doboot"))
  				updateCheckboxInput(session, "doboot", value = as.logical(inputvalues[[i]]))
  			for(i in which(inputIDs=="DoE"))
  				updateCheckboxInput(session, "DoE", value = as.logical(inputvalues[[i]]))
  			for(i in which(inputIDs=="LOODoE"))
  				updateRadioButtons(session, "LOODoE", selected =  as.logical(inputvalues[[i]]))



  		}else {
  			showNotification("No usable information detected in the provided file",type = "error")
  		}


  }

  }

  observeEvent(input$file1, {
  	if (is.null(input$file1))
  		return(NULL)
  	else{
  		fileReadData(input$file1$datapath)
  	}
    }
  )

  observe({
  	query <- parseQueryString(session$clientData$url_search)
  	if (!is.null(query[['example']])) {
  		fileExample =	paste0 ("example/",query[['example']],".ncb")
  		if (file.exists(fileExample))
  			fileReadData(fileExample )
  	}

  })

  output$save_inputs <- downloadHandler(
    filename = 'consensus.ncb',
    content = function(file) {
      inputList = reactiveValuesToList(input)
      inputList = inputList[names(inputList) != "file1"]
      inputList = inputList[names(inputList) != "validate"]
      inputList = inputList[names(inputList) != "dobayes"]
      inputList = inputList[names(inputList) != "dobayes_lap"]
      inputList = inputList[names(inputList) != "go"]
      inputList = inputList[names(inputList) != "dopool"]

      outString = "NICOB version=1.4"
      for (ii in c("lablabels","mean","se","df","units","coverage","DoE"))
      {
        if( !is.null(inputList[[ii]]))
        {
          outString = c(outString, paste0(ii, "=", inputList[[ii]]))
          inputList[[ii]] =NULL
        }

      }

      for (ii in 1:length(inputList))
      {
        outString = c(outString, paste0(names(inputList[ii]), "=", inputList[ii]))

      }

      write(outString, file = file, append = FALSE)
    }
  )

  #######################################
  ######## Validate
  #######################################



  val<- eventReactive(input$validate,{

    labmeans<-as.numeric(unlist(strsplit(gsub("\\s", "", input$mean),",")))

    labses<-as.numeric(unlist(strsplit(gsub("\\s", "", input$se),",")))

    labdfs<-as.numeric(unlist(strsplit(gsub("\\s", "", input$df),",")))

    lablabs<-strsplit(input$lablabels,",")[[1]]

    #################################
    ####### check inputs #########
    #################################
    validate(
      need(!any(is.na(labmeans)), 'Measured values must be numeric, separated by commas.'),
      need(!any(is.na(labses)) && !any(labses<=0), 'Standard uncertainties must be numeric (greater than zero), separated by commas.'),
      need(!any(is.na(labdfs)) && !any(labdfs<1), 'Degrees of freedom must be numeric (greater than or equal to 1), separated by commas.'),
      need(is.numeric(input$coverage) && input$coverage<=1 && input$coverage>=0, 'Invalid coverage probability, must be numeric between 0 and 1.'),
      need(length(labmeans)==length(labses), 'The number of measured values entered does not equal the number of standard uncertainties entered.'),
      need(length(labmeans)==length(lablabs) | length(lablabs)==0,'The number of study labels entered does not equal the number of measured values entered.'),
      need(!any(table(lablabs)>1) | length(lablabs)==0,"Laboratory labels should be unique."),
      need(length(labmeans)==length(labdfs) || length(labdfs)==0, "The number of degrees of freedom entered does not equal the number of measured values entered."),
      need(length(labmeans)>1, 'At least two measured values are required.'),
      need(length(labmeans)<=500, 'The NICOB will not combine more than 500 measured values.')


    )

    validate(
      need(is.na(input$DoEbootrep) | (is.numeric(input$DoEbootrep) && input$DoEbootrep>0 && input$DoEbootrep<=50000), 'Invalid number of bootstrap replicates for degrees of equivalence calculation, must be positive integer (less than 50000).')
    ) #not a perfect fix, still says it is ok if characters are entered...


    # icon("check","fa-2x")
    "Model inputs are valid."
    })

  output$val <- renderText({ paste("<font color=\"#54C571\"><b>", val(), "</b></font>") })


  #######################################
  ####### Write coverage probabilities on output side. Shiny apparently can only use each of these once, so I define 3 here.
  ####### Output for rma must be defined differently because I don't want the label to update until the results update.
  ####### Other three intervals update automatically when coverage probability changes.
  #######################################

  output$coverageProbabilityPercentBoot <- renderText({
    paste(input$coverage*100,"%",sep="")
  })
  output$coverageProbabilityPercentBayes <- renderText({
    paste(input$coverage*100,"%",sep="")
  })
  output$coverageProbabilityPercentPool <- renderText({
    paste(input$coverage*100,"%",sep="")
  })





  #######################################
  #######################################
  ####### rma
  #######################################
  #######################################


  labmeans<-eventReactive(input$go,{
    as.numeric(unlist(strsplit(gsub("\\s", "", input$mean),",")))
  }
  )

  labses<-eventReactive(input$go,{
    (as.numeric(unlist(strsplit(gsub("\\s", "", input$se),","))))
  }
  )

  labdfs<-eventReactive(input$go,{
    (as.numeric(unlist(strsplit(gsub("\\s", "", input$df),","))))
  }
  )

  lablabs<-eventReactive(input$go,{
    trimws( strsplit(input$lablabels,",")[[1]])
  }
  )

  seedrma<-eventReactive(input$go,{
    as.numeric(input$seed)
  })

  signifdl<-eventReactive(input$go,{
    as.numeric(input$signifDigits)
  })


  outrma<- eventReactive(input$go,{

    #################################
    ####### check inputs #########
    #################################
    validate(
      need(!any(is.na(labmeans())), 'Measured values must be numeric, separated by commas.'),
      need(!any(is.na(labses())) && !any(labses()<=0), 'Standard uncertainties must be numeric (greater than zero), separated by commas.'),
      need(!any(is.na(labdfs())) && !any(labdfs()<1), 'Degrees of freedom must be numeric (greater than or equal to 1), separated by commas.'),
      need(is.numeric(input$coverage) && input$coverage<=1 && input$coverage>=0, 'Invalid coverage probability, must be numeric between 0 and 1.'),
      need(length(labmeans())==length(labses()), 'The number of measured values entered does not equal the number of standard uncertainties entered.'),
      need(length(labmeans())==length(lablabs()) | length(lablabs())==0,'The number of study labels entered does not equal the number of measured values entered.'),
      need(!any(table(lablabs())>1) | length(lablabs())==0,"Laboratory labels should be unique."),
      need(length(labmeans())==length(labdfs()) || length(labdfs())==0, "The number of degrees of freedom entered does not equal the number of measured values entered."),
      need(length(labmeans())>1, 'At least two measured values are required.'),
      ####upper limits
      need(length(labmeans())<=500, 'The NICOB will not combine more than 500 measured values.'),
      need(!is.na(seedrma()), 'Random number generator seed must be numeric.')
    )


    set.seed(seedrma() )

    source("modifiedKHmethod.R")

  	if(length(lablabs())==0){
  	  lab=paste("L", 1:length(labmeans()), sep="")
  	}else{
  	  lab=lablabs()
  	}
  	#we remove the values of labs starting with "-"
    sanitize = !startsWith( lab,"-")
    nI=length(labmeans()[sanitize])
    ###Check number of included labs
    LOObool= as.logical(input$LOODoE)

    validate(
      need((nI>=3&&LOObool)||!LOObool, 'At least three labs are required to be included in the consensus for this method.')
    )



    if(input$knhacheck){
      modifiedKHmethod(x = labmeans()[sanitize],u = labses()[sanitize], coverageProb = input$coverage)
    }else{
      rma(labmeans()[sanitize],sei=labses()[sanitize],method="DL",level=input$coverage*100)
    }


  })


  rmacp<- eventReactive(input$go,{
    input$coverage

  })

  output$rmaest <- renderText({
    paste(signif(outrma()$b,input$signifDigits))

  })

  output$rmase <- renderText({
    paste(signif(outrma()$se,input$signifDigits))

  })

  output$rmaquant <- renderText({
    paste("The ",rmacp()*100,"% coverage interval ranges from: ",
          signif(outrma()$ci.lb,input$signifDigits),
          " to ",
          signif(outrma()$ci.ub,input$signifDigits),
          sep="")

  })

  output$rmatau <- renderText({
  	paste(signif(sqrt(outrma()$tau2),input$signifDigits))

  })



    rmaplot_func=function(){

    ################################
    ###### check inputs #########
    #################################
    validate(
      need(!any(is.na(labmeans())), 'Measured values must be numeric, separated by commas.'),
      need(!any(is.na(labses()))&& !any(labses()<=0), 'Standard uncertainties must be numeric (greater than zero), separated by commas.'),
      need(!any(is.na(labdfs())) && !any(labdfs()<1), 'Degrees of freedom must be numeric (greater than or equal to 1), separated by commas.'),
      need(is.numeric(input$coverage) && input$coverage<=1 && input$coverage>=0, 'Invalid coverage probability, must be numeric between 0 and 1.'),
      need(length(labmeans())==length(labses()), 'The number of measured values entered does not equal the number of standard uncertainties entered.'),
      need(length(labmeans())==length(labdfs()) || length(labdfs())==0, "The number of degrees of freedom entered does not equal the number of measured values entered."),
      # need(length(labmeans())==length(strsplit(input$lablabels,",")[[1]]) | length(strsplit(input$lablabels,",")[[1]])==0, 'The number of study labels entered does not equal the number of measured values entered.'),
      need(length(labmeans())==length(lablabs()) | length(lablabs())==0,'The number of study labels entered does not equal the number of measured values entered.'),
      need(!any(table(lablabs())>1) | length(lablabs())==0,"Laboratory labels should be unique."),
      need(length(labmeans())>1, 'At least two measured values are required.'),
      ####upper limits
      need(length(labmeans())<=500, 'The NICOB will not combine more than 500 measured values.')


    )

    par(fig=c(0,1,.6,1),lend=2)

    plot(c(1,1), c(1,1), type="n", axes=FALSE,xlab="", ylab="", bty="n")
    leg=legend("top",c(expression(paste("Consensus estimate, ", hat(mu),", and interval ",hat(mu)%+-%u(mu),sep="")), #")),#expression(paste("Coverage interval, ",hat(mu)%+-%u(mu),sep="")),
                       expression(paste("Measured value, ",x[j])),
                       expression(paste(x[j]%+-%u[j])),
                       expression(paste(x[j]%+-%sqrt(u[j]^2+hat(tau)^2)))
    ),
    lty=c(0,0,1,1),pch=c(15,19,-1,-1),col=c(col1,col3,col4,col5),pt.cex=c(3,1.5,3,1), lwd=c(1,1,11,5))

    legend("top",
           c(expression(paste("Consensus estimate, ", hat(mu),", and interval ",hat(mu)%+-%u(mu),sep="")), #")),#expression(paste("Coverage interval, ",hat(mu)%+-%u(mu),sep="")),
             expression(paste("Measured value, ",x[j])),
             expression(paste(x[j]%+-%u[j])),
             expression(paste(x[j]%+-%sqrt(u[j]^2+hat(tau)^2)))
           ),
           lty=c(1,0,1,1),pch=c(-1,19,-1,-1),col=c(col2,col3,col4,col5),pt.cex=c(3,1.5,3,1),lwd=c(1,1,11,5),bty="n")


    par(fig=c(0,1,0,.8), new=TRUE,lend=2)

    ##########new plot
    mu.x = outrma()$b
    mu.u = outrma()$se
    tau.x = sqrt(outrma()$tau2)
    # sigma.index = grep("sigma", dimnames(PCB28.Bayes$summary)[[1]])

	if(length(lablabs())==0){
      lab=paste("L", 1:length(labmeans()), sep="")
    }else{
      lab=lablabs()
    }

    sanitize = !startsWith( lab,"-")

	#we reorganize the lab to put those with a minus at the end

    lab = c(lab[sanitize],lab[!sanitize])
    x = c(labmeans()[sanitize],labmeans()[!sanitize])
    u = c(labses()[sanitize],labses()[!sanitize])

  	n = length(x)
  	#Number of labs taken into account for the consensus
    nI = sum(sanitize)

	  sigma.x=u

    xl = x-u
    xu = x+u

    xl.Sigma.Tau = x-sqrt(sigma.x^2 + tau.x^2)
    xu.Sigma.Tau = x+sqrt(sigma.x^2 + tau.x^2)


    if(input$units==""){
      label="Measured value"
    }else{
      label=paste("Measured value (",input$units,")",sep="")
    }

    plot(c(0.5,n+0.5), range(c(xl.Sigma.Tau, xu.Sigma.Tau)),
         axes=FALSE, xlab="",
         ylab=label, type="n", bty="n")
    axis(2)
    polygon(c(1-0.2, nI+0.2, nI+0.2, 1-0.2),
            c(mu.x-mu.u, mu.x-mu.u, mu.x+mu.u, mu.x+mu.u),
            col=col1, border=col1)
    abline(h=mu.x, col=col2)
    segments(1:n, xl, 1:n, xu, col=col4, lwd=11)
    segments(1:n, xl.Sigma.Tau, 1:n, xu.Sigma.Tau, col=col5, lwd=5)
    points(1:nI, x[1:nI], pch=19, cex=1.5, col=col3)
    points(nI:n, x[nI:n], pch=21, cex=1.5, col=col3)
    mtext(lab[seq(1, n, 2)], side=1, at=seq(1, n, 2), line=0, cex=1)
    mtext(lab[seq(2, n, 2)], side=1, at=seq(2, n, 2), line=1.2, cex=1)


    }

    output$rmaplot<- renderPlot({rmaplot_func()
    })

  output$downloadrmagraph <- downloadHandler(
    filename = function() { paste( 'DerSimonianLaird.pdf', sep='') },
    content = function(file) {
      pdf(file,height=9, width=8)
      rmaplot_func()
      dev.off()
    })
  #######################################
  #######################################
  ####### parametric bootstrap
  #######################################
  #######################################
  source("sampleFromTau2Dist.R")

  bootDL=function(K,thedat,themle){

    muB = numeric(K)
    n=nrow(thedat)
    indexInf = (thedat$nu == Inf)

    for (k in 1:K)
    {
      if (k %% 100 == 0) {incProgress(100/K, detail = paste("Iteration", k))}
      tau2B = sampleFromTau2Dist(thedat$x,thedat$u)
      xB = rnorm(n, mean=themle$mu, sd=sqrt(tau2B+thedat$u^2))
      uB=rep(NA,n)
      if (is.null(thedat$nu)) {uB = thedat$u} else {
        uB[indexInf] = thedat$u[indexInf]
        uB[!indexInf] = thedat$u[!indexInf] *
          sqrt(thedat$nu[!indexInf]/rchisq(sum(!indexInf), df=thedat$nu[!indexInf]))} #simulate new uncertainties


      muB[k] = rma(yi=xB, sei=uB, method="DL",knha=input$knhacheck)$b ### don't need modified version, only outputting b
    }

    return(muB)
  }

  outboot<- eventReactive(input$go,{
    if(is.numeric(input$bootrep)){
      #################################
      ####### check inputs #########
      #################################
      validate(
        need(!any(is.na(labmeans())), 'Measured values must be numeric, separated by commas.'),
        need(!any(is.na(labses()))&& !any(labses()<=0), 'Standard uncertainties must be numeric (greater than zero), separated by commas.'),
        need(!any(is.na(labdfs())) && !any(labdfs()<1), 'Degrees of freedom must be numeric (greater than or equal to 1), separated by commas.'),
        need(is.numeric(input$coverage) && input$coverage<=1 && input$coverage>=0, 'Invalid coverage probability, must be numeric between 0 and 1.'),
        need(length(labmeans())==length(labses()), 'The number of measured values entered does not equal the number of standard uncertainties entered.'),
        need(length(labmeans())==length(labdfs()) || length(labdfs())==0, "The number of degrees of freedom entered does not equal the number of measured values entered."),
        # need(length(labmeans())==length(strsplit(input$lablabels,",")[[1]]) | length(strsplit(input$lablabels,",")[[1]])==0, 'The number of study labels entered does not equal the number of measured values entered.'),
        need(length(labmeans())==length(lablabs()) | length(lablabs())==0,'The number of study labels entered does not equal the number of measured values entered.'),
        need(!any(table(lablabs())>1) | length(lablabs())==0,"Laboratory labels should be unique."),
        need(is.numeric(input$bootrep) && input$bootrep>0  && input$bootrep<=50000, 'Invalid number of bootstrap replicates, must be positive integer (less than 50000).'),
        need(length(labmeans())>1, 'At least two measured values are required.'),
        ####upper limits
        need(length(labmeans())<=500, 'The NICOB will not combine more than 500 measured values.')


      )

      ##################################

      withProgress(message = 'Running bootstrap', value = 0, {

        if(length(labdfs())==0){
          z=data.frame(x=labmeans(),u=labses())
        }else{
          z=data.frame(x=labmeans(),u=labses(),nu=labdfs())
        }

        bootDL(input$bootrep,z,list(mu=outrma()$b, tau2=outrma()$tau2))

      })
    }else{
      validate(
        need(is.numeric(input$bootrep) && input$bootrep>0 && input$bootrep<=50000,"Please enter the number of bootstrap replicates.")
      )
    }




  })


  output$textboot<-renderText({
    if(!is.null(outboot())){
      bootsd=sd(outboot(),na.rm=T)
      signifDigits = significator(outboot(),digits=TRUE)
      if (signifDigits==1)
      {
        signifString= paste(" (where 1 significant digit is believed to be reliable)")
      }else{
        signifString= paste(" (where", signifDigits, "significant digits are believed to be reliable)")
      }
      paste(signif(bootsd,input$signifDigits), signifString)
    }
  })

  output$bootquant <- renderText({
    bootcov=quantile(outboot(),
                     # c(.025,.975),
                     c((1-input$coverage)/2, (1+input$coverage)/2),
                     na.rm=T)

    paste(signif(bootcov[1],input$signifDigits),"to",signif(bootcov[2],input$signifDigits))
  })


  output$downloadDLbootout <- downloadHandler(
    filename = "DLout.csv",
    content = function(file) {
      write.csv(outboot(), file)
    }
  )



  #######################################
  #######################################
  ####### DerSimonian-Laird DoE
  #######################################
  #######################################

  source("symmetricalBootstrapCI.R")
  source("DoEUnilateralDerSimonianLaird.R")

  DLDoEUni=eventReactive(input$go   ,{
    validate(
      need(!any(is.na(labmeans())), 'Measured values must be numeric, separated by commas.'),
      need(!any(is.na(labses())) && !any(labses()<=0), 'Standard uncertainties must be numeric (greater than zero), separated by commas.'),
      need(!any(is.na(labdfs())) && !any(labdfs()<1), 'Degrees of freedom must be numeric (greater than or equal to 1), separated by commas.'),
      need(is.numeric(input$coverage) && input$coverage<=1 && input$coverage>=0, 'Invalid coverage probability, must be numeric between 0 and 1.'),
      need(length(labmeans())==length(labses()), 'The number of measured values entered does not equal the number of standard uncertainties entered.'),
      need(length(labmeans())==length(labdfs()) || length(labdfs())==0, "The number of degrees of freedom entered does not equal the number of measured values entered."),
      need(length(labmeans())==length(lablabs()) | length(lablabs())==0,'The number of study labels entered does not equal the number of measured values entered.'),
      need(!any(table(lablabs())>1) | length(lablabs())==0,"Laboratory labels should be unique."),
      need(length(labmeans())>1, 'At least two measured values are required.'),
      need(length(labmeans())<=500, 'The NICOB will not combine more than 500 measured values.')
    )

    validate(
      need(is.numeric(input$DoEbootrep) && input$DoEbootrep>0 && input$DoEbootrep<=50000, 'Invalid number of bootstrap replicates for degrees of equivalence calculation, must be positive integer (less than 50000).')
    )
    if(length(lablabs())==0){
      lab=paste("L", 1:length(labmeans()), sep="")
    }else{
      lab=lablabs()
    }
    #we remove the values of labs starting with "-"
    sanitize = !startsWith( lab,"-")
    nI=length(labmeans()[sanitize])
    ###Check number of included labs
    LOObool= as.logical(input$LOODoE)

    validate(
      need((nI>=3&&LOObool)||!LOObool, 'At least three labs are required to be included in the consensus for this method.')
    )


    withProgress(message = 'Calculating DoEs', value = 0, {



      if(length(lablabs())==0){
        lab=paste("L", 1:length(labmeans()), sep="")
      }else{
        lab=lablabs()
      }

      if(length(labdfs())!=0){
        nu=labdfs()
      }else{
        nu=NULL
        # nu=rep(NA,length(labmeans()))
      }

      DoEUnilateralDL(x=labmeans(), u=labses(), nu=nu,
                      lab=lab, K=input$DoEbootrep,LOO=as.logical(input$LOODoE),coverageProb=.95, DLRes=outrma()) # hard-coded coverageProb to match MRA definition and user's manual

    })



  })

  output$DLDoEConv <- renderUI({
    HTML(DLDoEUni()[["DoEwarn"]])
  })


  output$DLDoEUniTable <- renderTable(DLDoEUni()$DoE,
                                      digits=signifdl,
                                      display=rep("fg",6),
                                      caption = '<b> Unilateral Degrees of Equivalence',
                                      caption.placement='top',include.rownames=FALSE)


  DLDoEUniPlot_func=function(){
    par(fig=c(0,1,.6,1),lend=2)  #change this when DoE plot added
    plot(c(1,1), c(1,1), type="n", axes=FALSE,xlab="", ylab="", bty="n")
    legend("top",c("DoE estimate" ," 95% coverage interval cla"),
           lty=c(NA,1),pch=c(19,NA),col=c(col3,col4),pt.cex=c(1,1),lwd=c(5,5))

    par(fig=c(0,1,0,.8), new=TRUE,lend=2)

    plot(0,xlim=range(0.7,length(DLDoEUni()$DoE$DoE.x)+0.3),ylim=range(DLDoEUni()$DoE$DoE.Lwr,DLDoEUni()$DoE$DoE.Upr),type="n",     ylab="DoE",xlab="",xaxt="n", bty="n",axes=FALSE)
    abline(h=0,col="gray")

    if(length(lablabs())==0){
      lab=paste("L", 1:length(labmeans()), sep="")
    }else{
      lab=lablabs()
    }
    #we remove the values of labs starting with "-"
    sanitize = !startsWith( lab,"-")
    nI=length(labmeans()[sanitize])
    n=length(labmeans())

    lab=DLDoEUni()$DoE$Lab # use lab labels that were rearranged in DoE code, not original labels

    n = length(labmeans())
  	mtext(lab[seq(1, n, 2)], side=1, at=seq(1, n, 2), line=0, cex=1, col="Black")
    mtext(lab[seq(2, n, 2)], side=1, at=seq(2, n, 2), line=1.2, cex=1, col="Black")


    axis(2)

    for(i in 1:length(DLDoEUni()$DoE$DoE.x)){
      arrows(i,DLDoEUni()$DoE$DoE.Lwr[i],i,DLDoEUni()$DoE$DoE.Upr[i],length=0, col=col4, lwd=5)
    }
    points(x=1:length(DLDoEUni()$DoE$DoE.x),y=DLDoEUni()$DoE$DoE.x,pch=c(rep(19,nI),rep(21,n-nI)), cex=1.5, col=col3)


  }
    output$DLDoEUniPlot<- renderPlot({DLDoEUniPlot_func()})



  source("DoEBilateralDerSimonianLaird.R")


  DLDoEBi=eventReactive(input$go,{
    validate(
      need(!any(is.na(labmeans())), 'Measured values must be numeric, separated by commas.'),
      need(!any(is.na(labses())) && !any(labses()<=0), 'Standard uncertainties must be numeric (greater than zero), separated by commas.'),
      need(!any(is.na(labdfs())) && !any(labdfs()<1), 'Degrees of freedom must be numeric (greater than or equal to 1), separated by commas.'),
      need(is.numeric(input$coverage) && input$coverage<=1 && input$coverage>=0, 'Invalid coverage probability, must be numeric between 0 and 1.'),
      need(length(labmeans())==length(labses()), 'The number of measured values entered does not equal the number of standard uncertainties entered.'),
      need(length(labmeans())==length(labdfs()) || length(labdfs())==0, "The number of degrees of freedom entered does not equal the number of measured values entered."),
      # need(length(labmeans())==length(strsplit(input$lablabels,",")[[1]]) | length(strsplit(input$lablabels,",")[[1]])==0, 'The number of study labels entered does not equal the number of measured values entered.'),
      need(length(labmeans())==length(lablabs()) | length(lablabs())==0,'The number of study labels entered does not equal the number of measured values entered.'),
      need(!any(table(lablabs())>1) | length(lablabs())==0,"Laboratory labels should be unique."),
      need(length(labmeans())>1, 'At least two measured values are required.'),
      # need(!any(is.na(labmeans())), ''),
      # need(!any(is.na(labses())), ''),
      # need(is.numeric(input$coverage) && input$coverage<=1 && input$coverage>=0, '')
      need(length(labmeans())<=500, 'The NICOB will not combine more than 500 measured values.')
    )
    if(length(lablabs())==0){
      lab=paste("L", 1:length(labmeans()), sep="")
    }else{
      lab=lablabs()
    }
    #we remove the values of labs starting with "-"
    sanitize = !startsWith( lab,"-")
    nI=length(labmeans()[sanitize])
    ###Check number of included labs
    LOObool= as.logical(input$LOODoE)

    validate(
      need((nI>=3&&LOObool)||!LOObool, 'At least three labs are required to be included in the consensus for this method.')
    )

    withProgress(message = 'Calculating Bilateral DoEs', value = 0, {

      # removed lab input, don't use it in the function.

      if(length(labdfs())!=0){
        nu=labdfs()
      }else{
        nu=NULL
      }

      DoEBilateralDL(x=labmeans(), u=labses(), nu=nu,
                     K=input$DoEbootrep,coverageProb=.95,
                     DoEUnilateral=DLDoEUni())
    })

  })


  bold.allrows <- function(x) {
    h <- paste('<b>',x,'', sep ='')
    h
  }



  output$DLDoEBiTable <- renderTable(DLDoEBi()$B.x,
                                     digits = signifdl,
                                     rownames = T,
                                     caption = '<b> Bilateral Degrees of Equivalence (Estimates)',
                                     caption.placement='top',
                                     sanitize.rownames.function =  bold.allrows)

  output$DLDoE.UBiTable <- renderTable(DLDoEBi()$B.U,
                                       digits = signifdl,
                                       rownames = T,
                                       caption = '<b> Bilateral Degrees of Equivalence (Expanded Uncertainties U95)',
                                       caption.placement='top',
                                       sanitize.rownames.function =  bold.allrows)


  DLDoEBiPlot_func=function(){
    DoE.Bilateral.DL.x = DLDoEBi()$B.x
    DoE.Bilateral.DL.Lwr = DLDoEBi()$B.Lwr
    DoE.Bilateral.DL.Upr = DLDoEBi()$B.Upr
    n = dim(DoE.Bilateral.DL.x)[[1]]
    labNAMEs = dimnames(DoE.Bilateral.DL.x)[[1]]
    xx = (DoE.Bilateral.DL.Lwr * DoE.Bilateral.DL.Upr > 0)
    ij = which(xx, arr.ind=TRUE, useNames = TRUE)
    xx[xx == 1] = 0.9
    xx[xx == 0] = 0.1
    diag(xx) = rep(0.5, n)
    image(1:n, 1:n, xx, xlim=c(-1.5, n+1),breaks=c(-0.5, 0.4, 0.85, 1.5),
          col=c(col5,"grey36",col4),
          axes=FALSE, xlab="", ylab="", asp=1)
    K = nrow(ij)
    if (K > 0) {
      for (k in 1:K) {
        points(ij[k,1], ij[k,2], pch=8, font=3, col="Black") } }
    polygon(c(0.5,n+0.5,n+0.5,0.5), c(0.5,0.5,n+0.5,n+0.5),
            border="Black", lwd=2)

    segments(seq(1.5,n+0.5),0.5,seq(1.5,n+0.5),n+0.5)
    segments(n+0.5,seq(1.5,n+0.5),0.5,seq(1.5,n+0.5))

    mtext(labNAMEs, side=1, at=1:n, las=2, line=1.5, col="Blue", font=2,padj=0)
    text(rep(0, n), 1:n, row.names(DLDoEBi()$B.x), adj=1, col="Blue", font=2)

    }

  output$DLDoEBiPlot<- renderPlot({DLDoEBiPlot_func()})


  output$downloadDLUniDoEGraph <- downloadHandler(
    filename = function() { paste( 'DerSimonianLaird_Unilateral_DoE.pdf', sep='') },
    content = function(file) {
      pdf(file,height=8, width=8)
      DLDoEUniPlot_func()
      dev.off()
    })
  output$downloadDLBiDoEGraph <- downloadHandler(
    filename = function() { paste( 'DerSimonianLaird_Bilateral_DoE.pdf', sep='') },
    content = function(file) {
      pdf(file,height=8, width=8)
      DLDoEBiPlot_func()
      dev.off()
    })

  #######################################
  #######################################
  ########### bayesian
  #######################################
  #######################################
  source("bayesGelman.R")

  labmeansbayes<-eventReactive(input$dobayes,{
    (as.numeric(unlist(strsplit(gsub("\\s", "", input$mean),","))))
  }
  )
  labsesbayes<-eventReactive(input$dobayes,{
    (as.numeric(unlist(strsplit(gsub("\\s", "", input$se),","))))
  }
  )
  labdfsbayes<-eventReactive(input$dobayes,{
    (as.numeric(unlist(strsplit(gsub("\\s", "", input$df),","))))
  }
  )

  lablabsbayes<-eventReactive(input$dobayes,{
    trimws( strsplit(input$lablabels,",")[[1]])
  }
  )
  seedbayes<-eventReactive(input$dobayes,{
    as.numeric(input$seed)
  })
  signifbayes<-eventReactive(input$dobayes,{
    as.numeric(input$signifDigits)
  })



  outbayes<- eventReactive(input$dobayes,{


    #################################
    ####### check inputs #########
    #################################
    validate(
      need(!any(is.na(labmeansbayes())), 'Measured values must be numeric, separated by commas.'),
      need(!any(is.na(labsesbayes()))&& !any(labsesbayes()<=0), 'Standard uncertainties must be numeric (greater than zero), separated by commas.'),
      need(!any(is.na(labdfsbayes())) && !any(labdfsbayes()<1), 'Degrees of freedom must be numeric (greater than or equal to 1), separated by commas.'),
      need(is.numeric(input$coverage) && input$coverage<=1 && input$coverage>=0, 'Invalid coverage probability, must be numeric between 0 and 1.'),
      need(length(labmeansbayes())==length(labsesbayes()),  'The number of measured values entered does not equal the number of standard uncertainties entered.'),
      need(length(labmeansbayes())==length(labdfsbayes()) || length(labdfsbayes())==0, "The number of degrees of freedom entered does not equal the number of measured values entered."),

      need(length(labmeansbayes())==length(lablabsbayes()) | length(lablabsbayes())==0,'The number of study labels entered does not equal the number of measured values entered.'),
      need(!any(table(lablabsbayes())>1) | length(lablabsbayes())==0,"Laboratory labels should be unique."),

      need(input$niters>0 && input$niters<=1000000,'Total number of iterations must be a positive integer (less than 1000000).'),
      need(input$nburnin>=0 && input$nburnin<=500000,'Length of burn in must be a positive integer (less than 500000).'),
      need(input$nthin>0 && input$nthin<=1000,'Thinning rate must be a positive integer (less than 1000).'),
      need(input$nburnin<=input$niters/2,paste("Length of burn in too large (over half of the total number of iterations). Increase number of iterations to at least ",input$nburnin*2,".",sep="")),
      need(length(labmeansbayes())<=500, 'The NICOB will not combine more than 500 measured values.'),
      need(!is.na(seedbayes()), 'Random number generator seed must be numeric.'),
      ###Number of labs
      need(length(labmeansbayes())>1, 'At least two measured values are required.')


    )

    #################################
    set.seed(seedbayes())

    if(length(lablabsbayes())==0){
      lab=paste("L", 1:length(labmeansbayes()), sep="")
    }else{
      lab=lablabsbayes()
    }

    sanitize = !startsWith( lab,"-")

    #we reorganize the lab to remove those who starts with -

    lab = lab[sanitize]
    x = labmeansbayes()[sanitize]
    u = labsesbayes()[sanitize]

    nI = length(x)

    ###Check number of included labs

    LOObool= as.logical(input$LOODoE)

    validate(
      need((nI>=3&&LOObool)||!LOObool, 'At least three labs are required to be included in the consensus for this method.'),
      need((nI>=2&&!LOObool)||LOObool, 'At least two labs are required to be included in the consensus for this method.')
    )


    if(length(labdfsbayes())!=0){
      nu=labdfsbayes()[sanitize]
    }else{
      nu=NULL
    }

    if(is.na(input$halfCauchyScaleTau)){
      input_pscaletau=mad(x)#median(abs(labmeansbayes()-median(labmeansbayes()))) #tauPriorScale=mad(x)
      priorTauWarn=paste("Scale for half-Cauchy prior on between laboratory variance has either not
                      been entered or is not numeric. Instead, the median absolute deviation of the
                         measured values has been used: ",input_pscaletau,sep="")
    }else{
      input_pscaletau=input$halfCauchyScaleTau
      priorTauWarn=""
    }
    output$bayesTauWarn <- renderUI({
      HTML(priorTauWarn)
    })

    if(is.na(input$halfCauchyScaleSigma)){
      input_pscalesig=median(u)
      priorSigWarn=paste("Scale for half-Cauchy prior on within laboratory variances has either not
                      been entered or is not numeric. Instead, the median of the lab-specific standard
                      uncertainties has been used: ",input_pscalesig,sep="")
    }else{
      input_pscalesig=input$halfCauchyScaleSigma
      priorSigWarn=""
    }

    output$bayesSigWarn <- renderUI({
      HTML(priorSigWarn)
    })


    validate(
      need(input_pscaletau>0, 'Scale for half-Cauchy prior on between laboratory variance must be a positive number.'),
      need(input_pscalesig>0, 'Scale for half-Cauchy prior on within laboratory variances must be a positive number.')
    )

    bayesGelman(x,u,nu,tauPriorScale = input_pscaletau,sigmaPriorScale = input_pscalesig,ni = input$niters,nb = input$nburnin,nt = input$nthin)

  })


  output$downloadbayesout <- downloadHandler(
    filename = "MCMCout.csv",
    content = function(file) {
      write.table(outbayes()[["mcmcout"]][[1]], file,row.names = F)
    }
  )

  output$bayesConv <- renderUI({
    HTML(outbayes()[["warn"]])
  })


  output$bayesest <- renderText({
    bayesres=outbayes()[["mcmcout"]]
    signifDigits = significator(bayesres[[1]][,"mu"],digits=TRUE)
    if (signifDigits==1)
    {
      signifString= paste(" (where 1 significant digit is believed to be reliable)")
    }else{
     signifString= paste(" (where", signifDigits, "significant digits are believed to be reliable)")
    }
    paste(signif(mean(bayesres[[1]][,"mu"]),input$signifDigits), signifString)

  })

  output$bayesse <- renderText({
    bayesres=outbayes()[["mcmcout"]]
    paste(signif(sd(bayesres[[1]][,"mu"]),input$signifDigits))

  })

  output$bayesquant <- renderText({
    bayesres=outbayes()[["mcmcout"]]
    credible.interval=quantile(bayesres[[1]][,"mu"],c((1-input$coverage)/2, (1+input$coverage)/2))

    paste(signif(credible.interval[1],digits = input$signifDigits),
          "to",
          signif(credible.interval[2],digits = input$signifDigits))

  })

  output$bayestau<- renderText({
  	bayesres=outbayes()[["mcmcout"]]
  	postsumstats=summary(bayesres)$statistics

  	bayesres=outbayes()[["tau"]]
  	paste(signif(postsumstats["tau", "Mean"],input$signifDigits))

  })

  bayesplot_func=function(){
    bayesres=outbayes()[["mcmcout"]]

    postsumstats=summary(bayesres)$statistics
    if(length(lablabsbayes())==0){
      lab=paste("L", 1:length(labmeansbayes()), sep="")
    }else{
      lab=lablabsbayes()
    }

    sanitize = !startsWith( lab,"-")

    #we reorganize the lab to put those with a minus at the end

    lab = c(lab[sanitize],lab[!sanitize])
    x = c(labmeansbayes()[sanitize],labmeansbayes()[!sanitize])
    u = c(labsesbayes()[sanitize],labsesbayes()[!sanitize])

    n = length(x)
    #Number of labs taken into account for the consensus
    nI = length(labmeansbayes()[sanitize])

    postsumstats=summary(bayesres)$statistics
    postparameters=rownames(postsumstats)

    mu.x = postsumstats["mu","Mean"]
    mu.u = postsumstats["mu", "SD"]
    tau.x = postsumstats["tau", "Mean"]

    sigma.index = grep("sigma", dimnames(postsumstats)[[1]])
    sigma.names=rep(NA,nI)
    for(i in 1:nI){
      sigma.names[i] = paste("sigma[",i, "]",sep="")
    }

    if(length(sigma.index)==0){
      sigma.x = labsesbayes()[sanitize]
    }else{
      sigma.x = postsumstats[sigma.names, "Mean"]
    }


    if(input$units==""){
      label="Measured Value"
    }else{
      label=paste("Measured Value (",input$units,")",sep="")
    }


    par(fig=c(0,1,.6,1),lend=2)

    plot(c(1,1), c(1,1), type="n", axes=FALSE,xlab="", ylab="", bty="n")
    if(length(sigma.index)==0){
      leg=legend("top",c(expression(paste("Consensus estimate, ", hat(mu),", and interval ",hat(mu)%+-%u(mu),sep="")),
                         expression(paste("Measured value, ",x[j])),
                         expression(paste(x[j]%+-%u[j])),
                         expression(paste(x[j]%+-%sqrt(u[j]^2+hat(tau)^2)))
      ),
      lty=c(0,0,1,1),pch=c(15,19,-1,-1),col=c(col1,col3,col4,col5),pt.cex=c(3,1.5,3,1),lwd=c(1,1,2,.5))

      legend("top",
             c(expression(paste("Consensus estimate, ", hat(mu),", and interval ",hat(mu)%+-%u(mu),sep="")),
               expression(paste("Measured value, ",x[j])),
               expression(paste(x[j]%+-%u[j])),
               expression(paste(x[j]%+-%sqrt(u[j]^2+hat(tau)^2)))
             ),
             lty=c(1,0,1,1),pch=c(-1,19,-1,-1),col=c(col2,col3,col4,col5),pt.cex=c(3,1.5,3,1),lwd=c(1,1,11,5),bty="n")


    }else{
      leg=legend("top",c(expression(paste("Consensus estimate, ", hat(mu),", and interval ",hat(mu)%+-%u(mu),sep="")),
                         expression(paste("Measured value, ",x[j])),
                         expression(paste(x[j]%+-%hat(sigma)[j])),
                         expression(paste(x[j]%+-%sqrt(hat(sigma)[j]^2+hat(tau)^2)))
      ),
      lty=c(0,0,1,1),pch=c(15,19,-1,-1),col=c(col1,col3,col4,col5),pt.cex=c(3,1.5,3,1),lwd=c(1,1,11,5))

      legend("top",
             c(expression(paste("Consensus estimate, ", hat(mu),", and interval ",hat(mu)%+-%u(mu),sep="")),
               expression(paste("Measured value, ",x[j])),
               expression(paste(x[j]%+-%hat(sigma)[j])),
               expression(paste(x[j]%+-%sqrt(hat(sigma)[j]^2+hat(tau)^2)))
             ),
             lty=c(1,0,1,1),pch=c(-1,19,-1,-1),col=c(col2,col3,col4,col5),pt.cex=c(3,1.5,3,1),lwd=c(1,1,11,5),bty="n")


    }
    par(fig=c(0,1,0,.8), new=TRUE,lend=2)

    u.excluded=labsesbayes()[!sanitize]

    xl = x-c(sigma.x,u.excluded)
    xu = x+c(sigma.x,u.excluded)

    xl.Sigma.Tau.included = x[1:nI]-sqrt( sigma.x^2 +  tau.x^2)
    xu.Sigma.Tau.included = x[1:nI]+sqrt( sigma.x^2 +  tau.x^2)

    xl.Sigma.Tau.excluded = x[(nI+1):n]-sqrt( u.excluded^2 +  tau.x^2)
    xu.Sigma.Tau.excluded = x[(nI+1):n]+sqrt( u.excluded^2 +  tau.x^2)

    xl.Sigma.Tau=c(xl.Sigma.Tau.included,xl.Sigma.Tau.excluded)
    xu.Sigma.Tau=c(xu.Sigma.Tau.included,xu.Sigma.Tau.excluded)

    plot(c(0.5,n+0.5), range(c(xl-tau.x, xu+tau.x)),
         axes=FALSE, xlab="",
         ylab=label, type="n", bty="n")
    axis(2)
    polygon(c(1-0.2, nI+0.2, nI+0.2, 1-0.2),
            c(mu.x-mu.u, mu.x-mu.u, mu.x+mu.u, mu.x+mu.u),
            col=col1, border=col1)
    abline(h=mu.x, col=col2)
    segments(1:n, xl, 1:n, xu, col=col4, lwd=11)
    segments(1:n, xl.Sigma.Tau, 1:n, xu.Sigma.Tau, col=col5, lwd=5)

    points(1:nI, x[1:nI], pch=19, cex=1.5, col=col3)
    if(nI<n){
      points((nI+1):n, x[(nI+1):n], pch=21, cex=1.5, col=col3)
    }

    mtext(lab[seq(1, n, 2)], side=1, at=seq(1, n, 2), line=0, cex=1)
    mtext(lab[seq(2, n, 2)], side=1, at=seq(2, n, 2), line=1.2, cex=1)

    }

  output$bayesplot <- renderPlot({bayesplot_func()})

  output$bayesplotLegendWarning<- renderText({
    if(length(lablabsbayes())==0){
      lab=paste("L", 1:length(labmeansbayes()), sep="")
    }else{
      lab=lablabsbayes()
    }
    sanitize = !startsWith( lab,"-")

    bayesres=outbayes()[["mcmcout"]]
    postsumstats=summary(bayesres)$statistics
    sigma.index = grep("sigma", dimnames(postsumstats)[[1]])

    if(sum(!sanitize)!=0 && length(sigma.index)!=0){
      "* For the labs that have been excluded, the bars are calculated using the reported standard uncertainties."
    }else{
      ""
    }

  })


  output$downloadBayesGraph <- downloadHandler(
    filename = function() { paste( 'Bayesian.pdf', sep='') },
    content = function(file) {
      pdf(file,height=9, width=8)
      bayesplot_func()
      dev.off()
    })


  #######################################
  #######################################
  ####### Bayesian DoE
  #######################################
  #######################################

  source("DoEUnilateralBayes.R")

  BayesDoEUni=eventReactive(input$dobayes,{

    ### We need to remove the case when excluded labs mess with DoE computation

    if(length(lablabsbayes())==0){
      lab=paste("L", 1:length(labmeansbayes()), sep="")
    }else{
      lab=lablabsbayes()
    }
    sanitize = !startsWith( lab,"-")
    nlab = length(labmeansbayes())
    #Number of labs taken into account for the consensus
    nI = length(labmeansbayes()[sanitize])

    LOObool= as.logical(input$LOODoE)

    validate(
      need((nI>=3&&LOObool)||!LOObool, 'At least three labs are required to be included in the consensus for this method.'),
      need((nI>=2&&!LOObool)||LOObool, 'At least two labs are required to be included in the consensus for this method.')
    )

      withProgress(message = 'Calculating Unilateral DoEs', value = 0, style="old",{

      x=labmeansbayes()
      u=labsesbayes()
      if(length(labdfsbayes())!=0){
        nu=labdfsbayes()
      }else{
        nu=NULL
      }


      if(is.na(input$halfCauchyScaleTau)){
        input_pscaletau=mad(x)
      }else{
        input_pscaletau=input$halfCauchyScaleTau
      }

      if(is.na(input$halfCauchyScaleSigma)){
        input_pscalesig=median(u)
      }else{
        input_pscalesig=input$halfCauchyScaleSigma
      }

      DoEUnilateralBayes(x, u, nu, lab, LOO=as.logical(input$LOODoE),mcmc=outbayes()[["mcmcout"]],
                         ni = input$niters,nb = input$nburnin,nt = input$nthin,coverageProb = .95,#input$coverage) # hard-coded to match MRA definition and user's manual
                         UItauPriorScale=input_pscaletau, UIsigmaPriorScale=input_pscalesig)

    })


  })

  output$bayesDoEConv <- renderUI({
    HTML(BayesDoEUni()[["DoEwarn"]])
  })

  output$BayesDoEUniTable <- renderTable(BayesDoEUni()$DoE,
                                         digits=signifbayes,
                                         display=rep("fg",6),
                                         caption = '<b> Unilateral Degrees of Equivalence',
                                         caption.placement='top',
                                         include.rownames=FALSE)




  BayesDoEUniPlot_func=function(){

    par(fig=c(0,1,.6,1),lend=2)
    plot(c(1,1), c(1,1), type="n", axes=FALSE,xlab="", ylab="", bty="n")
    legend("top",c("DoE estimate" ," 95% coverage interval"),
           lty=c(NA,1),pch=c(19,NA),col=c(col3,col4),pt.cex=c(1,1),lwd=c(5,5))

    par(fig=c(0,1,0,.8), new=TRUE,lend=2)

    plot(0,xlim=range(0.7,length(BayesDoEUni()$DoE$DoE.x)+0.3),ylim=range(BayesDoEUni()$DoE$DoE.Lwr,BayesDoEUni()$DoE$DoE.Upr),type="n",     ylab="DoE",xlab="",xaxt="n", bty="n",axes=FALSE)
    abline(h=0,col="gray")

    if(length(lablabsbayes())==0){
      lab=paste("L", 1:length(labmeansbayes()), sep="")
    }else{
      lab=lablabsbayes()
    }
    #we remove the values of labs starting with "-"
    sanitize = !startsWith( lab,"-")
    nI=length(labmeansbayes()[sanitize])
    n=length(labmeansbayes())


    lab=BayesDoEUni()$DoE$Lab # order of the lab with excluded at the end


    n = length(labmeansbayes())
    mtext(lab[seq(1, n, 2)], side=1, at=seq(1, n, 2), line=0, cex=1, col="Black")
    mtext(lab[seq(2, n, 2)], side=1, at=seq(2, n, 2), line=1.2, cex=1, col="Black")


    axis(2)

    for(i in 1:length(BayesDoEUni()$DoE$DoE.x)){
      arrows(i,BayesDoEUni()$DoE$DoE.Lwr[i],i,BayesDoEUni()$DoE$DoE.Upr[i],length=0,col=col4, lwd=5)
    }
    points(x=1:length(BayesDoEUni()$DoE$DoE.x),y=BayesDoEUni()$DoE$DoE.x,pch=c(rep(19,nI),rep(21,n-nI)), cex=1.5, col=col3)


  }


  output$BayesDoEUniPlot<- renderPlot({BayesDoEUniPlot_func()})



  source("DoEBilateralBayes.R")


  BayesDoEBi=eventReactive(input$dobayes,{
    if(length(lablabsbayes())==0){
      lab=paste("L", 1:length(labmeansbayes()), sep="")
    }else{
      lab=lablabsbayes()
    }

    sanitize = !startsWith( lab,"-")

    #we reorganize the lab to remove those who starts with -

    lab = lab[sanitize]
    x = labmeansbayes()[sanitize]

    nI = length(x)

    ###Check number of included labs

    LOObool= as.logical(input$LOODoE)

    validate(
      need((nI>=3&&LOObool)||!LOObool, 'At least three labs are required to be included in the consensus for this method.'),
      need((nI>=2&&!LOObool)||LOObool, 'At least two labs are required to be included in the consensus for this method.')
    )

    withProgress(message = 'Calculating Bilateral DoEs', value = 0, style="old",{

      x=labmeansbayes()
      u=labsesbayes()
      if(length(labdfsbayes())!=0){
        nu=labdfsbayes()
      }else{
        nu=NULL
      }

      if(length(lablabsbayes())==0){
        lab=paste("L", 1:length(labmeansbayes()), sep="")
      }else{
        lab=lablabsbayes()
      }

      DoEBilateralBayes(x, u, nu,DoEUnilateral=BayesDoEUni(),coverageProb=.95)

    })

  })

  output$BayesDoEBiTable <- renderTable(BayesDoEBi()$B.x,
                                        digits = signifbayes,
                                        rownames = T,
                                        caption = '<b> Bilateral Degrees of Equivalence (Estimates)',
                                        caption.placement='top',
                                        sanitize.rownames.function =  bold.allrows)

  output$BayesDoE.UBiTable <- renderTable(BayesDoEBi()$B.U,
                                          digits = signifbayes,
                                          rownames = T,
                                          caption = '<b> Bilateral Degrees of Equivalence (Expanded Uncertainties U95)',
                                          caption.placement='top',
                                          sanitize.rownames.function =  bold.allrows)


  BayesDoEBiPlot_func=function(){

    DoE.Bilateral.Bayes.x = BayesDoEBi()$B.x
    DoE.Bilateral.Bayes.Lwr = BayesDoEBi()$B.Lwr
    DoE.Bilateral.Bayes.Upr = BayesDoEBi()$B.Upr

    n = dim(DoE.Bilateral.Bayes.x)[[1]]

    labNAMEs = row.names(DoE.Bilateral.Bayes.x)

    xx = (DoE.Bilateral.Bayes.Lwr * DoE.Bilateral.Bayes.Upr > 0)
    ij = which(xx, arr.ind=TRUE, useNames = TRUE)
    xx[xx == 1] = 0.9
    xx[xx == 0] = 0.1
    diag(xx) = rep(0.5, n)
    image(1:n, 1:n, xx, xlim=c(-1.5, n+1),breaks=c(-0.5, 0.4, 0.85, 1.5),
          col=c(col5,"grey36",col4),
          axes=FALSE, xlab="", ylab="", asp=1)
    K = nrow(ij)
    if (K > 0) {
      for (k in 1:K) {
        points(ij[k,1], ij[k,2], pch=8, font=3, col="Black") } }
    polygon(c(0.5,n+0.5,n+0.5,0.5), c(0.5,0.5,n+0.5,n+0.5),
            border="Black", lwd=2)

    segments(seq(1.5,n+0.5),0.5,seq(1.5,n+0.5),n+0.5)
    segments(n+0.5,seq(1.5,n+0.5),0.5,seq(1.5,n+0.5))

    mtext(labNAMEs, side=1, at=1:n, las=2, line=1.5, col="Blue", font=2,padj=0)
    text(rep(0, n), 1:n, labNAMEs, adj=1, col="Blue", font=2)

    }

  output$BayesDoEBiPlot<- renderPlot({BayesDoEBiPlot_func()})

  output$downloadBayesUniDoEGraph <- downloadHandler(
    filename = function() { paste( 'Bayesian_Unilateral_DoE.pdf', sep='') },
    content = function(file) {
      pdf(file,height=8, width=8)
      BayesDoEUniPlot_func()
      dev.off()
    })
  output$downloadBayesBiDoEGraph <- downloadHandler(
    filename = function() { paste( 'Bayesian_Bilateral_DoE.pdf', sep='') },
    content = function(file) {
      pdf(file,height=8, width=8)
      BayesDoEBiPlot_func()
      dev.off()
    })

  #######################################
  #######################################
  ########### bayesian Laplace
  #######################################
  #######################################
  source("bayesLaplace.R")

  labmeansbayes_lap<-eventReactive(input$dobayes_lap,{
    (as.numeric(unlist(strsplit(gsub("\\s", "", input$mean),","))))
  }
  )
  labsesbayes_lap<-eventReactive(input$dobayes_lap,{
    (as.numeric(unlist(strsplit(gsub("\\s", "", input$se),","))))
  }
  )
  labdfsbayes_lap<-eventReactive(input$dobayes_lap,{
    (as.numeric(unlist(strsplit(gsub("\\s", "", input$df),","))))
  }
  )

  lablabsbayes_lap<-eventReactive(input$dobayes_lap,{
    trimws( strsplit(input$lablabels,",")[[1]])
  }
  )
  seedbayes_lap<-eventReactive(input$dobayes_lap,{
    as.numeric(input$seed)
  })
  signifbayes_lap<-eventReactive(input$dobayes_lap,{
    as.numeric(input$signifDigits)
  })



  outbayes_lap<- eventReactive(input$dobayes_lap,{


    #################################
    ####### check inputs #########
    #################################
    validate(
      need(!any(is.na(labmeansbayes_lap())), 'Measured values must be numeric, separated by commas.'),
      need(!any(is.na(labsesbayes_lap()))&& !any(labsesbayes_lap()<=0), 'Standard uncertainties must be numeric (greater than zero), separated by commas.'),
      need(!any(is.na(labdfsbayes_lap())) && !any(labdfsbayes_lap()<1), 'Degrees of freedom must be numeric (greater than or equal to 1), separated by commas.'),
      need(is.numeric(input$coverage) && input$coverage<=1 && input$coverage>=0, 'Invalid coverage probability, must be numeric between 0 and 1.'),
      need(length(labmeansbayes_lap())==length(labsesbayes_lap()),  'The number of measured values entered does not equal the number of standard uncertainties entered.'),
      need(length(labmeansbayes_lap())==length(labdfsbayes_lap()) || length(labdfsbayes_lap())==0, "The number of degrees of freedom entered does not equal the number of measured values entered."),

      need(length(labmeansbayes_lap())==length(lablabsbayes_lap()) | length(lablabsbayes_lap())==0,'The number of study labels entered does not equal the number of measured values entered.'),
      need(!any(table(lablabsbayes_lap())>1) | length(lablabsbayes_lap())==0,"Laboratory labels should be unique."),

      need(input$niters_lap>0 && input$niters_lap<=1000000,'Total number of iterations must be a positive integer (less than 1000000).'),
      need(input$nburnin_lap>=0 && input$nburnin_lap<=500000,'Length of burn in must be a positive integer (less than 500000).'),
      need(input$nthin_lap>0 && input$nthin_lap<=1000,'Thinning rate must be a positive integer (less than 1000).'),
      need(input$nburnin_lap<=input$niters_lap/2,paste("Length of burn in too large (over half of the total number of iterations). Increase number of iterations to at least ",input$nburnin*2,".",sep="")),
      need(length(labmeansbayes_lap())<=500, 'The NICOB will not combine more than 500 measured values.'),
      need(!is.na(seedbayes_lap()), 'Random number generator seed must be numeric.'),
      ###Number of labs
      need(length(labmeansbayes_lap())>1, 'At least two measured values are required.')


    )

    #################################
    set.seed(seedbayes_lap())

    if(length(lablabsbayes_lap())==0){
      lab=paste("L", 1:length(labmeansbayes_lap()), sep="")
    }else{
      lab=lablabsbayes_lap()
    }

    sanitize = !startsWith( lab,"-")

    #we reorganize the lab to remove those who starts with -

    lab = lab[sanitize]
    x = labmeansbayes_lap()[sanitize]
    u = labsesbayes_lap()[sanitize]

    nI = length(x)

    ###Check number of included labs

    LOObool= as.logical(input$LOODoE)

    validate(
      need((nI>=3&&LOObool)||!LOObool, 'At least three labs are required to be included in the consensus for this method.'),
      need((nI>=2&&!LOObool)||LOObool, 'At least two labs are required to be included in the consensus for this method.')
    )


    if(length(labdfsbayes_lap())!=0){
      nu=labdfsbayes_lap()[sanitize]
    }else{
      nu=NULL
    }

    if(is.na(input$halfCauchyScaleTau_lap)){
      input_pscaletau=mad(x)#median(abs(labmeansbayes()-median(labmeansbayes()))) #tauPriorScale=mad(x)
      priorTauWarn_lap=paste("Scale for half-Cauchy prior on between laboratory variance has either not
                      been entered or is not numeric. Instead, the median absolute deviation of the
                         measured values has been used: ",input_pscaletau_lap,sep="")
    }else{
      input_pscaletau_lap=input$halfCauchyScaleTau_lap
      priorTauWarn_lap=""
    }
    output$bayesTauWarn_lap <- renderUI({
      HTML(priorTauWarn_lap)
    })

    if(is.na(input$halfCauchyScaleSigma_lap)){
      input_pscalesig_lap=median(u)
      priorSigWarn_lap=paste("Scale for half-Cauchy prior on within laboratory variances has either not
                      been entered or is not numeric. Instead, the median of the lab-specific standard
                      uncertainties has been used: ",input_pscalesig_lap,sep="")
    }else{
      input_pscalesig_lap=input$halfCauchyScaleSigma_lap
      priorSigWarn_lap=""
    }

    output$bayesSigWarn_lap <- renderUI({
      HTML(priorSigWarn_lap)
    })


    validate(
      need(input_pscaletau_lap>0, 'Scale for half-Cauchy prior on between laboratory variance must be a positive number.'),
      need(input_pscalesig_lap>0, 'Scale for half-Cauchy prior on within laboratory variances must be a positive number.')
    )

    bayesLaplace(x,u,nu,tauPriorScale = input_pscaletau_lap,sigmaPriorScale = input_pscalesig_lap,ni = input$niters_lap,nb = input$nburnin_lap,nt = input$nthin_lap)
  })


  output$downloadbayesout <- downloadHandler(
    filename = "MCMCout.csv",
    content = function(file) {
      write.table(outbayes_lap()[["mcmcout"]][[1]], file,row.names = F)
    }
  )

  output$bayesConv_lap <- renderUI({
    HTML(outbayes_lap()[["warn"]])
  })


  output$bayesest_lap <- renderText({
    bayesres=outbayes_lap()[["mcmcout"]]
    signifDigits = significator(bayesres[[1]][,"mu"],digits=TRUE)
    if (signifDigits==1)
    {
      signifString= paste(" (where 1 significant digit is believed to be reliable)")
    }else{
      signifString= paste(" (where", signifDigits, "significant digits are believed to be reliable)")
    }
    paste(signif(mean(bayesres[[1]][,"mu"]),input$signifDigits), signifString)

  })

  output$bayesse_lap <- renderText({
    bayesres=outbayes_lap()[["mcmcout"]]
    paste(signif(sd(bayesres[[1]][,"mu"]),input$signifDigits))

  })

  output$bayesquant_lap <- renderText({
    bayesres=outbayes_lap()[["mcmcout"]]
    credible.interval=quantile(bayesres[[1]][,"mu"],c((1-input$coverage)/2, (1+input$coverage)/2))

    paste(signif(credible.interval[1],digits = input$signifDigits),
          "to",
          signif(credible.interval[2],digits = input$signifDigits))

  })

  output$bayestau_lap<- renderText({
    bayesres=outbayes_lap()[["mcmcout"]]
    postsumstats=summary(bayesres)$statistics

    bayesres=outbayes_lap()[["tau"]]
    paste(signif(postsumstats["tau", "Mean"],input$signifDigits))

  })

  bayesplot_lap_func=function(){
    bayesres=outbayes_lap()[["mcmcout"]]

    postsumstats=summary(bayesres)$statistics
    if(length(lablabsbayes_lap())==0){
      lab=paste("L", 1:length(labmeansbayes_lap()), sep="")
    }else{
      lab=lablabsbayes_lap()
    }

    sanitize = !startsWith( lab,"-")

    #we reorganize the lab to put those with a minus at the end

    lab = c(lab[sanitize],lab[!sanitize])
    x = c(labmeansbayes_lap()[sanitize],labmeansbayes_lap()[!sanitize])
    u = c(labsesbayes_lap()[sanitize],labsesbayes_lap()[!sanitize])

    n = length(x)
    #Number of labs taken into account for the consensus
    nI = length(labmeansbayes_lap()[sanitize])


    postsumstats=summary(bayesres)$statistics
    postparameters=rownames(postsumstats)

    mu.x = postsumstats["mu","Mean"]
    mu.u = postsumstats["mu", "SD"]
    tau.x = postsumstats["tau", "Mean"]

    sigma.index = grep("sigma", dimnames(postsumstats)[[1]])

    sigma.names=rep(NA,nI)
    for(i in 1:nI){
      sigma.names[i] = paste("sigma[",i, "]",sep="")
    }

    if(length(sigma.index)==0){
      sigma.x = labsesbayes_lap()[sanitize]
    }else{
      sigma.x = postsumstats[sigma.names, "Mean"]
    }



    if(input$units==""){
      label="Measured Value"
    }else{
      label=paste("Measured Value (",input$units,")",sep="")
    }


    par(fig=c(0,1,.6,1),lend=2)

    plot(c(1,1), c(1,1), type="n", axes=FALSE,xlab="", ylab="", bty="n")
    if(length(sigma.index)==0){
      leg=legend("top",c(expression(paste("Consensus estimate, ", hat(mu),", and interval ",hat(mu)%+-%u(mu),sep="")),
                         expression(paste("Measured value, ",x[j])),
                         expression(paste(x[j]%+-%u[j])),
                         expression(paste(x[j]%+-%sqrt(u[j]^2+hat(tau)^2)))
      ),
      lty=c(0,0,1,1),pch=c(15,19,-1,-1),col=c(col1,col3,col4,col5),pt.cex=c(3,1.5,3,1),lwd=c(1,1,11,5))

      legend("top",
             c(expression(paste("Consensus estimate, ", hat(mu),", and interval ",hat(mu)%+-%u(mu),sep="")),
               expression(paste("Measured value, ",x[j])),
               expression(paste(x[j]%+-%u[j])),
               expression(paste(x[j]%+-%sqrt(u[j]^2+hat(tau)^2)))
             ),
             lty=c(1,0,1,1),pch=c(-1,19,-1,-1),col=c(col2,col3,col4,col5),pt.cex=c(3,1.5,3,1),lwd=c(1,1,11,5),bty="n")


    }else{
      leg=legend("top",c(expression(paste("Consensus estimate, ", hat(mu),", and interval ",hat(mu)%+-%u(mu),sep="")),
                         expression(paste("Measured value, ",x[j])),
                         expression(paste(x[j]%+-%hat(sigma)[j])),
                         expression(paste(x[j]%+-%sqrt(hat(sigma)[j]^2+hat(tau)^2)))
      ),
      lty=c(0,0,1,1),pch=c(15,19,-1,-1),col=c(col1,col3,col4,col5),pt.cex=c(3,1,3,1),lwd=c(1,1,11,5))

      legend("top",
             c(expression(paste("Consensus estimate, ", hat(mu),", and interval ",hat(mu)%+-%u(mu),sep="")),
               expression(paste("Measured value, ",x[j])),
               expression(paste(x[j]%+-%hat(sigma)[j])),
               expression(paste(x[j]%+-%sqrt(hat(sigma)[j]^2+hat(tau)^2)))
             ),
             lty=c(1,0,1,1),pch=c(-1,19,-1,-1),col=c(col2,col3,col4,col5),pt.cex=c(1,1,3,1),lwd=c(1,1,11,5),bty="n")


    }
    par(fig=c(0,1,0,.8), new=TRUE,lend=2)

    u.excluded=labsesbayes_lap()[!sanitize]

    xl = x-c(sigma.x,u.excluded)
    xu = x+c(sigma.x,u.excluded)

    xl.Sigma.Tau.included = x[1:nI]-sqrt( sigma.x^2 +  tau.x^2)
    xu.Sigma.Tau.included = x[1:nI]+sqrt( sigma.x^2 +  tau.x^2)

    xl.Sigma.Tau.excluded = x[(nI+1):n]-sqrt( u.excluded^2 +  tau.x^2)
    xu.Sigma.Tau.excluded = x[(nI+1):n]+sqrt( u.excluded^2 +  tau.x^2)

    xl.Sigma.Tau=c(xl.Sigma.Tau.included,xl.Sigma.Tau.excluded)
    xu.Sigma.Tau=c(xu.Sigma.Tau.included,xu.Sigma.Tau.excluded)

    plot(c(0.5,n+0.5), range(c(xl-tau.x, xu+tau.x)),
         axes=FALSE, xlab="",
         ylab=label, type="n", bty="n")
    axis(2)
    polygon(c(1-0.2, nI+0.2, nI+0.2, 1-0.2),
            c(mu.x-mu.u, mu.x-mu.u, mu.x+mu.u, mu.x+mu.u),
            col=col1, border=col1)
    abline(h=mu.x, col=col2)
    segments(1:n, xl, 1:n, xu, col=col4, lwd=11)
    segments(1:n, xl.Sigma.Tau, 1:n, xu.Sigma.Tau, col=col5, lwd=5)

    points(1:nI, x[1:nI], pch=19, cex=1.5, col=col3)
    if(nI<n){
      points((nI+1):n, x[(nI+1):n], pch=21, cex=1.5, col=col3)
    }

    mtext(lab[seq(1, n, 2)], side=1, at=seq(1, n, 2), line=0, cex=1)
    mtext(lab[seq(2, n, 2)], side=1, at=seq(2, n, 2), line=1.2, cex=1)

  }

  output$bayesplot_lap <- renderPlot({bayesplot_lap_func()})

  output$bayesplotLegendWarning_lap<- renderText({
    if(length(lablabsbayes_lap())==0){
      lab=paste("L", 1:length(labmeansbayes_lap()), sep="")
    }else{
      lab=lablabsbayes_lap()
    }
    sanitize = !startsWith( lab,"-")

    bayesres=outbayes_lap()[["mcmcout"]]
    postsumstats=summary(bayesres)$statistics
    sigma.index = grep("sigma", dimnames(postsumstats)[[1]])

    if(sum(!sanitize)!=0 && length(sigma.index)!=0){
      "* For the labs that have been excluded, the bars are calculated using the reported standard uncertainties."
    }else{
      ""
    }


  })

  output$downloadBayesGraph_lap <- downloadHandler(
    filename = function() { paste( 'Bayesian_lap.pdf', sep='') },
    content = function(file) {
      pdf(file,height=9, width=8)
      bayesplot_lap_func()
      dev.off()
    })

    #######################################
    #######################################
    ####### Bayesian Laplace DoE
    #######################################
    #######################################

    source("DoEUnilateralBayesLaplace.R")

    BayesDoEUni_lap=eventReactive(input$dobayes_lap,{

      ### We need to remove the case when excluded labs mess with DoE computation

      if(length(lablabsbayes_lap())==0){
        lab=paste("L", 1:length(labmeansbayes_lap()), sep="")
      }else{
        lab=lablabsbayes_lap()
      }
      sanitize = !startsWith( lab,"-")
      nlab = length(labmeansbayes_lap())
      #Number of labs taken into account for the consensus
      nI = length(labmeansbayes_lap()[sanitize])

      LOObool= as.logical(input$LOODoE)

      validate(
        need((nI>=3&&LOObool)||!LOObool, 'At least three labs are required to be included in the consensus for this method.'),
        need((nI>=2&&!LOObool)||LOObool, 'At least two labs are required to be included in the consensus for this method.')
      )

        withProgress(message = 'Calculating Unilateral DoEs', value = 0, style="old",{

        x=labmeansbayes_lap()
        u=labsesbayes_lap()
        if(length(labdfsbayes_lap())!=0){
          nu=labdfsbayes_lap()
        }else{
          nu=NULL
        }


        if(is.na(input$halfCauchyScaleTau_lap)){
          input_pscaletau=mad(x)
        }else{
          input_pscaletau=input$halfCauchyScaleTau_lap
        }

        if(is.na(input$halfCauchyScaleSigma_lap)){
          input_pscalesig=median(u)
        }else{
          input_pscalesig=input$halfCauchyScaleSigma_lap
        }

        DoEUnilateralBayesLaplace(x, u, nu, lab, LOO=as.logical(input$LOODoE),mcmc=outbayes_lap()[["mcmcout"]],
                           ni = input$niters_lap,nb = input$nburnin_lap,nt = input$nthin_lap,coverageProb = .95,#input$coverage) # hard-coded to match MRA definition and user's manual
                           UItauPriorScale=input_pscaletau, UIsigmaPriorScale=input_pscalesig)

      })


    })

    output$bayesDoEConv_lap <- renderUI({
      HTML(BayesDoEUni_lap()[["DoEwarn"]])
    })

    output$BayesDoEUniTable_lap <- renderTable(BayesDoEUni_lap()$DoE,
                                           digits=signifbayes_lap,
                                           display=rep("fg",6),
                                           caption = '<b> Unilateral Degrees of Equivalence',
                                           caption.placement='top',
                                           include.rownames=FALSE)




    BayesDoEUniPlot_lap_func=function(){

      par(fig=c(0,1,.6,1),lend=2)
      plot(c(1,1), c(1,1), type="n", axes=FALSE,xlab="", ylab="", bty="n")
      legend("top",c("DoE estimate" ," 95% coverage interval"),
             lty=c(NA,1),pch=c(19,NA),col=c(col3,col4),pt.cex=c(1.5,1.5),lwd=c(5,5))



      par(fig=c(0,1,0,.8), new=TRUE,lend=2)

      plot(0,xlim=range(0.7,length(BayesDoEUni_lap()$DoE$DoE.x)+0.3),ylim=range(BayesDoEUni_lap()$DoE$DoE.Lwr,BayesDoEUni_lap()$DoE$DoE.Upr),type="n",     ylab="DoE",xlab="",xaxt="n", bty="n",axes=FALSE)
      abline(h=0,col="gray")

      if(length(lablabsbayes_lap())==0){
        lab=paste("L", 1:length(labmeansbayes_lap()), sep="")
      }else{
        lab=lablabsbayes_lap()
      }
      #we remove the values of labs starting with "-"
      sanitize = !startsWith( lab,"-")
      nI=length(labmeansbayes_lap()[sanitize])
      n=length(labmeansbayes_lap())


      lab=BayesDoEUni_lap()$DoE$Lab # order of the lab with excluded at the end


      n = length(labmeansbayes_lap())
      mtext(lab[seq(1, n, 2)], side=1, at=seq(1, n, 2), line=0, cex=1, col="Black")
      mtext(lab[seq(2, n, 2)], side=1, at=seq(2, n, 2), line=1.2, cex=1, col="Black")


      axis(2)

      for(i in 1:length(BayesDoEUni_lap()$DoE$DoE.x)){
        arrows(i,BayesDoEUni_lap()$DoE$DoE.Lwr[i],i,BayesDoEUni_lap()$DoE$DoE.Upr[i],length=0, col=col4, lwd=5)
      }
      points(x=1:length(BayesDoEUni_lap()$DoE$DoE.x),y=BayesDoEUni_lap()$DoE$DoE.x,pch=c(rep(19,nI),rep(21,n-nI)), cex=1.5, col=col3)


    }


    output$BayesDoEUniPlot_lap<- renderPlot({BayesDoEUniPlot_lap_func()})



    source("DoEBilateralBayes.R")


    BayesDoEBi_lap=eventReactive(input$dobayes_lap,{
      if(length(lablabsbayes_lap())==0){
        lab=paste("L", 1:length(labmeansbayes_lap()), sep="")
      }else{
        lab=lablabsbayes_lap()
      }

      sanitize = !startsWith( lab,"-")

      #we reorganize the lab to remove those who starts with -

      lab = lab[sanitize]
      x = labmeansbayes_lap()[sanitize]

      nI = length(x)

      ###Check number of included labs

      LOObool= as.logical(input$LOODoE)

      validate(
        need((nI>=3&&LOObool)||!LOObool, 'At least three labs are required to be included in the consensus for this method.'),
        need((nI>=2&&!LOObool)||LOObool, 'At least two labs are required to be included in the consensus for this method.')
      )

      withProgress(message = 'Calculating Bilateral DoEs', value = 0, style="old",{

        x=labmeansbayes_lap()
        u=labsesbayes_lap()
        if(length(labdfsbayes_lap())!=0){
          nu=labdfsbayes_lap()
        }else{
          nu=NULL
        }

        if(length(lablabsbayes_lap())==0){
          lab=paste("L", 1:length(labmeansbayes_lap()), sep="")
        }else{
          lab=lablabsbayes_lap()
        }

        DoEBilateralBayes(x, u, nu,DoEUnilateral=BayesDoEUni_lap(),coverageProb=.95)

      })

    })

    output$BayesDoEBiTable_lap <- renderTable(BayesDoEBi_lap()$B.x,
                                          digits = signifbayes_lap,
                                          rownames = T,
                                          caption = '<b> Bilateral Degrees of Equivalence (Estimates)',
                                          caption.placement='top',
                                          sanitize.rownames.function =  bold.allrows)

    output$BayesDoE.UBiTable_lap <- renderTable(BayesDoEBi_lap()$B.U,
                                            digits = signifbayes_lap,
                                            rownames = T,
                                            caption = '<b> Bilateral Degrees of Equivalence (Expanded Uncertainties U95)',
                                            caption.placement='top',
                                            sanitize.rownames.function =  bold.allrows)


    BayesDoEBiPlot_lap_func=function(){

      DoE.Bilateral.Bayes.x = BayesDoEBi_lap()$B.x
      DoE.Bilateral.Bayes.Lwr = BayesDoEBi_lap()$B.Lwr
      DoE.Bilateral.Bayes.Upr = BayesDoEBi_lap()$B.Upr

      n = dim(DoE.Bilateral.Bayes.x)[[1]]

      labNAMEs = row.names(DoE.Bilateral.Bayes.x)

      xx = (DoE.Bilateral.Bayes.Lwr * DoE.Bilateral.Bayes.Upr > 0)
      ij = which(xx, arr.ind=TRUE, useNames = TRUE)
      xx[xx == 1] = 0.9
      xx[xx == 0] = 0.1
      diag(xx) = rep(0.5, n)
      image(1:n, 1:n, xx, xlim=c(-1.5, n+1),breaks=c(-0.5, 0.4, 0.85, 1.5),
            col=c(col5,"grey36",col4),
            axes=FALSE, xlab="", ylab="", asp=1)
      K = nrow(ij)
      if (K > 0) {
        for (k in 1:K) {
          points(ij[k,1], ij[k,2], pch=8, font=3, col="Black") } }
      polygon(c(0.5,n+0.5,n+0.5,0.5), c(0.5,0.5,n+0.5,n+0.5),
              border="Black", lwd=2)

      segments(seq(1.5,n+0.5),0.5,seq(1.5,n+0.5),n+0.5)
      segments(n+0.5,seq(1.5,n+0.5),0.5,seq(1.5,n+0.5))

      mtext(labNAMEs, side=1, at=1:n, las=2, line=1.5, col="Blue", font=2,padj=0)
      text(rep(0, n), 1:n, labNAMEs, adj=1, col="Blue", font=2)

      }

    output$BayesDoEBiPlot_lap<- renderPlot({BayesDoEBiPlot_lap_func()})

    output$downloadBayesUniDoEGraph_lap <- downloadHandler(
      filename = function() { paste( 'Bayesian_Laplace_Unilateral_DoE.pdf', sep='') },
      content = function(file) {
        pdf(file,height=8, width=8)
        BayesDoEUniPlot_lap_func()
        dev.off()
      })
    output$downloadBayesBiDoEGraph_lap <- downloadHandler(
      filename = function() { paste( 'Bayesian_Laplace_Bilateral_DoE.pdf', sep='') },
      content = function(file) {
        pdf(file,height=8, width=8)
        BayesDoEBiPlot_lap_func()
        dev.off()
      })

  ######################################################################
  ######################################################################
  ######################################################################
  ######################################################################
  ######################################################################
  ######################################################################
  ######################################################################
  ######## LINEAR POOLING ========================================
  ######################################################################
  # labmeanspool<-function(){
  #   return(as.numeric(unlist(strsplit(gsub("\\s", "", input$mean),","))))
  # }
  # labsespool<-function(){
  #   return(as.numeric(unlist(strsplit(gsub("\\s", "", input$se),","))))
  # }
  # labdfspool<-function(){
  #   return(as.numeric(unlist(strsplit(gsub("\\s", "", input$df),","))))
  # }
  # numericPoolWeights<-function(){
  #   return(as.numeric(unlist(strsplit(gsub("\\s", "", input$poolweights),","))))
  # }

  labmeanspool<-eventReactive(input$dopool,{
    (as.numeric(unlist(strsplit(gsub("\\s", "", input$mean),","))))
  })
  labsespool<-eventReactive(input$dopool,{
    (as.numeric(unlist(strsplit(gsub("\\s", "", input$se),","))))
  })
  labdfspool<-eventReactive(input$dopool,{
    (as.numeric(unlist(strsplit(gsub("\\s", "", input$df),","))))
  })
  numericPoolWeights<-eventReactive(input$dopool,{
    (as.numeric(unlist(strsplit(gsub("\\s", "", input$poolweights),","))))
  })
  lablabspool<-eventReactive(input$dopool,{
    trimws( strsplit(input$lablabels,",")[[1]])
  })
  seedpool<-eventReactive(input$dopool,{
        as.numeric(input$seed)
  })
  signifpool<-eventReactive(input$dopool,{
    as.numeric(input$signifDigits)
  })



  source("linearOP.R")

  outpool<- eventReactive(input$dopool,{

    #################################
    ####### check inputs #########
    #################################
    validate(
      need(!any(is.na(labmeanspool())), 'Measured values must be numeric, separated by commas.'),
      need(!any(is.na(labsespool()))&& !any(labsespool()<=0), 'Standard uncertainties must be numeric (greater than zero), separated by commas.'),
      need(!any(is.na(labdfspool())) && !any(labdfspool()<1), 'Degrees of freedom must be numeric (greater than or equal to 1), separated by commas.'),
      need(is.numeric(input$coverage) && input$coverage<=1 && input$coverage>=0, 'Invalid coverage probability, must be numeric between 0 and 1.'),
      need(!any(is.na(numericPoolWeights())) && !any(numericPoolWeights()<0), 'Weights must be numeric (greater than zero), separated by commas.'),
      need(sum(numericPoolWeights())>0 || length(numericPoolWeights())==0 , 'At least one of the weights should be greater than zero.'),
      need(is.numeric(input$linearOPrep), 'Sample size must be numeric.'),
      need(input$linearOPrep>0 && input$linearOPrep<=1000000, 'Sample size must be greater than zero and less than 1000000.'),
      need(length(labmeanspool())==length(numericPoolWeights()) || length(numericPoolWeights())==0, "The number of laboratory measured values does not equal the number of weights."),
      need(length(labmeanspool())==length(labsespool()),'The number of measured values entered does not equal the number of standard uncertainties entered.'),
      need(length(labmeanspool())==length(labdfspool()) || length(labdfspool())==0, "The number of degrees of freedom entered does not equal the number of measured values entered."),

      # need(length(labmeanspool())==length(strsplit(input$lablabels,",")[[1]]) | length(strsplit(input$lablabels,",")[[1]])==0, 'The number of study labels entered does not equal the number of measured values entered.'),
      need(length(labmeanspool())==length(lablabspool()) | length(lablabspool())==0,'The number of study labels entered does not equal the number of measured values entered.'),
      need(!any(table(lablabspool())>1) | length(lablabspool())==0,"Laboratory labels should be unique."),

      need(length(labmeanspool())>1, 'At least two measured values are required.'),
      need(length(labmeanspool())<=500, 'The NICOB will not combine more than 500 measured values.'),
      need(!is.na(seedpool()), 'Random number generator seed must be numeric.')

    )

    set.seed(seedpool() )


    if(length(lablabspool())==0){
      lab=paste("L", 1:length(labmeanspool()), sep="")
    }else{
      lab=lablabspool()
    }
    #we remove the values of labs starting with "-"
    sanitize = !startsWith( lab,"-")


    #we reorganize the lab to remove those who starts with -

    lab = lab[sanitize]
    x = labmeanspool()[sanitize]
    nI = length(x)

    ###Check number of included labs

    LOObool= as.logical(input$LOODoE)

    validate(
      need((nI>=3&&LOObool)||!LOObool, 'At least three labs are required to be included in the consensus for this method.'),
      need((nI>=2&&!LOObool)||LOObool, 'At least two labs are required to be included in the consensus for this method.')
    )
    if(length(numericPoolWeights())==0){
      LPweights=rep(1,length(labmeanspool()[sanitize]))
    }else{
      LPweights=numericPoolWeights()
    }

    linearOP(labmeanspool()[sanitize],labsespool()[sanitize],labdfspool()[sanitize],LPweights,m=input$linearOPrep)

  })


  output$poolest <- renderText({
    poolres=outpool()
    signifDigits = significator(poolres,digits=TRUE)
    if (signifDigits==1)
    {
      signifString= paste(" (where 1 significant digit is believed to be reliable)")
    }else{
      signifString= paste(" (where", signifDigits, "significant digits are believed to be reliable)")
    }
    paste(signif(mean(poolres),digits = input$signifDigits), signifString)

  })

  output$poolse <- renderText({
    poolres=outpool()
    paste(signif(sd(poolres),digits = input$signifDigits))
  })

  output$poolquant <- renderText({
    poolres=outpool()
    credible.interval=quantile(poolres,c((1-input$coverage)/2, (1+input$coverage)/2))

    paste(signif(credible.interval[1],digits = input$signifDigits),
          "to",
          signif(credible.interval[2],digits = input$signifDigits))
  })

  LPplot = function()
  {
  xVALUEs=outpool()

  z = density(xVALUEs,n=4096)

  L = quantile(xVALUEs, probs=(1-input$coverage)/2)
  U = quantile(xVALUEs, probs=(1+input$coverage)/2)
  iS = which.min(abs(z$x-L))
  iE = which.min(abs(z$x-U))

  ## NOTE: Need to pre-compute yRANGE so that none of the curves
  ## have their peaks clipped
  if (length(labdfspool())!=0){

    yRANGEmax=rep(NA,length(labdfspool()))
    for(j in 1:length(labdfspool())){
      if(labdfspool()[j]==Inf){
        yRANGEmax[j] = 1/(labsespool()[j]*sqrt(2*pi))
      }else if(labdfspool()[j]>2){
        # yRANGEmax[j] = gamma((labdfspool()[j]+1)/2)/gamma(labdfspool()[j]/2)/sqrt(labsespool()[j]^2*labdfspool()[j]*pi)*sqrt(labdfspool()[j]/(labdfspool()[j]-2))
        yRANGEmax[j] = exp(lgamma((labdfspool()[j]+1)/2)-lgamma(labdfspool()[j]/2))/sqrt(labsespool()[j]^2*labdfspool()[j]*pi)*sqrt(labdfspool()[j]/(labdfspool()[j]-2))
      }else{
        yRANGEmax[j] = gamma((labdfspool()[j]+1)/2)/gamma(labdfspool()[j]/2)/sqrt(labsespool()[j]^2*labdfspool()[j]*pi)
      }
    }
    yRANGE=c(0,max(yRANGEmax))


  }else{
    yRANGE = c(0, max(1/(labsespool()*sqrt(2*pi))))
  }

  xRANGE = c(quantile(xVALUEs, probs=0.00025),quantile(xVALUEs, probs=0.99975))


  if(input$units==""){
    label="Measured Value"
  }else{
    label=paste("Measured Value (",input$units,")",sep="")
  }

  if(length(lablabspool())==0){
    lab=paste("L", 1:length(labmeanspool()), sep="")
  }else{
    lab=lablabspool()
  }
  sanitize = !startsWith( lab,"-")

  plot(z$x, z$y, ylim=yRANGE,xlim=xRANGE,
       type="n", main="", bty="n",
       xlab=label, ylab="Probability Density")
  if (any(!sanitize)) # if any lab is excluded we add a legend for unuesd density
    legend("topright",
           c(expression(paste("Consensus estimate")),
             expression(paste("Used Density ")),
             expression(paste("Ignored Density")),
             expression(paste("Consensus Value")),
             expression(paste("95 % coeverage interval"))
           ),
           pch=c(-1,-1,-1,19,-1),col=c(col1,col4,col4,col3,col3),lty=c(1,1,2,-1,1),lwd=c(2,2,2,0,5),bty="o")
  else
    legend("topright",
           c(expression(paste("Consensus estimate")),
             expression(paste("Density ")),
             expression(paste("Consensus Value")),
             expression(paste("95 % coeverage interval"))
           ),
           pch=c(-1,-1,19,-1),col=c(col1,col4,col3,col3),lty=c(1,1,-1,1),lwd=c(2,2,0,5),bty="o")
  xp = c(L, U, z$x[iE:iS])
  yp = c(0, 0, z$y[iE:iS])
  polygon(xp, yp, border=F, col=col1)
  lines(z$x, z$y, type="l", col=col2, lwd=2)


  if (length(labdfspool())!=0){

    thenu=labdfspool()
    factorVar=numeric(length(thenu))
    indexVar = (thenu > 2)
    factorVar[indexVar] = labsespool()[indexVar]/
      sqrt(thenu[indexVar]/(thenu[indexVar]-2))
    factorVar[!indexVar] = labsespool()[!indexVar]

  }




  for (j in 1:length(labmeanspool())) {
    xx = seq(from=min(z$x), to=max(z$x), length=10000)

    if (length(labdfspool())!=0){

      if(labdfspool()[j]==Inf){
        yy = dnorm(xx, mean=labmeanspool()[j], sd=labsespool()[j])
      }else{
        yy=dt((xx-labmeanspool()[j])/factorVar[j],
              df=thenu[j])/
          factorVar[j]
      }
    }else{
      yy = dnorm(xx, mean=labmeanspool()[j], sd=labsespool()[j])
    }

    if(sanitize[j])
    {
      lines(xx, yy, col=col4, lwd=2)

    }else
    {
      lines(xx, yy, col=col4, lwd=2, lty=2)

    }

  }
  segments(L, 0, U, 0, lwd=5, col=col3)
  points(mean(xVALUEs),
         min(z$y)+0.4*par()$cxy[2], pch=19, col=col3)
  }



  output$poolplot <- renderPlot({LPplot()

  })

  output$downloadLPout <- downloadHandler(
    filename = "LPout.csv",
    content = function(file) {
      write.csv(outpool(), file)
    }
  )
  output$downloadLPgraph <- downloadHandler(
    filename = function() { paste( 'LinearPool.pdf', sep='') },
       content = function(file) {
         pdf(file,height=8, width=8)
         LPplot()
         dev.off()
  })



  #######################################
  #######################################
  ####### LP DoE
  #######################################
  #######################################

  source("DoEUnilateralLinearPool.R")
  source("DoEBilateralLinearPool.R")

  # LPcp<- eventReactive(input$dopool,{
  #   input$coverage
  #
  # })

  LPDoEUni=eventReactive(input$dopool,{

    validate(
      need(!any(is.na(labmeanspool())), 'Measured values must be numeric, separated by commas.'),
      need(!any(is.na(labsespool()))&& !any(labsespool()<=0), 'Standard uncertainties must be numeric (greater than zero), separated by commas.'),
      need(!any(is.na(labdfspool())) && !any(labdfspool()<1), 'Degrees of freedom must be numeric (greater than or equal to 1), separated by commas.'),
      need(is.numeric(input$coverage) && input$coverage<=1 && input$coverage>=0, 'Invalid coverage probability, must be numeric between 0 and 1.'),
      need(!any(is.na(numericPoolWeights())) && !any(numericPoolWeights()<0), 'Weights must be numeric (greater than zero), separated by commas.'),
      need(sum(numericPoolWeights())>0 || length(numericPoolWeights())==0, 'At least one of the weights should be greater than zero.'),
      need(is.numeric(input$linearOPrep), 'Sample size must be numeric.'),
      need(input$linearOPrep>0 && input$linearOPrep<=1000000, 'Sample size must be greater than zero and less than 1000000.'),
      need(length(labmeanspool())==length(numericPoolWeights()) || length(numericPoolWeights())==0, "The number of laboratory measured values does not equal the number of weights."),
      need(length(labmeanspool())==length(labsespool()),'The number of measured values entered does not equal the number of standard uncertainties entered.'),
      need(length(labmeanspool())==length(labdfspool()) || length(labdfspool())==0, "The number of degrees of freedom entered does not equal the number of measured values entered."),

      # need(length(labmeanspool())==length(strsplit(input$lablabels,",")[[1]]) | length(strsplit(input$lablabels,",")[[1]])==0, 'The number of study labels entered does not equal the number of measured values entered.'),
      need(length(labmeanspool())==length(lablabspool()) | length(lablabspool())==0,'The number of study labels entered does not equal the number of measured values entered.'),
      need(!any(table(lablabspool())>1) | length(lablabspool())==0,"Laboratory labels should be unique."),

      need(length(labmeanspool())>1, 'At least two measured values are required.'),
      need(length(labmeanspool())<=500, 'The NICOB will not combine more than 500 measured values.')

    )

    if(length(lablabspool())==0){
      lab=paste("L", 1:length(labmeanspool()), sep="")
    }else{
      lab=lablabspool()
    }

    #we remove the values of labs starting with "-"
    sanitize = !startsWith( lab,"-")
    nI = length(labmeanspool()[sanitize])

    ###Check number of included labs

    LOObool= as.logical(input$LOODoE)

    validate(
      need((nI>=3&&LOObool)||!LOObool, 'At least three labs are required to be included in the consensus for this method.'),
      need((nI>=2&&!LOObool)||LOObool, 'At least two labs are required to be included in the consensus for this method.')
    )

    withProgress(message = 'Calculating Unilateral DoEs', value = 0, style="old",{

      x=labmeanspool()
      u=labsespool()

      if(length(labdfspool())!=0){
        nu=labdfspool()
      }else{
        nu=NULL
      }

      if(length(numericPoolWeights())==0){
        LPweights=rep(1,length(labmeanspool()[sanitize]))
      }else{
        LPweights=numericPoolWeights()
      }
      DoEUnilateralLinearPool(x, u, nu, lab,weights=LPweights, K=input$linearOPrep,LOO=as.logical(input$LOODoE),coverageProb=.95, linearOpRes = outpool())#LPcp())


    })

  })
  output$LPDoEConv <- renderUI({
    HTML(LPDoEUni()[["DoEwarn"]])
  })

  output$LPDoEUniTable <- renderTable(LPDoEUni()$DoE,
                                      digits=signifpool,
                                      display=rep("fg",6),
                                      caption = '<b> Unilateral Degrees of Equivalence',
                                      caption.placement='top',include.rownames=FALSE)




    LPDoEUniPlot_func=function(){

      if(length(lablabspool())==0){
        lab=paste("L", 1:length(labmeanspool()), sep="")
      }else{
        lab=lablabspool()
      }

      #we remove the values of labs starting with "-"
      sanitize = !startsWith( lab,"-")
      nI = length(labmeanspool()[sanitize])
      n= length(labmeanspool())
    par(fig=c(0,1,.6,1),lend=2)  #change this when DoE plot added
    plot(c(1,1), c(1,1), type="n", axes=FALSE,xlab="", ylab="", bty="n")
    legend("top",c("DoE estimate" ," 95% coverage interval"),
           lty=c(NA,1),pch=c(19,NA),col=c(col3,col4),pt.cex=c(1,1) ,lwd=c(5,5))


    par(fig=c(0,1,0,.8), new=TRUE,lend=2)
    plot(0,xlim=range(0.7,length(LPDoEUni()$DoE$DoE.x)+0.3),ylim=range(LPDoEUni()$DoE$DoE.Lwr,LPDoEUni()$DoE$DoE.Upr),type="n",     ylab="DoE",xlab="",xaxt="n", bty="n",axes=FALSE)
    abline(h=0,col="gray")



    lab=LPDoEUni()$DoE$Lab

    n = length(labmeanspool())
    mtext(lab[seq(1, n, 2)], side=1, at=seq(1, n, 2), line=0, cex=1, col="Black")
    mtext(lab[seq(2, n, 2)], side=1, at=seq(2, n, 2), line=1.2, cex=1, col="Black")


    axis(2)


    for(i in 1:length(LPDoEUni()$DoE$DoE.x)){
      arrows(i,LPDoEUni()$DoE$DoE.Lwr[i],i,LPDoEUni()$DoE$DoE.Upr[i],length=0, col=col4,  lwd=5)
    }
    points(x=1:length(LPDoEUni()$DoE$DoE.x),y=LPDoEUni()$DoE$DoE.x,pch=c(rep(19,nI),rep(21,n-nI)), cex=1.5, col=col3)

    }

    output$LPDoEUniPlot<- renderPlot({LPDoEUniPlot_func()

    })

  ####################################################
  ####################################################
  ### Bilateral
  ####################################################
  ####################################################


  LPDoEBi=eventReactive(input$dopool,{

    validate(
      need(!any(is.na(labmeanspool())), 'Measured values must be numeric, separated by commas.'),
      need(!any(is.na(labsespool()))&& !any(labsespool()<=0), 'Standard uncertainties must be numeric (greater than zero), separated by commas.'),
      need(!any(is.na(labdfspool())) && !any(labdfspool()<1), 'Degrees of freedom must be numeric (greater than or equal to 1), separated by commas.'),
      need(is.numeric(input$coverage) && input$coverage<=1 && input$coverage>=0, 'Invalid coverage probability, must be numeric between 0 and 1.'),
      need(!any(is.na(numericPoolWeights())) && !any(numericPoolWeights()<0), 'Weights must be numeric (greater than zero), separated by commas.'),
      need(sum(numericPoolWeights())>0 || length(numericPoolWeights())==0, 'At least one of the weights should be greater than zero.'),
      need(is.numeric(input$linearOPrep), 'Sample size must be numeric.'),
      need(input$linearOPrep>0 && input$linearOPrep<=1000000, 'Sample size must be greater than zero and less than 1000000.'),
      need(length(labmeanspool())==length(numericPoolWeights()) || length(numericPoolWeights())==0, "The number of laboratory measured values does not equal the number of weights."),
      need(length(labmeanspool())==length(labsespool()),'The number of measured values entered does not equal the number of standard uncertainties entered.'),
      need(length(labmeanspool())==length(labdfspool()) || length(labdfspool())==0, "The number of degrees of freedom entered does not equal the number of measured values entered."),

      need(length(labmeanspool())==length(lablabspool()) | length(lablabspool())==0,'The number of study labels entered does not equal the number of measured values entered.'),
      need(!any(table(lablabspool())>1) | length(lablabspool())==0,"Laboratory labels should be unique."),

      need(length(labmeanspool())>1, 'At least two measured values are required.'),
      need(length(labmeanspool())<=500, 'The NICOB will not combine more than 500 measured values.')

    )

    if(length(lablabspool())==0){
      lab=paste("L", 1:length(labmeanspool()), sep="")
    }else{
      lab=lablabspool()
    }

    #we remove the values of labs starting with "-"
    sanitize = !startsWith( lab,"-")
    nI = length(labmeanspool()[sanitize])

    ###Check number of included labs

    LOObool= as.logical(input$LOODoE)

    validate(
      need((nI>=3&&LOObool)||!LOObool, 'At least three labs are required to be included in the consensus for this method.'),
      need((nI>=2&&!LOObool)||LOObool, 'At least two labs are required to be included in the consensus for this method.')
    )

    withProgress(message = 'Calculating Bilateral DoEs', value = 0, style="old",{

      x=labmeanspool()
      u=labsespool()

      if(length(labdfspool())!=0){
        nu=labdfspool()
      }else{
        nu=NULL
      }

      if(length(numericPoolWeights())==0){
        LPweights=rep(1,length(labmeanspool()[sanitize]))
      }else{
        LPweights=numericPoolWeights()
      }

      DoEBilateralPool(x, u, nu, lab,weights=LPweights, K=input$linearOPrep,coverageProb=.95)
    })

  })



  output$LPDoEBiTable <- renderTable(LPDoEBi()$B.x,
                                     digits=signifpool,
                                     rownames = T,
                                     caption = '<b> Bilateral Degrees of Equivalence (Estimates)',
                                     caption.placement='top',sanitize.rownames.function =  bold.allrows)

  output$LPDoE.UBiTable <- renderTable(LPDoEBi()$B.U,
                                     digits=signifpool,
                                     rownames = T,
                                     caption = '<b> Bilateral Degrees of Equivalence (Expanded Uncertainties U95)',
                                     caption.placement='top',sanitize.rownames.function =  bold.allrows)




  LPDoEBiPlot_func<- function(){

    DoE.Bilateral.LP.x = LPDoEBi()$B.x
    DoE.Bilateral.LP.Lwr = LPDoEBi()$B.Lwr
    DoE.Bilateral.LP.Upr = LPDoEBi()$B.Upr

    n = dim(DoE.Bilateral.LP.x)[[1]]

    labNAMEs = dimnames(DoE.Bilateral.LP.x)[[1]]

    xx = (DoE.Bilateral.LP.Lwr * DoE.Bilateral.LP.Upr > 0)
    ij = which(xx, arr.ind=TRUE, useNames = TRUE)
    xx[xx == 1] = 0.9
    xx[xx == 0] = 0.1
    diag(xx) = rep(0.5, n)
    image(1:n, 1:n, xx, xlim=c(-1.5, n+1),breaks=c(-0.5, 0.4, 0.85, 1.5),
          col=c(col5,"grey36",col4), #I also like "white" for the diagonal, rather than "lightblue4", to convey that this should be ignored
          axes=FALSE, xlab="", ylab="", asp=1)
    K = nrow(ij)
    if (K > 0) {
      for (k in 1:K) {
        points(ij[k,1], ij[k,2], pch=8, font=3, col="Black") } }
    polygon(c(0.5,n+0.5,n+0.5,0.5), c(0.5,0.5,n+0.5,n+0.5),
            border="Black", lwd=2)

    segments(seq(1.5,n+0.5),0.5,seq(1.5,n+0.5),n+0.5)
    segments(n+0.5,seq(1.5,n+0.5),0.5,seq(1.5,n+0.5))

    mtext(labNAMEs, side=1, at=1:n, las=2, line=1.5, col="Blue", font=2,padj=0)
    text(rep(0, n), 1:n, labNAMEs, adj=1, col="Blue", font=2)



  }

  output$LPDoEBiPlot<-renderPlot({LPDoEBiPlot_func()})



  output$downloadLPUniDoEGraph <- downloadHandler(
    filename = function() { paste( 'LinearPool_Unilateral_DoE.pdf', sep='') },
    content = function(file) {
      pdf(file,height=8, width=8)
      LPDoEUniPlot_func()
      dev.off()
    })
  output$downloadLPBiDoEGraph <- downloadHandler(
    filename = function() { paste( 'LinearPool_Bilateral_DoE.pdf', sep='') },
    content = function(file) {
      pdf(file,height=8, width=8)
      LPDoEBiPlot_func()
      dev.off()
    })





})
