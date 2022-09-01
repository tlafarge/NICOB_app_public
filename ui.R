######################################################################
######################################################################
####################### NIST Consensus Builder #######################
############################ Version 1.4 #############################
######################################################################
######################################################################

shinyUI(fixedPage(id = "fullpage",


    tags$head(
      tags$style(HTML("
                      .shiny-output-error-validation {
                      color: red;
                      }
                      "))
      ),

    tags$head(HTML("
      <!-- Google tag (gtag.js) -->
        <script async src='https://www.googletagmanager.com/gtag/js?id=G-5496WP7RNQ'></script>
        <script>
        window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      
      gtag('config', 'G-5496WP7RNQ');
      </script>
    ")),
    
    
    tags$head(tags$link(rel="shortcut icon", href="./favicon.ico")),
    tags$head(HTML("<title>NIST Consensus Builder</title>")),
    
    #NIST Header
    HTML('
      <link rel="stylesheet" href="css/nist-combined.css">
      <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
      <script src="js/nist-header-footer.js" type="text/javascript" defer="defer"></script>
      <html class="nist-footer-bottom">
    '),


  fluidRow(br(),column(6,h1("NIST Consensus Builder"),offset=0)),
  #fluidRow(),
  fluidRow(id = "content",
    column(12,offset=0,
           br(),
           navlistPanel(selected=h4("Enter data"),

             tabPanel(
               h3("About the NIST Consensus Builder"),
               h3(em(strong("Introduction"))),

               h5("The NIST Consensus Builder (NICOB) serves to combine measurement results
obtained by different laboratories or by application of different measurement
methods, into a consensus estimate of the value of a scalar measurand. The NICOB
qualifies the consensus estimate with an evaluation of measurement uncertainty that captures not only the stated uncertainties associated with the individual measured
values, but also any additional component of uncertainty that manifests itself
only when these measured values are inter-compared."),
               # br(),
               h5("In addition, the NICOB can also report the differences between the individual measured values and the consensus value, and the differences between pairs of values measured by different laboratories or methods, in both cases qualifying these differences
with evaluations of associated uncertainty. In the context of Key Comparisons,
these differences and associated uncertainties are called (unilateral, and bilateral, respectively) Degrees of Equivalence (DoE) (Comité International des Poids
et Mesures (CIPM), 1999)."),
               # br(),
               h5("When the reported measurement uncertainties associated with the
individual measured values are qualified with the numbers of degrees of freedom
that they are based on, these numbers are taken into account as well. In general, the
numbers of degrees of freedom convey the reliability of the evaluations of measurement
uncertainty, expressing the extent of the underlying evidentiary basis, be it
the size of the experimental data or the strength of other information used when
producing the evaluations."),

               # br(),
               #h5("For more information see the User's Manual at (web address)."),
               HTML("For more information see the <a  href='./NISTConsensusBuilder-UserManual.pdf'> User's Manual.</a>"),

               h3(em(strong("Quick Start"))),

               h4("Enter Data"),
               tags$ul(
                 tags$li(h5( "Labels designating the n participating laboratories
                             (required — character strings comprised of letters or
                             numbers, separated from one another by
                             commas)")),
                 tags$li(h5(HTML(paste("Measured values x",tags$sub("1"),",..., x",tags$sub("n")," produced by n different laboratories or measurement methods,
                                       independently of one another (required — numbers separated by commas)",sep="")))),
                 tags$li(h5( "Measurement units to qualify the numerical values of the
                             measured values which are used to label axes of plots (optional — character string)")),
                 tags$li(h5(HTML(paste("Standard uncertainties u",tags$sub("1"),",..., u",tags$sub("n")," associated with the measured values (required — positive numbers separated by commas)",sep="")))),
                 tags$li(h5(HTML(paste("Numbers of degrees of freedom ","&nu;",tags$sub("1"),",..., &nu;",tags$sub("n")," that the standard uncertainties are based on
                                       (optional — positive numbers separated by commas)",sep="")))),
                 tags$li(h5("Coverage probability (positive number between 0 and 1) desired for the coverage intervals (required, default: 0.95).")),
                 tags$li(h5("Indication, by means of a checkbox, of whether degrees of equivalence should be computed (default: Not computed).")),
                 h5("If this box is checked, the following additional inputs appear:"),
                 tags$ul(
                   tags$li(h6( "Indication, by means of a radio button, of whether degrees of equivalence should be
                               computed as defined in the MRA or based on leave-one-out estimates.")),
                   tags$li(h6( "Number of bootstrap replicates for degrees of equivalence uncertainty calculation. This is only used for the
                   DerSimonian-Laird procedure (default: 10000); the Hierarchical Bayes and Linear Pool procedures
                   use for this number the sample sizes of their method specific inputs."))
                 )

               ),
               h5(HTML(paste("Buttons at the bottom of the ", em("Enter data"), " page allow the user to load and save
configuration files with inputs for the NICOB. Clicking the button
labeled ", em("Save Configuration File"), " downloads a plain text
file named ", em("consensus.ncb")," to the local machine, which
specifies the current inputs for the NICOB. To use a previously saved
configuration file, search for and select the file using the ", em("Browse")," button.",sep=""))),

               h5(HTML(paste("Alternatively, the NICOB also accepts configuration files with inputs
                             specified as comma separated values and extensions ",em(".ncb"),", ",
                             em(".csv"),", or ",em(".txt"),". Each row of the file designates data
                             from a different laboratory or measurement method, and the file can
                             have two, three, or four comma separated columns. For each row, data
                             should by entered in the order: name (if available), measured value,
                             standard uncertainty, and number of degrees of freedom (if
                             available; missing or infinite degrees of freedom should be entered as Inf).",sep=""))),

               h4("Choose a method for analysis"),
               h5("The three methods of data reduction implemented in the NICOB are not meant
to be interchangeable, and the user should consider their characteristics, including advantages and disadvantages, to determine which may be best fit for the
purpose that the results of the data reduction are intended to serve."),
               HTML("For more information see the <a  href='./NISTConsensusBuilder-UserManual.pdf'> User's Manual.</a>"),

               tags$ul(
                 tags$li(h5("For the DerSimonian-Laird procedure:"),
                         h6("– If the modified Knapp-Hartung adjustment is desired (as described by Rover, Knapp, and Friede (2015, DOI 10.1186/s12874-015-0091-1)), check the corresponding box;"),
                         h6("– To apply the parametric bootstrap for uncertainty evaluation, check the corresponding box.
                            This reveals an input field for the desired number of bootstrap replicates (Suggested value: 10000).")
                 ),
                 tags$li(h5("For the hierarchical Bayesian procedure:"),
                         h6("– Positive numbers in the corresponding boxes which are used as the medians of
    the prior distributions for the between-laboratory and for
    the within-laboratory (or, between-method and within-method)
    variance components (default: robust indications of
    spread of the measured values for the between-laboratory
    variability, and for the laboratory-specific uncertainty)"),

                         h6("– Total number of iterations for the Markov Chain Monte Carlo
                         sampler (default: 250000)."),

                         h6("– Length of burn-in for the Markov chain, which is the number of
                            values discarded from the beginning of the realization of the chain (default: 50000)."),

                         h6("– Thinning rate for the Markov chain (default value is 25, meaning that only every 25th value generated in
                         the chain should be kept).")

                 ),
                 tags$li(h5("For the Linear Pool:"),
                         h6("– Weights (non-negative numbers separated by commas) to be associated with the different measurement results
                            (default: 1 for all)."),
                         h6("– Size of sample drawn from the mixture distribution of the measurand
                         (default: 100000).")

                         )
               ),
               h4("Output"),

               tags$ul(
                 tags$li(h5("Consensus estimate, associated standard uncertainty, and coverage interval for the true value of the measurand;")),
                 #tags$li(h5("Estimates, standard errors, and confidence intervals for the laboratory or method effects;")),
                 tags$li(h5("If degrees of equivalence  were requested, then estimates, standard uncertainties, and expanded uncertainties for
                            differences between measured values and the consensus value, and between pairs of measured values are reported and depicted graphically."))
               )





             ),
             HTML("<a  href='./NISTConsensusBuilder-UserManual.pdf' style='padding-left: 0px;'>User's Manual</a>"),


             "---------------------",

             tabPanel(
               title=h4("Enter data"),

               h4("List laboratory labels, measured values, standard uncertainties,
                        and (if available) numbers of degrees of freedom, separated by
                        commas." ),
               h5("Add a minus sign in front of the label of the labs to be left out of
                        the consensus value calculation,
                        but for which degrees of equivalence are desired" ),

               textInput("lablabels",
                         label = h4(span("Laboratories name"
                                         ,tags$span(title="Add a \"-\" in front of a lab label to extract it from the consensus computation but not the DoE"
                                                    ,tagList(icon("question-circle"))))),
               ),

               splitLayout(cellWidths = c("50%", "50%"),
                           textInput("mean",
                                     label = h4(HTML(paste("Measured values *",sep="")))#, ".8247,.8184,.8196,.8170,.8069,.8355,.8186,.8236"#,# (y",tags$sub("i"),")",sep=""))),
                                     #"34.3,32.9, 34.53, 32.42, 31.9, 35.8"
                           ),
                           textInput("units",
                                     label = h4("Measurement units, e.g. mg/kg ")
                                     #"mg/kg"
                           )
               ),


               textInput("se",
                         label = h4(HTML(paste("Standard uncertainties *",sep="")))#, ".0095,.0112,.0033,.007,.0072,.0130,.0038,.0058"#,# (u",tags$sub("i"),")",sep=""))),
                         #"1.03, 0.69, 0.83, 0.29, 0.4, 0.38"
               ),
               textInput("df",
                         label = h4("Numbers of Degrees of Freedom ")#,"60, 4, 18, 2, 13, 60"
               ),
               numericInput("coverage",
                            label = h4("Coverage probability *"),value = .95,step=.05, min = 0.01, max = 0.99
               ),
               numericInput("signifDigits",
                            label = h4("Number of significant digits desired for results"), value = 2,step=1,min=0
               ),
               numericInput("seed",
                            label = h4("Random number generator seed"),value = 5,step=1, min = 1, max = 100
               ),
              #HTML("<span style='font-size: 85%;'>* Required field</span>"),

               # ######## DoE
               br(),
               br(),
               h4("Degrees of equivalence"),
               #
               checkboxInput("DoE", label = "Compute degrees of equivalence", value = FALSE),

               conditionalPanel("input.DoE == 1",
                                radioButtons("LOODoE",label = strong("Type"),
                                             choices = c("DoEs conforming to MRA"=FALSE,"DoEs based on Leave-One-Out estimates"=TRUE),inline=T),
                                numericInput("DoEbootrep", label = strong("Number of bootstrap replicates (only used for DerSimonian-Laird procedure)"), value = 10000,#NA,
                                             step=1000,min=100)

               ),

              # ######## Validate
              br(),
              br(),

              actionButton("validate", "Validate model inputs", icon("check"),
                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
              htmlOutput("val",inline=T),
              br(),
              br(),
              #actionButton("load_inputs", "Load inputs"),
              fileInput('file1', 'Load configuration file',
                        accept=c('text/RDS','.ncb','.csv','.txt')),
              downloadButton('save_inputs', 'Save configuration file'),
              br(),
              br()
             ),
             "---------------------",
             "Choose a method for analysis",
             ######## rma with metafor
             tabPanel(

               title=h4("DerSimonian-Laird"),
               h4("Fit laboratory effects model using DerSimonian-Laird procedure"),

               checkboxInput("knhacheck", label = "Modified Knapp-Hartung adjustment", value = FALSE),

               checkboxInput("doboot", label = "Parametric bootstrap for uncertainty evaluation", value = FALSE),

               conditionalPanel("input.doboot == 1",
                                numericInput("bootrep", label = strong("Number of bootstrap replicates"), value = NA,step=1000,min=100)
               ),



               actionButton("go", "Fit the model"),
               br(),
               br(),


               conditionalPanel(
                 condition="input.go >= 1",

                                h2("DerSimonian-Laird"),
                                # numericInput("roundDec", label = h5("Number of decimal places to round to"), value = 2,step=1,min=0),

                                p("The consensus estimate is: ",textOutput("rmaest",inline=T)),
                                p("The standard uncertainty is: ",textOutput("rmase",inline=T)),
                                p(textOutput("rmaquant")),
                								p("The dark uncertainty (tau) is: ",textOutput("rmatau",inline=T)),
                                plotOutput("rmaplot", width = "auto", height = "700px"),
                                downloadButton('downloadrmagraph', 'Download plot'),

                 br(),
                 br()
               ),


               conditionalPanel("input.doboot ==true && input.go >= 1",
                                h3("Parametric bootstrap for uncertainty evaluation"),
                                p("The standard uncertainty is: ",textOutput("textboot",inline=T)),
                                p("The ",textOutput("coverageProbabilityPercentBoot",inline=T)," coverage interval ranges from: ",textOutput("bootquant",inline=T)),

                                br(),
                                br(),
                                downloadButton('downloadDLbootout', 'Download bootstrap output'),
                                br(),
                                br()
               ),




               conditionalPanel(
                 condition = "input.DoE == true && input.go >= 1",
                 h2("Unilateral degrees of equivalence"),
                 h5(htmlOutput("DLDoEConv")),
                 tags$head(tags$style("#DLDoEConv{color: red;
                                      font-style: italic;
                                      }")),
                 plotOutput("DLDoEUniPlot", width = "auto", height = "500px"),
                 downloadButton('downloadDLUniDoEGraph', 'Download unilateral DoE plot'),
                 br(),
                 br(),
                 tableOutput("DLDoEUniTable"),

                 h2("Bilateral degrees of equivalence"),
                 helpText(
                   "Yellow squares (with black asterisks in the center) indicate results that differ significantly from 0
                   at 95% coverage. Light blue squares indicate results that do not differ significantly at 95% coverage.
                   Dark squares are space-fillers for results when compared to themselves."
                 ),
                 plotOutput("DLDoEBiPlot", width = "auto", height = "500px"),
                 downloadButton('downloadDLBiDoEGraph', 'Download bilateral DoE plot'),
                 br(),
                 br(),

                 fluidRow(

                   splitLayout(cellWidths = c("50%", "50%"),tableOutput("DLDoEBiTable"),tableOutput("DLDoE.UBiTable"))

                 )

               )

             ),


             ###############Bayesian
             tabPanel(
               title=h4("Hierarchical Bayes (Gaussian)"),

               h4("Fit using Bayesian method with weakly informative priors"),
               uiOutput("dyn_input1"),
               helpText("Default is the median of the absolute values of the differences between the measured values and their median"),
               uiOutput("dyn_input2"),
               helpText("Default is the median of the lab-specific standard uncertainties"),

               numericInput("niters", label = h4("Total number of iterations"), value = "250000",step=1000,min=0),
               numericInput("nburnin", label = h4("Length of burn in"), value = "50000",step=1000,min=0),
               numericInput("nthin", label = h4("Thinning rate"), value = "25",min=0),


               actionButton("dobayes", "Fit the model"),
               br(),
               br(),



               conditionalPanel(
                 condition = "input.dobayes >= 1",


                 h2("Bayesian procedure "),
                 h4(div("Assuming weakly informative prior distributions and allowing for uncertainty in standard uncertainties"), style = "color:gray"),
                 h5(htmlOutput("bayesConv")),
                 tags$head(tags$style("#bayesConv{color: red;
                                      font-style: italic;
                                      }")),

                 h5(htmlOutput("bayesTauWarn")),
                 tags$head(tags$style("#bayesTauWarn{color: blue;
                                      font-style: italic;
                                      }")),
                 h5(htmlOutput("bayesSigWarn")),
                 tags$head(tags$style("#bayesSigWarn{color: blue;
                                      font-style: italic;
                                      }")),

                 p("The consensus estimate is: ",textOutput("bayesest",inline=T)),
                 p("The standard uncertainty is: ",textOutput("bayesse",inline=T)),
                 p("The ",textOutput("coverageProbabilityPercentBayes",inline=T)," credible interval ranges from: ",textOutput("bayesquant",inline=T)),
                 p("The dark uncertainty (tau) is: ",textOutput("bayestau",inline=T)),
                 plotOutput("bayesplot", width = "auto", height = "700px"),
                 p(textOutput("bayesplotLegendWarning",inline=T)),
                 
                 # br(),
                 # 
                 # br(),
                 # br(),

                 br(),
                 br(),
                 downloadButton('downloadbayesout', 'Download MCMC output'),
                 downloadButton('downloadBayesGraph', 'Download plot'),
                 br(),
                 br()

               ),




               conditionalPanel(
                 condition = "input.DoE == true && input.dobayes >= 1",
                 h2("Unilateral degrees of equivalence"),
                 h5(htmlOutput("bayesDoEConv")),
                 tags$head(tags$style("#bayesDoEConv{color: red;
                                      font-style: italic;
                                      }")),

                 plotOutput("BayesDoEUniPlot", width = "auto", height = "500px"),
                 downloadButton('downloadBayesUniDoEGraph', 'Download unilateral DoE plot'),
                 br(),
                 br(),
                 tableOutput("BayesDoEUniTable"),

                 h2("Bilateral degrees of equivalence"),
                 helpText(
                   "Yellow squares (with black asterisks in the center) indicate results that differ significantly from 0
                   at 95% coverage. Light blue squares indicate results that do not differ significantly at 95% coverage.
                   Dark squares are space-fillers for results when compared to themselves."
                 ),
                 plotOutput("BayesDoEBiPlot", width = "auto", height = "500px"),
                 downloadButton('downloadBayesBiDoEGraph', 'Download bilateral DoE plot'),

                 br(),
                 br(),


                 fluidRow(

                   splitLayout(cellWidths = c("50%", "50%"),tableOutput("BayesDoEBiTable"),tableOutput("BayesDoE.UBiTable"))

                 )
               )

             ),
             ###############Bayesian Laplace
             tabPanel(
               title=h4("Hierarchical Bayes (Laplace)"),

               h4("Fit using Bayesian method with weakly informative priors"),
               uiOutput("dyn_input1_lap"),
               helpText("Default is the median of the absolute values of the differences between the measured values and their median"),
               uiOutput("dyn_input2_lap"),
               helpText("Default is the median of the lab-specific standard uncertainties"),

               numericInput("niters_lap", label = h4("Total number of iterations"), value = "250000",step=1000,min=0),
               numericInput("nburnin_lap", label = h4("Length of burn in"), value = "50000",step=1000,min=0),
               numericInput("nthin_lap", label = h4("Thinning rate"), value = "25",min=0),


               actionButton("dobayes_lap", "Fit the model"),
               br(),
               br(),



               conditionalPanel(
                 condition = "input.dobayes_lap >= 1",


                 h2("Bayesian procedure"),
                 h4(div("Assuming weakly informative prior distributions and allowing for uncertainty in standard uncertainties"), style = "color:gray"),
                 h5(htmlOutput("bayesConv_lap")),
                 tags$head(tags$style("#bayesConv_lap{color: red;
                                       font-style: italic;
                                      }")),

                 h5(htmlOutput("bayesTauWarn_lap")),
                 tags$head(tags$style("#bayesTauWarn_lap{color: blue;
                                      font-style: italic;
                                      }")),
                 h5(htmlOutput("bayesSigWarn_lap")),
                 tags$head(tags$style("#bayesSigWarn_lap{color: blue;
                                      font-style: italic;
                                      }")),

                 p("The consensus estimate is: ",textOutput("bayesest_lap",inline=T)),
                 p("The standard uncertainty is: ",textOutput("bayesse_lap",inline=T)),
                 p("The ",textOutput("coverageProbabilityPercentBayes_lap",inline=T)," credible interval ranges from: ",textOutput("bayesquant_lap",inline=T)),
                 p("The dark uncertainty (tau) is: ",textOutput("bayestau_lap",inline=T)),
                 plotOutput("bayesplot_lap", width = "auto", height = "700px"),
                 p(textOutput("bayesplotLegendWarning_lap",inline=T)),

                 br(),
                 br(),
                 downloadButton('downloadbayesout_lap', 'Download MCMC output'),
                 downloadButton('downloadBayesGraph_lap', 'Download plot'),
                 br()

               ),




               conditionalPanel(
                 condition = "input.DoE == true && input.dobayes_lap >= 1",
                 h2("Unilateral degrees of equivalence"),
                 h5(htmlOutput("bayesDoEConv_lap")),
                 tags$head(tags$style("#bayesDoEConv_lap{color: red;
                                      font-style: italic;
                                      }")),

                 plotOutput("BayesDoEUniPlot_lap", width = "auto", height = "500px"),
                 downloadButton('downloadBayesUniDoEGraph_lap', 'Download unilateral DoE plot'),
                 br(),
                 br(),
                 tableOutput("BayesDoEUniTable_lap"),

                 h2("Bilateral degrees of equivalence"),
                 helpText(
                   "Yellow squares (with black asterisks in the center) indicate results that differ significantly from 0
                   at 95% coverage. Light blue squares indicate results that do not differ significantly at 95% coverage.
                   Dark squares are space-fillers for results when compared to themselves."
                 ),
                 plotOutput("BayesDoEBiPlot_lap", width = "auto", height = "500px"),
                 downloadButton('downloadBayesBiDoEGraph_lap', 'Download bilateral DoE plot'),

                 br(),
                 br(),


                 fluidRow(

                   splitLayout(cellWidths = c("50%", "50%"),tableOutput("BayesDoEBiTable_lap"),tableOutput("BayesDoE.UBiTable_lap"))

                 )
               )

             ),



             ###############opinion pooling
             tabPanel(
               title=h4("Linear Pool"),
               h4("Linear opinion pooling"),
               textInput("poolweights",
                         label = h4("Weights") #,"1,1,1,1,1,1"
                         ),
               helpText("If no weights are entered, default weights will all be equal to 1."),
               helpText("Only include weights for labs included in the consensus value calculation."),
               numericInput("linearOPrep", label = h4("Sample size"), value = 1e5,step=1000,min=0),

               actionButton("dopool", "Fit the model"),

               conditionalPanel(
                 condition = "input.dopool >= 1",

                 h2("Linear opinion pooling"),
                 p("The consensus estimate is: ",textOutput("poolest",inline=T)),
                 p("The standard uncertainty is: ",textOutput("poolse",inline=T)),
                 p("The ",textOutput("coverageProbabilityPercentPool",inline=T)," coverage interval ranges from: ",textOutput("poolquant",inline=T)),
                 plotOutput("poolplot"),

                 br(),
                 br(),
                 downloadButton('downloadLPout', 'Download linear pool output'),
                 downloadButton('downloadLPgraph', 'Download plot'),
                 br(),
                 br()
               ),





               conditionalPanel(
                 condition = "input.DoE == true && input.dopool >= 1",
                 h2("Unilateral degrees of equivalence"),
                 h5(htmlOutput("LPDoEConv")),
                 tags$head(tags$style("#LPDoEConv{color: red;
                                      font-style: italic;
                                      }")),
                 plotOutput("LPDoEUniPlot", width = "auto", height = "500px"),

                 downloadButton('downloadLPUniDoEGraph', 'Download unilateral DoE plot'),
                 br(),
                 br(),

                 tableOutput("LPDoEUniTable"),

                 h2("Bilateral degrees of equivalence"),

                 helpText(
                   "Yellow squares (with black asterisks in the center) indicate results that differ significantly from 0
                   at 95% coverage. Light blue squares indicate results that do not differ significantly at 95% coverage.
                   Dark squares are space-fillers for results when compared to themselves."
                 ),
                 plotOutput("LPDoEBiPlot", width = "auto", height = "500px"),
                 downloadButton('downloadLPBiDoEGraph', 'Download bilateral DoE plot'),

                 br(),
                 br(),

                 fluidRow(

                   splitLayout(cellWidths = c("50%", "50%"),tableOutput("LPDoEBiTable"),tableOutput("LPDoE.UBiTable"))

                 )
               )
             )
           )
         )
)





  ))
