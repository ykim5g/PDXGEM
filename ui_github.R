ui <- fluidPage(
  
  ## App title
  titlePanel("PDXGEM:"),

  #  shiny::runApp('PDX')
  sidebarLayout(
      sidebarPanel(
        
        ##### Input: Selector for choosing drug name
        selectInput( inputId = "MarkerSet", 
                     label="Select Biomarker Discovery Set", 
                     choices=c("Novartis", 
                          
                               "else"), 
                     selected = 'Novartis'
        ),

                        
        conditionalPanel(
          condition = "1==1", #input.MarkerSet %in% c('EuroPDX-Breast')",
          ## User Interface for selecting a drug or a drug combination
          uiOutput("uiDrug") ## inputID='drug'
        ),
        
        
        conditionalPanel(
          condition="input.drug != 'else'",
          uiOutput("uiTissue") ## inputID= tissueSite
        ), ## end of condition Paenl
        
        ## Input: Selector for choosing Original Cancer Type
        selectInput( inputId = "origin", 
                     label="Select Human Cancer Type for drug response prediction", 
                     choices=c("Breast", 
                               "Pancreatic",
                               "Colorectal", #"CRC-MEEC",
                                "NSCLC")
        )
        
      ), ## end of 'sidebarPanel'
    

          
      ## Main panel for displaying outputs  
      mainPanel(
        #h3(textOutput("Step 1: Drug Sensitivity Biomarker Discovery", container=span)),
        h4('Step 1: In vivo PDX Drug Sensitivity Biomarkers'),
        plotOutput(outputId="step1"),
        #highcharter::highchartOutput("hcontainerMarkerInfo",height = "300px"),

                
        conditionalPanel(
          condition="input.tissueSite != 'else'", 
          selectInput( inputId = "BiomarkerDiscoveryMethod", 
                       label="Step 1: Select a method for initial drug sensitivity biomarkers", 
                       choices=list("correlation"="correlation", 
                                    "correlation-rank"="correlation-rank", 
                                    "t-test"="t-test", 
                                    "t-test2"="t-test2",
                                   "else"="else"),
                       selected="else"
          )
        ), ## end of condition Paenl

               
        ## Biomarker discovery: test statistics vs p-values
        plotOutput(outputId="markerinfo"),
        textOutput(outputId="initMarker"), 
        fluidRow(
          ## buttone
          column(width=5, downloadButton(outputId="downInitMarker", label="Download") )
        ), 

        
        h4('Step 2: CCEA analysis of drug sensitivity biomarkers'),
        
        
        fluidRow(
          ## CCC Cutoff slide
          sliderInput('sliderCCC', label=h5('Adjust a threhold of CCEC'), min=0, max=1, value=0.2),
                      
          ## Download button
          column(width=5, downloadButton(outputId="downCCCMarker", label="Download") )
        ), 
        
        
        ## Input: Selector for choosing drug name
        selectInput( inputId = "useCCC", 
                     label="Do you want to filter initial drug sensitivity biomarkers by Concordant Co-Expression Analysis?", 
                     choices=c("Yes", "No"), 
                     selected = 'Yes'
        ),
        
        textOutput(outputId="outputCCCMarker"),
        ## textOutput(outputId="summary.drug.marker")
        
        

        h4('Step 3: Builing a multi-gene expression model for predictin a clinical response to the drug of interest'),
        ## Input: Selector for choosing Original Cancer Type
        selectInput( inputId = "ranked.response", 
                     label="Select response transformation", 
                     choices=c("None", "rank-transform"), ##, 'box-cox'),
                     selected='None'
        ),

        
        h4('Step 3-1: Training Random Forest-based model (PDXGEM) training'),

        plotOutput(outputId="step4RFtrainingset", height="400px"),
        fluidRow(
          ## Download button
          textOutput(outputId="outputRandomForestVariableImportance"), 
          column(width=5, downloadButton(outputId="downRandomForest", label="Download RF importance") ),
          column(width=5, downloadButton(outputId="downRFprediction", label="Download RF predictor") )
        ), 
        
                

        h4('Step 3-2: Validation of PDXGEM on an indepenent patient cohort'),
        fluidRow(
          # Input: Select a file ----
          # Multipurpose Internet Mail Extensions (MIME) type 
          fileInput("file_userTestData", "Choose CSV or Tab delimited File",
                    multiple = FALSE,
                    accept = c("text/csv", 
                               "text/comma-separated-values,text/plain",
                               ".csv", 
                               'text/tsv',
                               "text/tab-separated-values, text/plain",
                               ".tsv",
                               ".txt"))
        ),  
        
        
        conditionalPanel(
          condition= "1==1", #length(output.outputCCCMarker) > 0", ## the number of Concord biomarkers
          selectInput( inputId = "testData", 
                       label="Select a Model Test set", 
                       choices=c(
                         "BR:GSE20194(MAQCII)",
                         "BR:GSE41998:Horak(paclitaxel)", 
                         "BR:GSE20271:Tabchy(TFAC)", 
                         "BR:GSE20271:Tabchy(FAC)",
                         "BR:GSE25055:Hatzis", 
                         "BR:GSE32646:Miyake", 
                         "BR:USO(FAC+TX)", 
                         "BR:USO(FAC+TX+Trastuzumab)",
                         "BR:GSE22226:I-SPY(ACT+Trastuzumab):Agilent",
                         "BR:GSE22358:Agilent:Trastuzumab",
                         
                         
                         "CRC:Khambata(Cetuximab)",
                         "CRC:Khambata:KRASwild(Cetuximab)",
                         "CRC:Khambata:KRASmutation(Cetuximab)",
                         
                         "CRC_GSE62322_FOLFIRI_Primary", ## Liver met diseases and gene expression profiled on primary tissue
                         
                         "CRC_GSE39582_5FU+FOLFIRI",
                         "CRC_GSE39582_5FU+FOLFIRI_Metastasis",
                         "CRC_GSE39582_5FU+FOLFIRI_NoMetastasis",
                         
                         "CRC_GSE39582_5FU_NoMetastasis",
                         
                         "CRC_GSE39582_FOLFIRI_Metastasis",
                         "CRC_GSE39582_FOLFIRI_NoMetastasis", # FOLFIRI is not effective at primary disease setting

                         "CRC_GSE52735_Fluoropyrimidine",
                         
                         "CRC_GSE14095_LiverMetastasis", ## Japan data on Liver-metastasis CRC; No row data. 
                         "CRC_GSE19860_Metastasis_FOLFOX+Bevacizumab", ## Japan (Watanabe et al); No raw data, metastatic CRC
                                 "NSCLC:GSE31625CellLines:Erlotinib", 
                                 "NSCLC:GSE37138:Bevacizumab+Erlotinib",
                                 "NSCLC:BATTLE(Sorafenib|Erlotinib)",
                                 "NSCLC:BATTLE(All)",
                         
                                 "NSCLC:GSE68793:TCGA:NoErlotinib",
                         
                                 "HNSCC_GSE65021(Illumina)_Cetuximab",
                                 
                                 "Pancreatic:GSE57495:Moffitt",
                                 "Pancreatic:EMEXP2780",
                                 "Pancreatic:GSE17891:Patient",
                                 "Pancreatic_GSE28735_HuEXst1",
                                 "Pancreatic:ICGC_RNAseq"
                         
                       ), 
                       selected='BR:GSE20194(MAQCII)'
          )
        ), ## end of condition Paenl
        
        
        fluidRow(
          ## Download button
          column(width=5, downloadButton(outputId="downRFvalid", 
                                         label="Download a validation results (GEM score and clinical outcomes)") )
        ), 
        
        plotOutput(outputId="step4RFtestset", height="880px"),

        h4('Step 3: Model validation')
      ) ## end of 'main panel'

  ) ## end of 'sidebarLayout'
) ## end of 'fluidPage
  

#  shiny::runApp('Deconv')