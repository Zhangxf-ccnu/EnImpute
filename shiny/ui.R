# ui.R
library(shiny)
library(shinyBS)
library(EnImpute)
library(Seurat)

source("tsne_cluster.R")
shinyUI(navbarPage("EnImpute-Shiny", id="master",
                   fluidPage(
                     sidebarLayout(
                       sidebarPanel(

                         fileInput("count",
                                   label = "Count matrix",
                                   multiple = FALSE),

                         numericInput("scale.factor",
                                      label = "scale.factor",
                                      value = 10000, min =  1, step = 1000),

                         numericInput("trim",
                                      label = "trim",
                                      value = 0.3, min =  0, max = 0.5, step = 0.1),

                         checkboxInput("ALRA", "ALRA", value = TRUE),
                         checkboxInput("DCA", "DCA", value = TRUE),
                         checkboxInput("DrImpute", "DrImpute", value = TRUE),
                         checkboxInput("MAGIC", "MAGIC", value = TRUE),
                         checkboxInput("SAVER", "SAVER", value = TRUE),
                         checkboxInput("scImpute", "scImpute", value = TRUE),
                         checkboxInput("scRMD", "scRMD", value = TRUE),
                         checkboxInput("Seurat", "Seurat", value = TRUE),

                         numericInput("ALRA.k",
                                      label = "ALRA.k",
                                      value = 0, min =  0, step = 5),

                         numericInput("ALRA.q",
                                      label = "ALRA.q",
                                      value = 10, min =  1,  step = 5),


                         selectInput('DCA.normtype', 'DCA.normtype ', c("deseq", "zheng"), selected = "zheng"),
                         selectInput('DCA.type', 'DCA.type', c("normal", "poisson","nb", "nb-shared",
                                                           "nb-conddisp", "nb-fork", "zinb", "zinb-shared",
                                                           "zinbconddisp","zinb-fork"),
                                     selected = "zinb-conddisp"),
                         numericInput("DCA.l2",
                                      label = "DCA.l2",
                                      value = 0, min =  0, step = 1),
                         numericInput("DCA.l1",
                                      label = "DCA.l1",
                                      value = 0, min =  0, step = 1),
                         numericInput("DCA.l2enc",
                                      label = "DCA.l2enc",
                                      value = 0, min =  0, step = 1),
                         numericInput("DCA.l1enc",
                                      label = "DCA.l1enc",
                                      value = 0, min =  0, step = 1),
                         numericInput("DCA.ridge",
                                      label = "DCA.ridge",
                                      value = 0, min =  0, step = 1),
                         numericInput("DCA.gradclip",
                                      label = "DCA.gradclip",
                                      value = 5, step = 1),
                         selectInput('DCA.activation', 'DCA.activation',
                                     c("sigmoid", "tanh", "elu", "selu", "softplus", "softsign",
                                       "relu", "relu6", "crelu", "relu_x", "dropout"),
                                     selected = "relu"),
                         textInput("DCA.hiddensize", label = "DCA.hiddensize", value = "64,32,64"),
                         checkboxInput("DCA.hyper", "DCA.hyper", value = FALSE),
                         numericInput("DCA.hypern", label = "DCA.hypern", value = 1000),


                         textInput('DrImpute.ks', 'DrImpute.ks', "10,11,12,13,14,15"),
                         selectInput("DrImpute.dists", "DrImpute.dists",
                                     c("spearman", "pearson"), selected =  c("spearman", "pearson"), multiple = TRUE),
                         selectInput('DrImpute.method', 'DrImpute.method', c("mean", "med"), selected = "mean"),


                         numericInput("MAGIC.k",
                                      label = "MAGIC.k",
                                      value = 10, step = 1,  min=1),
                         numericInput("MAGIC.alpha",
                                      label = "MAGIC.alpha",
                                      value = 15, step = 1,  min=1),
                         textInput('MAGIC.t', 'MAGIC.t', "auto", placeholder="int, optional, default: 'auto'"),
                         numericInput("MAGIC.npca",
                                      label = "MAGIC.npca",
                                      value = 20, step = 1,  min=1),
                         numericInput("MAGIC.t.max",
                                      label = "MAGIC.t.max",
                                      value = 20, step = 1,  min=1),
                         selectInput("MAGIC.knn.dist.method", "MAGIC.knn.dist.method",
                                     c("euclidean", "cosine"), selected =  "euclidean"),
                         numericInput("MAGIC.n.jobs",
                                      label = "MAGIC.n.jobs",
                                      value = 1),


                         checkboxInput("SAVER.do.fast", "SAVER.do.fast", value = TRUE),
                         numericInput("SAVER.ncores", label = "SAVER.ncores", value = 2),


                         numericInput("scImpute.drop_thre",
                                      label = "scImpute.drop_thre",
                                      value = 0.5, step = 0.1,  min=0, max =1),
                         numericInput("scImpute.Kcluster",
                                      label = "scImpute.Kcluster",
                                      value = 5, step = 1,  min=1),
                         numericInput("scImpute.ncores",
                                      label = "scImpute.ncores",
                                      value = 1, step = 1,  min=1),

                         numericInput("scRMD.candidate",
                                      label = "scRMD.tau",
                                      value = 0.05, step = 0.05,  min=0, max =1),
                         numericInput("scRMD.lambda",
                                      label = "scRMD.lambda",
                                      value = NULL, min = 0, step = 1),
                         numericInput("scRMD.tau",
                                      label = "scRMD.tau",
                                      value = NULL, min = 0, step = 1),

                         checkboxInput("Seurat.gram", "Seurat.gram", value = TRUE),

                         actionButton("Run_EnImpute",
                                      label = "Run EnImpute",
                                      style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),

                         width = 2
                       ),
                       mainPanel(

                         fluidRow(
                           column(9,

                                  plotOutput("plot"),
                                  br(),
                                  textOutput("summary1"),
                                  br(),
                                  textOutput("summary2")
                           ),
                           column(3,
                                  radioButtons("scale",
                                               label = "Output File Scale",
                                               choices = c("log (logarithmic scale)" = "exp",
                                                           "exp (Exponential  scale)" = "log"),
                                               selected = "exp"),
                                  radioButtons("fileformat",
                                               label = "Output File Format",
                                               choices = c(".txt (tab-delimited text)" = "txt",
                                                           ".csv (comma-separated values)" = "csv"),
                                               selected = "csv"),
                                  textInput("dlname",
                                            label = "Output File Name (Do not include file extension)"),
                                  downloadButton("download",
                                                 label = "Download")
                           )

                         )
                     )
                   )
                       )
)
)
