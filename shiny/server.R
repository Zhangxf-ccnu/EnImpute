# server.R

options(shiny.maxRequestSize=100*1024^2)

function(input, output, session) {
  observeEvent(input$Run_EnImpute, {
    # Read in expression dataset
    count = read.csv(file = input$count$datapath, header = TRUE, row.names = 1)
    count = as.matrix(count)
    
    time.used = system.time({Impute.count = try(EnImpute(count, scale.factor = input$scale.factor, trim = input$trim, 
                                                         ALRA = input$ALRA, ALRA.k = input$ALRA.k, ALRA.q = input$ALRA.q, 
                                                         DCA = input$DCA, DCA.normtype = input$DCA.normtype,
                                                         DCA.type = input$DCA.type, DCA.l2 = input$DCA.l2, DCA.l1 = input$DCA.l1, 
                                                         DCA.l2enc = input$DCA.l2enc, DCA.l1enc = input$DCA.l1enc, DCA.ridge = input$DCA.ridge, 
                                                         DCA.gradclip = input$DCA.gradclip, DCA.activation = input$DCA.activation, 
                                                         DCA.hiddensize = input$DCA.hiddensize, DCA.hyper = input$DCA.hyper,
                                                         DCA.hypern = input$DCA.hypern, 
                                                         DrImpute = input$DrImpute, DrImpute.ks = as.numeric(unlist(strsplit(input$DrImpute.ks,","))),
                                                         DrImpute.dists = input$DrImpute.dists, DrImpute.method = input$DrImpute.method,
                                                         DrImpute.cls = NULL, 
                                                         MAGIC = input$MAGIC, MAGIC.k = input$MAGIC.k, MAGIC.alpha = input$MAGIC.alpha,
                                                         MAGIC.t = input$MAGIC.t, MAGIC.npca = input$MAGIC.npca, MAGIC.t.max = input$MAGIC.t.max,
                                                         MAGIC.knn.dist.method =  input$MAGIC.knn.dist.method, MAGIC.n.jobs = input$MAGIC.n.jobs, 
                                                         SAVER = input$SAVER,
                                                         SAVER.do.fast = input$SAVER.do.fast, SAVER.ncores = input$SAVER.ncores, 
                                                         SAVER.size.factor = NULL,
                                                         SAVER.npred = NULL, SAVER.null.model = FALSE, SAVER.mu = NULL,
                                                         scImpute = input$scImpute, scImpute.drop_thre = input$scImpute.drop_thre, 
                                                         scImpute.Kcluster = input$scImpute.Kcluster,
                                                         scImpute.labeled = FALSE, scImpute.labels = NULL,
                                                         scImpute.genelen = NULL, scImpute.ncores = input$scImpute.ncores, 
                                                         Seurat = input$Seurat, Seurat.genes.use = NULL, 
                                                         Seurat.genes.fit = NULL, Seurat.gram = input$Seurat.gram))})
    
    output$summary1 = renderText("Imputation is finished (EnImpute)")
    output$summary2 = renderText({
      paste("Time used (seconds):", round(as.numeric(time.used[3]),2))
    })
    
   
    
    output$plot = renderPlot({
      
      count.tsne = tsne_cluster(count)
      count.EnImpute.tsne = tsne_cluster(Impute.count$count.EnImpute.exp)
      par(mfrow = c(1, 2)) 
      plot(count.tsne@dr$tsne@cell.embeddings, col = as.numeric(count.EnImpute.tsne@ident), pch = 19, axes = FALSE, 
           frame.plot = TRUE,  cex = 0.6, main="Before EnImpute")
      plot(count.EnImpute.tsne@dr$tsne@cell.embeddings, col = as.numeric(count.EnImpute.tsne@ident), pch = 19, axes = FALSE, 
           frame.plot = TRUE, cex = 0.6, main="After EnImpute")
      cat("Imputation finished. You can see the tsne visualization and download the imputed data now.")
    })
    
    
    
    output$download = downloadHandler(
      
      filename = function(){ paste(input$dlname, ".", input$fileformat, sep="") },
      
      content = function(file) {
        if (input$fileformat == "txt"){
          if (input$scale=="exp")
            write.table(Impute.count$count.EnImpute.exp, file, sep = "\t", quote = FALSE)
          else
            write.table(Impute.count$count.EnImpute.log, file, sep = "\t", quote = FALSE)
        }
        else{
          if (input$scale=="exp")
            write.csv(Impute.count$count.EnImpute.exp, file, quote = FALSE)
          else
            write.csv(Impute.count$count.EnImpute.log, file,  quote = FALSE)
        }
      }
    )
    
  })
  
}







