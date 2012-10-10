setClass("OptimizedHierarchy",
         representation(nodes="matrix", edges="matrix", pathScores="vector"))


setMethod("summary", signature(object="OptimizedHierarchy"), function(object){
  cat(sprintf("An optimized hierarchy with %d nodes, %d edges, and %d paths\n", dim(object@nodes)[2], dim(object@edges)[2], length(object@pathScores)))
})


setMethod("plot", signature(x="OptimizedHierarchy", y="ANY"),
          function(x, 
                   phenotypeScores, 
                   uniformColors=FALSE, 
                   ylab=NULL, 
                   xlab=NULL, 
                   colors=c('blue','cyan','yellow','red'), 
                   edgeWeights=TRUE, 
                   edgeLabels=TRUE, 
                   nodeLabels=TRUE, 
                   min.score=NA, 
                   max.score=NA, 
                   cell.proportions=NULL, 
                   min.proportion=NA, 
                   max.proportion=NA, 
                   proportion.colors=c("black", "white"),
                   node.lwd=5,
                   root.name='All Cells',
                   legend.size=1.25,
                   plot.legend=TRUE){
            if (!is.vector(phenotypeScores))
              stop("phenptypeScores must be a numeric vector.")
            if (!is.logical(uniformColors))
              stop("uniformColors must be a logical value.")
            if (!(is.null(ylab) | is.character(ylab)))
              stop("ylab must be a character string or NULL.")
            if (!(is.null(xlab) | is.character(xlab)))
              stop("xlab must be a character string or NULL.")
            if (!is.vector(colors))
              stop("colors must be a color vector.")
            if (!is.vector(proportion.colors))
              stop("proportion.colors must be a color vector.")
            if (!is.logical(edgeWeights))
              stop("edgeWeights must be a logical value.")
            if (!is.logical(edgeLabels))
              stop("edgeLabels must be a logical value.")
            if (!is.logical(nodeLabels))
              stop("nodeLabels must be a logical value.")
            if (!(is.na(min.score) | is.numeric(min.score)))
              stop("min.score must be a numeric value or NA.")
            if (!(is.na(max.score) | is.numeric(max.score)))
              stop("max.score must be a numeric value or NA.")
            if (!(is.vector(cell.proportions) | is.null(cell.proportions)))
              stop("cell.proportions must be a numeric vector.")
            if (!(is.na(min.proportion) | is.numeric(min.proportion)))
              stop("min.proportion must be a numeric value or NA.")
            if (!(is.na(max.proportion) | is.numeric(max.proportion)))
              stop("max.proportion must be a numeric value or NA.")
            if (!is.numeric(node.lwd))
              stop("node.lwd must be a numeric value.")
            if (!is.character(root.name))
              stop("root.name must be a character string.")
            
            
            change.adrin.nima.base3.type<-function(phenotype)
              {
                subChar <- function(the.str, old.val, new.val)
                  {
                    split.string <- strsplit(the.str, '')[[1]]
                    split.string[which(split.string==old.val)] <- new.val
                    paste(split.string, sep='', collapse="")
                  }
                phenotype=subChar(phenotype,'0','3')
                phenotype=subChar(phenotype,'1','0')
                phenotype=subChar(phenotype,'3','1')
                return(phenotype)
              }
            
            SetTextContrastColor <- function(color){
              ifelse( mean(col2rgb(color)) > 127, "black", "white")
            }
            
            
            OptH=x
            ab=OptH
            Scores=phenotypeScores
            ed <- vector("list", length=dim(ab@nodes)[2])
            V=vector();
            for (i in 1:dim(ab@nodes)[2]){  
              V[i]=ab@nodes[2,i]
            }
            
            V[1]=root.name
            names(ed) <- V
            for (i in 1:dim(ab@nodes)[2]){  
              ed[[i]] <- list(edges=c(), weights=c(), labels=c())
            }
            for (i in 1:dim(ab@edges)[2]){
              temp=RchyOptimyx:::getNodeNumber(ab@edges[1,i], ab@nodes[1,])
              ed[[temp[1]]]$edges=c(ed[[temp[1]]]$edges,temp[2])
              if(edgeWeights){
                ab@edges[2,i] <- max(0, as.numeric(ab@edges[2,i]))
                ed[[temp[1]]]$weights=c(ed[[temp[1]]]$weights,as.numeric(ab@edges[2,i]))
                ##     Ves <- unlist(strsplit(ab@edges[1,i],'~'))
                ##      Ves[1]=ab@nodes[2,which(ab@nodes[1,]==Ves[1])]
                ##      Ves[2]=ab@nodes[2,which(ab@nodes[1,]==Ves[2])]
                ed[[temp[1]]]$labels=c(ed[[temp[1]]]$labels,ab@edges[3,i])      
              }
            }
            g <- new("graphNEL", nodes=V, edgeL=ed, edgemode='directed')
            attrs <- list(node = list(shape = "ellipse", fixedsize = FALSE))
            
            colorFunc <- colorRampPalette(colors)
            colorFunc2 <- colorRampPalette(proportion.colors)
            scores=unlist(lapply(1:length(ab@nodes[1,]), 
              function(i) {return(Scores[as.intBase(change.adrin.nima.base3.type(ab@nodes[1,i]), base=3)])}))
            
            
            if (is.na(min.score))
              min.score<-min(scores)
            if (is.na(max.score))
              max.score<-max(scores)
            
            proportions=vector();
            if(!is.null(cell.proportions)){
              proportions=unlist(lapply(1:length(ab@nodes[1,]), 
                function(i) {return(cell.proportions[as.intBase(change.adrin.nima.base3.type(ab@nodes[1,i]), base=3)])}))
            }
            
            if(!is.null(cell.proportions)){
              if (is.na(min.proportion))
                min.proportion<-min(proportions)
              if (is.na(max.proportion))
                max.proportion<-max(proportions)
            }
            
            if(!uniformColors){
              z=pretty(c(min.score, max.score))
              colinds=unlist(lapply(1:length(scores), function(i){max(1,ceiling((scores-min(z))*1000/(max(z)-min(z)))[i])}))
              cols=colorFunc(1000)[colinds]
              
              if (!is.null(cell.proportions)){
                z2=pretty(c(min.proportion, max.proportion))
                colinds2=unlist(lapply(1:length(proportions), function(i){max(1,ceiling((proportions-min(z2))*1000/(max(z2)-min(z2)))[i])}))
                cols2=colorFunc2(1000)[colinds2]
              }
            }
            
            nAttrs <- list()
            nAttrs$fillcolor=rep('white', length(ab@nodes[3,]))
            nAttrs$width=rep(0.75, length(ab@nodes[3,]))
            nAttrs$style=rep('bold', length(ab@nodes[3,]))
            ##nAttrs$shape=rep('ellipse', length(ab@nodes[3,]))
            if (uniformColors==FALSE){
              if (!is.null(cell.proportions)){
                nAttrs$color=rep(cols2, length(ab@nodes[3,]))
              }
              nAttrs$fillcolor=cols
            }
            nAttrs$label=rep('',length(V))
            if(nodeLabels)
              nAttrs$label=V
            ##nAttrs$fillcolor=ab$nodes[3,]
            
            names(nAttrs$fillcolor)=V
            names(nAttrs$width)=V
            names(nAttrs$label)=V
            names(nAttrs$style)=V
            ##names(nAttrs$shape)=V
            if (!is.null(cell.proportions)){
              names(nAttrs$color)=V
            }  
            eAttrs <- list()
            edgevalues <- as.numeric(ab@edges[2,])
            eAttrs$label=rep('',length(ab@edges))
            if(edgeLabels)
              eAttrs$label=unlist(lapply(1:length(ed), function(x) { return(ed[[x]]$labels)}))
            eAttrs$color=rep('gray', length(ab@edges[3,]))
            #eAttrs$arrowhead=rep('empty', length(ab@edges[3,]))
            #names(eAttrs$arrowhead)=edgeNames(g)
            names(eAttrs$color)=edgeNames(g)
            names(eAttrs$label)=edgeNames(g)
            
            if ((!uniformColors)){
              if (plot.legend){
                delta.sizex=legend.size/dev.size()[1]
                delta.sizey=legend.size/dev.size()[2]
                ##delta.sizex=0.2
                ##delta.sizey=0.2   
                if (is.null(cell.proportions))
                  split.screen(t(cbind(c(0,1-delta.sizex,0,1),c(1-delta.sizex,1,0,1))))
                if (!is.null(cell.proportions))
                  split.screen(t(cbind(c(0,1-delta.sizex,delta.sizey,1),c(1-delta.sizex,1,delta.sizey,1),c(0,1-delta.sizex,0,delta.sizey))))
                screen(1)
                par(mar=c(0,0,0,0))
              }
            }
            lag <- layoutGraph(g, nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs)

            edgeRenderInfo(lag) <- list(lwd = rep(1,length(ab@edges[2,])))
            if (uniformColors==FALSE) {
                if (!is.null(cell.proportions)) {
                    lwds <- rep(node.lwd, length(ab@nodes[3,]))
                    names(lwds) <- V
                    nodeRenderInfo(lag) <- list(lwd = lwds)
                }
                if (edgeWeights) {
                    lwds <-  as.numeric(unlist(lapply(1:length(ed), function(x) { return(0.75+ 15*(ed[[x]]$weights - min(edgevalues))/(max(edgevalues) - min(edgevalues)))})))
                    names(lwds) <- edgeNames(g)
                    edgeRenderInfo(lag) <- list(lwd = lwds)
                }
                
                textCols <- unlist(sapply(cols, SetTextContrastColor))
                names(textCols) <- V
                nodeRenderInfo(lag) <- list(textCol = textCols)
                cexs <- rep(0.5, length(ab@edges[2,]))
                names(cexs) <- edgeNames(g)
                edgeRenderInfo(lag) <- list(cex = cexs)
            }
            renderGraph(lag)
                    
            if (!uniformColors){
              if (plot.legend){
                screen(2)
                par(mar=c(1,4,1,0.2))
                image(matrix(1:2500, 50), col = colorFunc(50), xaxt='n', ylab='', yaxt='n')
                par(mgp=c(2.5,1,0))    
                title(ylab=ylab);
                ##z=z[which(z<max(scores))]
                axis(2, at=c((1/(length(z)-1))*(0:(length(z)-1))), labels=z)
                ##axis(2, at=z*max(scores), labels=z)
                if (!is.null(cell.proportions)){
                  screen(3);
                  par(mar=c(4,1,0.2,1))
                  image(t(matrix(1:2500,50)), col = colorFunc2(50), yaxt='n', xlab='', xaxt='n')
                  par(mgp=c(2.5,1,0))
                  title(xlab=xlab);
                  axis(1, at=c((1/(length(z2)-1))*(0:(length(z2)-1))), labels=z2)
                }
                close.screen(all.screens = TRUE)
              }              
            }
          })
