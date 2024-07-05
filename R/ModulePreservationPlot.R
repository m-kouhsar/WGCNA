module.preservation.plot <- function(MP.Data , Modules,title.medrank="",title.zsummary=""){
  suppressMessages(library(WGCNA))
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggrepel))
  ref = 1 
  test = 2 
  
  PlotData=cbind.data.frame(modColors = rownames(MP.Data$preservation$observed[[ref]][[test]]) ,
                            moduleSizes = MP.Data$preservation$Z[[ref]][[test]][, 1],
                            Zsummary = MP.Data$preservation$Z[[ref]][[test]][, 2],
                            MedianRank=MP.Data$preservation$observed[[ref]][[test]][, 2])
  PlotData <- PlotData[PlotData$modColors %in% Modules,]
  
  # plot median rank
  p1 <- ggplot(PlotData, aes(moduleSizes, MedianRank,label=modColors,colour = modColors)) +         
    geom_point(size = 10)+ 
    scale_colour_manual(values=PlotData$modColors)+
    geom_text_repel(colour = "black",max.overlaps = 10,size=4, vjust = 0, nudge_y = 0.5,hjust = -0.3,nudge_x = 0)+
    #xlim(-50,1200)+
    xlab("Module size")+ylab("Median rank")+ggtitle(title.medrank)+ theme_bw() +
    theme(legend.position = "none")+scale_y_reverse() #(n.breaks =6,limits = c(12,0))
  
  # plot Z summary
  ylim_ <- c(-1 , round(max(PlotData$Zsummary),digits = 0)+5)
  xlim_ <- c(-50 , max(PlotData$moduleSizes)+100)
  p2 <- ggplot(PlotData, aes(moduleSizes, Zsummary,label=modColors,colour = modColors)) +         
    geom_point(size = 10)+
    scale_colour_manual(values=PlotData$modColors)+
    geom_text_repel(colour = "black",max.overlaps = 10,size=4, vjust = 0, nudge_y = 0.5,hjust = -0.3,nudge_x = 0)+
    xlab("Module size")+ylab("Z summary")+ggtitle(title.zsummary)+theme_bw() +
    theme(legend.position = "none")+geom_hline(yintercept=0, color = "black")+
    geom_hline(yintercept=5, linetype="dashed", color = "blue")+
    geom_hline(yintercept=10, linetype="dashed", color = "red") +
    scale_y_continuous(n.breaks = 8,limits = ylim_)+
    scale_x_continuous(n.breaks = 6,limits = xlim_)
  
  return(list(Plot.Med.Rank=p1 , Plot.Z.Summary=p2))
}

#####################################################################################################################################

#setwd("C:/Users/mk693/OneDrive - University of Exeter/Desktop/2021/NIH/Brain EWAS Paper")
#library(stringr)

#ModulePres_file="WGCNA/mad.0.5.lm.v2.pow3/Pitts.BDR.All.Braak.max.1000.gold.1000_ModulePreservation.Rdata"
#module.names <- "seashell4,mediumpurple3,thistle2,green4,firebrick3,sienna2,lightcoral,antiquewhite,magenta1,lightcyan,burlywood,linen,thistle3,lavenderblush,darkgoldenrod1,tan2,greenyellow,whitesmoke,orange,midnightblue,mistyrose,darkred,lightblue1"
#out_pref = "Figures/wgcna.mad.0.5.v2.pow3/Module.Preservation/Pitts.All.BDR.All3"

#load(ModulePres_file)
#module.names = str_split(module.names,pattern = ',',simplify = T)[1,]

#result <- module.preservation.plot(MP.Data = mp , Modules = module.names,title.medrank = "Median rank for module preservation of PITT-ADRC in BDR",
#                                   title.zsummary = "Z summary for module preservation of PITT-ADRC in BDR")

#tiff(file=paste0(out_pref,".ModulePreservation.MedianRank.tiff"), wi=10, h=8,units = "in",res = 500) 
#result$Plot.Med.Rank
#graphics.off()


#tiff(file=paste0(out_pref,".ModulePreservation.ZSummary.tiff"), wi=10, h=8,units = "in",res = 500) 
#result$Plot.Z.Summary
#graphics.off()


