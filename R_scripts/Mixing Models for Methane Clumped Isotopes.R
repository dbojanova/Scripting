##Script for mixing model

library(ggpubr)
library(ggplot2)
library(ggforce)
library(gridExtra)
library(reshape)

setwd("")

#dataframe of possible endmembers
col <- c("darkgreen","blue","red")
sites <- c("U1545","U1546","U1547")
em1 <- list()
em2 <- list()
em3 <- list()

 #U1545 
      em1[[1]] <- data.frame(describe = c("U1545-methano-lowG-30C","U1545-methano-lowG-30C"),
                        dD = c(-225,-225), 
                        d13C = c(-80,-80), 
                        D13CD = c(5.3,5.3), 
                        DDD = c(5,5))
      
      em2[[1]] <- data.frame(describe = c("U1545-thermo-300C","U1545-abiotic"),
                        dD = c(-150,-150), 
                        d13C = c(-40,-40), 
                        D13CD = c(1.6,0.5), 
                        DDD = c(-2,-3))

 #U1546

           em1[[2]] <- data.frame(describe = c("U1546-aom","U1546-methano"),
                             dD = c(-225,-225), 
                             d13C = c(-80,-80), 
                             D13CD = c(1,1), 
                             DDD = c(-30,-30))
           
           em2[[2]] <- data.frame(describe = c("U1546-thermo","U1546-abiotic"),
                             dD = c(-150,-100), 
                             d13C = c(-40,-20), 
                             D13CD = c(1.6,0.5), 
                             DDD = c(-2,-3))
# #U1547  
    
          em1[[3]] <- data.frame(describe = c("U1547-methano-lowG-70C","U1547-methano-70C"),
                            dD = c(-225,-225), 
                            d13C = c(-80,-80), 
                            D13CD = c(4.3,4.3), 
                            DDD = c(5,5))
          
          em2[[3]] <- data.frame(describe = c("U1547-thermo","U1547-abiotic"),
                            dD = c(-150,-100), 
                            d13C = c(-40,-20), 
                            D13CD = c(0.5,0.5), 
                            DDD = c(-3,-3))



################################################################################
  
for(i in 1:length(col)){
  num = nrow(em1[[i]])
  #calculate d13CD and dDD from deltas
  em1[[i]]$d13CD <- (((1+em1[[i]]$D13CD/1000)*(1+em1[[i]]$dD/1000)*(1+em1[[i]]$d13C/1000))-1) * 1000
  em1[[i]]$dDD <- (((1+em1[[i]]$DDD/1000)*(1+em1[[i]]$dD/1000)*(1+em1[[i]]$dD/1000))-1) * 1000
  
  em2[[i]]$d13CD <- (((1+em2[[i]]$D13CD/1000)*(1+em2[[i]]$dD/1000)*(1+em2[[i]]$d13C/1000))-1) * 1000
  em2[[i]]$dDD <- (((1+em2[[i]]$DDD/1000)*(1+em2[[i]]$dD/1000)*(1+em2[[i]]$dD/1000))-1) * 1000
  
  #for loop to go through all endmember combos
  plots <- list()
  plots_b <- list()
  plot_num <- c("a","b","c")
  
  for(end in 1:num){
    #endmember dD, d13C, d13CD, dDD
    e1 <- as.numeric(em1[[i]][end,c(2,3,6,7)])
    e2 <- as.numeric(em2[[i]][end,c(2,3,6,7)])
    
    #fractions
    f1 <- seq(0,1,0.1)
    f2 <- rev(f1)
    
    #multiple input by fraction -- then add the two endmembers
    fracs <- outer(f1,e1) + outer(f2,e2)
    colnames(fracs) <- colnames(em1[[i]][,c(2,3,6,7)])
    rownames(fracs) <- f1
    
    #calculate Deltas
    delta_13CD <- (((1+fracs[,3]/1000)/((1+fracs[,1]/1000)*(1+fracs[,2]/1000)))-1) *1000
    delta_12CDD <- (((1+fracs[,4]/1000)/((1+fracs[,1]/1000)*(1+fracs[,1]/1000)))-1) *1000
    mix_data <- data.frame(delta_13CD,delta_12CDD)
    
    clump_spec <-clump_spec_eq + 
      geom_errorbar(data= c_isotopes[sepa[[i]],],mapping = aes(D13CH3D, D12CH2D2, ymin=D12CH2D2-D12CH2D2.error, ymax=D12CH2D2+D12CH2D2.error, color = col[i]))+
      geom_errorbarh(data=c_isotopes[sepa[[i]],],mapping = aes(D13CH3D, D12CH2D2,xmin=D13CH3D-D13CH3D.error, xmax=D13CH3D+D13CH3D.error, color = col[i]))+
      geom_point(data = c_isotopes[sepa[[i]],], mapping = aes(D13CH3D, D12CH2D2, fill=col[i]),colour= "black",pch=23,size =5,stroke = 1.5)+
      theme_classic()+
      expand_limits(x=c(-5,7), y=c(-70, 30))+
      scale_x_continuous(expand = c(0,0))+
      ggtitle("Clumped  Isotopes")+
      ylab(expression(Delta^12*"CH"[2]*"D"[2]))+
      xlab(expression(Delta^13*"CH"[3]*"D"))+
      scale_fill_manual(values=col[i])+
      scale_color_manual(values=col[i])+
      geom_text(data = c_isotopes[sepa[[i]],],mapping = aes(D13CH3D, D12CH2D2, label = c(1,2,3)),color = "white",hjust = 0.5,vjust=0.5,
                size = 4)+
      theme(
        legend.position = "none"
      )
    
    trad_spec <- ggplot() + 
      geom_errorbar(data = isotopes[sepa[[i]],],mapping =aes(d13C, dD,ymin=dD-dD.error, ymax=dD+dD.error),color = col[i])+
      geom_errorbarh(data = isotopes[sepa[[i]],],mapping = aes(d13C, dD,xmin=d13C-d13C.error, xmax=d13C+d13C.error),color = col[i])+
      geom_point(data = isotopes[sepa[[i]],], mapping = aes(d13C, dD, fill = site), color = "black", pch=23,size =5,stroke = 1.5)+
      #geom_ellipse(data = isotopes,aes(x0=d13C, y0=dD, a = d13C.error, b = dD.error, angle = 0,fill=site, alpha = 0.05))+
      theme_classic()+
      #expand_limits(x=c(-120,0), y=c(-310,-100))+
      ggtitle("Bulk Isotopes")+
      ylab(expression(delta*"D"))+
      xlab(expression(delta^13*"C"))+
      scale_fill_manual(values=col[i])+
      scale_color_manual(values=col[i])+
      scale_shape_manual(values=c(23,23))+
      geom_text(data = isotopes[sepa[[i]],],mapping = aes(d13C, dD, label = c(1,2,3)),colour = "white",hjust = 0.5,vjust=0.5,
                size = 4)+
      theme(
        legend.position = "none"
      )+
      geom_polygon(regions, mapping= aes(biox,bioy), color = NA,fill = "orange",alpha = 0.2)+
      scale_x_continuous(expand = c(0,0))+
      scale_y_continuous(expand = c(0,0))+
      geom_polygon(regions, mapping= aes(thermx,thermy), color = NA,fill = "purple",alpha = 0.2)+
      geom_polygon(regions, mapping= aes(abiox,abioy), color = NA,fill = "grey",alpha = 0.2)
    
    
    #plot alongside data -- see Isotopes script
    labs <- paste0(names(em1[2])," = ",em2$dD[end], " / ", em1$dD[end]," \n ",
                   names(em1[3])," = ",em2$d13C[end], " / ",em1$d13C[end]," \n ",
                   names(em1[4])," = ",em2$D13CD[end], " / ",  em1$D13CD[end]," \n ",
                   names(em1[5])," = ", em2$DDD[end], " / ", em1$DDD[end])
    
    plots[[end]] <- clump_spec +  geom_path(data = mix_data, mapping = aes(delta_13CD, delta_12CDD),color="black",fill = "orange",size =1)+
      geom_point(data = mix_data, mapping = aes(delta_13CD, delta_12CDD),color="black",fill = "orange",size =2)+
      ggtitle(plot_num[end])
    #geom_rect(mapping = aes(xmin = 2.5, xmax = 7, ymin = -70, ymax = -40), fill = "grey",color = "black", alpha = 0.5)+
    #geom_text(aes(x=4.5,y=-55),label = labs,size=2)
    
    ##BULK##
    #make fracs matrix into a dataframe
    fracs_b <- as.data.frame(fracs)
    
    plots_b[[end]] <- trad_spec+  geom_path(data = fracs_b, mapping = aes(d13C,dD),color = "black",size =1)+
      geom_point(data = fracs_b, mapping = aes(d13C,dD),color = "black",size =2)+
      ggtitle(plot_num[end])
    #geom_rect(mapping = aes(xmin = -110, xmax = -90, ymin = -120, ymax = -90), fill = "grey",color = "black", alpha = 0.5)
    #geom_text(aes(x=-100,y=-105),label = labs_b)
  }
  
  #plot all on same page
  together <- ggarrange(plots_b[[1]],plots[[1]],
                        plots_b[[2]],plots[[2]],
                        ncol = 2,nrow = 2)
  
  final <- annotate_figure(together, top = text_grob(sites[i],color = col[i], face = "bold", size = 14 ))
  final
  #export as pdf
  #ggsave(paste("Mixing of", sites[i],".png"), final,width = 11, height = 9)
}
     

