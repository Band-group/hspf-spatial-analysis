#Predictive scores

#basic packages and parallel computing packages (add more if needed)
list.of.packages <- c("ROCR","dplyr","tibble")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

results <- list()
for (l in 1:length(Pfalleles)){
#l=1
# Use lapply to get a list of file names for each pattern
predmodels <- lapply(Pfalleles[l], function(p) {
  list.files(
    path = "output/csv/",
    pattern = paste0("Pfoutput.*", p, ".*\\.csv$"),
    full.names = TRUE
  )
})

# Flatten the list into a single vector
predmodels <- unlist(predmodels)
predall <- lapply(predmodels, FUN = read.csv)
predall <- as.data.frame(dplyr::bind_rows(predall))
rownames(predall)<-1:nrow(predall)
predall$model <-as.factor(predall$model)
predall$region <-as.factor(predall$region)
predall$country <-as.factor(predall$country)
predall$pf <- as.factor(Pfalleles[l])
#levels(predall$model)<-c("aspatial","spatial",levels(predall$model)[-c(1,length(levels(predall$model)))]) 

# Compute metrics for each model
#compute for regional model

# Conditional computation for 'alt_hat.mean' and 'alt_hat.sd'
if ('alt_hat.mean' %in% colnames(predall)) {
  regionalresults <- predall[predall$model=='regional',] %>% 
    group_by(model,region,pf) %>%
    summarise(
      r0 = round(mean(r0,na.rm=TRUE),1),
      sigma0 = round(mean(sigma0),1),
      mean_intercept = round(mean(intercept.mean,na.rm=TRUE),2),
      sd_intercept = round(mean(intercept.sd,na.rm=TRUE),2),
      mean_beta_HbS = round(mean(HbS_hat.mean,na.rm=TRUE),2),
      sd_beta_HbS = round(mean(HbS_hat.sd,na.rm=TRUE),2),
      mean_beta_alt = round(mean(alt_hat.mean, na.rm = TRUE), 2),
      sd_beta_alt = round(mean(alt_hat.sd, na.rm = TRUE), 2),
      MAE = round(mean(abs(obs - pred)) * 100, 3),
      RMSE = round(sqrt(mean((obs - pred)^2)) * 100, 3),
      LogLoss = round(-mean(obs * log(pred) + (1 - obs) * log(1 - pred)) * 100, 3),
      BSref = round( (mean(obs) * (1 - mean(obs)))*100,3),
      BSS = {
        1 - ((mean((obs - pred)^2)) / (mean(obs) * (1 - mean(obs))) )
      }* 100 %>% round(., 3)
    )} else {
      regionalresults <- predall[predall$model=='regional',] %>% 
        group_by(model,region,pf) %>%
        summarise(
      r0 = round(mean(r0,na.rm=TRUE),1),
      sigma0 = round(mean(sigma0),1),
      mean_intercept = round(mean(intercept.mean,na.rm=TRUE),2),
      sd_intercept = round(mean(intercept.sd,na.rm=TRUE),2),
      mean_beta_HbS = round(mean(HbS_hat.mean,na.rm=TRUE),2),
      sd_beta_HbS = round(mean(HbS_hat.sd,na.rm=TRUE),2),
      MAE = round(mean(abs(obs - pred)) * 100, 3),
      RMSE = round(sqrt(mean((obs - pred)^2)) * 100, 3),
      LogLoss = round(-mean(obs * log(pred) + (1 - obs) * log(1 - pred)) * 100, 3),
      BSref = round( (mean(obs) * (1 - mean(obs)))*100,3),
      BSS = {
        1 - ((mean((obs - pred)^2)) / (mean(obs) * (1 - mean(obs))) )
      }* 100 %>% round(., 3)
      ) 
    }
 
#compute for country models
if ('alt_hat.mean' %in% colnames(predall)) {
  countryresults <- predall[predall$model=='country',] %>% 
    group_by(model,country,pf) %>%
    summarise(
      r0 = round(mean(r0,na.rm=TRUE),1),
      sigma0 = round(mean(sigma0),1),
      mean_intercept = round(mean(intercept.mean,na.rm=TRUE),2),
      sd_intercept = round(mean(intercept.sd,na.rm=TRUE),2),
      mean_beta_HbS = round(mean(HbS_hat.mean,na.rm=TRUE),2),
      sd_beta_HbS = round(mean(HbS_hat.sd,na.rm=TRUE),2),
      mean_beta_alt = round(mean(alt_hat.mean, na.rm = TRUE), 2),
      sd_beta_alt = round(mean(alt_hat.sd, na.rm = TRUE), 2),
      MAE = round(mean(abs(obs - pred)) * 100, 3),
      RMSE = round(sqrt(mean((obs - pred)^2)) * 100, 3),
      LogLoss = round(-mean(obs * log(pred) + (1 - obs) * log(1 - pred)) * 100, 3),
      BSref = round( (mean(obs) * (1 - mean(obs)))*100,3),
      BSS = {
        1 - ((mean((obs - pred)^2)) / (mean(obs) * (1 - mean(obs))) )
      }* 100 %>% round(., 3)
    )
  
} else {
  countryresults <- predall[predall$model=='country',] %>% 
    group_by(model,country,pf) %>%
    summarise(
      r0 = round(mean(r0,na.rm=TRUE),1),
      sigma0 = round(mean(sigma0),1),
      mean_intercept = round(mean(intercept.mean,na.rm=TRUE),2),
      sd_intercept = round(mean(intercept.sd,na.rm=TRUE),2),
      mean_beta_HbS = round(mean(HbS_hat.mean,na.rm=TRUE),2),
      sd_beta_HbS = round(mean(HbS_hat.sd,na.rm=TRUE),2),
      MAE = round(mean(abs(obs - pred)) * 100, 3),
      RMSE = round(sqrt(mean((obs - pred)^2)) * 100, 3),
      LogLoss = round(-mean(obs * log(pred) + (1 - obs) * log(1 - pred)) * 100, 3),
      BSref = round( (mean(obs) * (1 - mean(obs)))*100,3),
      BSS = {
        1 - ((mean((obs - pred)^2)) / (mean(obs) * (1 - mean(obs))) )
      }* 100 %>% round(., 3)
    )
}

names(countryresults)[names(countryresults) == "country"] <- "region"
results[[l]] <- rbind(regionalresults,countryresults)
}

myresult <- do.call(rbind, results)
# Export to latex and csv
write.table(x = myresult, 
            file = "output/predscores.tex", 
            sep = "&",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            eol = " \\\\ \n",
            na = "",
            dec = ".")
write.csv(x = myresult, 
            file = "output/predscores.csv", 
            row.names = FALSE,
            quote = FALSE,
            na = "")

#save proportion of HbS mean estimates > 0 
library(dplyr)
HbSpos <- myresult[c("model","region","pf","mean_beta_HbS")] %>%
  group_by(pf,model,region) %>%
  summarize(HbSpos = round((sum(mean_beta_HbS > 0)/n()),2))
HbSposlatex <- xtable::xtable(data.frame(HbSpos))
print(HbSposlatex, file=paste0("output/csv/HbSpos_allmodels_latex.txt"))



for (l in 1:length(Pfalleles)){

#plot HbS estimate by model
library(ggplot2)
library(ggsci)
HbSp2 <- ggplot(myresult[myresult$pf==Pfalleles[l],], aes(model, mean_beta_HbS,color=region)) +
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6,
    .width = 0,
    justification = -.3,
    point_colour = NA) +
  #scale_color_npg()+
  geom_point(
    size = 0.8,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )) +
      geom_boxplot(
        width = .25,
        outlier.shape = NA,
        color='black',fill='transparent'
      ) +
   guides(color = "none")+
  coord_cartesian(xlim = c(1.2, NA), clip = "off")+
  ggthemes::theme_few(22)+ theme(panel.border=element_rect(linewidth = 0.3))+
  ylab(bquote(hat(beta)[HbS]))+
  xlab("Model")
ggsave(HbSp2,file=paste0("output/pdf/HbScoef_distribution",Pfalleles[l],".pdf"),width=17,height=4)

message(paste0("\nEND Pf_predscores.R for ", Pfalleles[l]))
}
#END predictive scores

#Make graph for Figure 2
# Define the color palette
country_colors <- c("Gambia" = "#0000cd", "Mali" = "#42426f", "Ghana" = "#03b4cd", "DRC" = "#2E8B57", "Tanzania" = "#ee5c42", "Ethiopia" = "#ee5500")
region_colors <- c(
  "West Africa" = "#0E4C92",   #Yale Blue; Royal Blue: "#4169E1"
  "East Africa" = "#DA680F",    #Burgundyred#8D021F, Orangered: #D9534F
  "Africa" = "grey35"           # Dark grey
)
#define region and country levels for wrap plots
alevels <- c("East Africa","West Africa","Gambia","Mali","Ghana","DRC","Tanzania")

modnames <- list('country','regional')
allpalettes <- list(country_colors,region_colors)
for (l in 1:length(Pfalleles)){
  #plot hbs effects for both country and regional models
  hbsdata <- myresult[(myresult$pf==Pfalleles[l]),]
  hbsdata$region <- factor(hbsdata$region,levels=rev(alevels))
  # Create a ggplot graph
  hbshatplot <- ggplot(hbsdata, aes(x = mean_beta_HbS, y = region)) +
    geom_pointrange(aes(xmin = mean_beta_HbS - 1.96*sd_beta_HbS, xmax = mean_beta_HbS + 1.96*sd_beta_HbS),
                    color = "black", alpha = 1,linewidth=1.25) +
    geom_point(aes(color = region), size = 6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey35",linewidth=1)+
    scale_color_manual(values = unlist(allpalettes)) +
    labs(title = paste0("Estimated HbS effects on ",Pfalleles[l]," (log odd-ratios)")) +
    theme_bw(16) +
    theme(legend.position = "none",axis.title=element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  ggsave(file=paste0("output/pdf/HbShat_allmodels","_", Pfalleles[l],".pdf"),hbshatplot,dpi = 100,width = 14,height = 5)
  ggsave(file=paste0("output/svg/HbShat_allmodels","_", Pfalleles[l],".svg"),hbshatplot,width = 14,height = 5)
  
  
  #plot hbs effects for country models
 for (i in 1:length(modnames)){
  hbsdata <- myresult[(myresult$model==modnames[i] & myresult$pf==Pfalleles[l]),]
  if (modnames[i]=='regional'){
  rlevels <- c("All","West Africa","East Africa")
  hbsdata$region <- factor(hbsdata$region,levels=rev(rlevels))
 } 
  if (modnames[i]=='country'){
  clevels <- c("All","Gambia","Mali","Ghana","DRC","Tanzania")
  hbsdata$region <- factor(hbsdata$region,levels=rev(clevels))
 }
# Create a ggplot graph
hbshatplot <- ggplot(hbsdata, aes(x = mean_beta_HbS, y = region)) +
  geom_pointrange(aes(xmin = mean_beta_HbS - 1.96*sd_beta_HbS, xmax = mean_beta_HbS + 1.96*sd_beta_HbS),
                  color = "black", alpha = 1,linewidth=1.25) +
  geom_point(aes(color = region), size = 6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey35",linewidth=1)+
   scale_color_manual(values = unlist(allpalettes[i])) +
  labs(title = paste0("Estimated HbS effects on ",Pfalleles[l]," (log odd-ratios)")) +
  theme_bw(16) +
  theme(legend.position = "none",axis.title=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
ggsave(file=paste0("output/pdf/HbShat_",modnames[i],"_", Pfalleles[l],".pdf"),hbshatplot,dpi = 100,width = 14,height = 4)
ggsave(file=paste0("output/svg/HbShat_",modnames[i],"_", Pfalleles[l],".svg"),hbshatplot,width = 14,height = 4)
  }
}
# #plot hbs effects for regional models
# hbsdata <- predall[(predall$model=='regional' & predall$pf==Pfalleles[l]),]
# # Create a ggplot graph
# hbshatplot <- ggplot(hbsdata, aes(x = HbS_hat.mean, y = region)) +
#   geom_pointrange(aes(xmin = HbS_hat.mean - HbS_hat.sd, xmax = HbS_hat.mean + HbS_hat.sd), color = "black", alpha = 0.5) +
#   geom_point(aes(color = region), size = 3) +
#   scale_color_manual(values = region_colors) +
#   labs(title = paste0("Estimated HbS effects on ",Pfalleles[l]," (log odd-ratios)")) +
#   theme_minimal(22) +
#   theme(legend.position = "none",axis.title=element_blank())
# ggsave(file=paste0("output/pdf/HbShatregion_",Pfalleles[l],".pdf"),hbshatplot,dpi = 100,width = 16,height = 7)
# ggsave(file=paste0("output/svg/HbShatregion_",Pfalleles[l],".svg"),hbshatplot,width = 16,height = 7)
# }
message(paste0("\nEND Pf_predscores.R for all Pf alleles"))
gc()  
