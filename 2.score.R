library(tidyverse)
rm(list=ls())
setwd("workingDirectory/")
fgsea.clue.order <- readRDS("./result/cp/fgsea.clue.order_with_id.RDS")
#head(fgsea.clue.order )
#screening strategy1
# fgsea.clue.screen1 <- purrr::map(fgsea.clue.order,~ dplyr::filter(.x,padj<0.05,NES>0))
# fgsea.clue.screen2 <- purrr::map(fgsea.clue.screen1,~  .x$id)
# fgsea.clue.screen3 <- Reduce(intersect,fgsea.clue.screen2)
# good.id <- fgsea.clue.screen3

#screening strategy2
#Considering that certain drugs significantly enhance specific pathways but show slightly insufficient significance in one or two other pathways, to prevent the oversight of such drugs, we have implemented an alternative threshold screening allowing for some flexibility. Specifically, we allow an adjusted p-value between 0.05 and 0.2 for enrichment results in two pathways, while the adjusted p-value for other enriched pathways must be less than 0.05.
fgsea.clue.screen1 <- purrr::map(fgsea.clue.order,~ dplyr::filter(.x,padj<0.2,NES>0))
fgsea.clue.screen2 <- purrr::map(fgsea.clue.screen1,~  .x$id)
fgsea.clue.screen3 <- Reduce(intersect,fgsea.clue.screen2)
fgsea.clue.screen4 <- purrr::map(fgsea.clue.order,~ .x[as.integer(fgsea.clue.screen3),])
fgsea.clue.screen5 <- purrr::map(fgsea.clue.screen4 ,~ .x[,3])
fgsea.clue.screen6 <- Reduce(cbind,fgsea.clue.screen5)
fgsea.clue.screen6 <- cbind(fgsea.clue.screen6,fgsea.clue.screen4[[1]]$id)
fgsea.clue.screen6 <- fgsea.clue.screen6 %>% as.data.frame()
fgsea.clue.screen6.row <- purrr::map(as.data.frame(t(fgsea.clue.screen6[,1:5])), ~ .x)
fgsea.clue.screen6.row <- purrr::map(fgsea.clue.screen6.row,~ as.numeric(.))
fgsea.clue.screen6.count <- purrr::map(fgsea.clue.screen6.row,function(x){count <- 0;
for(i in 1:5){if(x[i]<0.05){count <- count+1}};return(count);count <- 0})
unlist(fgsea.clue.screen6.count) %>% as.data.frame() -> temp
colnames(fgsea.clue.screen6)[1:5] <- paste(colnames(fgsea.clue.screen6)[1],1:5,sep = "");
colnames(fgsea.clue.screen6)[6] <- "id"
fgsea.clue.screen6 <- fgsea.clue.screen6  %>% mutate(num_of_0.05=temp$.)
fgsea.clue.screen.end <- dplyr::filter(fgsea.clue.screen6,num_of_0.05>=4)
fgsea.clue.screen.id <- fgsea.clue.screen.end$id
good.id <- fgsea.clue.screen.id

##screen data
load("./template.RData")
list.good <- list()
for(i in 1:5){list.good[[i]] <- fgsea.clue.order[[i]][as.numeric(good.id),]}
good.df <- as.data.frame(matrix(data = NA,nrow =dim(list.good[[1]])[1],ncol = 5));
{for(i in 1:5){good.df[,i] <- list.good[[i]][,5]}}
colnames(good.df) <- names(fgsea.clue.order)
colnames(good.df)=rev(c('Reactive Oxygen Species', 'Lipid peroxidation', 'Iron metabolism', 'Glutathione','Fatty acid'))
good.score <- vector();

################################################################################
#Fedrug-Score
for(i in 1:dim(list.good[[1]])[1]){
  good.score[i] <-  -0.1718+
   0.3418*good.df[,'Reactive Oxygen Species'][i]+
   0.3776*good.df[,"Lipid peroxidation"][i]+
   0.3841*good.df[,"Iron metabolism"][i]+
   0.3314*good.df[,"Glutathione"][i]+
   0.354*good.df[,'Fatty acid'][i]}
FeDrug.condition <- readRDS("./Data/Compound_Data/FeDrug.condition")
good.df.withscore <- data.frame(good.score,good.df)
good.df.withscore <- data.frame(rownames(good.df.withscore),good.df.withscore)
colnames(good.df.withscore)[1] <- "good_id_index"
good.df.withscore.order <- good.df.withscore %>% arrange(desc(good.score))
good.id.order <- good.id[as.numeric(good.df.withscore.order$good_id_index)]
good.id.order.trt_number <- trt_cp_number[as.numeric(good.id.order)]
sig_info.choose <- FeDrug.condition[good.id.order.trt_number,]
sig_info.choose <- data.frame(1:dim(list.good[[1]])[1],sig_info.choose)
colnames(sig_info.choose)[1] <- "id"
sig_info.choose <- sig_info.choose %>% mutate(Pert_Score=good.df.withscore.order$good.score)
sig_info.choose.all=sig_info.choose 
sig_info.choose.fgsea <- cbind(sig_info.choose.all,good.df.withscore.order)

################################################################################
sig_info.choose<- sig_info.choose.all 
drug.choose.detail <- table(sig_info.choose$trt_iname) %>% as.data.frame()%>% arrange(desc(Freq))
colnames(drug.choose.detail)[1] <- "Compound"
colnames(drug.choose.detail)[2] <- "Times in cell"


library(dplyr)
pert_stats <- sig_info.choose %>%
  group_by(trt_iname) %>%
  summarise(
    pert_score_max = round(max(Pert_Score, na.rm = TRUE), 2),
  ) %>%
  rename(Compound = trt_iname)  
drug.choose.detail <- drug.choose.detail %>%
  left_join(pert_stats, by = "Compound")

fun_paste <- function(x,y){
  paste(x,y,sep="_")
}
drug.choose.detail <- purrr::map_df(drug.choose.detail,~ unlist(.)) 
drug.choose.detail[is.na( drug.choose.detail)] <- -666 #in original metadata -666 means NA

#save in excel file-type
library(openxlsx)
cpound=data.frame(table(FeDrug.condition$trt_iname))
colnames(cpound)=c("Compound","treat")
drug.choose.final=merge(drug.choose.detail,cpound,by="Compound")
openxlsx::write.xlsx(x = drug.choose.final , file = "results_with_Fedrug.xlsx",
                     sheetName = "screenResult", rownames = FALSE)


