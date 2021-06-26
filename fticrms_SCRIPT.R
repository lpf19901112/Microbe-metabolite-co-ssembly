rm(list = ls())
setwd("H:/2019连续标记/高分辨质谱结果/assemble/rawdata")
library('xlsx')
library('dplyr')
library('Hmisc')
#library('plyr')
library('reshape2')
# compound category
FTICR_category = function(df){
  
  df = mutate(df,
              O_C = O/C,
              H_C = H/C,
              N_C = N/C,
              DBE_C = DBE/C,
              # https://www.nature.com/articles/s41467-018-02922-9
              # https://www.sciencedirect.com/science/article/pii/S0016703711000378
              Gibbs = 60.3 - 28.5*(-((1+4*C+H-3*N-2*O+5*P-2*S)/C)+4), 
              # https://onlinelibrary.wiley.com/doi/full/10.1002/rcm.2386
              AI = ifelse((1+C-O-S-0.5*H) >= 0 & (C-O-N-S-P) > 0, (1+C-O-S-0.5*H)/(C-O-N-S-P),0))
  
  
  df = mutate(df,
              # https://pubs.acs.org/doi/10.1021/es405246a
              # https://pubs.acs.org/doi/10.1021/acs.est.6b03971
              group = ifelse(0    <= O_C & O_C <= 0.3  & 1.5 <= H_C & H_C <= 2,     'lipids',
                             ifelse(0.3  <= O_C & O_C <= 2/3  & 1.5 <= H_C & H_C <= 2.2,   'protein_amino_sugar',
                                    ifelse(2/3  <= O_C & O_C <= 1.2  & 1.5 <= H_C & H_C <= 2,     'carbohydrates',
                                           ifelse(0    <= O_C & O_C <= 0.1  & 0.7 <= H_C & H_C <= 1.5,   'unsatured_hydrocarbons',
                                                  ifelse(0.1  <= O_C & O_C <= 2/3  & 0.7 <= H_C & H_C <= 1.5,   'ligins',
                                                         ifelse(2/3  <= O_C & O_C <= 1.2  & 0.5 <= H_C & H_C <= 1.5,   'tannins',
                                                                ifelse(0    <= O_C & O_C <= 2/3  & 0.2 <= H_C & H_C <= 0.7,   'condensed_aromatics','others'))))))))
  df
}

# a blank data.frame to save result of mass and compound formula
result = data.frame()

library('stringr')

### now batch process
for (smp in list.files(path = "../rawdata")){
  
  #(1) input data and adjust some variables
  
  tempdat = read.xlsx(paste0('../rawdata/',smp),sheetIndex = 1) %>%
    # tempdat = read.xlsx(paste0('../01_data/1','.xlsx'),sheetIndex = 1) %>%
    
    # rename some variables, rename
    rename(iso = "isotope.kind", mz = "m.z",  sn = "s.n")  %>% 
    
    # generate a new variable to filter unwanted peaks
    mutate(OP = O-P) %>% 
    
    # if the number of H is odd, then plus 1, and add 1.007825 to the corresponding mz
    mutate(mz = ifelse((H %% 2) == 0, mz, mz + 1.007825), 
           H = ifelse((H %% 2) == 0,H,H+1)) %>%
    
    select(category, iso, C, H, N, O, S, P, mz, sn, ppm, DBE, OP,int) 
  
  
  # (2) assign formula to peaks
  # https://doi.org/10.1002/2017JG003967
  # Chemical formulae were assigned based on the following criteria: 
  # (1) S/N > 7, 
  # (2) mass measurement error < 1 ppm, taking into consideration the presence of C, H, O, N, S, and P and excluding other elements. 
  # (3) assignment of 1 P atoms requires the presence of at least 4 O atoms.
  
  # replace NA with 4
  tempdat$OP = impute(tempdat$OP, 4)
  
  
  assigned = tempdat %>% 
    subset(sn>7 & abs(ppm)<1 & OP>3) %>%  
    select(-c(sn,ppm,OP))
  
  
  for (j in 1:dim(assigned)[1]) {
    
    ###提取每一行不含NA的列
    temp = assigned[j,c('iso','C','H','N','O','S','P')] %>% .[,colSums(is.na(.)) == 0]
    
    ###将原子和对应原子数粘合作为chemical formula
    assigned$ID[j] = paste0(colnames(temp), as.matrix(temp), collapse = '') 
  }
  
  # those peaks that are not assigned chemical formula
  unassigned = tempdat[!(rownames(tempdat) %in% rownames(assigned)),] %>% select(-c(sn,ppm,OP))
  
  unassigned$ID = paste0('unassigned_',seq(nrow(unassigned)))
  
  tempdat = rbind(assigned, unassigned) %>% select(-(iso))
  
  # (3) compound category
  
  # replace NA with 0
  tempdat$N <- impute(tempdat$N,0)
  tempdat$S <- impute(tempdat$S,0)
  tempdat$P <- impute(tempdat$P,0)
  
  
  # compound category
  smp = gsub('.xlsx','',smp)
  
  
  # (4) generate some other indices
  
  tempdat = FTICR_category(tempdat) %>% 
    mutate(Sample = paste('S',smp,sep = '')) %>%
    subset(select = c(category,AI,DBE_C,Gibbs,mz,O_C,H_C,N_C,ID,group,Sample,int))
  
  result = rbind(result,tempdat)
  
  
  # remove those temperal variables
  remove(tempdat,temp_category,assigned, unassigned)
  
  
  # monitor the batch process
  print(c(date(),smp))
}
# varify if there were NA
library('mice')
md.pattern(result)
result = na.omit(result)
# discard Cl-containing peaks
filter_Cl_result = subset(result, !(category %in% grep('Cl',result$category,value = T)))

write.csv(filter_Cl_result,'formulae.csv')
# for component classification
classification = dcast(data = filter_Cl_result,Sample ~ group, length)
md.pattern(classification)
write.csv(classification, file = 'classification.csv')
library('plyr')
# calculate gibbs free energy
gibbs = filter_Cl_result %>% 
  select(Sample,Gibbs) %>% 
  ddply(.(Sample), 
        summarize,
        mean = mean(Gibbs),
        median = median(Gibbs)) %>% melt()

write.csv(gibbs,'gibbs.csv')

library('dplyr')
unassigned = subset(filter_Cl_result, ID %in% grep(pattern = 'unassigned', filter_Cl_result$ID, value = T))
assigned = subset(filter_Cl_result, !(ID %in% grep(pattern = 'unassigned', filter_Cl_result$ID, value = T)))
unassigned_new = data.frame(unique(unassigned$mz)) %>% rename(c("unique.unassigned.mz." = "mz"))
unassigned_new$newID = paste0('unassigned_',1:nrow(unassigned_new))
head(unassigned_new)
unassigned = merge(unassigned, unassigned_new, by = 'mz') %>% select(-ID) %>% rename(c("newID" = "ID"))

new_result0 = rbind(unassigned, assigned)
write.csv(new_result0,"final_formulae.csv")
df = subset(new_result0, select = c(int,ID,Sample,group,O_C,H_C))
aq_dat=melt(df, id=c("ID","Sample","group","O_C","H_C"), na.rm=TRUE)
#write.csv(aq_dat,"aq_dat.csv")


k <- read.csv("formula1.csv")
x=dcast(aq_dat,  group + ID + O_C + H_C ~ Sample, fill = 0)
x=dcast(k,  mz + category + AI + KMD + DBE + DBE_C + DBE_O + Gibbs + intensity + C + H + N + O + S + P + O_C + H_C + N_C + ID + group1 + group2 ~ Sample, value.var="intensity",fun.aggregate = toString)
rownames(x)=x$mz
write.csv(x,file = "final_data1.csv")
data1=x[,-1:-4]
data1[data1>0]=1
data2=data1[rowSums(data1)>1,]
idx=rownames(data2) %in% rownames(x)
sub_data2=data2[idx,]
final_data=x[rownames(sub_data2),]
write.csv(final_data,file = "final_data3.csv")

datak1 <- read.csv("otu_forest_20000.csv", row.names = 1, header = T)
datak2 <- datak1[,-1:-4]
datak2[datak2>0]=1
datak3=datak2[rowSums(datak2)>1,]
idxK=rownames(datak3) %in% rownames(datak1)
sub_datak3=datak3[idxK,]
final_dataL=datak1[rownames(sub_datak3),]
write.csv(datak2,file = "binary.csv")

write.csv(datak3, "binary_data2.csv")

library(vegan)
library(ape)
ds <- vegdist(t(datak3), method = "binomial")
pc <- pcoa(ds)
ve <- pc$vectors
plot(ve[,1:2])
write.csv(ve, "vectors_binomial.csv")


variable <- read.csv("group.csv", row.names = 1, header = T)
adonis(ds ~ Sampletime*Year, data = variable, method = "binomial")#bray,euclidean
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'euclidean', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 鈥?***鈥? 0.001 鈥?**鈥? 0.01 鈥?*鈥? 0.05 鈥?.鈥? 0.1 鈥? 鈥? 1")
  return(pairw.res)
  
} 

pairwise.adonis(t(datak3),variable$detail)
