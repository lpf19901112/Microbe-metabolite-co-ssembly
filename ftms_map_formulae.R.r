#为了顺利加载'xlsx'包，需先确认电脑上已成功安装java，然后依次下载如下包:
rm(list = ls())
library('readxl')
library('dplyr')
library('Hmisc')

setwd("H:/ftms/长期施肥数据/rawdata")

# compound category
FTICR_category = function(df){
  
  df = mutate(df,
              O_C = O/C,
              H_C = H/C,
              N_C = N/C,
              DBE_C = DBE/C,
              DBE_O = DBE/O,
              P_C = P/C,
              P_N = P/N,
              # https://www.nature.com/articles/s41467-018-02922-9
              # https://www.sciencedirect.com/science/article/pii/S0016703711000378
              # the net charge of OC speciation was assumed to be 0
              Gibbs = 60.3 - 28.5*(-((0+4*C+H-3*N-2*O+5*P-2*S)/C)+4), 
              # https://onlinelibrary.wiley.com/doi/full/10.1002/rcm.2386
              AI = ifelse((1+C-O-S-0.5*H) >= 0 & (C-O-N-S-P) > 0, (1+C-O-S-0.5*H)/(C-O-N-S-P),0))


  df = mutate(df,
              # https://pubs.acs.org/doi/10.1021/es405246a
              # https://pubs.acs.org/doi/10.1021/acs.est.6b03971
    group1 = ifelse(2/3  <  AI,                                                         'polycyclic aromatics originating from combustion',
             ifelse(0.5  <= AI  & AI <= 0.67,                                           'polyphenols derived from plants',
             ifelse(0.50 >= AI  & H_C < 1.5,                                            'unsaturated and phenolic compounds',
             ifelse(1.5  <= H_C & H_C <= 2.0,                                           'aliphatic compounds','others')))),
    
    group2 = ifelse(0    <= O_C & O_C <= 0.3  & 1.5 <= H_C & H_C <= 2,                  'lipids',
             ifelse(0.3  <= O_C & O_C <= 0.67  & 1.5 <= H_C & H_C <= 2.2 & N_C >= 0.05, 'protein_amino_sugar',
             ifelse(0.67 <= O_C & O_C <= 1.2  & 1.5 <= H_C & H_C <= 2,                  'carbohydrates',
             ifelse(0    <= O_C & O_C <= 0.1  & 0.7 <= H_C & H_C <= 1.5,                'unsatured_hydrocarbons',
             ifelse(0.1  <= O_C & O_C <= 0.67 & 0.7 <= H_C & H_C <= 1.5 & AI <= 0.67,   'ligins',
             ifelse(0.67 <= O_C & O_C <= 1.2  & 0.5 <= H_C & H_C <= 1.5 & AI <= 0.67,   'tannins',
             ifelse(0    <= O_C & O_C <= 0.67 & 0.2 <= H_C & H_C <= 0.7 & AI >= 0.67,   'condensed_aromatics','others'))))))))
  df
}

# a blank data.frame to save result of mass and compound formula
result0 = data.frame()
#tempdat0 = read_excel("1.xlsx", sheet = 1)
#for (smp in seq(1,71,3)){
#list.files(path = "./00_data")
### now batch process
for (smp in seq(54)){

  #(1) input data and adjust some variables
  
   tempdat = read_excel(paste0(smp,'.xlsx'), sheet = 1, 
                        col_types = c(rep("text",2),rep("numeric",12),rep("text",1),rep("numeric",3))) %>%
     
     subset(is.na(category) == FALSE) %>%
    
    # rename some variables, rename
    rename(iso = 'isotope kind',mz = "m/z",  sn = "s/n") %>%
    
    # generate a new variable to filter unwanted peaks
    mutate(OP = O - P,P_C = P/C, P_N = P/N) %>% 
    
    # if the number of H is odd, then plus 1, and add 1.007825 to the corresponding mz
    mutate(mz = ifelse((H %% 2) == 0, mz, mz + 1.007825), 
           H = ifelse((H %% 2) == 0, H, H+1)) %>%
     
    # we do remove all isotope peaks
    subset(is.na(iso) == TRUE) %>%
    
     # we do remove peaks that containing 13C atom
    select(category,iso, C, H, N, O, S, P, mz, intensity, sn, ppm, DBE, OP, P_C,P_N, KMD)  
  
  
  # (2) assign formula to peaks
  # https://doi.org/10.1002/2017JG003967
  # Chemical formulae were assigned based on the following criteria: 
  # (1) S/N > 7, 
  # (2) mass measurement error < 1 ppm, taking into consideration the presence of C, H, O, N, S, and P and excluding other elements. 
  # (3) assignment of 1 P atoms requires the presence of at least 4 O atoms.
  
  # replace NA with 4
  tempdat$OP = impute(tempdat$OP, 4)
  
  
  assigned = tempdat %>% 
    subset(sn>=6 & abs(ppm)<1 & OP >=4) %>%  
    select(-c(sn,ppm,OP))
  

  for (j in 1:dim(assigned)[1]) {
    
    ###提取每一行不含NA的列
    temp = assigned[j,c('C','H','N','O','S','P')] %>% .[,colSums(is.na(.)) == 0]
    
    ###将原子和对应原子数粘合作为chemical formula
    assigned$ID[j] = paste0(colnames(temp), as.matrix(temp), collapse = '') 
  }
  
  # those peaks that are not assigned chemical formula
  unassigned = tempdat[!(rownames(tempdat) %in% rownames(assigned)),] %>% select(-c(sn,ppm,OP))
  
  unassigned$ID = paste0('unassigned_',seq(nrow(unassigned)))
  
  tempdat = rbind(assigned, unassigned)
  
  
  # replace NA with 0
  tempdat$N <- impute(tempdat$N,0)
  tempdat$S <- impute(tempdat$S,0)
  tempdat$P <- impute(tempdat$P,0)
  
 
  # (3) generate some other indices
  
  tempdat = FTICR_category(tempdat) %>% 
    mutate(Sample = paste('S',smp,sep = '')) %>%
    subset(select = c(category,iso,AI,KMD,DBE,DBE_C,DBE_O,Gibbs,mz,intensity,C, H, N, O, S, P,O_C,H_C,N_C,P_C,P_N,ID,group1,group2,Sample))

  result0 = rbind(result0,tempdat)
  
  
  # remove those temperal variables
  remove(tempdat,assigned, unassigned)
  
  
  # monitor the batch process
  print(c(date(),smp))
}

write.csv(result0,"formula0.csv")
# rename their unassigned peaks that has same mz value but different name among samples
unassigned = result0 %>% subset(ID %in% grep(pattern = 'unassigned', result0$ID, value = T))
assigned = result0 %>% subset(!(ID %in% grep(pattern = 'unassigned', result0$ID, value = T)))
unassigned_new = data.frame(unique(unassigned$mz)) %>% rename(mz = "unique.unassigned.mz.")
unassigned_new$newID = paste0('unassigned_',1:nrow(unassigned_new))
head(unassigned_new)
unassigned = merge(unassigned, unassigned_new, by = 'mz') %>% 
  select(-ID) %>% 
  rename(ID = newID) %>%
  unique.data.frame()

result0 = rbind(unassigned, assigned)

max(result0$O_C)
hist(result0$O_C)
result0 = subset(result0, O_C <=1.2)
write.csv(result0,'formulae.csv')


# for component classification
library(reshape2)
(group1 = dcast(data = result0, Sample ~ group1, length))
(group2 = dcast(data = result0, Sample ~ group2, length))
write.csv(group1, 'group1.csv')
write.csv(group2, 'group2.csv')



(peaks_per_smp = table(result0$Sample) %>% as.data.frame() %>% rename(sample = Var1, peaks = Freq) %>% arrange(peaks))
write.csv(peaks_per_smp, 'peaks_per_sample.csv', row.names = F)

sink('peaks_sta.txt')

"summary(unique(result0$mz))"
summary(unique(result0$mz))
cat('\n\n\n')


"total quality peaks"
length(unique(result0$ID))

"assigned peaks"
length(unique(assigned$ID))
cat('\n\n\n')

subset(result0, select = c(ID, category)) %>% 
  unique.data.frame() %>%
  with(table(category))/length(unique(result0$ID))
cat('\n\n\n')

subset(result0, select = c(ID, group1)) %>% 
  unique.data.frame() %>% 
  with(table(group1))/length(unique(result0$ID))
cat('\n\n\n')

subset(result0, select = c(ID, group2)) %>% 
  unique.data.frame() %>% 
  with(table(group2))/length(unique(result0$ID))

sink()



# calculate gibbs free energy
(gibbs = aggregate(Gibbs ~ Sample, median, data = result0) %>% arrange(Gibbs))
write.csv(gibbs,'gibbs.csv')




# (optional) discard those peaks that only observed only once within a Treatment
map = read.table('map.txt', row.names = 1, header = T)
result6 = data.frame()

for (trt in unique(map$Treatment)){
  
  submap = subset(map, Treatment == trt)
  
  temp = subset(result0, Sample %in% rownames(submap)) %>%
    select(ID) %>% 
    table() %>%
    data.frame() %>%
    rename(ID = '.') %>% 
    subset(Freq >= 6)
  
  result6 = rbind(result6,subset(result0, Sample %in% rownames(submap) & ID %in% temp$ID))
  remove(submap, temp)
  
}

library(mice)
md.pattern(result6)

write.csv(result6,file = './peaks exist at least 6 within treatment/formulae.csv')

(peaks_per_smp_at_least_6_within_treatment = table(result6$Sample) %>% as.data.frame() %>% rename(sample = Var1, peaks = Freq) %>% arrange(peaks))
write.csv(peaks_per_smp_at_least_6_within_treatment, 
          './peaks exist at least 6 within treatment/peaks_per_sample.csv', row.names = F)


sink('./peaks exist at least 6 within treatment/peaks_sta.txt')

"summary(unique(result6$mz))"
summary(unique(result6$mz))
cat('\n\n\n')


"total quality peaks"
length(unique(result6$ID))
"assigned peaks"
length(unique(subset(result6, ID %in% assigned$ID)$ID))
cat('\n\n\n')

subset(result6, select = c(ID, category)) %>% 
  unique.data.frame() %>%
  with(table(category))/length(unique(result6$ID))
cat('\n\n\n')

subset(result6, select = c(ID, group1)) %>% 
  unique.data.frame() %>% 
  with(table(group1))/length(unique(result6$ID))
cat('\n\n\n')

subset(result6, select = c(ID, group2)) %>% 
  unique.data.frame() %>% 
  with(table(group2))/length(unique(result6$ID))

sink()


# calculate gibbs free energy
(gibbs_6 = aggregate(Gibbs ~ Sample, median, data = result6) %>% arrange(Gibbs))
write.csv(gibbs_6,'./peaks exist at least 6 within treatment/gibbs.csv')

