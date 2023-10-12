#Skeletal Muscle Proteomic Analysis
#Written by: Sean Perez

#Import Libraries
library(dplyr)
library(tidyverse)
library(gplots)
library(ggpubr)

########################################################## Data Read In ########################################################## 
########Control Exercise Data####################################
#FEvsME
names(DeqMS_FEvsME)[4] <- 'FDR'
DeqMS_FEvsME <- filter(DeqMS_FEvsME, FDR<0.05)
write(DeqMS_FEvsMEC$`Gene Name`, file = '[YOUR PATH HERE]/FEvsMEgenes.txt')
#Go to https://www.uniprot.org/id-mapping and load txt file in box. From database: Gene Name, To database: UniProtKB/Swiss-Port, Taxonomy: Homo Sapien
#Download results as uncompressed Excel file with only "Gene Names" and "Protein Names" selected
#Import Data set as Uniprot_FEvsME
Uniprot_FEvsME <- cbind(DeqMS_FEvsME[match(Uniprot_FEvsME$...1, DeqMS_FEvsME$...1),], Uniprot_FEvsME$Accession, Uniprot_FEvsME$`Protein Name`)
rownames(Uniprot_FEvsME) <- 1:nrow(Uniprot_FEvsME)
names(Uniprot_FEvsME)[6] <- 'Accession'
names(Uniprot_FEvsME)[7] <- 'Protein.names'
Uniprot_FEvsME <- Uniprot_FEvsME[,c("...1", "Accession", "Gene Name", "LogFC", "FDR", "Direction", "Protein.names")]


#FCvsMC
names(DeqMS_FCvsMC)[4] <- 'FDR'
DeqMS_FCvsMC <- filter(DeqMS_FCvsMC, FDR<0.05)
write(DeqMS_FCvsMC$`Gene Name`, file = '[YOUR PATH HERE]/FCvsMCgenes.txt')
#Go to https://www.uniprot.org/id-mapping and load txt file in box. From database: Gene Name, To database: UniProtKB/Swiss-Port, Taxonomy: Homo Sapien
#Download results as uncompressed Excel file with only "Gene Names (primary)" and "Protein Names" selected
#Import Data set as Uniprot_FCvsMC
Uniprot_FCvsMC <- cbind(DeqMS_FCvsMC[match(Uniprot_FCvsMC$...1, DeqMS_FCvsMC$...1),], Uniprot_FCvsMC$Accession, Uniprot_FCvsMC$`Protein Name`)
rownames(Uniprot_FCvsMC) <- 1:nrow(Uniprot_FCvsMC)
names(Uniprot_FCvsMC)[6] <- 'Accession'
names(Uniprot_FCvsMC)[7] <- 'Protein.names'
Uniprot_FCvsMC <- Uniprot_FCvsMC[,c("...1", "Accession", "Gene Name", "LogFC", "FDR", "Direction", "Protein.names")]

#MEvsMS
names(DeqMS_MEvsMS)[4] <- 'FDR'
DeqMS_MEvsMS <- filter(DeqMS_MEvsMS, FDR<0.05)
write(DeqMS_MEvsMS$`Gene Name`, file = '[YOUR PATH HERE]/MEvsMSgenes.txt')
#Go to https://www.uniprot.org/id-mapping and load txt file in box. From database: Gene Name, To database: UniProtKB/Swiss-Port, Taxonomy: Homo Sapien
#Download results as uncompressed Excel file with only "Gene Names (primary)" and "Protein Names" selected
#Import Data set as Uniprot_MEvsMS
Uniprot_MEvsMS <- cbind(DeqMS_MEvsMS[match(Uniprot_MEvsMS$...1, DeqMS_MEvsMS$...1),], Uniprot_MEvsMS$Accession, Uniprot_MEvsMS$`Protein Name`)
rownames(Uniprot_MEvsMS) <- 1:nrow(Uniprot_MEvsMS)
names(Uniprot_MEvsMS)[6] <- 'Accession'
names(Uniprot_MEvsMS)[7] <- 'Protein.names'
Uniprot_MEvsMS <- Uniprot_MEvsMS[,c("...1", "Accession", "Gene Name", "LogFC", "FDR", "Direction", "Protein.names")]

#FEvsFC
names(DeqMS_FEvsFC)[4] <- 'FDR'
DeqMS_FEvsFC <- filter(DeqMS_FEvsFC, FDR<0.05)
write(DeqMS_FEvsFC$`Gene Name`, file = '[YOUR PATH HERE]/FEvsFCgenes.txt')
#Go to https://www.uniprot.org/id-mapping and load txt file in box. From database: Gene Name, To database: UniProtKB/Swiss-Port, Taxonomy: Homo Sapien
#Download results as uncompressed Excel file with only "Gene Names (primary)" and "Protein Names" selected
#Import Data set as Uniprot_FEvsFC
Uniprot_FEvsFC <- cbind(DeqMS_FEvsFC[match(Uniprot_FEvsFC$...1, DeqMS_FEvsFC$...1),], Uniprot_FEvsFC$Accession, Uniprot_FEvsFC$`Protein Name`)
rownames(Uniprot_FEvsFC) <- 1:nrow(Uniprot_FEvsFC)
names(Uniprot_FEvsFC)[6] <- 'Accession'
names(Uniprot_FEvsFC)[7] <- 'Protein.names'
Uniprot_FEvsFC <- Uniprot_FEvsFC[,c("...1", "Accession", "Gene Name", "LogFC", "FDR", "Direction", "Protein.names")]

#MEvsMC
names(DeqMS_MEvsMC)[4] <- 'FDR'
DeqMS_MEvsMC <- filter(DeqMS_MEvsMC, FDR<0.05)
write(DeqMS_MEvsMC$`Gene Name`, file = '[YOUR PATH HERE]/MEvsMCgenes.txt')
#Go to https://www.uniprot.org/id-mapping and load txt file in box. From database: Gene Name, To database: UniProtKB/Swiss-Port, Taxonomy: Homo Sapien
#Download results as uncompressed Excel file with only "Gene Names (primary)" and "Protein Names" selected
#Import Data set as Uniprot_MEvsMC
Uniprot_MEvsMC <- cbind(DeqMS_MEvsMC[match(Uniprot_MEvsMC$...1, DeqMS_MEvsMC$...1),], Uniprot_MEvsMC$Accession, Uniprot_MEvsMC$`Protein Name`)
rownames(Uniprot_MEvsMC) <- 1:nrow(Uniprot_MEvsMC)
names(Uniprot_MEvsMC)[6] <- 'Accession'
names(Uniprot_MEvsMC)[7] <- 'Protein.names'
Uniprot_MEvsMC <- Uniprot_MEvsMC[,c("...1", "Accession", "Gene Name", "LogFC", "FDR", "Direction", "Protein.names")]


#MSvsMC
names(DeqMS_MSvsMC)[4] <- 'FDR'
DeqMS_MSvsMC <- filter(DeqMS_MSvsMC, FDR<0.05)
write(DeqMS_MSvsMC$`Gene Name`, file = '[YOUR PATH HERE]/MSvsMCgenes.txt')
#Go to https://www.uniprot.org/id-mapping and load txt file in box. From database: Gene Name, To database: UniProtKB/Swiss-Port, Taxonomy: Homo Sapien
#Download results as uncompressed Excel file with only "Gene Names (primary)" and "Protein Names" selected
#Import Data set as Uniprot_MSvsMC
Uniprot_MSvsMC <- cbind(DeqMS_MSvsMC[match(Uniprot_MSvsMC$...1, DeqMS_MSvsMC$...1),], Uniprot_MSvsMC$Accession, Uniprot_MSvsMC$`Protein Name`)
rownames(Uniprot_MSvsMC) <- 1:nrow(Uniprot_MSvsMC)
names(Uniprot_MSvsMC)[6] <- 'Accession'
names(Uniprot_MSvsMC)[7] <- 'Protein.names'
Uniprot_MSvsMC <- Uniprot_MSvsMC[,c("...1", "Accession", "Gene Name", "LogFC", "FDR", "Direction", "Protein.names")]

########Creation of All Endurance Data (AllEnd)############
AllEnd <- rbind(Uniprot_MEvsMC, Uniprot_FEvsFC)
AllEnd <- AllEnd[,c("Gene Name","LogFC", "Direction", "Protein.names")]
hold <- AllEnd[!duplicated(AllEnd$Protein.names),]
AllEnd <- AllEnd[,c("Protein.names", "LogFC")]
AllEnd <- aggregate(.~Protein.names, AllEnd, mean)
AllEnd <- cbind(AllEnd, hold[match(AllEnd$Protein.names, hold$Protein.names), c("Gene Name", "Direction")])
#AllEnd <- AllEnd[!duplicated(AllEnd$Protein.names),]
#COL4A1 was only gene oppositely regulated between ME/MC vs FE/FC so we omit it for AllEnd data comparisons
AllEnd <- AllEnd[-c(which(AllEnd$`Gene Name`=="COL4A1")),]

########Aging Data Read In (Simplified)#######
#https://elifesciences.org/articles/49874/figures#content (Figure 1--Source Data 5)
#https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNDk4NzQvZWxpZmUtNDk4NzQtZmlnMS1kYXRhNS12MS54bHN4/elife-49874-fig1-data5-v1.xlsx?_hash=AGpktMQtPvuP28t82swLrP4OSXULHjjaejn5%2F3REzT4%3D

Age_DEPs <- read_excel("Downloads/elife-49874-fig1-data5-v1.xlsx", 
                       +     sheet = "Sheet1")
names(Age_DEPs)[2] <- "LogFC"
Age_DEPs$Direction[Age_DEPs$LogFC > 0 ] = "UP"
Age_DEPs$Direction[Age_DEPs$LogFC < 0 ] = "DOWN"
Uniprot_Age <- cbind(Age_DEPs[match(Uniprot_Age$From, Age_DEPs$`Acc No.`),],
                     Uniprot_Age$`Gene Name`,
                     Uniprot_Age$Accession,
                     Uniprot_Age$`Protein names`)
rownames(Uniprot_Age) <- 1:nrow(Uniprot_Age)
names(Uniprot_Age)[1] <- "Entry"
names(Uniprot_Age)[6] <- "Gene Name"
names(Uniprot_Age)[7] <- "Accession"
names(Uniprot_Age)[8] <- "Protein.names"
Uniprot_Age <- Uniprot_Age[,c("Entry", "Accession", "Gene Name", "LogFC", "Direction", "Protein.names")]


########HIIT Data Read In########
#https://elifesciences.org/articles/69802/figures#content (Figure 1--Source Data 1)
#https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNjk4MDIvZWxpZmUtNjk4MDItZmlnMS1kYXRhMS12MS54bHN4/elife-69802-fig1-data1-v1.xlsx?_hash=y0yHJuD3YqCnZMEIDbmcp5s46rfwftirgIAjztM9%2BkE%3D
HIIT_DEPs <- read_excel("Downloads/elife-69802-supp2-v1.xlsx", 
                        +     sheet = "Sheet2")
names(HIIT_DEPs)[1] <- "Accession"
names(HIIT_DEPs)[2] <- "LogFC"
HIIT_DEPs$Direction[HIIT_DEPs$LogFC > 0 ] = "UP"
HIIT_DEPs$Direction[HIIT_DEPs$LogFC < 0 ] = "DOWN"
names(Uniprot_HIIT)[2] <- "Accession"
names(Uniprot_HIIT)[3] <- "Gene Name"
names(Uniprot_HIIT)[4] <- "Protein.names"
Uniprot_HIIT <- cbind(HIIT_DEPs[match(Uniprot_HIIT$From, HIIT_DEPs$Accession),],
                      Uniprot_HIIT$`Gene Name`,
                      Uniprot_HIIT$Accession,
                      Uniprot_HIIT$Protein.names)
rownames(Uniprot_HIIT) <- 1:nrow(Uniprot_HIIT)
names(Uniprot_HIIT)[4] <- "Gene Name"
names(Uniprot_HIIT)[5] <- "Accession"
names(Uniprot_HIIT)[6] <- "Protein.names"
Uniprot_HIIT <- Uniprot_HIIT[,c("Entry", "Accession", "Gene Name", "LogFC", "Direction", "Protein.names")]

########################################################## Main Comparisons ###################################################### 
########ME/MC vs HIIT#######
MEMCintHIIT <- intersect(Uniprot_HIIT$`Accession`, Uniprot_MEvsMC$Accession)
MEMCvsHIIT <- cbind(Uniprot_HIIT[match(MEMCintHIIT, Uniprot_HIIT$Protein.names),c("Accession", "Gene Name")],
                    MEMCintHIIT, Uniprot_HIIT[match(MEMCintHIIT, Uniprot_HIIT$Protein.names),c("LogFC", "Direction")],
                    Uniprot_MEvsMC[match(MEMCintHIIT, Uniprot_MEvsMC$Protein.names),c("LogFC", "Direction")])
colnames(MEMCvsHIIT) <- c("Accession", "Gene.name", "Protein.names", "HIITLogFC", "HIIT", "MEMCLogFC", "MEvsMC")
MEMCvsHIIT <- data.frame(MEMCvsHIIT)
MEMCvsHIIT$Direction <- ifelse(MEMCvsHIIT$HIIT == MEMCvsHIIT$MEvsMC, "Same", "Opposite")


HIITb <- ggplot(MEMCvsHIIT, aes(Direction, fill = Direction)) +
  geom_bar(show.legend = F) + 
  ggtitle('MEvs.MC Compared to HIIT \n (n=126)') + 
  xlab("Regulation Direction") + ylab("Number of Proteins") +
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.title = element_text(size=20), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), 
        legend.title = element_text(size=20), legend.text = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  scale_fill_manual(values= c( "Opposite" = "darkgrey","Same" = "lightgray")) + 
  geom_text(aes(label=..count..), stat = "count", vjust=-0.5, size=7) 


########All Endurance vs Aging##########
AllEndintAge <- intersect(Uniprot_Age$Protein.names, AllEnd$Protein.names)
AllEndvsAge <- cbind(Uniprot_Age[match(AllEndintAge, Uniprot_Age$Protein.names),c("Accession", "Gene Name")],
                     AllEndintAge, Uniprot_Age[match(AllEndintAge, Uniprot_Age$Protein.names),c("LogFC", "Direction")],
                     AllEnd[match(AllEndintAge, AllEnd$Protein.names),c("LogFC", "Direction")])
colnames(AllEndvsAge) <- c("Accession", "Gene.name", "Protein.names", "AgeLogFC", "Age", "AllEndLogFC", "AllEnd")
AllEndvsAge <- data.frame(AllEndvsAge)
AllEndvsAge$Direction <- ifelse(AllEndvsAge$Age == AllEndvsAge$AllEnd, "Same", "Opposite")


Ab <- ggplot(AllEndvsAge, aes(Direction, fill = Direction)) +
  geom_bar(show.legend=F) + 
  ggtitle('All Endurance Compared to Age \n (n=1265)') + 
  xlab("Regulated Direction") + ylab("Number of Proteins") +
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.title = element_text(size=20),
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), 
        legend.title = element_text(size=20), legend.text = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  scale_fill_manual(values= c( "Opposite" = "darkgrey","Same" = "lightgray")) + 
  geom_text(aes(label=..count..), stat = "count", vjust=-0.5, size=7)


########ME/MC vs T2DM##########

Uniprot_T2DM2$Direction <-ifelse(Uniprot_T2DM2$`FC (log2)` < 0, "DOWN", "UP")
MEMCintT2DM2 <- intersect(Uniprot_T2DM2$UniProtID, Uniprot_MEvsMC$Accession)
MEMCvsT2DM2 <- cbind(Uniprot_T2DM2[match(MEMCintT2DM2, Uniprot_T2DM2$UniProtID),c("UniProtID")],
                     Uniprot_MEvsMC[match(MEMCintT2DM2, Uniprot_MEvsMC$Accession),c("Gene Name", "Protein.names")],
                     Uniprot_T2DM2[match(MEMCintT2DM2, Uniprot_T2DM2$UniProtID),c("FC (log2)", "Direction")],
                     Uniprot_MEvsMC[match(MEMCintT2DM2, Uniprot_MEvsMC$Accession),c("LogFC", "Direction")])
colnames(MEMCvsT2DM2) <- c("Accession", "Gene.name", "Protein.names","T2DM2LogFC", "T2DM2","MEMCLogFC", "MEvsMC")
MEMCvsT2DM2 <- data.frame(MEMCvsT2DM2)
MEMCvsT2DM2$Direction <- ifelse(MEMCvsT2DM2$T2DM2 == MEMCvsT2DM2$MEvsMC, "Same", "Opposite")


T2DM2b <- ggplot(MEMCvsT2DM2, aes(Direction, fill = Direction)) +
  geom_bar(show.legend = F) + 
  ggtitle('MEvs.MC Compared to T2DM \n (n=420)') + 
  xlab("Regulation Direction") + ylab("Number of Proteins") +
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.title = element_text(size=20),
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12),
        legend.title = element_text(size=20), legend.text = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  scale_fill_manual(values= c( "Opposite" = "darkgrey","Same" = "lightgray")) + 
  geom_text(aes(label=..count..), stat = "count", vjust=-0.5, size=7)

########ME/MC vs IFG##########
names(Uniprot_IFG)[4] <- "pvalue"
Uniprot_IFG <- filter(Uniprot_IFG, pvalue<0.05)
Uniprot_IFG$Direction <-ifelse(Uniprot_IFG$`FC (log2)` < 0, "DOWN", "UP")
MEMCintIFG <- intersect(Uniprot_IFG$UniProtID, Uniprot_MEvsMC$Accession)
MEMCvsIFG <- cbind(Uniprot_IFG[match(MEMCintIFG, Uniprot_IFG$UniProtID),c("UniProtID")],
                   Uniprot_MEvsMC[match(MEMCintIFG, Uniprot_MEvsMC$Accession),c("Gene Name", "Protein.names")],
                   Uniprot_IFG[match(MEMCintIFG, Uniprot_IFG$UniProtID),c("FC (log2)", "Direction")],
                   Uniprot_MEvsMC[match(MEMCintIFG, Uniprot_MEvsMC$Accession),c("LogFC", "Direction")])
colnames(MEMCvsIFG) <- c("Accession", "Gene.name", "Protein.names","IFGLogFC", "IFG","MEMCLogFC", "MEvsMC")
MEMCvsIFG <- data.frame(MEMCvsIFG)
MEMCvsIFG$Direction <- ifelse(MEMCvsIFG$IFG == MEMCvsIFG$MEvsMC, "Same", "Opposite")


IFGb <- ggplot(MEMCvsIFG, aes(Direction, fill = Direction)) +
  geom_bar(show.legend = F) + 
  ggtitle('MEvs.MC Compared to IFG \n (n=95)') + 
  xlab("Regulation Direction") + ylab("Number of Proteins") +
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.title = element_text(size=20),
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12),
        legend.title = element_text(size=20), legend.text = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  scale_fill_manual(values= c( "Opposite" = "darkgrey","Same" = "lightgray")) + 
  geom_text(aes(label=..count..), stat = "count", vjust=-0.5, size=7)

########ME/MC vs IGT##########
names(Uniprot_IGT)[4] <- "pvalue"
Uniprot_IGT <- filter(Uniprot_IGT, pvalue<0.05)
Uniprot_IGT$Direction <-ifelse(Uniprot_IGT$`FC (log2)` < 0, "DOWN", "UP")
MEMCintIGT <- intersect(Uniprot_IGT$UniProtID, Uniprot_MEvsMC$Accession)
MEMCvsIGT <- cbind(Uniprot_IGT[match(MEMCintIGT, Uniprot_IGT$UniProtID),c("UniProtID")],
                   Uniprot_MEvsMC[match(MEMCintIGT, Uniprot_MEvsMC$Accession),c("Gene Name", "Protein.names")],
                   Uniprot_IGT[match(MEMCintIGT, Uniprot_IGT$UniProtID),c("FC (log2)", "Direction")],
                   Uniprot_MEvsMC[match(MEMCintIGT, Uniprot_MEvsMC$Accession),c("LogFC", "Direction")])
colnames(MEMCvsIGT) <- c("Accession", "Gene.name", "Protein.names","IGTLogFC", "IGT","MEMCLogFC", "MEvsMC")
MEMCvsIGT <- data.frame(MEMCvsIGT)
MEMCvsIGT$Direction <- ifelse(MEMCvsIGT$IGT == MEMCvsIGT$MEvsMC, "Same", "Opposite")


IGTb <- ggplot(MEMCvsIGT, aes(Direction, fill = Direction)) +
  geom_bar(show.legend = F) + 
  ggtitle('MEvs.MC Compared to IGT \n (n=67)') + 
  xlab("Regulation Direction") + ylab("Number of Proteins") +
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.title = element_text(size=20),
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12),
        legend.title = element_text(size=20), legend.text = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  scale_fill_manual(values= c( "Opposite" = "darkgrey","Same" = "lightgray")) + 
  geom_text(aes(label=..count..), stat = "count", vjust=-0.5, size=7)

########ME/MC vs Prediabetic & T2D##########
MPT <- data.frame(Direction=rep(c("Same", "Opposite","Opposite","Opposite")),
                   Condition=rep(c("T2D", "IFG", "IGT", "T2D")),
                   Count=c(33, 66, 57, 196))
head(MPT)
MPT$Condition <- as.factor(MPT$Condition)
ggplot(data=MPT, aes(x=Condition, y=Count, fill=Direction)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c('darkgray', 'lightgrey'))+
  ggtitle("ME vs. Prediabetics and T2D \n (n=95, n=67, & n=420)") +
  theme(plot.title = element_text(hjust = 0.5, size=30),
        axis.title = element_text(size=20),
        axis.text.y = element_text(size=15), axis.text.x = element_text(size=20),
        legend.title = element_text(size=20), legend.text = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  ylab("Number of Proteins") +
  xlab("")+
  geom_text(aes(label=Count),position = position_dodge(0.9), vjust=-0.5,size=7)


MPT2 <- data.frame(Condition=c(rep("IFG",95), rep("IGT", 67), rep("T2D", 420)),
                      Direction=c(rep("Opposite", 66), rep("Same", 0), rep("None",29), 
                                  rep("Opposite", 57), rep("Same", 0), rep("None",10),
                                  rep("Opposite", 196), rep("Same", 33), rep("None",191)))
MPT2 %>% table()  
MPT2 %>%
  ggplot(aes(x=Condition)) +
  geom_bar(aes(fill=Direction), position = "fill") +
  scale_fill_manual(values=c('white', 'darkgray', 'lightgrey'))+
  ggtitle("MEMC vs IFG, IGT, & T2D") +
  theme(plot.title = element_text(hjust = 0.5, size=30),
        axis.title = element_text(size=20),
        axis.text.y = element_text(size=15), axis.text.x = element_text(size=20),
        legend.title = element_text(size=20), legend.text = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  ylab("Percentage Shared DRPs")
