setwd("~/Dropbox/Duke Statistics Courses/Spring 2016 Duke Course/STA 613 COMSCI BIO/Final Project/Data")

dat = read.table("dat.txt", header = TRUE)

names(dat)
colnames(dat)
dat.df = data.frame(dat)   # 61468   288

View(dat.df)
breed = c("wildCanids", "wildCanids", 
          "wildCanids", "wildCanids",
          "wildCanids", "wildCanids",
          "acientBreeds", "acientBreeds",
          "acientBreeds", "acientBreeds",
          "acientBreeds", "acientBreeds",
          "acientBreeds", "acientBreeds", 
          "acientBreeds", "acientBreeds",
          "acientBreeds", "acientBreeds",
          "spitz", "spitz",
          "spitz", "spitz",
          "spitz", "spitz",
          "spitz", "spitz",
          "spitz", "spitz",
          "spitz", "spitz",
          "spitz", "spitz",
          "toy", "toy", 
          "toy", "toy", 
          "herding", "herding",
          "herding", "herding",
          "herding", "herding",
          "herding", "herding",
          "herding", "herding",
          "mastiff", "mastiff",
          "mastiff", "mastiff",
          "mastiff", "mastiff",
          "mastiff", "mastiff",
          "mastiff", "mastiff",
          "mastiff", "mastiff",
          "mastiff", "mastiff",
          "mastiff", "mastiff",
          "mastiff", "mastiff",
          "mastiff", "mastiff",
          "mastiff", "mastiff",
          "mastiff", "mastiff",
          "mastiff", "mastiff",
          "mastiff", "mastiff",
          "mastiff", "mastiff")

# extract the breeds that will be used for the analysis
df.pop = cbind(dat.df[,1:6],
               dat.df$CoyoteAllele1, dat.df$CoyoteAllele2, 
               dat.df$RedWolfAllele1, dat.df$RedWolfAllele2, 
               dat.df$WolfAllele1, dat.df$WolfAllele2, 
               dat.df$BasenjiAllele1, dat.df$BasenjiAllele2, 
               dat.df$SalukiAllele1, dat.df$SalukiAllele2, 
               dat.df$AfghanHoundAllele1, dat.df$AfghanHoundAllele2,
               dat.df$KuvaszAllele1, dat.df$KuvaszAllele2, 
               dat.df$ChihuahuaAllele1, dat.df$ChihuahuaAllele2, 
               dat.df$IbizanHoundAllele1, dat.df$IbizanHoundAllele2,
               ## spitz
               dat.df$ChineseSharPeiAllele1, dat.df$ChineseSharPeiAllele2,
               dat.df$ChowChowAllele1, dat.df$ChowChowAllele2,
               dat.df$AkitaAllele1, dat.df$AkitaAllele2,
               dat.df$AlaskanMalamuteAllele1, dat.df$AlaskanMalamuteAllele2,
               dat.df$SiberianHuskyAllele1, dat.df$SiberianHuskyAllele2,
               dat.df$SamoyedAllele1, dat.df$SamoyedAllele2,
               dat.df$AmericanEskimoDogAllele1, dat.df$AmericanEskimoDogAllele2,
               ## toy
               dat.df$PekingeseAllele1, dat.df$PekingeseAllele2,
               dat.df$ShihTzuAllele1, dat.df$ShihTzuAllele2,
               ## herding
               dat.df$GermanShepherdDogAllele1, dat.df$GermanShepherdDogAllele2,
               dat.df$StandardSchnauzerAllele1, dat.df$StandardSchnauzerAllele2,
               dat.df$BriardAllele1, dat.df$BriardAllele1,
               dat.df$GiantSchnauzerAllele1, dat.df$GiantSchnauzerAllele2,
               dat.df$OldEnglishSheepdogAllele1,dat.df$OldEnglishSheepdogAllele2,
               ## mastiff
               dat.df$PugAllele1, dat.df$PugAllele2,
               dat.df$BostonTerrierAllele1, dat.df$BostonTerrierAllele1,
               dat.df$FrenchBulldogAllele1, dat.df$FrenchBulldogAllele2,
               dat.df$StaffordshireBullTerrierAllele1, dat.df$StaffordshireBullTerrierAllele2,
               dat.df$MiniatureBullTerrierAllele1, dat.df$MiniatureBullTerrierAllele2,
               dat.df$BullTerrierAllele1, dat.df$BullTerrierAllele2,
               dat.df$BulldogAllele1, dat.df$BulldogAllele2,
               dat.df$BoxerAllele1, dat.df$BoxerAllele2,
               dat.df$MastiffAllele1, dat.df$MastiffAllele2,
               dat.df$BullMastiffAllele1, dat.df$BullMastiffAllele2,
               dat.df$GreatDaneAllele1, dat.df$GreatDaneAllele2,
               dat.df$RottweilerAllele1, dat.df$RottweilerAllele2,
               dat.df$SaintBernardAllele1, dat.df$SaintBernardAllele2,
               dat.df$BerneseMountainDogAllele1, dat.df$BerneseMountainDogAllele2,
               dat.df$NewfoundlandAllele1, dat.df$NewfoundlandAllele2)

# dim(df.pop) = 61468    82

# rename the dataframe
names(df.pop) = c(names(dat.df)[1:6],substring(names(df.pop)[7:82], 8))


# build the allele frequency table
freq.Table = matrix(0, nrow = 61468, ncol = 38)
for (i in 1:38) {
  #i =38
  freq.Table[,i] = df.pop[, 5 + 2*i] /(df.pop[, 5 + 2*i] + df.pop[, 6 + 2*i])
}

# save the freq table as a data frame
freq.Table= as.data.frame(freq.Table)
# 61468    38

freq.TableName = c()
for (i in 1:38) {
  dummyName = paste(names(df.pop)[5 + 2*i], "Freq", sep = '')
  freq.TableName = c(freq.TableName, dummyName)
}
  
# rename the frequency data
names(freq.Table) = freq.TableName
dim(freq.Table) # 61468    38
# bind the two table
df.pop.freq = cbind(df.pop, freq.Table)

# if value is missing
# build the allele frequency table
missing.Table = matrix(0, nrow = 61468, ncol = 38)
for (i in 1:38) {
  #i =38
  missing.Table[,i] = (df.pop[, 5 + 2*i] + df.pop[, 6 + 2*i])
}

missing.Table = as.data.frame(missing.Table)

# delete the rows that have 0
dummy.freq.Table = freq.Table

for (i in 1:38) {
  control = missing.Table[,i] != 0 
  dummy.freq.Table = dummy.freq.Table[control, ]
  missing.Table = missing.Table[control, ]
}

#> dim(dummy.freq.Table)
#[1] 38312    38
#> dim(freq.Table)
#[1] 61468    38

# dim(missing.Table) 38312    38
# dim(dummy.freq.Table) 38312    38
################
# dummy.freq.Table is the table for final analyzing
save(file = "dummy.freq.Table.Rdata", dummy.freq.Table)
