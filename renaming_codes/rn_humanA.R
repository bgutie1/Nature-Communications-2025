#This code replaces the gene symbols in data file with approved gene symbols from the HGNC database
#Input file: Column 1: ENS code, Column 2: Gene Symbol
#If there isn't a straightforward match based on ENS, looks at alternate ENS (i.e. ENs.supplied.by.emsembl),
#symbol names, alias symbols, and previous symbols for a match.


library(stringr)
require(data.table)

#genenames.tsv downloaded from https://www.genenames.org/download/statistics-and-files/
#custom download (with columns specified) of all approved symbols. 
genenames<-as.data.frame(fread("genenames.tsv")) 

#from rnaseq output
bgenes<-as.data.frame(fread("features.tsv",header=FALSE))

#Set up working matrix
bgenes2<-bgenes;
bgenes2$result <- 0
bgenes2$new_name <- "";


#Loops over all genes and replaces symbols
#Result 1: Matches ENS, matches Approved Symbol
#Result 2: Matches ENS, doesn't match approved symbol
#Result 3: No Match ENS, but match ENS.supplied.by.ensembl
#Result 4: No Match ENS, Match Approved Symbol
#Result 5: No Match ENS, No Match Approved Symbol, Unique Match Previous Symbol
#Result 6: No Match ENS, No Match Approved Symbol, Non-unique Match Previous Symbol
#Result 7: No Match ENS, No Match Approved Symbol, No Match Previous Symbol, Match Unique Alias
#Result 8: No Match ENS, No Match Approved Symbol, No Match Previous Symbol, Match Non-Unique Alias
#Result 9: No Match at all, no symbol assigned

for (ii in 1:dim(bgenes)[1]){
  name<-""
  matchesENS<-sum(genenames$Ensembl.gene.ID==bgenes[ii,1]);
  if(matchesENS==1){
    idx<-which(genenames$Ensembl.gene.ID==bgenes[ii,1]);
    matchesSymbol<-sum(genenames$Approved.symbol==bgenes[ii,2]);
    if(matchesSymbol==1){
      #Result 1: Matches ENS, matches Approved Symbol
      bgenes2$result[ii]<-1;
      name<-genenames$Approved.symbol[idx];
    } else {
      #Result 2: Matches ENS, doesn't match approved symbol
      bgenes2$result[ii]<-2;
      name<-genenames$Approved.symbol[idx];
    }
  } else if (sum(genenames$Ensembl.ID.supplied.by.Ensembl.==bgenes[ii,1])==1){
    #Result 3: No Match ENS, but match ENS.supplied.by.ensembl
    idx<-which(genenames$Ensembl.ID.supplied.by.Ensembl.==bgenes[ii,1]);
    bgenes2$result[ii]<-3;
    name<-genenames$Approved.symbol[idx];
  } else if (sum(genenames$Approved.symbol==bgenes[ii,2])==1) { 
    #Result 4: No Match ENS, Match Approved Symbol
    idx<-which(genenames$Approved.symbol==bgenes[ii,2]);
    bgenes2$result[ii]<-4;
    name<-genenames$Approved.symbol[idx];
  }else if (sum(str_detect(genenames$Previous.symbols,sprintf("(?<!\\w)(?<!-)%s(?!\\w)(?!-)",bgenes[ii,2])))==1){
    #Result 5: No Match ENS, No Match Approved Symbol, Unique Match Previous Symbol
    bgenes2$result[ii]<-5;
    idx<-which(str_detect(genenames$Previous.symbols,sprintf("(?<!\\w)(?<!-)%s(?!\\w)(?!-)",bgenes[ii,2])));
    name<-genenames$Approved.symbol[idx];
  }  else if (sum(str_detect(genenames$Previous.symbols,sprintf("(?<!\\w)(?<!-)%s(?!\\w)(?!-)",bgenes[ii,2])))>1){
    #Result 6: No Match ENS, No Match Approved Symbol, Non-unique Match Previous Symbol
    bgenes2$result[ii]<-6;
    idx<-which(str_detect(genenames$Previous.symbols,sprintf("(?<!\\w)(?<!-)%s(?!\\w)(?!-)",bgenes[ii,2])));
  } else if (sum(str_detect(genenames$Alias.symbols,sprintf("(?<!\\w)(?<!-)%s(?!\\w)(?!-)",bgenes[ii,2])))==1){
    #Result 7: No Match ENS, No Match Approved Symbol, No Match Previous Symbol, Match Unique Alias
    bgenes2$result[ii]<-7;
    idx<-which(str_detect(genenames$Alias.symbols,sprintf("(?<!\\w)(?<!-)%s(?!\\w)(?!-)",bgenes[ii,2])));
    name<-genenames$Approved.symbol[idx];
  } else if (sum(str_detect(genenames$Alias.symbols,sprintf("(?<!\\w)(?<!-)%s(?!\\w)(?!-)",bgenes[ii,2])))>1){
    #Result 8: No Match ENS, No Match Approved Symbol, No Match Previous Symbol, Match Non-Unique Alias
    bgenes2$result[ii]<-8;
    idx<-which(str_detect(genenames$Alias.symbols,sprintf("(?<!\\w)(?<!-)%s(?!\\w)(?!-)",bgenes[ii,2])));
  }else {
    #Result 9:no matching at all
    bgenes2$result[ii]<-9;
  }
  
  #check for duplicates, add .n to name if there is a duplicate found
  numDups<-sum(bgenes2$new_name==name);
  if(numDups==0){
    bgenes2$new_name[ii]<-name;
  } else{
    if(name!=""){
      bgenes2$new_name[ii]<-sprintf("%s.%d",name,numDups)
    }
  }
  if(ii%%1000==0){
    print(ii);
  }
}

#If no new name found, set the name as old name, also appends .n to the name for duplicates
bgenes3<-bgenes2;
bgenes3$new_name[bgenes2$new_name==""]<-bgenes3$V2[bgenes2$new_name==""]
for (ii in 2:dim(bgenes3)[1]){
  name<-bgenes3$new_name[ii];
  numDups<-sum(bgenes3$new_name[1:ii-1]==name);
  if(numDups>1){
   bgenes3$new_name[ii]<-sprintf("%s.%d",name,numDups-1)
  }
}

#Exports new names as a feature file
featuresb_r<-bgenes;
featuresb_r$V1<-bgenes3$V1;
featuresb_r$V2<-bgenes3$new_name;
featuresb_r$V3<-bgenes3$V3;
write.table(featuresb_r,file="featuresb_r.tsv",quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)

