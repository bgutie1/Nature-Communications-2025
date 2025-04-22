#This code replaces the gene symbols in data file with approved gene symbols from the HGNC database
#Input file: expression matrix with first column as gene symbols


library(stringr)
require(data.table)

#genenames.tsv downloaded from https://www.genenames.org/download/statistics-and-files/
#custom download (with columns specified) of all approved symbols. 
genenames<-as.data.frame(fread("genenames.tsv")) 
exprMatrix.kriegstein<-as.data.frame(fread("exprMatrix.kriegstein.tsv.gz",header=TRUE))
kgenes<-as.data.frame(exprMatrix.kriegstein$gene);

#set uf working matrix
kgenes2<-kgenes;
kgenes2$result <- 0;
kgenes2$new_name <- "";
colnames(kgenes2)<-c("gene","result","new_name");

#Result 1: Unique Matches Approved Symbol
#Result 2: Unique Match Previous Symbol
#Result 4: Unique match alias
#Result 6: No match

for (ii in 1:dim(kgenes)[1]){
  
  #gene symbol matches Approved.symbol
  matchesSymbol<-sum(genenames$Approved.symbol==kgenes[ii,1]);

  if(matchesSymbol==1){
    #Result 1: Unique Matches Approved Symbol
    idx<-which(genenames$Approved.symbol==kgenes[ii,1]);
    kgenes2$result[ii]<-1;
    kgenes2$new_name[ii]<-genenames$Approved.symbol[idx];
  } else if (sum(str_detect(genenames$Previous.symbols,sprintf("(?<!\\w)(?<!-)%s(?!\\w)(?!-)",kgenes[ii,1])))==1){
    #Result 2: Unique Match Previous Symbol
    kgenes2$result[ii]<-2;
    idx<-which(str_detect(genenames$Previous.symbols,sprintf("(?<!\\w)(?<!-)%s(?!\\w)(?!-)",kgenes[ii,1])));
    kgenes2$new_name[ii]<-genenames$Approved.symbol[idx];
  } else if (sum(str_detect(genenames$Alias.symbols,sprintf("(?<!\\w)(?<!-)%s(?!\\w)(?!-)",kgenes[ii,1])))==1){
    #Result 4: Unique match alias
    kgenes2$result[ii]<-4;
    idx<-which(str_detect(genenames$Alias.symbols,sprintf("(?<!\\w)(?<!-)%s(?!\\w)(?!-)",kgenes[ii,1])));
    kgenes2$new_name[ii]<-genenames$Approved.symbol[idx];
  } 
  else {
    #Result 6: No match
    kgenes2$result[ii]<-6
  }
  if(ii%%1000==0){print(ii);}
}


#If no new name found, set the name as old name, also appends to name for duplicates
kgenes3<-kgenes2;
kgenes3$new_name[kgenes2$new_name==""]<-kgenes3$gene[kgenes2$new_name==""]

for (ii in 2:dim(kgenes3)[1]){
  name<-kgenes3$new_name[ii];
  if(name!=""){
    numDups<-sum(kgenes3$new_name[1:ii-1]==name);
    if(numDups>0){
      kgenes3$new_name[ii]<-sprintf("%s.%d",name,numDups)
    }
  }
  if(ii%%1000==0){print(ii);}
}

#Replaces name of genes in expression matrix, exports it.
exprMatrix.kriegstein.r<-exprMatrix.kriegstein;
exprMatrix.kriegstein.r$gene<-kgenes3$new_name;
write.table(exprMatrix.kriegstein.r,file="k_r.tsv",quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)

