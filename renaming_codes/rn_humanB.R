#This code replaces the gene symbols in data file with approved gene symbols from the HGNC database
#Input file: expression matrix with first column as gene symbols, multiple files
#If there isn't a straightforward match based on ENS, looks at alternate ENS (i.e. ENs.supplied.by.emsembl),
#symbol names, alias symbols, and previous symbols for a match.


library(stringr)
require(data.table)
install.packages("janitor", dependencies=TRUE)
library(janitor)

#genenames.tsv downloaded from https://www.genenames.org/download/statistics-and-files/
#custom download (with columns specified) of all approved symbols. 
genenames<-as.data.frame(fread("genenames.tsv")) 

GSM3587923_AML1012.D0.dem <- read.delim("GSM3587923_AML1012-D0.dem.txt")
sgenes<-data.frame(GSM3587923_AML1012.D0.dem$Gene);
sgenes2<-sgenes;
sgenes2$result <- 0;
sgenes2$new_name <- "";
colnames(sgenes2)<-c("gene","result","new_name");


#Result 1: Unique Matches Approved Symbol
#Result 2: Unique Match Previous Symbol
#Result 4: Unique match alias
#Result 6: No match


for (ii in 1:dim(sgenes)[1]){
  
  matchesSymbol<-sum(genenames$Approved.symbol==sgenes[ii,1]);

  if(matchesSymbol==1){
    #Result 1: Unique Matches Approved Symbol
    idx<-which(genenames$Approved.symbol==sgenes[ii,1]);
    sgenes2$result[ii]<-1;
    sgenes2$new_name[ii]<-genenames$Approved.symbol[idx];
  } else if (sum(str_detect(genenames$Previous.symbols,sprintf("(?<!\\w)(?<!-)%s(?!\\w)(?!-)",sgenes[ii,1])))==1){
    #Result 2: Unique Match Previous Symbol
    sgenes2$result[ii]<-2;
    idx<-which(str_detect(genenames$Previous.symbols,sprintf("(?<!\\w)(?<!-)%s(?!\\w)(?!-)",sgenes[ii,1])));
    sgenes2$new_name[ii]<-genenames$Approved.symbol[idx];
  }  
  else if (sum(str_detect(genenames$Alias.symbols,sprintf("(?<!\\w)(?<!-)%s(?!\\w)(?!-)",sgenes[ii,1])))==1){
    #Result 4: Unique match alias
    sgenes2$result[ii]<-4;
    idx<-which(str_detect(genenames$Alias.symbols,sprintf("(?<!\\w)(?<!-)%s(?!\\w)(?!-)",sgenes[ii,1])));
    sgenes2$new_name[ii]<-genenames$Approved.symbol[idx];
  } 
  else {
    #Result 6: No match
    sgenes2$result[ii]<-6
  }
  
  if(ii%%1000==0){
    print(ii);
  }
  
}

#If no new name found, set the name as old name
sgenes3<-sgenes2;
sgenes3$new_name[sgenes3$new_name==""]<-sgenes3$gene[sgenes2$new_name==""]


#load each expression matrix in each file, replace gene names, save matrix

counts_files<- list.files(path = ".", full.names = TRUE, pattern = "*dem.txt.gz")
for (rr in 1:length(counts_files)){
  print(sprintf("Processing No. %d: %s", rr, counts_files[rr]))
  dataname<-sub(pattern = "(.*?)\\..*$", replacement = "\\1", basename(counts_files[rr]))
  y<-read.table(counts_files[rr], sep="\t");
  y<-y %>% row_to_names(row_number=1);
  y$Gene<-sgenes3$new_name;
  write.table(y,file=sprintf("%s_r.txt",dataname),quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)
}

