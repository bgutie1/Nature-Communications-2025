library(stringr)
require(data.table)

#mouse gene names downloaded from https://www.informatics.jax.org/downloads/reports/index.html#marker
#MRK_List1.rpt.gz (including withdrawn marker symbols) 
mgenenames2<-as.data.frame(fread("MRK_List1.rpt.gz")) 
colnames(mgenenames2)<-c("MGI.Accession.ID","Chr"                                    
                         ,"cM.Position"                            
                         ,"genome.coordinate.start"                
                         ,"genome.coordinate.end"                  
                         ,"strand"                                 
                         ,"Marker.Symbol"                          
                         ,"Status"                                 
                         ,"Marker.Name"                            
                         ,"Marker.Type"                            
                         ,"Feature.Type"                           
                         ,"Marker.Synonyms..pipe.separated."       
                         ,"Current.MGI.Accession.ID..if.withdrawn."
                         ,"Current.Marker.Symbol..if.withdrawn.");


#Read features file
GSM5039270_G1_CellRanger_outs_features.tsv <- read.delim("GSM5039270_G1_CellRanger_outs_features.tsv.gz", header=FALSE)

#Set up working matrix
m3genes<-data.frame(GSM5039270_G1_CellRanger_outs_features.tsv);
lastidx<-dim(m3genes)[1];
m3genes2<-data.frame(GSM5039270_G1_CellRanger_outs_features.tsv)[1:lastidx,];
m3genes2$result <- 0;
m3genes2$new_name <- "";


#Loops over all genes and replaces symbol if the old symbol was withdrawn
#Result 1: Matches symbol, not withdrawn
#Result 2: Matches symbol, withdrawn
#Result 3: No Match

for (ii in 1:lastidx){
  idx<-match(m3genes2[ii,2],mgenenames2$Marker.Symbol);
  name<-"";
  if(!is.na(idx)){
    
    if(mgenenames2$Status[idx]=="O"){
      m3genes2$result[ii]<-1;
      name<-mgenenames2$Marker.Symbol[idx];
    }
    if(mgenenames2$Status[idx]=="W"){
      m3genes2$result[ii]<-2;
      name<-mgenenames2$Current.Marker.Symbol..if.withdrawn.[idx];
    }

  }
  if(is.na(idx)){
    m3genes2$result[ii]<-3
  }
  m3genes2$new_name[ii]<-name;
  if(ii%%1000==0) print(ii);
}

#If no new name found, set the name as old name
m3genes3<-m3genes2;
m3genes3$new_name[m3genes2$new_name==""]<-m3genes2$V2[m3genes2$new_name==""]

#check for duplicates, add .n to name if there is a duplicate found
for (ii in 2:dim(m3genes3)[1]){
  name<-m3genes3$new_name[ii];
  numDups<-sum(m3genes3$new_name[1:ii-1]==name);
  if(numDups>0){
    m3genes3$new_name[ii]<-sprintf("%s.%d",name,numDups)
  }
}

#Exports new names as a feature file, same features file in other data set
featuresm3_1.tsv<-  read.delim("GSM5039270_G1_CellRanger_outs_features.tsv.gz", header=FALSE);
featuresm3_1_r<-featuresm3_1.tsv
featuresm3_1_r$V2<-m3genes3$new_name;
write.table(featuresm3_1_r,file="featuresm3_1_r.tsv",quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)
write.table(featuresm3_1_r,file="featuresm3_2_r.tsv",quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)

