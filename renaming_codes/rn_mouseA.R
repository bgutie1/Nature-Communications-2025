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


#Set up working matrix
mgenes<-data.frame(mRG_Li_exprMatrix$gene);
mgenes2<-mgenes;
mgenes2$result <- 0;
mgenes2$new_name <- "";
colnames(mgenes2)<-c("gene","result","new_name");


#Loops over all genes and replaces symbol if the old symbol was withdrawn
#Result 1: Matches symbol, not withdrawn
#Result 2: Matches symbol, withdrawn
#Result 3: No Match

for (ii in 1:dim(mgenes2)[1]){
  idx<-match(mgenes2[ii,1],mgenenames2$Marker.Symbol);
  name<-"";
  if(!is.na(idx)){
    if(mgenenames2$Status[idx]=="O"){
      mgenes2$result[ii]<-1;
      name<-mgenenames2$Marker.Symbol[idx];
    }
    if(mgenenames2$Status[idx]=="W"){
      mgenes2$result[ii]<-2;
      name<-mgenenames2$Current.Marker.Symbol..if.withdrawn.[idx];
    }
  }
  if(is.na(idx)){
    mgenes2$result[ii]<-3
  }
  mgenes2$new_name[ii]<-name;
  if(ii%%1000==0) print(ii);
}

#If no new name found, set the name as old name
mgenes3<-mgenes2;
mgenes3$new_name[mgenes3$new_name==""]<-mgenes3$gene[mgenes3$new_name==""]

#check for duplicates, add .n to name if there is a duplicate found
for (ii in 2:dim(mgenes3)[1]){
  name<-mgenes3$new_name[ii];
  numDups<-sum(m3genes3$new_name[1:ii-1]==name);
  if(numDups>0){
    mgenes3$new_name[ii]<-sprintf("%s.%d",name,numDups)
  }
}

#Replaces names , exports into expression matrix
mRG_Li_exprMatrix <- read.delim("mRG_Li_exprMatrix.tsv")
mRG_Li_exprMatrix_r<-mRG_Li_exprMatrix
mRG_Li_exprMatrix_r$gene<-mgenes3$new_name;
write.table(mRG_Li_exprMatrix_r,file="mRG_Li_r.tsv",quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)

