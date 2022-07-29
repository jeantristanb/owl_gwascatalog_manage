#gwas-catalog-v1.0.3-ancestries-r2022-07-09.tsv  gwas-catalog-v1.0.3-studies-r2022-07-09.tsv  gwascat.r
DataAncestry<-read.csv('gwas-catalog-v1.0.3-ancestries-r2022-07-09.tsv', sep='\t')
DataPuplication<-read.csv('gwas-catalog-v1.0.3-studies-r2022-07-09.tsv', sep='\t')
#DataPubCancer<-DataPuplication[grep('Cancer',DataPuplication$DISEASE.TRAIT,value=F, ignore.case=T),]

listuri<-gsub(' ','',unique(unlist(strsplit(as.character(DataPuplication$MAPPED_TRAIT_URI),split=','))))
writeLines(listuri, file='uri_gc')
Data<-read.table('list.key')
table(listuri %in% Data$V1)


DataUriCancer<-read.table('MONDO_0045024.descent')
DataUriCancer$V1<-as.character(DataUriCancer$V1)

DataPuplication$IsCancer<-sapply(strsplit(as.character(DataPuplication$MAPPED_TRAIT_URI),split=','),function(x){
					 x<-gsub(' ','',x);
                                         return(any(x %in% DataUriCancer$V1))
})

DataPuplication$IsCancer2<-F 
DataPuplication$IsCancer2[grep('cancer', DataPuplication$DISEASE.TRAIT,ignore.case=T)]<-T
DataPuplication$IsCancer2[grep('cancer', DataPuplication$MAPPED_TRAIT,ignore.case=T)]<-T
DataPuplication[DataPuplication$IsCancer2 & DataPuplication$IsCancer==F,'DISEASE.TRAIT']
write.csv(DataPuplication[DataPuplication$IsCancer2 & DataPuplication$IsCancer==F,c('DISEASE.TRAIT','MAPPED_TRAIT_URI')], file='error_cancer.csv', row.names=F)





