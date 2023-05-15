#gwas-catalog-v1.0.3-ancestries-r2022-07-09.tsv  gwas-catalog-v1.0.3-studies-r2022-07-09.tsv  gwascat.r
library(ggplot2)
DataAncestry<-read.csv('../datai/gwas-catalog-v1.0.3-ancestries-r2022-07-09.tsv', sep='\t')
DataAncestry<-DataAncestry[DataAncestry$STAGE=='initial',]
DataPublication<-read.csv('../datai/gwas-catalog-v1.0.3-studies-r2022-07-09.tsv', sep='\t')
#DataPubCancer<-DataPublication[grep('Cancer',DataPublication$DISEASE.TRAIT,value=F, ignore.case=T),]

listuri<-gsub(' ','',unique(unlist(strsplit(as.character(DataPublication$MAPPED_TRAIT_URI),split=','))))
writeLines(listuri, con='../result/uri_gc')
Data<-read.table('list.key')
table(listuri %in% Data$V1)


DataUriCancer<-read.table('../result/MONDO_0045024.descent')
DataUriCancer$V1<-as.character(DataUriCancer$V1)

DataPublication$IsCancer<-sapply(strsplit(as.character(DataPublication$MAPPED_TRAIT_URI),split=','),function(x){
					 x<-gsub(' ','',x);
                                         return(any(x %in% DataUriCancer$V1))
})

DataPublicationCancer<-DataPublication[DataPublication$IsCancer,]
DataPublicationCancer2<-merge(DataPublicationCancer,DataAncestry,all.x=T,by=c('STUDY.ACCESSION','PUBMED.ID'))
DataPublicationCancer3<-unique(DataPublicationCancer2[,c('STUDY.ACCESSION','PUBMED.ID', 'DISEASE.TRAIT','INITIAL.SAMPLE.SIZE','MAPPED_TRAIT','NUMBER.OF.INDIVIDUALS', 'BROAD.ANCESTRAL.CATEGORY')])

DataPublicationCancer3_big<-DataPublicationCancer3[DataPublicationCancer3$NUMBER.OF.INDIVIDUALS>10000,]
DataPublicationCancer3_big$PUBMED.ID<-as.character(DataPublicationCancer3_big$PUBMED.ID)
DataPublicationCancer3_big$ANCESTRAL.CATEGORY<-as.character(DataPublicationCancer3_big$BROAD.ANCESTRAL.CATEGORY)
DataPublicationCancer3_big$MAPPED_TRAIT<-as.character(DataPublicationCancer3_big$MAPPED_TRAIT)
DataPublicationCancer3_big$BROAD.ANCESTRAL.CATEGORY<-as.character(DataPublicationCancer3_big$BROAD.ANCESTRAL.CATEGORY)
DataPublicationCancer3_big$NUMBER.OF.INDIVIDUALS<-as.character(DataPublicationCancer3_big$NUMBER.OF.INDIVIDUALS)
DataPublicationCancer3_big$INITIAL.SAMPLE.SIZE<-as.character(DataPublicationCancer3_big$INITIAL.SAMPLE.SIZE)
mergepubid<-aggregate(cbind(MAPPED_TRAIT,BROAD.ANCESTRAL.CATEGORY, NUMBER.OF.INDIVIDUALS,INITIAL.SAMPLE.SIZE)~PUBMED.ID, data=DataPublicationCancer3_big, paste, collapse=';')
write.csv(mergepubid, file='../result/cancerstudies_more10000.csv', row.names=F)


write.csv(DataPublicationCancer3, file='../result/cancer_ethnicity.csv',row.names=F)
DataPublicationCancer3_big$BROAD.ANCESTRAL.CATEGORY<-as.character(DataPublicationCancer3_big$BROAD.ANCESTRAL.CATEGORY)


DataPublicationCancer3$BROAD.ANCESTRAL.CATEGORY<-gsub('NR',DataPublicationCancer3$BROAD.ANCESTRAL.CATEGORY)
DataPublicationCancer3$BROAD.ANCESTRAL.CATEGORY2<-sapply(strsplit(as.character(DataPublicationCancer3$BROAD.ANCESTRAL.CATEGORY), split=','), function(x){
x2<-x[nchar(gsub(' ','', x))>1]
DataPublicationCancer3_big$PUBMED.ID<-as.character(DataPublicationCancer3_big$PUBMED.ID)
if(length(x)==0)return('NR')
if(length(x)==1)return(x2)
return ('TransEthnic')
})

DataPublication$IsCancer2<-F 
DataPublication$IsCancer2[grep('cancer', DataPublication$DISEASE.TRAIT,ignore.case=T)]<-T
DataPublication$IsCancer2[grep('cancer', DataPublication$MAPPED_TRAIT,ignore.case=T)]<-T
DataPublication[DataPublication$IsCancer2 & DataPublication$IsCancer==F,'DISEASE.TRAIT']

DataPublicationCancer3$BROAD.ANCESTRAL.CATEGORY2
DataPublicationCancer3$Ancestry<-DataPublicationCancer3$BROAD.ANCESTRAL.CATEGORY2
DataPublicationCancer3$Ancestry[DataPublicationCancer3$BROAD.ANCESTRAL.CATEGORY2 %in% c('Asian unspecified', 'East Asian', 'South East Asian', 'South Asian')]<-'Asian'
DataPublicationCancer3$Ancestry[DataPublicationCancer3$BROAD.ANCESTRAL.CATEGORY2 %in% c('Other', 'NR')]<-'No specified / other'
DataPublicationCancer3$Ancestry[DataPublicationCancer3$BROAD.ANCESTRAL.CATEGORY2 %in% c('African unspecified', 'African American or Afro-Caribbean')]<-'African ancestry (other)'
#ggplot(data, aes(x = factor(x), fill = factor(x))) +
#  geom_bar() +
#  geom_text(aes(label = ..count..), stat = "count", vjust = 1.5, colour = "white")
DataPublicationCancer4<-unique(DataPublicationCancer3[,c('STUDY.ACCESSION', 'BROAD.ANCESTRAL.CATEGORY')])
ggplot(data=DataPublicationCancer3, aes(x=factor(Ancestry), fill=factor(Ancestry))) +geom_bar() +   coord_flip() + ylab('Studies number') + xlab('') + geom_text(aes(label = ..count..), stat = "count", vjust = 1.5, colour = "black")
ggsave('../result/dist_cancer_studies_gc.jpeg')

DataPublicationCancer4<-unique(DataPublicationCancer3[,c('STUDY.ACCESSION', 'NUMBER.OF.INDIVIDUALS','Ancestry')])
DataPublicationCancer4<-aggregate(NUMBER.OF.INDIVIDUALS~STUDY.ACCESSION+Ancestry,data=DataPublicationCancer4,mean)
#DataPublicationCancer4$NUMBER.OF.INDIVIDUALS[DataPublicationCancer4$NUMBER.OF.INDIVIDUALS>10000]<-10000
stat_box_data <- function(y, upper_limit = max(DataPublicationCancer4$NUMBER.OF.INDIVIDUALS) * 1.15) {
  return(
    data.frame(
      y = 0.75 * upper_limit,
      label = paste('N=', length(y),sep='')
    )
  )
}
ggplot(data=DataPublicationCancer4, aes(x=factor(Ancestry),y=NUMBER.OF.INDIVIDUALS , fill=factor(Ancestry))) +geom_boxplot() + 
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.9
  ) + coord_flip()+ theme(legend.position="none") + ylab('n')
ggsave('../result/dist2_cancer_studies_gc.jpeg')


#write.csv(DataPublication[DataPublication$IsCancer2 & DataPublication$IsCancer==F,c('DISEASE.TRAIT','MAPPED_TRAIT_URI')], file='error_cancer.csv', row.names=F)


