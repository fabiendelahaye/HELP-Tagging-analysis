### R script to analyze HELP-Tagging data ###

HpaII annotation http://wasp.einstein.yu.edu/downloads/anno_hpaii/anno_hpaii_hg19.txt.bz2

MspI reference http://wasp.einstein.yu.edu/downloads/anno_hpaii/mspi_hg19.hcount.bz2

### From HPAII count to methylation score ###

options(stringsAsFactors=F) #avoid to have to deal with factors type data later on


msp1=read.table("directory for msp1 file",skip=1) #read in msp1 reference
msp1=cbind(msp1[,1:3],msp1[,6]) #select column of interest

name=c("sample_1","sample_2",...,"sample_n") #create a list of name of file we need to upload

for (i in c(1:length(names))){
hcount=read.table(paste("directory of your file",name[i],".hcount",sep=""))
#confidence score
c=sum(hcount[,6])/sum(msp1[,4])



m=merge(msp1[,1:3],cbind(hcount[,1],hcount[,6]),by=1,all.x=T)
assign(name[i],m)} #merge all the file with the reference file msp1

final=cbind(msp1[,1:3],msp1[,4]/sum(na.omit(msp1[,4])),"sample_1"[,4]/sum(na.omit("sample_1"[,4])),"sample_2"[,4]/sum(na.omit("sample_2"[,4])),...,"sample_n"[,4]/sum(na.omit("sample_n"[,4]))) #create a matrix with all the sample as column and cpg as row

angle=final
for(i in c(1:length(name)){
print(i)
angle[,i+4]=atan(final[,i+4]/final[,4])*2/pi*100
} #calculate the angle/methylation score


###linear model analysis ###

###we are comparing 3 groups CTRL/IUGR/LGA with 20 samples in each group

options(stringsAsFactors=F)
mcore=T
ncores=6
outputfile="file_name_1"
outputfile2="file_name_2"
outputdir="output_directory"
require(plyr)
require(doMC)
registerDoMC(cores=ncores)
data=read.table("file_with_angle_value")
batch=read.table("matrix_with_batch_info_by_samples") #table with each column=batch and each row=sample

split_data<-split(data,as.character(data[,2])) #split the data basd on chromosome

llply(split_data, function(eachchrom){
chrnum=as.character(eachchrom[1,2])
table<-t(c("ids","x","IUGRvsC","LGAvsC","LGAvsIUGR","diffIUGRvsC","diffLGAvsC","diffLGAvsIUGR")) #"x" score to know if at least one group have 5 samples with score > 0

write.table(table,file=paste(outputdir,chrnum,"_",outputfile,sep=""),row.names=FALSE,col.names=FALSE,append=F,quote=F, sep="\t") #to get the header of our table

CTRL=as.matrix(eachchrom[,4:23]) #selection of angle value for each group
IUGR=as.matrix(eachchrom[,24:43])
LGA=as.matrix(eachchrom[,44:63])

len=nrow(eachchrom)

sapply(1:len,function(i) {
print(c(paste(eachchrom[1,2],"_",i,sep=""),"out of",len))
#selection of batch of interest
f1=c(batch[which(!is.na(CTRL[i,])),2],batch[(20+which(!is.na(IUGR[i,]))),2],batch[(40+which(!is.na(LGA[i,]))),2])
f2=c(batch[which(!is.na(CTRL[i,])),3],batch[20+which(!is.na(IUGR[i,])),3],batch[40+which(!is.na(LGA[i,])),3])
if(all((table(factor(table(factor(f1,levels=c(1,2,3)))>1,levels=c(FALSE,TRUE)))[2]>1),(length(unique(f2))>1))) {

C=as.vector(na.omit(CTRL[i,]))
I=as.vector(na.omit(IUGR[i,]))
L=as.vector(na.omit(LGA[i,]))
id=as.numeric(eachchrom[i,1])
q=c(C,I,L)
new=data.frame(x=q,treatment=f1,sequence=f2)
new1=new[new$x!="NULL",]
a=try(summary(aov(lm(as.numeric(new1$x)~factor(new1$treatment)*factor(new1$sequence)))),new1)
b=try(TukeyHSD(aov(lm(as.numeric(new1$x)~factor(new1$treatment)*factor(new1$sequence)))),new1)
e=any(c((length(which(((C)>0)))>5),(length(which(((I)>0)))>5),(length(which(((L)>0)))>5)))
 if (length(c(b[[1]][,4],b[[1]][,1]))<2){
d=t(c((id),rep("NaN", 7)))

write.table(d,file=paste(outputdir,chrnum,"_",outputfile,sep=""),row.names=FALSE,col.names=FALSE,append=T,quote=F, sep="\t")

    }else{
    
d=t(c(id,e,b[[1]][,4],b[[1]][,1]))
     names(d)<-NULL
write.table(d, file=paste(outputdir,chrnum,"_",outputfile,sep=""),row.names=FALSE,col.names=FALSE,append=T,quote=F, sep="\t")}}

##matrix exporting number of sample for each comparison

f1factor<-factor(f1,levels=c("1","2","3"))
f2factor<-factor(f2,levels=c("1","2"))
testout<-cbind(table(f1factor,f2factor)[,1], table(f1factor,f2factor)[,2])
rownames(testout)<-paste(rownames(testout),"_",eachchrom[i,1],sep="")     
write.table(testout,file=paste(outputdir,chrnum,"_",outputfile2,sep=""),row.names=T,col.names=FALSE,append=T,quote=F, sep="\t")      
})}
 , .progress=progress_text(char = "-"), .parallel=mcore) 


