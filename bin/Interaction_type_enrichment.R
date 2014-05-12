
list.of.packages <- c("argparse","ggplot2","scales")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages,repos="http://cran.us.r-project.org");
}
suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(scales))

# creat parser object
parser <- ArgumentParser(description="plot the statistical significance for enrichment of different interaction types")

# specify options
parser$add_argument("interaction",help="the strong interaction file,[required]")
parser$add_argument("-n","--num",type="integer", nargs='+', default=c(4,11), help="Column numbers for the type name in two part,[default: %(default)s]")
parser$add_argument("-o","--output",default="interaction_type.pdf",help="output pdf figure file, [default: %(default)s]")


# get the input
args <- parser$parse_args()
interaction <- args$interaction
col_n <- args$num
output <- args$output

## this is just for test.
#name="ACCT"
#folder ="/data2/sysbio/UCSD-sequencing/2013-11-27-Bharat_Tri_Shu/Undetermined_indices/Sample_lane8/" 
#interaction = read.table(paste(folder,name,"_combine/",name,"_fragment_paired_align_rmSingleFragment.txt",sep='')
#                     ,sep='\t',header=F)
#col_n=c(6,15)

interaction = read.table(interaction,sep='\t',header=F)

interaction[as.character(interaction[,col_n[1]])==as.character(interaction[,col_n[1]+2]),col_n[1]+2]="."
interaction[as.character(interaction[,col_n[2]])==as.character(interaction[,col_n[2]+2]),col_n[2]+2]="."
# type:subtype
part1 = apply(interaction,1,function(x) 
  paste(as.character(x[c(col_n[1],col_n[1]+2)]),collapse=":"))
part2 = apply(interaction,1,function(x) 
  paste(as.character(x[c(col_n[2],col_n[2]+2)]),collapse=":"))

part1=gsub("protein_coding","mRNA",part1)
part2=gsub("protein_coding","mRNA",part2)

pairs = apply(cbind(part1,part2),1,function(x) paste(x,collapse="--"))

count = as.data.frame(table(pairs))

temp=apply(count,1,function(x) strsplit(as.character(x[1]),"--")[[1]])

count=cbind(count,t(temp))
colnames(count)[3:4] = c("part1","part2")



# remove those with non and count less than 0.1% of total interaction
name = t(temp)[(!grepl("non",count$pairs))&(count$Freq>0.001*length(pairs)),]
count = count[(count$part1 %in% c(name[,1],"miRNA:."))&(count$part2 %in% c(name[,2],"miRNA:.")),]
count = count[count$Freq>0,]

log_ratio <- function(count,part1,part2) {
  temp = strsplit(as.character(count[1]), "--")[[1]]
  count1 = sum(part1==temp[1])
  count2 = sum(part2==temp[2])
  #ratio = as.numeric(count[2])*length(part1)/count1/count2
  p_value = phyper(as.integer(count[2])-1,count1,length(part1)-count1,count2,lower.tail=F)
  return(-log(p_value))
} 

logRatio = apply(count,1,function(x) log_ratio(x,part1,part2))

count = cbind(count,logRatio)

min = min(count$logRatio)
count$logRatio[count$logRatio>300]=300
max = max(count$logRatio)

#b=c(seq(min,0,length.out=50),seq(0,max,length.out=50))
b=10^(seq(0,log(max,10),length.out=100))
my.colors = colorRampPalette(c("white",'pink','red'))(100)
pdf(output,width=8,height=6)
po.nopanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank()))
plot = ggplot(count)+geom_tile(aes(x=part1,y=part2,fill=logRatio),colour='black')+
  scale_fill_gradientn(values=rescale(b),colours=my.colors,name='-ln(p_value))')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(plot)
dev.off()
