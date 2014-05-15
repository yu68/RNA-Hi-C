list.of.packages <- c("argparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages,repos="http://cran.us.r-project.org");
}
suppressPackageStartupMessages(require(argparse))

# creat parser object
parser <- ArgumentParser(description="Call intersections based on gene pairs")

# specify options
parser$add_argument("interactionA",help="the interaction file a,[required]")
parser$add_argument("interactionB",help="the interaction file b,[required]")
parser$add_argument("-n","--num",type="integer", nargs='+', default=c(5,12), help="Column numbers for the gene name in two part,[default: %(default)s]")
parser$add_argument("-p","--pvalue",action="store_true", help="set to do 100 permutations for p-value of overlap")
parser$add_argument("-r","--release",action="store_true",help="set to only require match of chromosome and RNA name, but not subtype")
parser$add_argument("-o","--output",default="intersect.txt",help="output intersection file name, pairs in A that overlap with B, [default: %(default)s]")

# get the input
args <- parser$parse_args()
interactionA <- args$interactionA
interactionB <- args$interactionB
col_n <- args$num
output <- args$output
release <- args$release

# read interaction data
interactionA = read.table(interactionA,sep='\t',header=F)
interactionB = read.table(interactionB,sep='\t',header=F)

# get gene pair string for each interaction
Pair <- function(interaction,col_n,release) {
  # release option only require chromosome and name to be the same but not subtype
  if release {
    part1 = paste(as.character(interaction[c(1,col_n[1])]),collapse=":")
    part2 = paste(as.character(interaction[c(1+col_n[2]-col_n[1],col_n[2])]),collapse=":")
  } else {
    part1 = paste(as.character(interaction[c(1,col_n[1],col_n[1]+1)]),collapse=":")
    part2 = paste(as.character(interaction[c(1+col_n[2]-col_n[1],col_n[2],col_n[2]+1)]),collapse=":")
  }
  return(paste(sort(c(part1,part2)),collapse="--"))
}

# names of gene pairs in A and B
NamesA=apply(interactionA,1, function(x) Pair(x,col_n,release))
NamesB=apply(interactionB,1, function(x) Pair(x,col_n,release))


write.table(interactionA[NamesA %in% NamesB,],output,quote=F,sep='\t',col.names=F,row.names=F)

real_intersection_n = sum(NamesA %in% NamesB)
cat("Intersection of gene-pairs:",real_intersection_n,"\n")


if (args$pvalue==T) {
  perm_N=100
  perm_num = rep(0,perm_N)
  for (i in 1:perm_N) {
    interactionB[,c(1+col_n[2]-col_n[1],col_n[2],col_n[2]+1)]=
      interactionB[sample(1:dim(interactionB)[1]),c(1+col_n[2]-col_n[1],col_n[2],col_n[2]+1)]
    interactionA[,c(1+col_n[2]-col_n[1],col_n[2],col_n[2]+1)]=
      interactionA[sample(1:dim(interactionA)[1]),c(1+col_n[2]-col_n[1],col_n[2],col_n[2]+1)]
    NamesA=apply(interactionA,1, function(x) Pair(x,col_n,release))
    NamesB=apply(interactionB,1, function(x) Pair(x,col_n,release))
    perm_num[i]=sum(NamesA %in% NamesB)
    cat("Perm:",i,"Intersections:",perm_num[i],"\r")
  }
  Mean = mean(perm_num)
  SD = sd(perm_num)
  cat("Permutation p-value:",mean(perm_num>real_intersection_n),"mean:",Mean,"sd:",SD,"\n")
  cat("Permutation distribution p-value",1-pnorm(real_intersection_n,Mean,SD),"\n")
}

