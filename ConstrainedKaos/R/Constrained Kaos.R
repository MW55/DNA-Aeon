############################################################################################
# Constrained Kaos
# By Hannah Franziska Loechel
# 03.08.2020
############################################################################################


library(kaos)
library(viridis)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

############################################################
# Color plots
############################################################
color.plot = function(data, color) {
  
  matrixplot = ggplot(melt(data), aes(x = Var1, y = Var2)) +
    
    geom_raster(aes(fill = value)) +
    
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none")+
    
    scale_fill_gradient(low = "white", high=color) +
    
    labs(x = "", y = "", title = "")
  
  matrixplot
  
  
}


#############################################################
# Homopolymers and motifs
#############################################################


hp<-function(sequence,l){
  sequence=strsplit(sequence,"")[[1]]
  x=cgr(sequence,seq.base = c("C","G","A","T"),res=2)$x[length(sequence)]
  y=cgr(sequence,seq.base = c("C","G","A","T"),res=2)$y[length(sequence)]
  init= cgr2(x,y,2^length(sequence))$matrix
  if(l>length(sequence)){
    for (i in (length(sequence)+1):l) {
      init<- matrix(1, nrow =2, ncol = 2) %x% init +init %x% matrix(1, nrow =2, ncol = 2)
      }
    }
  init
}

codewords<-function(matrix,n){
  m=matrix(0,ncol=ncol(matrix),nrow=nrow(matrix))
  m[which(matrix==n,arr.ind = T)]=1
  m
}

######################
# Helper function to obtain cgr from last coordinates
######################
cgr2=function (x,y,res) 
{
  r=1
  
  data.length = length(x)
  
  A = matrix(data = 0, ncol = res, nrow = res)
  
  for (i in 1:data.length) {
    
    
    x.matrix = ceiling((x[i] + r) * res/(2 * r))
    y.matrix = ceiling((y[i] + r) * res/(2 * r))
    A[x.matrix, y.matrix] = A[x.matrix, y.matrix] + 1
  }
  #A=A+matrix(1,ncol=res,nrow=res)
  return(list(matrix = A, x = x, y = y))
}


# Can be freely chosen
sequence="AA"

hp1=hp(sequence, 6)
# Sequences without the substring "AA"
color.plot(codewords(hp1,0),"black")
# Sequences containing the substring "AA", colored based on the occurrence of the substring
color.plot(hp1,rainbow(4))

hp2=hp("GG", 6)
hp3=hp("CC", 6)
hp4=hp("TT", 6)
# Plot for hp>=2, colored based on the occurrence of the substrings, white are codewords with no hp >=2
color.plot(hp1+hp2+hp3+hp4,viridis(20))
# Sequences without hp>=2
color.plot(codewords(hp1+hp2+hp3+hp4,0),"black")


#############################################################
# GC Content and hamming distance
#############################################################


distance<-function(init,l){
  if(l==1){
    return (init)
  }
  
  else{
    a<-matrix(1, nrow =nrow(init), ncol = nrow(init)) %x% init +init %x% matrix(1, nrow =nrow(init), ncol = nrow(init))
    if(l>=3){
      for (i in 3:l) {
        a<- matrix(1, nrow =nrow(a), ncol =nrow(a)) %x% init + a %x% matrix(1, nrow =nrow(init), ncol = nrow(init))
        
      }
    }
    a
    
  }
}


# Flower
flower=matrix(c(1,0,0,1),nrow = 2,ncol = 2)
# Lines
lines=matrix(c(1,1,0,0),nrow = 2,ncol = 2)
# Hamming
hamming=matrix(c(1,1,1,0, 1,1,0,1, 1,0,1,1, 0,1,1,1), nrow =4, ncol = 4)

# GC plot for diagonal order of A and T
f<-distance(flower,4)
# Amount of A and T per sequence
color.plot(f,brewer.pal(4,"PiYG"))
# Sequences with 50 % GC content
color.plot(codewords(f,2),"black")

# GC plot for vertical order of A and T
l<-distance(lines,6)
# Amount of A and T per sequence
color.plot(l,brewer.pal(4,"PiYG"))
# Sequences with 50 % GC content
color.plot(codewords(l,3),"black")

# Hamming distance for sequences of length = 4
h<-distance(hamming,4)
# Hamming distance by colors 
color.plot(h,brewer.pal(4,"YlGnBu"))
# Sequences with hamming distance of exactly 2
color.plot(codewords(h,1),"black")

#############################################################
# Hamming distance for single sequences
#############################################################

hamming=function(sequence){
  
  sequence=strsplit(sequence,"")[[1]]
  
  if(length(sequence)==1){
    d=getM(sequence[1])
  }
  else if(length(sequence)==2){
    start=getM(sequence[1])
    n=getM(sequence[2])
    
    d<-matrix(1, nrow =2, ncol = 2) %x% start +n %x% matrix(1, nrow =2, ncol = 2)
  }
  
  else{
    
    start=getM(sequence[1])
    n=getM(sequence[2])
    
    d<-matrix(1, nrow =2, ncol = 2) %x% start +n %x% matrix(1, nrow =2, ncol = 2)
    
    
    for (i  in c(3:length(sequence))) {
      n=getM(sequence[i])
      d<- matrix(1, nrow =2, ncol = 2) %x% d + n %x% matrix(1, nrow =nrow(d), ncol =nrow(d))
    }}
  d
  
}


######################
# Helper function to obtain matrices for each base
######################
getM=function(s){
  g=matrix(c(0,1,1,1),nrow = 2,ncol = 2)
  c=matrix(c(1,0,1,1),nrow = 2,ncol = 2)
  a=matrix(c(1,1,0,1),nrow = 2,ncol = 2)
  t=matrix(c(1,1,1,0),nrow = 2,ncol = 2)
  if(s=="A"){
    return(a)
  }
  else if(s=="T"){
    return(t)
    
  }
  else if(s=="G"){
    return(g)
    
  } 
  
  else return(c)
}

# Hamming distance to a specific sequence
ha<-hamming("GGATGGA")
# Colored plot of different hamming distances
color.plot(ha,brewer.pal(6,"Spectral"))
# Maximum hamming distance (corresponds to the wordlength)
color.plot(codewords(ha,7),"black")

