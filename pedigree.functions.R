

#################
### NEP.calc()
#################

#Description
#This funciton calculations non-exclusion power for any given set of loci with any given set of alleles
#originally from Jamieson, A. 1997

#Input Parameters:
#tmp - a data.frame with only two columns, the first is the name of the loci and the second the frequencies
#      of the alleles at that locus - in long form 

NEP.calc <- function(tmp, opt=1){
  loci <- unique(tmp[,1])
  
  i <- NULL
  OUT <- NULL
  for(i in 1:length(loci)){
    tmp1 <- tmp[tmp[,1] == loci[i],2]
    p1 <- -2*sum(tmp1^2)
    p1
    
    p2 <- sum(tmp1^3)
    p2
    
    p3 <- 2*sum(tmp1^4)
    p3
    
    p4 <- -3*sum(tmp1^5)
    p4
    
    p5 <- -2*(sum(tmp1^2)^2)
    p5
    
    p6 <- 3*sum(tmp1^2)*sum(tmp1^3)
    p6
    
    P <- 1+(sum(p1,p2,p3,p4,p5,p6))
    
    OUT <- c(OUT,P)
  }
  
  OUT
  
  i<-1
  i <- NULL
  OUT1 <- NULL
  
  for(i in 1:length(loci)){
    tmp1 <- tmp[tmp[,1] == loci[i],2]
    p1 <- -4*sum(tmp1^2)
    p1
    
    p2 <- 2*(sum(tmp1^2)^2)
    p2
    
    p3 <- 4*sum(tmp1^3)
    p3
    
    p4 <- -3*sum(tmp1^4)
    p4
    
    
    P <- 1+(sum(p1,p2,p3,p4))
    
    OUT1 <- c(OUT1,P)
  }
  
  OUT1
  
  
  
  i<-1
  i <- NULL
  OUT2 <- NULL
  
  for(i in 1:length(loci)){
    tmp1 <- tmp[tmp[,1] == loci[i],2]
    
    p1 <- 4*sum(tmp1^4)
    p1
    
    p2 <- -4*sum(tmp1^5)
    p2
    
    p3 <- -3*sum(tmp1^6)
    p3
    
    p4 <- -8*(sum(tmp1^2)^2)
    p4
    
    p5 <- 8*sum(tmp1^2)*sum(tmp1^3)
    p5
    
    p6 <- 2*(sum(tmp1^3)^2)
    p6
    
    P <- 1+(sum(p1,p2,p3,p4,p5,p6))
    
    OUT2 <- c(OUT2,P)
  }
  
  
  #parparing for first for loops
  i <-  NULL
  PId <- NULL
  
  for(i in 1:length(loci)){
    
    #grabbing each locus individually and ordering some smallest to largest
    tmp1 <- tmp[tmp[,1] == loci[i],]
    tmp1 <- tmp1[order(tmp1[,2]),]
    tmp1
    
    #raising the allele frequencies to the fourth power
    p14 <-  sum(tmp1[,2]^4)
    p14 #makes test in excel
    
    #parparing for second for loops
    j <- NULL
    calc <- NULL
    
    for(j in 1:nrow(tmp1)){
      
      #take each allele frequencies one at a time
      tmp1.f <- tmp1[j,2]
      tmp1.f
      
      #and taking all allele frequencies greater than that and (2*Pi*Pj)^2 and summing
      tmp2 <- tmp1[-1:-j,2]    
      calc1 <- sum((2*tmp1.f*tmp2)^2)
      
      #adding to calc
      calc <- c(calc,calc1)
    }
    
    #summing all calcs
    calc <- sum(calc)
    calc
    
    #adding p14 and calc
    PId1 <- p14 +calc
    
    #adding to PId
    PId <- c(PId,PId1)
  }
  
  
  
  #parparing for first for loops
  i <-  NULL
  PrSib <- NULL
  
  for(i in 1:length(loci)){
    
    #grabbing each locus individually and ordering some smallest to largest
    tmp1 <- tmp[tmp[,1] == loci[i],]
    tmp1 <- tmp1[order(tmp1[,2]),]
    tmp1
    
    tmp1$p2 <- tmp1[,2]^2
    tmp1$p4 <- tmp1[,2]^4
    
    tmp1
    
    0.5*sum(tmp1$p2)
    (0.5*sum(tmp1$p2))^2
    -0.25*sum(tmp1$p4)
    
    p1 <- 0.25
    p2 <- 0.5*sum(tmp1$p2)
    p3 <- 0.5*(sum(tmp1$p2)^2)
    p4 <- -0.25*(sum(tmp1$p4))
    c(p1,p2,p3,p4)
    ps1 <- sum(p1,p2,p3,p4)
    PrSib <- c(PrSib,ps1)
  }
  
  
  #getting the power calcs accross all loci and adding it to the individual locus calcs
  Locus <- c(loci,"Combined")
  Locus
  
  NE_2P <- c(OUT , prod(1-OUT))
  NE_2P
  
  NE_1P <- c(OUT1 , prod(1-OUT1))
  NE_1P
  
  NE_PP <- c(OUT2 , prod(1-OUT2))
  NE_PP
  
  PrId <- c(PId,prod(PId))
  PrId
  
  PrIdSib <- c(PrSib,prod(PrSib))
  PrIdSib
  
  Output <- data.frame(Locus,NE_1P,NE_2P,NE_PP, PrId,PrIdSib, stringsAsFactors=F)
  
  if(opt==1){
    return(Output)
  }
  
  if(opt==2){
    Output <- Output[Output$Locus == "Combined",]
    return(Output)
  }
}


#################
### sim.gt.data2()
#################

#Description
# This fuction creates genotypic data for parents and offspring using a breeding matrix

# Input parameters
# 1) alfs - a simple data frame first column 1:number of loci and 2 the frequecies of alleles at each locus
# 2) breeding matrix - three column data frame with mom id, dad id, and RS for each pair
# 3) error - defining the amount of random error in the dataset
# 4) Nunrelated - defining the number of unrelated individuals wanted in the dataset
# 5) write2file - defualt is set to FALSE and returns data to console, if TRUE writes data to file in current working directory

#Defining function
sim.gt.data2 <-  function(alfs, breeding.matrix, error= 0, Nunrelated= 0, write2file = F){
  
  #breeding.matrix names
  breeding.matrix [,1] <- paste("Mom",breeding.matrix[,1])
  breeding.matrix [,2] <- paste("Dad",breeding.matrix[,2])
  
  #figuring our how many parents are here
  moms <- unique(breeding.matrix[,1])
  dads <- unique(breeding.matrix[,2])
  par.names <- c(moms,dads)
  Nadults <- sum(length(moms),length(dads))
  
  #getting afreqs input
  afreqs <-alfs
  
  #making a function to simulate genotypes for parents
  
  OUT <- NULL
  sims <- function(sims){
    
    #table allele frequencies
    #getting just the frequencies changing generic name to three character genotypes
    alleles2 <- afreqs[afreqs[,1] == loci[i],]
    alleles3 <- cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])                  
    
    #create homozygote allele frequencies
    homos <- (alleles3[,2])^2                                                          
    homos2 <- cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)
    
    #create heterozygote allele frequencies
    #first create all possible combinations of freqs now make expected freqs
    hets <- t(combn(alleles3[,2],2))                                
    hetfreq <- 2*(hets[,1]*hets[,2])
    
    #create heterozygote allele names
    #now getting the het. genotype combinations and combining them
    hetvals <- t(combn(as.character(alleles3[,1]),2))                                  
    hets2 <- cbind(hetvals,hetfreq)
    
    #combine hets and homos and create genotypes
    gfreqs <- rbind(hets2,homos2)                                                      
    
    #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
    n <- 1000000                                    
    
    #create genotypes(by coloumn, 1 for each allele)
    #for the first and second allele
    gfreqs1 <- rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))                            
    gfreqs2 <- rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
    
    #combining them and mixing them up
    gtypes <- cbind(gfreqs1,gfreqs2)
    gtypes <- gtypes[sample(1:length(gtypes[,1]),replace= F),]
    
    #now getting a random sample of the gentoypes for the parents
    sg1 <- gtypes[sample(1:length(gtypes[,1]),Nadults),]
    OUT<<-cbind(OUT,sg1)
    
  } # end of sim function
  
  #defining the number of loci there are
  loci.length <- length(unique(afreqs[,1]))
  loci <- unique(afreqs[,1])
  loci
  
  #doing the simulation process for each locus
  for(i in 1:length(loci)){
    lapply(1,sims)
  } # end for loop
  OUT
  #saving OUT into object called parents
  parents <- OUT
  
  parents <- as.data.frame(parents,stringsAsFactors = F)
  for(i in 1:ncol(parents)){parents[,i] <- as.numeric(parents[,i])}
  
  #adding names of parents 
  parents$names <- par.names
  
  i <- NULL
  out <- NULL
  for(i in 1:nrow(breeding.matrix)){
    
    #getting the parent ids of interest
    mom <- breeding.matrix[i,1]
    mom.id <-as.numeric(gsub("Mom ","",breeding.matrix[i,1]))
    
    dad <- breeding.matrix[i,2]
    dad.id <-as.numeric(gsub("Dad ","",breeding.matrix[i,2]))
    
    #now picking those parents gts
    mom.gts <- parents[parents$names == mom,-ncol(parents)]
    dad.gts <- parents[parents$names == dad,-ncol(parents)]
    
    #combining parents gts
    parent.gts <- rbind(mom.gts,dad.gts)
    
    #identifying where a new locus starts in the DF
    loci.locs <- seq(from= 1, to=ncol(mom.gts),by = 2)
    
    #making kids genotype
    j<- NULL
    out1 <- NULL 
    for(j in loci.locs){
      
      #mom gts sampled
      alleles <- parent.gts[1,j]
      alleles <- c(alleles, parent.gts[1,j+1])
      
      mom.allele <- sample(alleles,size = breeding.matrix[i,3],replace = T)
      
      #dad gts sampled
      alleles <- parent.gts[2,j]
      alleles <- c(alleles, parent.gts[2,j+1])
      dad.allele <- sample(alleles,size = breeding.matrix[i,3],replace = T)
      
      #combining parents alleles
      alleles <- cbind(mom.allele,dad.allele)
      
      #saving to out1
      out1 <- cbind(out1,alleles)
    }
    
    off.id <-paste("Offspring",paste(mom.id,dad.id, sep="."))
    off.id <- paste(off.id,1:breeding.matrix[i,3],sep=".")  
    out1 <- cbind(off.id,out1)  
    out <- rbind(out,out1)
  }
  
  Offs <- data.frame(out,stringsAsFactors = F)
  for(i in 2:ncol(Offs))(Offs[,i] <- as.numeric(Offs[,i]))
  
  #seperating the parents
  parents <- parents[,c(ncol(parents),1:(ncol(parents)-1))]
  
  Dads <- parents[grepl("Dad",parents[,1])==T,]
  Moms <- parents[grepl("Mom",parents[,1])==T,]
  
  #add error======================================================================#
  #first for the dads
  #error=0
  #getting the gts
  Dadsg <- Dads[,-1]
  
  #figuring out how many gts to make errors in
  ldad <- length(Dadsg)*error
  
  #randomly selecting where to put the errors
  pdad1 <- sample(1:length(Dadsg),ldad,replace= FALSE)
  
  #randomly selecting some gets to replace in the those locations
  pdad2 <- Dadsg[sample(1:length(Dadsg),ldad,replace=FALSE)]
  
  #actually putting in the error
  Dadsg2 <- replace(Dadsg,pdad1,pdad2)
  
  #replacing the old Dads gts
  Dads <- cbind(Dads[,1],Dadsg2)
  
  #doing the same thing for the kids
  Offsg=Offs[,-1]
  loff=length(Offsg)*error
  poff1=sample(1:length(Offsg),loff,replace=FALSE)
  poff2=Offsg[sample(1:length(Offsg),loff,replace=FALSE)]
  Offsg2=replace(Offsg,poff1,poff2)
  Offs=cbind(Offs[,1],Offsg2)
  
  
  #############################
  ### unrelated individuals ###
  #############################
  
  #Add unrealated individuals <- ====================================================#
  
  
  if (Nunrelated>0){
    
    #here to be used as the number of breeders (2* the total number of pairs and number of offspring)
    Nadults <- Nunrelated*3                                                            
    afreqs <- alfs
    OUT <- NULL
    
    #making the allele frequencies again, this way they are completely random and unrelated
    #doing the same thing as lines 34-113
    sims <- function(sims) {
      alleles2 <- afreqs[which(afreqs[,1] == z),]
      #table allele frequencies
      alleles3 <- cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])
      
      #create homozygote allele frequencies
      homos <- (alleles3[,2])^2                                                         
      homos2 <- cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)
      
      #create heterozygote allele frequencies
      hets <- t(combn(alleles3[,2],2))                                                   
      hetfreq <- 2*(hets[,1]*hets[,2])
      
      #create heterozygote allele names
      hetvals <- t(combn(as.character(alleles3[,1]),2))                                  
      hets2 <- cbind(hetvals,hetfreq)
      
      #combine hets and homos and create genotypes
      gfreqs <- rbind(hets2,homos2)                                                      
      
      #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
      n <- 1000000         
      
      #create genotypes(by coloumn, 1 for each allele)
      gfreqs1 <- rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))                            
      gfreqs2 <- rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
      gtypes <- cbind(gfreqs1,gfreqs2)
      gtypes <- gtypes[sample(1:length(gtypes[,1]),replace=F),]
      sg1 <- gtypes[sample(1:length(gtypes[,1]),Nadults),]
      
      OUT <<-cbind(OUT,sg1)
      
    } #end of sims function
    
    #now applying function
    z1 <- length(unique(afreqs[,1]))
    for(z in 1:z1) {lapply(z,sims)}
    
    #repeating lines 122-204
    parents <- OUT
    c <- c(1:(ncol(OUT)))
    odd <- 2*(unique(round(((c-2))/2)))+1
    l <- length(odd) * 1000
    codes <- seq(from=1,to=l,by=1000)
    cols <- sort(rep(codes,2))-1
    Anumbs <- matrix(cols,Nadults,ncol(OUT),byrow=T)
    parents <- as.numeric(parents)+Anumbs
    parents <- as.numeric(parents)-Anumbs
    
    #to called this thing unrelated
    unrelated <- cbind(paste("Individual",1:length(parents[,1])),parents)
    Lid <- length(unrelated[,1])/3
    
    #splitting the unrelated individuals among dads, moms, and kids
    m1 <- unrelated[1:Lid,]
    d1 <- unrelated[(Lid+1):(Lid*2),]
    j1 <- unrelated[((Lid*2)+1):(Lid*3),]
    
    #adding the unrelated individuals to the moms, dads, and kids
    m1 <- data.frame(m1)
    colnames(m1) <- colnames(Moms)
    m1[,1] <- gsub(pattern = "Individual", replacement = "Ind.mom",x = m1[,1])
    Moms <- rbind(Moms,m1)
    for(i in 2:ncol(Moms)){Moms[,i] <- as.numeric(Moms[,i])}
    
    d1 <- data.frame(d1)
    colnames(d1) <- colnames(Dads)
    d1[,1] <- gsub(pattern = "Individual", replacement = "Ind.dad",x = d1[,1])
    Dads <- rbind(Dads,d1)
    for(i in 2:ncol(Dads)){Dads[,i] <- as.numeric(Dads[,i])}
    
    j1 <- data.frame(j1)
    colnames(j1) <- colnames(Offs)
    j1[,1] <- gsub(pattern = "Individual", replacement = "Ind.off",x = j1[,1])
    Offs <- rbind(Offs,j1)
    for(i in 2:ncol(Offs)){Offs[,i] <- as.numeric(Offs[,i])}
  }
  ### final touches
  
  #making generic locus names
  colsnmz <- rep("Locus",ncol(Offs)-1)
  colsnmz2 <- sort(rep(1:(length(colsnmz)/2),2))
  colsnmz3 <- paste(colsnmz,colsnmz2)
  colsnmz3 <- c("IDs",colsnmz3)
  
  #usig them to make the colnames for the Offs, Dads, and Moms
  colnames(Offs)<- colsnmz3
  colnames(Dads)<- colsnmz3
  colnames(Moms)<- colsnmz3
  
  df <- rbind(Moms,Dads,Offs)
  
  p.count <- nrow(Moms)+nrow(Dads)
  k.count <- nrow(Offs)
  Type <- c(rep("Parents",times = p.count),rep("Kids",times = k.count))
  
  df <- data.frame(Type,df,stringsAsFactors = F)
  
  if(write2file==F){
    
    return(df)
    
  } else{
    
    write.table(Moms, "moms.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
    write.table(Dads, "dads.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
    write.table(Offs, "kids.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
  }
}


################
### RRS.test ###
################

#Description
#This funciton calculates relative reproductive success (RRS) between two different factors using a randomization procedure

#Input Parameters:
#df - need two column dataframe with the type you want to look at in column one (can only be two different types) 
#n- sets the number of random RRS you want to calculate
#type- defines what you want to type you want the RS to be relative too

#type <- "W"
#df <- tmp
#head(tmp)
#n <- 1000

rrs.test <- function(df, n=1000, type="Wild"){
  
  #making sure the type column is as.character
  df[,1] <- as.character(df[,1])
  
  #getting the overall and group sample sizes
  Nadults <- nrow(df)
  Ntype1 <-nrow(df[df[,1] == type,])
  Ntype2 <-nrow(df[df[,1] != type,])
  
  Ntype1
  Ntype2
  
  #getting the groupe means
  Ntype1.mean <-mean(df[df[,1] == type,2])
  Ntype1.mean
  Ntype2.mean <-mean(df[df[,1] != type,2])
  Ntype2.mean
  
  #calculating the actual relative reproducive success
  rrs <- Ntype2.mean/Ntype1.mean
  rrs
  
  #preparing for the randomization test
  
  #making generic type names, combining them together, and randomizing them
  type1 <- rep("t1", times= Ntype1)
  type2 <- rep("t2", times= Ntype2)
  types <- c(type1,type2)
  types <- sample(types,size=length(types),replace=T)
  
  
  head(df)
  #making a for loop to calculate the RRS for two randomly selection groups n number of times
  OUT <- NULL
  for(i in 1:n){
    #making a generic types column randomly selecting from types
    df$types <- sample(types,size=length(types),replace=T)
    df
    #getting means
    mean.t1 <- mean(df[df$types == "t1",2])
    mean.t1
    mean.t2 <- mean(df[df$types == "t2",2])
    mean.t2
    
    RRS <- mean.t1/mean.t2
    
    OUT <- c(OUT,RRS)
  } # end of random RRS for loop
  
  OUT <- data.frame(OUT)
  head(OUT)
  
  #preparing data for output
  p.val <- length(OUT[OUT[,1]<rrs,])/nrow(OUT)
  p.val  
  
  rrs <- round(rrs,2)
  output <- data.frame(rrs,p.val,stringsAsFactors=F)
  
  write.table(OUT,"Random.RRS.txt",append=F,quote=F,sep="\t",row.names=F,col.names=T)
  print("Wrote randomized RRS output to current working directory: RRS.output.txt ")
  return(output)
  
  
} #end of RRS permutation

##################
### RRS.test.u ###
##################

#Description
#This funciton calculates relative reproductive success (RRS) between two different factors using a randomization procedure in an
#unbiased way - see Araki and Blouin 2005

#Input Parameters:
#df - need two column dataframe with the type you want to look at in column one (can only be two different types) 
#n- sets the number of random RRS you want to calculate
#type- defines what you want to type you want the RS to be relative too
#Noff - Defualt 0, defines the number of offspring that were used
#Nassign - Defualt 0, defines the number of offspring that assigned to parents
#Npar - Default 0, defines the number of parents used to assign offspring
#b - Default 0, defines the Type B error, see Araki and Blouin 

rrs.test.u <- function(df, n=1000, type="Wild", Noff=0, Nassign=0, Npar=0, b=0){
  
  #calculating the correction to make RRS unbiased
  cor1 <- (Noff-Nassign)/Npar
  cor2 <- b/(1-b)
  cor <- cor1*cor2
  
  #making sure the type column is as.character
  df[,1] <- as.character(df[,1])
  
  #getting the overall and group sample sizes
  Nadults <- nrow(df)
  Ntype1 <-nrow(df[df[,1] == type,])
  Ntype2 <-nrow(df[df[,1] != type,])
  
  #getting the groupe means
  Ntype1.mean <-mean(df[df[,1] == type,2])
  Ntype2.mean <-mean(df[df[,1] != type,2])
  
  #adding in the correction
  Ntype1.mean <- Ntype1.mean-cor
  Ntype2.mean <- Ntype2.mean-cor
  
  #calculating the actual relative reproducive success
  rrs <- Ntype2.mean/Ntype1.mean
  
  #preparing for the randomization test
  
  #making generic type names, combining them together, and randomizing them
  type1 <- rep("t1", times= Ntype1)
  type2 <- rep("t2", times= Ntype2)
  types <- c(type1,type2)
  types <- sample(types,size=length(types),replace=T)
  
  #making a for loop to calculate the RRS for two randomly selection groups n number of times
  OUT <- NULL
  for(i in 1:n){
    #making a generic types column randomly selecting from types
    df$types <- sample(types,size=length(types),replace=T)
    
    #getting means
    mean.t1 <- mean(df[df$types == "t1",2])
    mean.t2 <- mean(df[df$types == "t2",2])
    
    #calculating RRS
    RRS <- mean.t1/mean.t2
    
    OUT <- c(OUT,RRS)
  } # end of random RRS for loop
  
  OUT <- data.frame(OUT)
  head(OUT)
  
  #preparing data for output
  p.val <- length(OUT[OUT[,1]<rrs,])/nrow(OUT)
  p.val  
  
  rrs <- round(rrs,2)
  output <- data.frame(rrs,p.val,stringsAsFactors=F)
  
  write.table(OUT,"Random.RRS.txt",append=F,quote=F,sep="\t",row.names=F,col.names=T)
  print("Wrote randomized RRS output to current working directory: RRS.output.txt ")
  return(output)
  
  
} #end of RRS permutation


#################
### rs.counter
#################

#Description
#This fuction counts the reproductive success (RS) for each individual in that was used in parentage program to assign parents to offspring.

#Input Parameters:
#ped - a dataframe with three columns in this order offspring IDs, mother IDs,and father IDs 
#sup.info - a dataframe with the first column as sample name and the second it's sex (Female ("F") or Male ("M")) and any number of columns after it 

rs.counter <- function(sup.info, ped){
  
  #making a column for RS
  sup.info$rs <- 0
  
  #making sure the necessary columns are in character form
  ped[,2] <-as.character(ped[,2]) 
  ped[,3] <-as.character(ped[,3]) 
  sup.info[,1] <- as.character(sup.info[,1])
  
  j <- 10
  j <- NULL
  #here is the for loop
  for(j in 1:nrow(sup.info)){
    print(j)
    
    if(sup.info[j,2] == "F"){
      
      sup.info$rs[j] <- nrow(ped[ped[,2] == sup.info[j,1],])  
    }
    
    if(sup.info[j,2] == "M"){
      
      sup.info$rs[j] <- nrow(ped[ped[,3] == sup.info[j,1],]
      )  
    }  
  }
  return(sup.info)
}


#################
### ms.counter
#################

#Description
#This fuction counts the mating success (MS) for each individual in that was used in parentage program to assign parents to offspring.

#Input Parameters:
#ped - a dataframe with three columns in this order offspring IDs, mother IDs,and father IDs 
#sup.info - a dataframe with the first column as sample name and the second it's sex (Female ("F") or Male ("M")) and any number of columns after it 

ms.counter <- function(sup.info, ped, missing.parent="UK"){
  
  mp <- missing.parent
  sup.info$mates <- 0
  sup.info$mates.uk <- 0
  
  for(i in 1:nrow(sup.info)){
    tmp1 <- NULL
    print(i)
    if(sup.info$sex[i] == "F"){
      
      tmp1 <- ped[ped[,2] == sup.info[i,1],]
      sup.info$mates[i] <- length(unique(tmp1[,3]))
      sup.info$mates.uk[i] <- length(tmp1[tmp1[,3] == mp, 3])
    }
    
    if(sup.info$sex[i] == "M"){
      tmp1 <- ped[ped[,3] == sup.info[i,1],]
      sup.info$mates[i] <- length(unique(tmp1[,2]))
      sup.info$mates.uk[i] <- length(tmp1[tmp1[,2] == mp, 2])
    } 
  }
  return(sup.info)
}