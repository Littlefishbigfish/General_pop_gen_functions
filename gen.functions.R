#some general functions for working with genetic data in R.

#################
### alf.freqs()
#################

#Description
# This fuction was originally written by Mark Christie, but has been modified to add additional
#parameters.

#Input Parameters:

# population - a table with sample.name in column one, and loci after that
# three.char - Default is F, which means the missing data is 0
#   if T then missing data is coded as something else, which is defined by missing.data
# missing.data - Default is "000", only envoked when three.char is true

alf.freq <- function(population, missing.data = "000", one.pop=F, opt = 1) {
  
  if(one.pop==T){
    
    x <- colnames(population)
    population <- data.frame("Pop",population,stringsAsFactors=F)
    colnames(population) <- c("Pops",x)
    
    pops <- unique(population[,1])
    pops <- as.character(pops)
    population <- population[,c(-2)]
    
  } else {
   
    pops <- unique(population[,1])
    population <- population[,c(-2)]
    
  }
  
  OUT1 <- NULL
  i <- 1
  i <- NULL
  for(i in 1:length(pops)){
    
    population1 <- population[population[,1] == pops[i], ]
    head(population1)
  
    #removing the sample names
    population1 <- population1[,-1]
  
  #if the data are in three character form then it replaces whatever missing data is defined as 
  if(missing.data != "000"){  
    population1[population1 == missing.data] <- "0"
    population1 <- data.frame(lapply(population1, as.numeric))
  } else {
    population1 <- data.frame(lapply(population1, as.numeric))
  }
  
  #this bit grabs the locus names for the population1
  L <- 1:ncol(population1)
  locus.locs <- L[seq(1,ncol(population1),2)]
  locus.names <- colnames(population1)
  
  OUT <- NULL 
  
  #this for loop grabs the two columns of alleles for each locus
  #and uses table to count the number of times each allele is observed
  for (x in locus.locs) {
    
    alleles <- c(population1[,x],population1[,x+1])
    alleles2 <- as.data.frame(table(alleles))
    missing.alleles <- alleles2[which(alleles2[,1]==0),2]
    
    #removes the missing values  
    if(length(which(alleles2[,1]==0))>0) {
      alleles2=alleles2[-which(alleles2[,1]==0),]
    }
    
    #this calculates the frequency and puts it all in a table
    alleles4 <- cbind(alleles2,alleles2[,2]/sum(alleles2[,2])) 
    output <- cbind(pops[i],locus.names[x],alleles4)                      
    OUT <- as.data.frame(rbind(OUT,output),stringsAsFactors=F)
  }

    OUT1 <- as.data.frame(rbind(OUT1,OUT),stringsAsFactors=F)
    OUT1  
  }
  
  #adds the column names
  colnames(OUT1) <- c("Pops","Locus","allele","count","frequency")
  Allelefreqs <- OUT1
  
  if(opt == 1){ 
    Allelefreqs$Pops <- as.character(Allelefreqs$Pops)
    Allelefreqs$Locus <- as.character(Allelefreqs$Locus)
    Allelefreqs$allele <- as.character(Allelefreqs$allele)
    return(Allelefreqs) }

  if(opt == 2){
    
    alfs <- Allelefreqs[,c("Locus","frequency")]
    alfs$Locus <- as.character(alfs$Locus)
    locs <- unique(alfs$Locus)
    
    OUT <- NULL
    for(i in 1:length(locs)){
      num <- nrow(alfs[alfs$Locus == locs[i],])
      OUT1 <- rep(as.numeric(paste(i)),times=num)
      OUT <- c(OUT,OUT1)
    }
    alfs$Locus <- OUT
    return(alfs)
  }

}


#################
### gt.summary()
#################

#Description
# This fuction was originally written by Mark Christie, but has been modified to add additional
# lines to calculated observed and expected heterozygosity. In addition it calculates the
# percent genotyped

#Input Parameters:

# population - a table with sample.name in column one, and loci after that
# three.char - Default is F, which means the missing data is 0
#   if T then missing data is coded as something else, which is defined by missing.data
# missing.data - Default is "000", only envoked when three.char is true
# pop.sum - Default is F, which means it will not summarize pop gene stats accross loci
# one.pop - Default is T, which means that it will automatically assume you only have one population
# round.num - Default set to 4, meaning it will round the population summary to 4 decimals



gt.summary <- function(population, missing.data = "000", pop.sum=F, one.pop=T, round.num=4) {
  
  if(one.pop==T){
    
    x <- colnames(population)
    population <- data.frame("Pop",population,stringsAsFactors=F)
    colnames(population) <- c("Pops",x)
    
    pops <- unique(population[,1])
    pops <- as.character(pops)
    population <- population[,c(-2)]
    
  } else{
    
    pops <- unique(population[,1])
    pops <- as.character(pops)
    population <- population[,c(-2)]
    
  }
  

  i <- NULL
  OUT1 <- NULL
  for(i in 1:length(pops)){
    
    population1 <- population[population[,1] == pops[i], ]

    #removing the sample names
    population1 <- population1[,-1]
    
    #if the data are in three character form then it replaces whatever missing data is defined as 

    population1[population1 == missing.data] <- "0"
    population1 <- data.frame(lapply(population1, as.numeric))
    
    #this bit grabs the locus names for the population1
    L <- 1:ncol(population1)
    locus.locs <- L[seq(1,ncol(population1),2)]
    locus.names <- colnames(population1)
    
    OUT <- NULL 

    #counting the number of individuals in the sample
    x <- NULL  
    for (x in locus.locs) {
      tot <- nrow(population1)
      
      #calculating percent genotyped
      missing <- nrow(population1[population1[,x]==0,])
      missing1 <- nrow(population1[population1[,x]==0,])
      Percent.gt <- round((((tot-missing)/tot)*100),2)
      
      #calculating expected heterozygosity
      alleles <- c(population1[,x],population1[,x+1])
      num.alleles <- length(unique(alleles))
      alleles2 <- as.data.frame(table(alleles))
      missing <- alleles2[which(alleles2[,1]==0),2]
      
      #removing zeros
      if(length(which(alleles2[,1]==0))>0) {
        alleles2 <- alleles2[-which(alleles2[,1]==0),]
      }
      
      #calculating p for all loci and then squaring them for all alleles, and calculating He using Nei 1987 eq. 8.4
      N <- tot- as.numeric(missing1)
      alleles4 <- data.frame(alleles2,alleles2[,2]/sum(alleles2[,2])) 
      alleles4$p2 <- alleles4[,3]^2
      bias.he <- 1-sum(alleles4$p2)
      bias.he
      correction <- (2*N)/(2*N-1)
      correction
      exp.het <- round(correction*bias.he,2)
      exp.het
      
      #calculating observed heterozygosity
      alleles <- data.frame(cbind(population1[,x],population1[,x+1]))
      alleles$truth <- ifelse(alleles$X1 == alleles$X2,T,F)
      missing <- nrow(subset(alleles, X1 == 0 & X2 == 0))
      tot2 <- nrow(alleles)-missing
      het <- nrow(subset(alleles, truth ==  F))
      obs.het <- round(het/tot2,2)
      fis <- round(1- (obs.het/exp.het),2)
      
      output <- data.frame(pops[i],locus.names[x], num.alleles, exp.het, obs.het, fis, tot, missing, Percent.gt)
      OUT <- as.data.frame(rbind(OUT,output),stringsAsFactors=F)
    }
    
    OUT1 <- rbind(OUT1,OUT)
    
  }
  
  colnames(OUT1) <- c("Pop","Locus","A","He","Ho","Fis","N","N.Missing","P.GT")
  OUT1 <- as.data.frame(OUT1,stringsAsFactors=F)  
  OUT1$Pop <- as.character(OUT1$Pop)
  OUT1$Locus <- as.character(OUT1$Locus)
  
  if(pop.sum==F){
    return(OUT1)
  } else {
    
    #getting unique pop identifiers
    pops <- as.character(unique(OUT1[,1]))
    
    #
    df <- OUT1
    i <- 1
    i <- NULL
    OUT <- NULL
    
    for(i in 1:length(pops)){
      #getting just the pop of interest
      df1 <- df[df[,1] == pops[i],]
      df1
      
      #getting mean and standard error for number of alleles
      avg.a <- round(mean(df1[,"A"]),round.num)
      se.a <- round(st.err(df1[,"A"]),round.num)
      
      #getting mean and standard error for expected heterozygosity
      avg.he <- round(mean(df1[,"He"]),round.num)
      se.he <- round(st.err(df1[,"He"]),round.num)  
      
      #getting mean and standard error for observed heterozygosity
      avg.ho <- round(mean(df1[,"Ho"]),round.num)
      se.ho <- round(st.err(df1[,"Ho"]),round.num)  
      
      #getting mean and standard error for Fis
      avg.fis <- round(mean(df1[,"Fis"]),round.num)
      se.fis <- round(st.err(df1[,"Fis"]),round.num)  
      
      #getting mean and standard error for n.missing
      avg.nm <- round(mean(df1[,"N.Missing"]),round.num)
      se.nm <- round(st.err(df1[,"N.Missing"]),round.num)  
      
      #getting mean and standard error for percent genotyped
      avg.pgt <- round(mean(df1[,"P.GT"]),round.num)
      se.pgt <- round(st.err(df1[,"P.GT"]),round.num)  
      
      #putting it all together and making a nice neat DF
      OUT1 <- cbind(df1[1,"N"], avg.a, avg.he, avg.ho, avg.fis, avg.nm, avg.pgt, se.a, se.he, se.ho, se.fis, se.nm, se.pgt)
      OUT1 <- as.data.frame(OUT1, stringsAsFactors=F)
      colnames(OUT1) <- c("N","A","He","Ho","Fis","N.Missing","P.GT","A.se","He.se","Ho.se","Fis.se","N.Missing.se","P.GT.se")
      OUT1$Pop <- pops[i]
      OUT1 <- move.me(OUT1,"Pop","first")
      OUT <- rbind(OUT,OUT1)
    }
    return(OUT)
  }
}


#################
### six.char()
#################

#Description
# This fuction takes genotypes in three character form and puts them into six character form

#Input Parameters:

# population - a table with sample.name in column one, and loci after that. All loci must be in three character form. Examples 000, 111, 999
six.char <- function(population, missing.data = "0", one.pop=T) {
  
  missing.data <- as.numeric(missing.data)
  missing.data <- paste0(missing.data,missing.data)
  missing.data1 <- paste0(missing.data,missing.data,missing.data) 
  cnames <- colnames(population)
  
  if(one.pop==T){  
    pop.names <- population[,1]
    population=population[,-1]  
  } else{
    pop.names <- population[,1]
    id.names <- population[,2]
    population=population[,-c(1,2)]
  }

  L<-1:ncol(population)
  locus.locs <- L[seq(1,ncol(population),2)]
  locus.names <- colnames(population)
  
  head(population)
  OUT <- NULL
  for (x in locus.locs) {
    output=paste(population[,x],population[,x+1],sep="")
    OUT <- cbind(OUT,output)
  }
  
  locus.names <- locus.names[seq(1,length(locus.names),2)]
  
  head(OUT)
  if(one.pop==T){
    OUT <- data.frame(pop.names,OUT)
    colnames(OUT) <- c(cnames[1],locus.names)
  } else {
    OUT <- data.frame(pop.names,id.names,OUT)
    head(OUT)
    colnames(OUT) <- c(cnames[1:2],locus.names)
  }

  OUT <- data.frame(lapply(OUT, as.character), stringsAsFactors=FALSE)
  OUT[OUT == missing.data] <- missing.data1
  head(OUT)
  return(OUT)
}

###################
### dup.identifer()
###################

#Description
# This fuction identifies any genotypes that are the same and reports which ones are

#Input Parameters:

# tmp - a table with sample.name in column one, and loci after that
# If the remove.dups = T, then the first genotype in the file is kept and all others that
# have the same genotype are removed.


dup.identifer <- function(tmp, one.pop = T){
  
  if(one.pop == T){
    #this pastes all the genotypes for each locus together  
    tmp$gt<- apply(tmp[,-1],1,paste,collapse="")  
  } else {
    #this pastes all the genotypes for each locus together  
    tmp$gt<- apply(tmp[,-c(1,2)],1,paste,collapse="")      
  }

  
  #the duplicated function is used to duplicated() is used to identify any duplicated genotypes
  tmp$logic <- duplicated(tmp$gt)
  
  if(one.pop == T){
    #this grabs the duplicated ids
    dup.ids <- tmp[tmp$logic == T, 1]
    dup.ids  
  } else {
    #this grabs the duplicated ids
    dup.ids <- tmp[tmp$logic == T, 2]
    dup.ids      
  }

  head(tmp)
  if(length(dup.ids) == 0){
    
    #cleaning up df
    tmp$gt <- NULL
    tmp$logic <- NULL
    print("There were no duplicated genotypes present in data")
    return(tmp)
    
  } else {
   
    #grabbing the genotypes that were genotypes
    dup.gts <- tmp[tmp$logic == T,]$gt
    dup.gts <- unique(dup.gts)
    ids <- paste0(rep("dup.group."),1:length(dup.gts))
    
    tmp$ids <- "unique"
    for(i in 1:length(dup.gts)){
      
      #now changing the logic to T for all genotypes that were duplicated, which is how to find all
      #genotypes that were duplicated
      tmp$ids[tmp$gt ==  dup.gts[i]]<- ids[i]
    }
    
    #cleaning up df
    tmp$gt <- NULL
    tmp$logic <- NULL
    return(tmp)
  } 
}


#################
### three.char()
#################

#Description
# This fuction takes genotypes in six character form and puts them into three character form
#Input Parameters:

# population - a table with sample.name in column one, and loci after that. All loci must be in six character form. Examples 000000, 111111, 999999

three.char <- function(tmp, one.pop=T){
  
  if(one.pop==T){
    
    #saving sample names into a seperate file
    OUT <- data.frame(tmp[,1])
    colnames(OUT) <- colnames(tmp)[1]
    
    #just grabbing all the genotypes with no sample names
    gts <- tmp[,-1]  
  } else{
    
    #saving sample names into a seperate file
    OUT <- data.frame(tmp[,1:2])
    colnames(OUT) <- colnames(tmp)[1:2]
    
    #just grabbing all the genotypes with no sample names
    gts <- tmp[,-1:-2]
  }

  
  #counting the number of loci and saving their names
  L <- ncol(gts)
  Locus.names <- colnames(gts)
  
  #using a small for loop to take each column and split it into two columns and cbinding that to OUT, which originally just had the sample names in it
  for(i in 1:ncol(gts)){
    locus <- data.frame(gts[,i])
    colnames(locus) <- Locus.names[i]
    
    locus$x1 <- substr(locus[,1],1,3)
    locus$x2 <- substr(locus[,1],4,6)
    locus[,1] <- NULL
    
    colnames(locus) <- c(paste(Locus.names[i],".1", sep=""),paste(Locus.names[i],".2",sep=""))
    OUT <- cbind(OUT,locus)
  }
  OUT[,1] <- as.character(OUT[,1])
  return(OUT)
}

####################
### m.data.counter()
####################

#Description
# This fuction counts how many loci are missing data for each sample
# It requires data that are in three character form

# Input parameters
# 1) tmp - table with genotypes in three character form and sample names in column one
# 2) opt - options that can be set to 1 or 2
#             Option 1 - returns a summary of the number of missing loci in dataset
#             Option 2 - returns entire dataset with a new column that counts the number of loci missing data
# 3) missing.data - default "000", but if its different such as "999", this can be changed


m.data.counter <- function(tmp, opt = 1, missing.data = "000", one.pop=T){
  
  #saving original genotypes
  tmp1 <- tmp
  
  if(one.pop==T){
    
    #replacing missing data with 0 and making entire data frame numeric
    tmp[tmp == missing.data] <- "0"
    tmp <- data.frame(tmp[,1],lapply(tmp[,-1], as.numeric))
    colnames(tmp) <- c("Sample.Names",colnames(tmp)[-1])
    
    #saving the names for each sample
    ids <- setNames(data.frame(tmp[,1],stringsAsFactors=F),"Sample.Name")
    
    #removing the names from the original data frame
    tmp <- tmp[,-1]
  
  } else {
    
    #replacing missing data with 0 and making entire data frame numeric
    tmp[tmp == missing.data] <- "0"
    tmp <- data.frame(tmp[,c(1,2)],lapply(tmp[,-c(1,2)], as.numeric))
    colnames(tmp) <- colnames(tmp1)
    
    #saving the names for each sample
    ids <- data.frame(tmp[,c(1,2)],stringsAsFactors=F)
    colnames(ids) <- colnames(tmp1)[1:2]
    
    #removing the names from the original data frame
    tmp <- tmp[,-c(1,2)]
  }
  
  #identifying the missing data
  tmp[tmp == 0] <- -10000 
  tmp[tmp != -10000] <- 0 
  tmp[tmp == -10000] <- 1
  
  #calculating the number of missing loci
  tmp$sum <- rowSums(tmp[,c(1:ncol(tmp))])
  tmp$sum <- tmp$sum/2
  
  tmp1$m.data <- tmp$sum
  tmp1[,1] <- as.character(tmp1[,1])
    
  return(tmp1)
}


#################
### uniq.alleles
#################

#Description
# This fuction takes output from alf.freqs() and The program compares
# the alleles at each locus between the two programs and outputs only the alleles
# that are unique to each population along with their counts and frequencies.

#Input Parameters:

# alfs - requires output from alf.freqs() and the data.frame used in that function
# had only two populations


uniq.alleles <- function(alfs) {
  
  x1 <- NULL
  y1 <- NULL
  
  loci <- as.character(unique(alfs$Locus))
  
  for(i in 1:length(loci)){
    alfs1 <- alfs[alfs$Locus == loci[i],]
    
    pops <- as.character(unique(alfs1$Pops))
    
    pop1 <- alfs1[alfs1$Pops == pops[1],]
    pop1 <- as.numeric(as.character(pop1$allele))
    
    pop2 <- alfs1[alfs1$Pops == pops[2],]
    pop2 <- as.numeric(as.character(pop2$allele))
    
    x <- data.frame(pop1,match(pop1,pop2))
    colnames(x) <- c("alleles","matches")
    x <- x[is.na(x$matches) == T,1]
    x <- alfs[alfs$Pops == pops[1] & alfs$Locus == loci[i] & alfs$allele %in% x, ]
    x1 <- rbind(x1,x)
    
    y <- data.frame(pop2,match(pop2,pop1))
    colnames(y) <- c("alleles","matches")
    y <- y[is.na(y$matches) == T,1]
    y <- alfs[alfs$Pops == pops[2] & alfs$Locus == loci[i] & alfs$allele %in% y, ]
    y1 <- rbind(y1,y)
  }
  
  OUT <- rbind(x1,y1)
  return(OUT)
}

#################
### st.err()
#################

#Description
# This fuction calculates standard error

st.err <- function(x){ sd(x)/sqrt(length(x))}

#################
### move.me()
#################

#Description
#This fuction reorders columns in your dataframe in anyway you want. It was orignally developed by Ananda Mahto and found it
#on stackoverflow at this link: http://stackoverflow.com/questions/18339370/reordering-columns-in-a-large-dataframe. By default
#it takes the column(s) you want and moves them to the end of the file, but if designated can move them to the front or infront of
#any specified column

#Input Parameters:
#data - a dataframe
#tomove - can be single or multiple character string
#ba - tells the program what columns to stick tomove columns in front of

move.me <- function(data, tomove, where = "last", ba = NULL) {
  temp <- setdiff(names(data), tomove)
  x <- switch(
    where,
    first = data[c(tomove, temp)],
    last = data[c(temp, tomove)],
    before = {
      if (is.null(ba)) stop("must specify ba column")
      if (length(ba) > 1) stop("ba must be a single character string")
      data[append(temp, values = tomove, after = (match(ba, temp)-1))]
    },
    after = {
      if (is.null(ba)) stop("must specify ba column")
      if (length(ba) > 1) stop("ba must be a single character string")
      data[append(temp, values = tomove, after = (match(ba, temp)))]
    })
  x
}

#################
### sim.pop.gts()
#################

#Description
# This fuction creates possion distributed random data to represent microsatellite alleles. 
#  It will repeat that process for however many individuals, loci, and populations you with.

# Input parameters
# 1) n - the number of individuals per population
# 2) loci - the number of loci desired
# 3) missing.data - can put some percentage of missing data randomly distributed thoroughout the data
# 4) pops - the number of populations desired
# 5) k - mean number of alleles per loci
# 6) k.sd - standard devation of number of alleles per loci


#making the function
sim.pop.gts <- function(n=1000,loci=10,missing.data=0, pops=2, k=20, k.sd=3){
  
  #Pops
  Pops <- rep("POP",times=n*pops)
  pop.ids <- rep(1:pops,each=n)
  Pops <- paste(Pops, pop.ids, sep= ".")
  
  output <- NULL
  for(j in 1:pops){
    
    #IDs
    nums <-1:n
    IDs <- rep("ID.",times=n)
    IDs <- paste0(IDs,nums)
    

    i <- NULL
    OUT <- NULL
    for(i in 1:loci){
      
      #getting the alleles for each locus
      x <- rpois(n=100000, lambda=115)
      x <- sort(x)
      
      x <- x[x>100 & x < 200]
      xs <- unique(x)
      y <- sample(xs,size=round(rnorm(n=1,mean=k,sd=k.sd)),replace=F)
      x <- x[x %in% y]
      
      
      #sampling those alleles the number of times required to
      #make complete genotypes (allele 1 and allele 2)
      z.1 <- sample(x,size=n,replace=T)
      z.2 <- sample(x,size=n,replace=T)
      
      #combining to make a complete genotype in three character form
      OUT1 <- cbind(z.1,z.2)
      
      #making sure the smaller allele is first
      OUT1 <- as.data.frame(OUT1, stringsAsFactors=F)
      OUT1$logic <- OUT1[,1] > OUT1[,2]
      OUT1.1 <- OUT1[OUT1$logic == F,-ncol(OUT1)]
      OUT1.2 <- OUT1[OUT1$logic != F,-ncol(OUT1)]
      OUT1.2 <- OUT1.2[,c(2,1)]
      colnames(OUT1.2) <- colnames(OUT1.1)
      OUT1 <- rbind(OUT1.1,OUT1.2)
      
      if(missing.data != 0){
        #making missing data
        num.md <- missing.data*n
        md.locs <- sample(1:n,size=num.md,replace=F)
        OUT1[md.locs,1] <- "000" 
        OUT1[md.locs,2] <- "000"
      }
      OUT1 <- as.matrix(OUT1)
      OUT <- cbind(OUT,OUT1)
    } # end of for loci loop 
    
    df <- data.frame(IDs,OUT,stringsAsFactors=F)
    loc.names <- rep("Locus",times=2*loci)
    loc.nums <- rep(1:loci, each=2)
    loc.alleles <- rep(1:2,times=loci)
    loc.names <- paste(loc.names,loc.nums,sep="_")
    loc.names <- paste(loc.names,loc.alleles,sep=".")
    colnames(df) <- c("Sample.name",loc.names)
    
    output <- rbind(output,df)
  } #end of Pops loop
  
  #getting rid of any factors!
  output <- data.frame(Pops,output,stringsAsFactors=F)
  return(output)
} #end of sim.pop.gts



####################
### genepop.format()
####################

#Description
#This fuction takes three character genotypic data with a pop and sample id column in the beginning and puts it into 
#a format acceptable to genepop

#Input Parameters:
#genos - a dataframe with two columns - "pop" and "sample.name", then the three character genotypic data
#title - the name of the data frame you are putting into genepop - Default: "Genotypic data"

genepop.format <- function(genos, title="Genotypic data"){
  
  #gts names for genepop
  genos[,1] <- paste(genos[,1], ",")
  head(genos)
  
  #now for the genepop prep
  my.sep <- genos
  my.sep[1,] <- c("Pop",rep("",(ncol(genos)-1)))
  my.sep <- my.sep[1,]
  my.sep
  
  #getting the unique pop identifiers
  my.ids <- unique(genos[,1])
  
  tmp2 <- NULL
  for(i in 1:length(my.ids)){
    x <- genos[genos[,1] == my.ids[i],]
    tmp2 <- rbind(tmp2,my.sep,x)
  }
  
  #just a check to make sure I have eight pops
  nrow(tmp2[tmp2[,1] == "Pop",])
  
  #making into six.char()
  tmp2 <- six.char(tmp2,missing.data="000",one.pop=F)
  head(tmp2)
  
  #removing the Pop ID
  tmp2 <- tmp2[,-2]
  head(tmp2)
  
  
  #prepping for the header for the genepop file
  loci <- c(paste0(colnames(tmp2)[-c(1, length(colnames(tmp2)))],", "),paste0(colnames(tmp2)[length(colnames(tmp2))]), "")
  loci
  
  head(tmp2)
  ncol(tmp2)
  my.sep <- tmp2
  my.sep[1,] <- c(title,rep("",(ncol(tmp2)-1)))
  my.sep <- my.sep[1,]
  my.sep
  
  tmp2 <- rbind(my.sep,loci,tmp2)
  
  
  print("Writing Genepop formatted data to current working directory")
  #writing to file now
  write.table(tmp2,"genepop.gts.txt",append=F,quote=F,sep="\t",row.names=F,col.names=F)
}

