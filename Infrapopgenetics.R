#################################################Estimating population genetics parameters for infrapopulations#############################################################

This function allow to simulate 100 (or whatever) sampling scenarios. This is a way to use all your data when you have some infrapopulations (i.e. number of symbiont individuals per individual host) better sampled than others, by using one random individual per infrapopulation each time for iteration calculations. Later you will able to compare iteration with general results. An example, you have some infrapopulations with 5 individuals sampled and other with only one. Here, in each iteration, only one individual is used per infrapopulation (always the same for the infrapopulation with one individual). Then, you can use the "Infraglm function" (see below) to calculate Genearalized linear models for over your variables of interest.
*In some circumstances, can be realistic work with uncertainty in your variables (e.g. abundance). This function allows to calculate a new value into a range per iteration. In this example, Minimum, Median, Quartile 95 and Maximun are simulated between the 50% and the 100% of the "real" value.


#Preparing data

Take a look over two examples files here (https://www.dropbox.com/sh/kctmuslqg57cjwj/AACwSplUtAUnhV9YSjdXJ17ga?dl=0)
Aligment_example.csv and Pop_par_example.csv

The alignmentcsv (DNA sequences as characters) and the pop_parameterscsv (Population parameters ordered as alignment-i.e. Raw order have to be the same as species order in alignment-)


#Output file

The output.csv includes the results for 100 iterations.

In this example:

Genetic parameters estimated: Tajimas D, R2, nucleotide diversity (pi), haplotype number and haplotype proportion.
Population parameters simulated: Minimum, Median, Quartile 95 and Maximun.

In addition, sequencesused.csv include all the sequences used in the simulation.


###################################################################################################################

require (ape) 
require (pegas)
require (dplyr)
##

Arguments:

iteration #The total number of iterations
alignmentcsv #Input DNA file (see above)
pop_parameterscsv# Input population parameters file (see above)
species# Number of species in the input files

###

Infrapopgenetics <- function(alignmentcsv, pop_parameterscsv, iteration, species) { 

  data=read.csv("alignmentcsv") 
  abun=read.csv("pop_parameterscsv") 
  indvectorit=numeric()
  seq=numeric()
  pi=numeric()
  pisp=numeric()
  h=numeric()
  prop=numeric()
  d=numeric()
  pvd=numeric()
  ramon=numeric()
  orsin=numeric()
  aligned= numeric()
  prop= numeric()
  median= numeric()
  mean= numeric()
  q95= numeric()
  max= numeric()
  
  for (it in 1:iteration){   
          
                    for(i in 1:species)
                      { 
                            
                            a=data.frame(filter(data, sp ==i)) 
                            indvector = unique (a$ind,comparable = TRUE) 
                            indvector= as.matrix(indvector)
                            N=length(indvector) 
                            birdsp=numeric()
                                    
                                for (b in 1:N){ 
                          
                                    r= subset(a, ind %in% c(b)) 
                                    m=as.data.frame(r)
                                    bird=sample_n(m,1,replace=TRUE) 
                                    birdsp=rbind(birdsp, bird) 
                                   
                                   }
                           
                            birdm=as.matrix(birdsp)
                            birdmm= birdm[,-1]  
                            birdalignment= birdmm[,-1]         
                            align=as.DNAbin(birdalignment)
                            pi= nuc.div(align)                          
                            haplo=haplotype(align)
                            h=nrow(haplo)
                            prop=h/nrow(align)
                            t= tajima.test(align)
                            d= t$D
                            pvd=t$Pval.normal
                            r=R2.test (align, plot=FALSE, quiet = TRUE)
                            ramon= r$R2
                            orsin= r$P.val
                            m=  runif(1, 0.5, 1) 
                            mediancolum=select(abun,QUANTILE50)
                            median1=slice(mediancolum,i)
                            median= median1*m
                            meancolum= select(abun,MEAN)
                            mean1=slice(meancolum,i)
                            mean=mean1*m
                            q95colum= select(abun,QUANTILE95)
                            q951=slice(q95colum,i)
                            q95=q951*m
                            maxcolum=select(abun,MAX)
                            max1=slice(maxcolum,i)
                            max=max1*m
                            pisp= rbind(pisp,c(it,i,pi, d,pvd,ramon,orsin,h,prop,median,mean,q95,max))
                            alignment=  as.character(align)
                            al=as.data.frame(alignment)                          
                            v= 1:N
                            al$new.col <- v
                            repet= rep(i,N)
                            al$new.cole <- repet
                            aligned= rbind(aligned, al) 
                            
                                 
                      }
   write.csv(pisp, "output.csv")
   write.csv (aligned, "sequencesused.csv")               
  }
#################################################################################################################

#Infraglm function

This function calculate a Genearalized linear model over the data generated from the previous simulation (i.e. one model per iteration).
*In this example, is wrote for a Gaussian model (link = "identity") accounting for the differences in sample size (weights option over an sample size column manually added to the output.csv file)

#Input data

data=read.csv("output.csv") #From the previous function
*If you are interested in use the weight option of glm, you have to nclude a column named "n" with the sample size in your data.

#Output data

The output file contains:

iteration:the number of iteration
pvalue: the model pvalue
intcp: the intercept
slp: the slope
f: the F statistic
es: the Nagelkerke R2
df: the degree of freedom
dev: the deviance explained by the model.


Infraglm <- function(input,iteration) {

data=read.csv("input")
pvalue=numeric()
intcp=numeric()
slp=numeric()
iteration=numeric()
f=numeric()
es=numeric()
df=numeric()
dev=numeric()

      for (it in 1:iteration){ 
      
              sub1 = filter(data, V1==it)
              sub=as.data.frame(sub1)
              iteration=c(it,iteration)
              glm=glm(V3 ~ QUANTILE50, data = sub,  family = gaussian(link = "identity"), weights = n) 
              int = coef(glm)[1] 
              intcp= c(intcp,as.numeric(int))
              sl = coef(glm)[2]
              slp= c(slp,as.numeric(sl))
              pval <- summary(glm)$coef[, "Pr(>|t|)"]
              pvalmat=as.matrix(pval)
              pv=pvalmat[2,]
              pvalue=c(pvalue,pv)
              deviance=glm$deviance
              null=glm$null.deviance
              D2= ((null-deviance)/(null))*100 
              dev=c(dev,D2)
              anova=anova(glm,test="F")
              f1=anova$F[2]
              f=c(f,f1)
              df1=anova$"Resid. Df"[2]
              df=c(df,df1)
              r2=NagelkerkeR2(glm)
              es1=r2$R2
              es=c(es,es1)

             
                          
            }  
  results= cbind(iteration, pvalue, intcp, slp,f,es,df,dev)
  write.csv(results, "glmresults.csv")
  

}
