#################################################Estimating population genetics parameters for infrapopulations#############################################################

#This function allow to simulate 100 (or whatever) sampling scenarios. For example, you have some infrapopulations (i.e. number of symbiont individuals per individual host) with more than one individual sampled and other with only one. Here, in each iteration, only use one individual per infrapopulation. Then, is prepared to calculate Genearalized linear models (Infraglm function) for your variables of interest.
#Sometimes, could be realistic work with uncertainty in your population parameters (e.g. abundance). This function allows to calculate a new value into a range per iteration. In this example, Minimum, Median, Quartile 95 and Maximun are simulated between the 50% and the 100% of the "real" value.


#Preparing data

You need two input files. The alignmentcsv (DNA sequences as characters) and the pop_parameterscsv (Population parameters ordered as alignment-i.e. Column order have to be the same as species order in alignment-)
You can find two examples files in this same repository; Aligment_example.csv and Pop_par_example.csv


#Output file

The output.csv includes the results for 100 iterations.

Genetic parameters estimated: Tajima's D, R2, nucleotide diversity (pi), haplotype number and haplotype proportion.
Population parameters simulated: Minimum, Median, Quartile 95 and Maximun.

In addition, sequencesused.csv includes all the sequences used in the simulation.

require (ape) 
require (pegas)
require (dplyr)

###################################################################################################################

Infrapopgenetics <- function(alignmentcsv, pop_parameterscsv, iteration) {

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
  
  for (it in 1:iteration){   #The total number of iterations
          
                    for(i in 1:17)#  The total number of species in the file data
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

This function calculate one Genearalized linear model over the data generated from the previous simulation (i.e. one per iteration).
In this example, is wrote for do a Gaussian model (link = "identity") accounting for the differences in sample size (weights option over an n column)
#Input data

data=read.csv("output.csv") #From the previous function
*If you are interested in use the weight option of glm, you have to include a column named "n" with the sample size in your data.

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
              glm=glm(V3 ~ QUANTILE50, data = sub,  family = gaussian(link = "identity"), weights = n) #PREGUNTAR A ROGER
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
  write.csv(results, "pi_glm_final.csv")
  

}
