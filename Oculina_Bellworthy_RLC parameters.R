library(tidyverse)
library(lubridate)
library(ggplot2)  
library(phytotools)
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)
library(car)
library(plotrix)
library(cowplot)


#calculating RLC parameters using phytotools (pratt 1980)
#setwd
d <- read.csv('rETR_OculinaProject.csv')

# Set column order PAR - ETR - ID
rlc.data <- d[c(3,4,2)]
rlc.data$colony <- as.factor(rlc.data$colony)

rlc.data<- setNames(rlc.data, c("par","etr", "id"))
head(rlc.data)

rlc.data$etr <- na_if(rlc.data$etr, 0)
View(rlc.data)

ncurves <- length(unique(rlc.data$id)) # number of unique ids in the data 
ids <- unique(rlc.data$id) # store the unique ids 

str(rlc.data)

rlc.parameters <- data.frame(
  id = ids, 
  alpha = 0, 
  beta = 0, 
  ETRmax = 0, 
  Ek = 0, 
  ps = 0
)


for (i in 1:ncurves){
  
  temp.id = ids[i] # extract the id of the curve to be fitted
  
  print(paste("Now fitting curve ", as.character(temp.id))) # to keep track what's happening if the data has many curves
  
  temp.rlc.data <- rlc.data[rlc.data$id==temp.id,] # extract the the data of a single curve into a temporary variable
  PAR = temp.rlc.data$par 
  ETR = temp.rlc.data$etr
  
  fit = fitPGH(PAR, ETR, fitmethod = "Port") # for more options and explanation see package phytotools manual
  
  # store the fitted RLC values into temporary variables
  alpha.rlc = fit$alpha[1]
  beta.rlc = fit$beta[1]
  ps.rlc = fit$ps[1]
  
  # store the parameters
  rlc.parameters$id[i] <- temp.id
  rlc.parameters$alpha[i] <- alpha.rlc
  rlc.parameters$beta[i] <- beta.rlc
  rlc.parameters$ps[i] <- ps.rlc
  
  # calculate ETRmax and Ek for the PGH model (see e.g.Ralph & Gademann 2005 Aquatic Botany 82 (3): 222 - 237). 
  # Note that the equation depends on the model fitted, the code below applies only to the PGH model! 
  # Model equations are documented in the phytotools package code examples (and in the original papers): https://cran.r-project.org/web/packages/phytotools/phytotools.pdf
  
  ETRmax = ps.rlc*(alpha.rlc/(alpha.rlc + beta.rlc))*(beta.rlc/(alpha.rlc+beta.rlc))^(beta.rlc/alpha.rlc)
  Ek = ETRmax/alpha.rlc 
  
  # store the variables
  rlc.parameters$ETRmax[i] <- ETRmax
  rlc.parameters$Ek[i] <- Ek
  
  #plotting the curve and fitted model into a tiff file. By default the file name is the id of the curve. 
  tiff(file=paste0(temp.id, ".tiff"), compression="lzw")
  
  # plot the data, 
  plot(x=PAR, y=ETR, main=temp.id) 
  
  # plot the model fit
  with(fit, {
    P <- ps.rlc*(1-exp(-1*alpha.rlc*PAR/ps.rlc))*exp(-1*beta.rlc*PAR/ps.rlc) # the PGH model equation
    lines(PAR,P)
  }
  ) # end of with
  dev.off() #close the plotting devide. if this is not done, the next run of the loop will override the plot. 
  
}

# now the data frame rlc.parameters contains the fitted values for each curve. Tiff plots should be in current working directory. 
rlc.parameters
View(rlc.parameters)

warnings()

write.csv(file="Oculina_rlc.parameters.csv",rlc.parameters)


####check statistical differences
### Manually copy the FVFM to the rlc parameter data  
d <- read.csv('Oculina_rlc.parameters_NEW.csv')

head(d)
View(d)
d$depth <- as.factor(d$depth)
d$id <- as.factor(d$id)

sh = subset(d, depth == 'shallow')
deep = subset(d, depth == 'deep')

#etrMAX 
shapiro.test(sh$ETRmax)
shapiro.test(deep$ETRmax)
leveneTest(ETRmax~depth,d=d)
t.test(d$ETRmax~d$depth)

#alpha
shapiro.test(sh$alpha)
shapiro.test(deep$alpha)
leveneTest(alpha~depth,d=d)
t.test(d$alpha~d$depth)

#eK
shapiro.test(sh$Ek) # not normal
shapiro.test(deep$Ek)
leveneTest(Ek~depth,d=d)
kruskal.test(d$Ek~d$depth)

#FV/FM
shapiro.test(sh$fvfm)
shapiro.test(deep$fvfm) 
leveneTest(fvfm~depth,d=d)
t.test(d$fvfm~d$depth)


###beta
shapiro.test(sh$beta) #not normal
shapiro.test(deep$beta) #not normal
leveneTest(beta~depth,d=d)
kruskal.test(d$fvfm~d$depth)


## summary of results 
Sum_all <- d %>% 
  group_by(depth) %>% 
  summarise_each(funs(mean(., na.rm=TRUE), n = sum(!is.na(.)),
                      se = sd(., na.rm=TRUE)/sqrt(sum(!is.na(.)))), alpha:fvfm)

View(Sum_all)
write.csv(file="Oculina_rlc.summary.csv", Sum_all)





#### Graphics in Box plots  #####
mytheme = theme_classic()+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 6), axis.title = element_text(size = 7))

mynames <- c( "Shallow","Deep")
d$depth = factor(d$depth, levels = c("shallow", "deep"))

#FVFM
fvfm = ggplot(d, aes(y = fvfm, x = depth)) +
  geom_boxplot(aes(fill = depth), outlier.shape = "", fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.5)+
  labs(y= ~F[V] ~F[M])+
  scale_x_discrete(labels = mynames)+
  scale_fill_manual(values = c("#D55E00","#0072B2"))+
  scale_y_continuous(limits = c(0.5, 0.7), breaks = seq(0.5,0.7,0.1))+
  mytheme+
  guides(fill = FALSE)

fvfm


#ETRmax
ETRmax = ggplot(d, aes(y = ETRmax, x = depth)) +
  geom_boxplot(aes(fill = depth), outlier.shape = "", fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.5)+
  labs(y= ~ETR[MAX])+
  scale_x_discrete(labels = mynames)+
  scale_fill_manual(values = c("#D55E00","#0072B2"))+
  scale_y_continuous(limits = c(20, 110), breaks = seq(20,100,20))+
  mytheme+
  guides(fill = FALSE)

ETRmax


#Ek
Ek = ggplot(d, aes(y = Ek, x = depth)) +
  geom_boxplot(aes(fill = depth), outlier.shape = "", fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.5)+
  labs(y= ~E[K] ~(µmol ~m^-2 ~s^-1))+
  scale_x_discrete(labels = mynames)+
  scale_fill_manual(values = c("#D55E00","#0072B2"))+
  scale_y_continuous(limits = c(50, 150), breaks = seq(50,150,20))+
  mytheme+
  guides(fill = FALSE)

Ek

#ALPHA
alpha = ggplot(d, aes(y = alpha, x = depth)) +
  geom_boxplot(aes(fill = depth), outlier.shape = "", fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.5)+
  labs(y= "alpha")+
  scale_x_discrete(labels = mynames)+
  scale_fill_manual(values = c("#D55E00","#0072B2"))+
  scale_y_continuous(limits = c(0.4, 0.8), breaks = seq(0.4,0.8,0.20))+
  mytheme+
  guides(fill = FALSE)

alpha


#BETA
beta =ggplot(d, aes(y = beta, x = depth)) +
  geom_boxplot(aes(fill = depth), outlier.shape = "", fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.5)+
  labs(y= "beta")+
  scale_x_discrete(labels = mynames)+
  scale_fill_manual(values = c("#D55E00","#0072B2"))+
  scale_y_continuous(limits = c(0, 5), breaks = seq(0,5,1))+
  mytheme+
  guides(fill = FALSE)

beta


# combine plots and save
rlc.plots <- plot_grid(fvfm,ETRmax,Ek,alpha, beta, labels = c('A', 'B', 'C', 'D', 'E'),
                       label_y = 0.985, label_size = 11,ncol = 5, align = "h", byrow = F)

ggsave("PAM_OculinaProject2020.jpeg", plot = rlc.plots, width = 16, height = 4,dpi=300, 
       units = "cm")
ggsave("PAM_OculinaProject2020.pdf", plot = rlc.plots, width = 16, height = 4,dpi=300, 
       units = "cm")



# add the physiology plots to the PAM
all.plots <- plot_grid(physio.plots, rlc.plots2, ncol = 1, align = "h", byrow = F)
ggsave("allplots_OculinaProject2020.pdf", plot = all.plots, width = 16, height = 8,dpi=300, 
       units = "cm")
ggsave("allplots_OculinaProject2020.jpeg", plot = all.plots, width = 16, height = 8,dpi=300, 
       units = "cm")



