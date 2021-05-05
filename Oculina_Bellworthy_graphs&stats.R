# setwd
physio <- read.csv('Oculina_physio.csv')
View(physio)
physio$depth = factor(physio$depth, levels = c("shallow", "deep"))


library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(car)
library(cowplot)

mytheme = theme_classic()+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 6), axis.title = element_text(size = 7))


#### Graphics in Box plots
p.protein = ggplot(physio, aes(y = protein.cm, x = depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = depth), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.5)+
  labs(y= µg ~protein ~cm^-2)+
  mytheme +
  scale_fill_manual(values = c("#D55E00","#0072B2"))+
  scale_x_discrete(labels = c("Shallow", "Deep"))+
  guides(fill = FALSE)
p.protein

p.cal = ggplot(physio, aes(y = calcification_cm, x = depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = depth), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.5)+
  labs(y= µmol ~CaCO[3] ~cm^-2 ~h^-1)+
  mytheme +
  scale_fill_manual(values = c("#D55E00","#0072B2"))+
  scale_x_discrete(labels = c("Shallow", "Deep"))+
  guides(fill = FALSE)
p.cal

p.cellcm = ggplot(physio, aes(y = cell.cm, x = depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = depth), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.5)+
  labs(y = ~cells ~x10^6 ~cm^-2)+
  mytheme +
  scale_fill_manual(values = c("#D55E00","#0072B2"))+
  scale_x_discrete(labels = c("Shallow", "Deep"))+
  guides(fill = FALSE)
p.cellcm

p.chl = ggplot(physio, aes(y = chl.cm, x = depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = depth), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.5)+
  labs(y= ug ~chlorophyll[a] ~cm^-2)+
  mytheme +
  scale_fill_manual(values = c("#D55E00","#0072B2"))+
  scale_x_discrete(labels = c("Shallow", "Deep"))+
  guides(fill = FALSE)
p.chl

p.chlcell = ggplot(physio, aes(y = chl.cell, x = depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = depth), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.5)+
  labs(y= pg ~chlorophyll[a] ~cell^-1)+
  mytheme +
  scale_fill_manual(values = c("#D55E00","#0072B2"))+
  scale_x_discrete(labels = c("Shallow", "Deep"))+
  guides(fill = FALSE)+
  scale_y_continuous(limits = c(0,6))
p.chlcell


# Combine plots to grid and save
physio.plots <- plot_grid(p.protein, p.cal, p.cellcm, p.chl, p.chlcell, labels = c('A', 'B', 'C', 'D', 'E'),
                   label_y = 0.985, label_size = 11,ncol = 5, align = "h", byrow = F)
physio.plots
ggsave("Physio_OculinaProject2020.jpeg", plot = physio.plots, width = 16, height = 4,dpi=300, 
       units = "cm")
ggsave("Physio_OculinaProject2020.pdf", plot = physio.plots, width = 16, height = 4,dpi=300, 
       units = "cm")



##### Stats ########
# Calcification rate
shapiro.test(physio$calcification_cm)  
leveneTest(calcification_cm~depth,d=physio) 
t.test(calcification_cm~depth, data=physio)

# chlorophyll cm
shapiro.test(physio$chl.cm) # not normal
leveneTest(chl.cm~depth,d=physio) 
kruskal.test(chl.cm~depth, data=physio)

# cellcm
shapiro.test(physio$cell.cm)  
leveneTest(cell.cm~depth,d=physio) 
t.test(cell.cm~depth, data=physio)

#chl/ cell
shapiro.test(physio$chl.cell) # not normal
leveneTest(chl.cell~depth,d=physio) 
kruskal.test(chl.cell~depth, data=physio)

# Protein.cm
shapiro.test(physio$protein.cm)  
leveneTest(protein.cm~depth,d=physio) 
t.test(protein.cm~depth, data=physio)



############# Isotopes ###############
iso = read.csv("Oculina_isotopes.csv")
View(iso)
str(iso)
head(iso)

iso$Part_Dis = as.factor(iso$Part_Dis)
iso$Host_Sym_Tot = as.factor(iso$Host_Sym_Tot)
iso$depth = as.factor(iso$depth)

iso_names = c('d' = "Dissolved", 'p' = "Particulate")
iso$depth = factor(iso$depth, levels = c("shallow", "deep"))

iso = subset(iso, Host_Sym_Tot != 't') # remove the total component from the data

# NITROGEN
# conduct t.test to add to graphs
stat.test <- iso %>%
  group_by(Part_Dis, Host_Sym_Tot) %>%
  t_test(N_cm ~ depth) %>%
  add_significance() %>%
  p_round(digits = 3)
stat.test = stat.test%>% add_xy_position(x = "Host_Sym_Tot")
stat.test

#plot in boxplot
p.Nitrogen = ggplot(iso, aes(y = N_cm, x = Host_Sym_Tot)) +
  geom_boxplot(outlier.colour = "NA", fatten = 0.5, alpha = 0.7, lwd = 0.2, (aes(fill = depth)))+
 geom_jitter(position=position_dodge(width=0.75),aes(group=depth), size = 0.7)+
  labs(y= nmol ~N ~cm^-2 ~h^-1)+
  theme_classic()+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black", size =  6), axis.title = element_text(size = 7), 
        legend.title = element_blank(), strip.text = element_text(size = 8),
        strip.background = element_rect(size = 0.5, colour = NA, fill = NA))+
  scale_fill_manual(labels = c("Shallow","Deep"), values = c("#D55E00","#0072B2"))+
  scale_x_discrete(labels = c("Host", "Symbiont"))+
  facet_wrap("Part_Dis", scales = "free_y", labeller = as_labeller(iso_names))+
  stat_pvalue_manual(stat.test, label = "{p.signif} {p}", tip.length = 0.005, hide.ns = TRUE, size = 1.7)+
  guides(fill= FALSE)
p.Nitrogen


# CARBON
# conduct t.test to add to graphs
stat.test.C <- iso %>%
  group_by(Part_Dis, Host_Sym_Tot) %>%
  t_test(C_cm ~ depth) %>%
  add_significance() %>%
  p_round(digits = 3)
stat.test.C = stat.test.C%>% add_xy_position(x = "Host_Sym_Tot")
stat.test

#plot in boxplot
p.Carbon = ggplot(iso, aes(y = C_cm, x = Host_Sym_Tot)) +
  geom_boxplot(outlier.colour = "NA", fatten = 0.5, alpha = 0.7, lwd = 0.2, aes(fill = depth))+
  geom_jitter(position=position_dodge(width=0.75),aes(group=depth), size = 0.7)+
  labs(y= nmol ~C ~cm^-2 ~h^-1)+
  theme_classic()+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black", size =  6), 
        axis.title = element_text(size = 7), legend.title = element_blank(), 
        strip.text = element_text(size = 8),
        strip.background = element_rect(size = 0.5, colour = NA, fill = NA))+
  scale_fill_manual(labels = c("Shallow","Deep"), values = c("#D55E00","#0072B2"))+
  scale_x_discrete(labels = c("Host", "Symbiont"))+
  facet_wrap("Part_Dis", scales = "free_y", labeller = as_labeller(iso_names))+
  stat_pvalue_manual(stat.test.C, label = "{p.signif} {p}", tip.length = 0.005, hide.ns = TRUE, size = 1.7)+
  guides(fill = FALSE)
p.Carbon

#Print and save plots
iso.plots <- plot_grid(p.Carbon, p.Nitrogen, labels = c('A', 'B'),
                          label_x = 0.06, label_y = 0.985, label_size = 11,ncol = 1, align = "v", byrow = F, hjust =2)
iso.plots
ggsave("Isotopes_OculinaProject2020.jpeg", plot = iso.plots, width = 10, height = 10,dpi=300, units = "cm")
ggsave("Isotopes_OculinaProject2020.pdf", plot = iso.plots, width = 10, height = 10,dpi=300, units = "cm")



##### Additional Isotope Stats ####
# Compare between shallow and deep for N and C in each compartment

# Subset by Part_Dis and Host or symbiont compartment
Dis_H = subset(iso, Part_Dis == "d" & Host_Sym_Tot == "h") 
View(Dis_H)
Dis_S = subset(iso, Part_Dis == "d" & Host_Sym_Tot == "s") 
View(Dis_S)

Part_H = subset(iso, Part_Dis == "p" & Host_Sym_Tot == "h") 
View(Part_H)
Part_S = subset(iso, Part_Dis == "p" & Host_Sym_Tot == "s") 
View(Part_S)



### CARBON  #####
# Host
# Further subset on depth
Dis_H_s = subset(Dis_H, depth == 'shallow')
Dis_H_d = subset(Dis_H, depth == 'deep')
Part_H_s = subset(Part_H, depth == 'shallow')
Part_H_d = subset(Part_H, depth == 'deep')

## Dissolved Carbon
shapiro.test(Dis_H_s$C_cm) 
shapiro.test(Dis_H_d$C_cm) 
leveneTest(C_cm~depth,d=Dis_H)  
t.test(Dis_H$C_cm~Dis_H$depth)

## Particulate Carbon
shapiro.test(Part_H_s$C_cm) 
shapiro.test(Part_H_d$C_cm)
leveneTest(C_cm~depth,d=Part_H) 
t.test(Part_H$C_cm~Part_H$depth)


# Symbiont
# Further subset on depth
Dis_S_s = subset(Dis_S, depth == 'shallow')
Dis_S_d = subset(Dis_S, depth == 'deep')
Part_S_s = subset(Part_S, depth == 'shallow')
Part_S_d = subset(Part_S, depth == 'deep')

## Dissolved Carbon
shapiro.test(Dis_S_s$C_cm) 
shapiro.test(Dis_S_d$C_cm) 
leveneTest(C_cm~depth,d=Dis_S)
t.test(Dis_S$C_cm~Dis_S$depth)

## Particulate Carbon
shapiro.test(Part_S_s$C_cm) 
shapiro.test(Part_S_d$C_cm) 
leveneTest(C_cm~depth,d=Part_S)  
t.test(Part_S$C_cm~Part_S$depth)


#### NITROGEN  ######
# Host
## Dissolved Nitrogen
shapiro.test(Dis_H_s$N_cm)
shapiro.test(Dis_H_d$N_cm) 
leveneTest(N_cm~depth,d=Dis_H)  
t.test(Dis_H$N_cm~Dis_H$depth)

## Particulate Nitrogen
#normality assumptions 
shapiro.test(Part_H_s$N_cm)
shapiro.test(Part_H_d$N_cm) 
leveneTest(N_cm~depth,d=Part_H) 
t.test(Part_H$N_cm~Part_H$depth)


# Symbiont
## Dissolved Nitrogen
#normality assumptions 
shapiro.test(Dis_S_s$N_cm) # p = 0.0224 - NOT NORMAL
shapiro.test(Dis_S_d$N_cm) 
leveneTest(N_cm~depth,d=Dis_S)  
# Kruskal - Wallis non parametric test
kruskal.test(Dis_S$N_cm~Dis_S$depth)

## Particulate Carbon
shapiro.test(Part_S_s$N_cm)
shapiro.test(Part_S_d$N_cm) # p = 0.002681 -  NOT NORMAL
leveneTest(N_cm~depth,d=Part_S) 
# Kruskal - Wallis non parametric test
kruskal.test(Part_S$N_cm~Part_S$depth)


### Test for differences in between Host and Symbiont compartments within each depth
### SHALLOW ###
# Dissolved Carbon and Nitrogen (Host_Sym_Total)
Dis_shallow = subset(iso, Part_Dis == "d" & depth == "shallow")

# Carbon - already know all groups meet parametric assumptions from above
t.test(Dis_shallow$C_cm~Dis_shallow$Host_Sym_Tot)
# Nitrogen - 
t.test(Dis_shallow$N_cm~Dis_shallow$Host_Sym_Tot)

# Particulate Carbon and Nitrogen (Host_Sym_Total)
Part_shallow = subset(iso, Part_Dis == "p" & depth == "shallow")

# Carbon - already know all groups meet parametric assumptions from above
t.test(Part_shallow$C_cm~Part_shallow$Host_Sym_Tot)
# Nitrogen - 
t.test(Part_shallow$N_cm~Part_shallow$Host_Sym_Tot)



### DEEP ###
# Dissolved Carbon and Nitrogen (Host_Sym_Total)
Dis_deep = subset(iso, Part_Dis == "d" & depth == "deep")

# Carbon - already know all groups meet parametric assumptions from above
t.test(Dis_deep$C_cm~Dis_deep$Host_Sym_Tot)
# Nitrogen - 
t.test(Dis_deep$N_cm~Dis_deep$Host_Sym_Tot)

# Particulate Carbon and Nitrogen (Host_Sym_Total)
Part_deep = subset(iso, Part_Dis == "p" & depth == "deep")

# Carbon - already know all groups meet parametric assumptions from above
t.test(Part_deep$C_cm~Part_deep$Host_Sym_Tot)
# Nitrogen - 
t.test(Part_deep$N_cm~Part_deep$Host_Sym_Tot)



############ Trophic Position ###############
TP_PCA = read.csv("Oculina_TP_PCA.csv")
TP_PCA$Depth = factor(TP_PCA$Depth, levels = c("Shallow", "Deep"))

deep = subset(TP_PCA, Depth == "Deep")
shallow = subset(TP_PCA, Depth == "Shallow")
shapiro.test(deep$Trophic.Position) 
shapiro.test(shallow$Trophic.Position)
t_test(data = deep, Trophic.Position ~ H_S)
t_test(data = shallow, Trophic.Position ~ H_S) 

# Check for depth dependant differences within a compartment to add t.test to the plots
stat.test.TP <- TP_PCA %>%
  group_by(H_S) %>%
  t_test(Trophic.Position ~ Depth) %>%
  add_significance() %>%
  p_round(digits = 3)
stat.test.TP = stat.test.TP%>% add_xy_position(x = "H_S") # this test is not significant and therefore will not be shown on graph

# Plot
mytheme = theme_classic()+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 10), axis.title = element_text(size = 12))

p.TP = ggplot(TP_PCA, aes(y = Trophic.Position, x = H_S)) +
  geom_boxplot(outlier.colour = "NA", fatten = 0.5, alpha = 0.7, lwd = 0.2, aes(fill = Depth))+
  geom_point(position=position_dodge(width=0.75),aes(group= Depth), size = 0.7)+
  labs(y= "Trophic Position")+
mytheme+
  scale_fill_manual(values = c("#D55E00","#0072B2"))+
  scale_x_discrete(labels = c("Host", "Symbiont"))+
stat_pvalue_manual(stat.test.TP, label = "{p.signif} {p}", tip.length = 0.005, y.position = "y.position", hide.ns = TRUE, size = 2.5)+
theme(legend.position = c(0.8, 0.85), legend.title = element_blank())
p.TP

ggsave("TP_OculinaProject2020.pdf", plot = p.TP, width = 8, height = 6,dpi=300, units = "cm")
ggsave("TP_OculinaProject2020.jpeg", plot = p.TP, width = 8, height = 6,dpi=300, units = "cm")
ggsave("TP_OculinaProject2020.svg", plot = p.TP, width = 8, height = 6,dpi=300, units = "cm")




###########  MORPHOLOGY ###########
data<- read.csv("Oculina_Binocular.csv")
View(data)

data$depth = as.factor(data$depth)
data$depth = factor(data$depth, levels = c("Shallow", "Deep"))

mytheme = theme_classic()+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 10), axis.title = element_text(size = 10))


#Distance between polyps
a = ggplot(data, aes(y = polyp_dist_mm, x = depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = depth), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.25)+
  labs(y= "Calyx spacing (mm)")+
  mytheme +
  scale_fill_manual(values = c("#D55E00","#0072B2"))+
  guides(fill = FALSE)
a

#Polyp Diameter
b = ggplot(data, aes(y = calyx_width_mm, x = depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = depth), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.25)+
  mytheme +
  scale_fill_manual(values = c("#D55E00","#0072B2"))+
  labs(y= "Calyx diameter (mm)")+
  scale_y_continuous(breaks = seq(0, 4, by= 1), limits = c(0.5, 3.5))+
  guides(fill = FALSE)
b


# STATS For Binocular data: subset
s = subset(data, depth=="Shallow")
View(s)
d = subset(data, depth=="Deep")
View(d)
 
# Distance between polyps
#normality assumptions 
hist(data$polyp_dist_mm)
hist(s$polyp_dist_mm)
hist(d$polyp_dist_mm)
shapiro.test(s$polyp_dist_mm) 
shapiro.test(d$polyp_dist_mm) # p = 0.0002672 Not normal
leveneTest(polyp_dist_mm~depth,d=data)  
# Kruskal - Wallis non parametric test
kruskal.test(data$polyp_dist_mm~sub$depth)


# Calyx diameter
#normality assumptions 
# Removing one high outlier from shallow to maybe improve normality
sub = subset(sub, calyx_width_mm < 3)
s = subset(sub, depth=="Shallow")
View(s)
d = subset(sub, depth=="Deep")
View(d)

hist(sub$calyx_width_mm)
hist(s$calyx_width_mm)
hist(d$calyx_width_mm)
shapiro.test(s$calyx_width_mm) # p = 0.01526 - STILL NOT NORMAL -
shapiro.test(d$calyx_width_mm) 
leveneTest(calyx_width_mm~depth,d=sub)  
# Kruskal - Wallis non parametric test
kruskal.test(sub$calyx_width_mm~sub$depth)

# Data summary
Sum_all <- sub %>% 
  group_by(depth) %>% 
  summarise_at(c("polyp_dist_mm", "calyx_width_mm"), list(min = min, max = max, mean = mean, sd = sd), na.rm = TRUE)

write.csv(file="Bino.averages.csv", Sum_all)



### SEM
SEM<- read.csv("Oculina_SEM.csv")
SEM$depth = factor(SEM$depth, levels = c("Shallow", "Deep"))

#separate shallow and deep data to look at distributions
s = subset(SEM, depth=="Shallow")
View(s)
d = subset(SEM, depth=="Deep")
View(d)

#normality assumptions 
hist(s$septa_width)
hist(d$septa_width)
shapiro.test(s$septa_width) 
shapiro.test(d$septa_width) 
leveneTest(septa_width~depth,d=SEM)
t.test(SEM$septa_width~SEM$depth) 

# Width of septa in a boxplot
e = ggplot(SEM, aes(y = septa_width, x = depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = depth), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.25)+
  labs(y= "Septa width (µm)")+
  mytheme+
  scale_fill_manual(values = c("#D55E00","#0072B2"))+
  guides(fill = FALSE)
e


#  Summary
Sum_SEM <- SEM %>% 
  group_by(depth) %>% 
  summarise_at(vars(septa_width:centre_diameter), list(min = min, max = max, mean = mean, sd = sd), na.rm = TRUE)

write.csv(file="SEM.averages.csv", Sum_SEM)


#Print figure
plots <- plot_grid(a, b, e, ncol = 1, align = "v")

ggsave("morphology_OculinaProject2020.jpeg", plot = plots, width = 5, height = 14, dpi=600,units = "cm")
ggsave("morphology_OculinaProject2020.pdf", plot = plots, width = 5, height = 14, dpi=600,units = "cm")




#######  Amino Acids Stats  #######
aa = read.csv("Oculina_AA.csv")


# Subset to H_S
aa_H = subset(aa, H_S == 'Host')
View(aa_H)
aa_S = subset(aa, H_S == 'Symbiont')

# Further subset on depth
aa_H_s = subset(aa_H, depth == 'Shallow')
aa_H_d = subset(aa_H, depth == 'Deep')
aa_S_s = subset(aa_S, depth == 'Shallow')
aa_S_d = subset(aa_S, depth == 'Deep')

### HOST ####
#Carbon
## Val
shapiro.test(aa_H_s$Val) 
shapiro.test(aa_H_d$Val) 
leveneTest(Val~depth,d=aa_H)  
t.test(aa_H$Val~aa_H$depth)

## Leu
shapiro.test(aa_H_s$Leu) 
shapiro.test(aa_H_d$Leu) 
leveneTest(Leu~depth,d=aa_H)  
t.test(aa_H$Leu~aa_H$depth)

## Iso
shapiro.test(aa_H_s$Iso) 
shapiro.test(aa_H_d$Iso) 
leveneTest(Iso~depth,d=aa_H)  
t.test(aa_H$Iso~aa_H$depth)

## Met
shapiro.test(aa_H_s$Met) 
shapiro.test(aa_H_d$Met) 
leveneTest(Met~depth,d=aa_H)  
t.test(aa_H$Met~aa_H$depth)

## Phen
shapiro.test(aa_H_s$Phen) 
shapiro.test(aa_H_d$Phen) 
leveneTest(Phen~depth,d=aa_H)  
t.test(aa_H$Phen~aa_H$depth)

#Nitrogen
## Glut_Nitrogen
shapiro.test(aa_H_s$Glut_N) 
shapiro.test(aa_H_d$Glut_N) # not normal
leveneTest(Glut_N~depth,d=aa_H)  
kruskal.test(aa_H$Glut_N~aa_H$depth)

## Phen_Nitrogen
shapiro.test(aa_H_s$Phen_N) 
shapiro.test(aa_H_d$Phen_N) 
leveneTest(Phen_N~depth,d=aa_H)  
t.test(aa_H$Phen_N~aa_H$depth)


### SYMBIONT ####
#Carbon
## Val
shapiro.test(aa_S_s$Val) 
shapiro.test(aa_S_d$Val) 
leveneTest(Val~depth,d=aa_S)  
t.test(aa_S$Val~aa_S$depth)

## Leu
shapiro.test(aa_S_s$Leu) 
shapiro.test(aa_S_d$Leu) 
leveneTest(Leu~depth,d=aa_S)  
t.test(aa_S$Leu~aa_S$depth)

## Iso
shapiro.test(aa_S_s$Iso) 
shapiro.test(aa_S_d$Iso) 
leveneTest(Iso~depth,d=aa_S)  
t.test(aa_S$Iso~aa_S$depth)

## Met
shapiro.test(aa_S_s$Met) 
shapiro.test(aa_S_d$Met) 
leveneTest(Met~depth,d=aa_S)  
t.test(aa_S$Met~aa_S$depth)

## Phen
shapiro.test(aa_S_s$Phen) 
shapiro.test(aa_S_d$Phen) 
leveneTest(Phen~depth,d=aa_S)  
t.test(aa_S$Phen~aa_S$depth)

#Nitrogen
## Glut_Nitrogen
shapiro.test(aa_S_s$Glut_N) 
shapiro.test(aa_S_d$Glut_N) 
leveneTest(Glut_N~depth,d=aa_S)  
t.test(aa_S$Glut_N~aa_S$depth)

## Phen_Nitrogen
shapiro.test(aa_S_s$Phen_N) 
shapiro.test(aa_S_d$Phen_N) # Not normal
leveneTest(Phen_N~depth,d=aa_S)  
kruskal.test(aa_S$Phen_N~aa_S$depth)
