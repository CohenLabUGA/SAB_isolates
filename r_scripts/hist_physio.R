#### plot historical physiology data ####
library(tidyverse)
library(patchwork)
library(gridExtra)
library(rstatix)
library(ggpubr)

#### growth rate ####
setwd("/Users/lquirk/Library/CloudStorage/OneDrive-UniversityofGeorgia/Thesis_g/raw.data")
path1 = "/Users/lquirk/Library/CloudStorage/OneDrive-UniversityofGeorgia/Thesis_g/Figures/Physio"
path2= "/Users/lquirk/Library/CloudStorage/OneDrive-UniversityofGeorgia/Thesis_g/output/Physio/"

##### read in data #####
hist.mu <- read.csv('./Growth-Rate/hist.growth.csv')
hist.mu$Organism <- str_replace(hist.mu$Organism, 'E. huxleyi','G. huxleyi')
hist.av.mu <- hist.mu %>% group_by(Organism,Treatment,Taxa)%>%
  summarize('sd'=sd(GrowthRate),'GrowthRate'=mean(GrowthRate))
colnames(hist.av.mu)[1:3] <- c('culture', 'treatment','taxa')
hist.av.mu <- hist.av.mu %>% group_by(treatment, culture, taxa)

##### test growth and plot #####

remove.outliers <- function(data,var){
  
  quartiles <- quantile(data[[var]], probs=c(.25, .75), na.rm = FALSE)
  IQR <- IQR( data[[var]])
  
  Lower <- quartiles[1] - 1.5*IQR
  Upper <- quartiles[2] + 1.5*IQR 
  
  data_no_outlier <- subset(data,  data[[var]] > Lower &  data[[var]] < Upper)
  data_no_outlier
}


high.vs.low <- function(mu, organism){
  df <- filter(mu,Organism==eval(quote(organism)))
  df$Treatment <- factor(df$Treatment)
  high.fe <- filter(df,Treatment=="High Fe")
  print(str(high.fe))
  low.fe <- filter(df,Treatment=="Low Fe")
  #remove outliers
  high.fe <- remove.outliers(high.fe, 'GrowthRate')
  low.fe <- remove.outliers(low.fe, 'GrowthRate')
  if (nrow(low.fe) >= 3){ #run normality check if there are enough rows
    normal.h <-shapiro_test(high.fe$GrowthRate)
    normal.l <- shapiro_test(low.fe$GrowthRate)
    
    if (normal.h$p.value < 0.05 || normal.l$p.value < 0.05){ #if not normal,...
      print("Sample not normal, run a wilcoxon rank test")
      stat.test <- wilcox_test(df,GrowthRate~Treatment) %>% 
        adjust_pvalue(method='fdr') %>%
        add_significance('p.adj')
      stat.test
    }else{
      print("Sample is normal, run an f test")
      f.test <- var.test(GrowthRate~Treatment,df)  # to test equal variance between treatments, we run an F test
      print(f.test)
      
      if(f.test$p.value >0.05){
        # the H0: sample variances are equal, so if p-value > 0.05, run a t test with
        # equal variances. 
        print("run a t.test with equal variance")
        stat.test <- t_test(df, GrowthRate~Treatment) %>% 
          adjust_pvalue(method='fdr') %>%
          add_significance('p.adj')
        stat.test
      }else{
        #if p-value < 0.05, run a t test with unequal variances 
        print("run a t test with unequal variance")  
        stat.test <- t_test(df, GrowthRate~Treatment, var.equal = FALSE) %>% 
          adjust_pvalue(method='fdr') %>%
          add_significance('p.adj')
        print(stat.test)
      }
    }
  }else{
    print("Sample not long enough, run a single sample t test")
    stat.test <- high.fe %>% t_test(x~1, mu=0.01743318, alternative = "two.sided") %>% 
      adjust_pvalue(method='fdr') %>%
      add_significance('p.adj')
  }
  stat.test$Organism <- organism
  stat.test 
}


hist04 <- high.vs.low(hist.mu,'C. closterium 4')
hist08 <- high.vs.low(hist.mu,'C. closterium 8')
hist06 <- high.vs.low(hist.mu,'G. oceanica')
hist13 <- high.vs.low(hist.mu,'G. huxleyi')

hist.mu.sig <- data.frame(
  'culture'=c('C. closterium 4','C. closterium 8','G. oceanica','G. huxleyi'),
  'sig'=c('**','***','***','NS'),'treatment'='High Fe')
hist.av.mu <- full_join(hist.av.mu,hist.mu.sig)

hist.av.mu$culture <- factor(hist.av.mu$culture, levels=c('G. oceanica','G. huxleyi',
                                                          'C. closterium 8','C. closterium 4'))
hist.mu.plot <- ggplot(hist.av.mu, aes(culture, GrowthRate, fill=treatment))+
  geom_bar(size=2, color='black',stat='identity', position = position_dodge2(0.8))+
  geom_errorbar(aes(ymin=GrowthRate-sd, ymax=GrowthRate+sd),
                position = position_dodge(width=0.9),
                width=0.2, size=2) +
  geom_text(aes(label=sig, y=GrowthRate+sd, 
                vjust=-1.6),
            size=c(40/.pt),
            show.legend=FALSE) +
  facet_grid(~taxa, scales="free_x", space="free_x", switch="x")+
  labs(y='Growth Rate (day^1)')+
  theme_pubr()+
  theme(text = element_text(size=40),
        axis.title.x=element_blank(), 
        axis.text.x = element_text(face="italic"),
        legend.position = 'none',
        strip.placement = 'outside',
        strip.text = element_text(size=40, face='bold'),
        strip.background = element_rect(fill='white', color='white'))  +
  scale_fill_manual(values=c("black","white")) +
  scale_x_discrete(labels=c(
    "C. closterium 8" = "C. closterium UGA8\n(Inner)",
    "G. oceanica" = "G. oceanica\n(Inner)",
    "C. closterium 4" = "C. closterium UGA4\n(Outer)",
    "G. huxleyi" = "G. huxleyi\n(Outer)")) +
  scale_y_continuous(expand=expansion(mul=c(0,0.2)))

hist.mu.plot

#calculate mu/muMax for historical growth rates
hist.mumax <- {
  ag <- hist.mu %>% group_by(Organism, Treatment) %>%
    summarise(m=mean(GrowthRate),sd=sd(GrowthRate))
  highm <- filter(ag, Treatment=="High Fe") 
  lowm <- filter(ag, Treatment=="Low Fe" ) 
  hist.muMax <- data.frame('µ/µMAX'=signif((lowm$m/highm$m),2))
  hist.muMax <- bind_cols(Organism =lowm$Organism, hist.muMax)
  hist.muMax <- hist.muMax[c(4,3,2,1),]
  print(hist.muMax)
  write.csv(hist.muMax, file="/Users/lquirk/Library/CloudStorage/OneDrive-UniversityofGeorgia/Thesis_g/output/Physio/hist_mu.muMax.csv")
  hist.muMax
}

hist_muMax_df <- mutate(hist.muMax, taxa=c('Coccolithophore','Coccolithophore','Diatom','Diatom'))

tbs <- lapply(split(hist_muMax_df, hist_muMax_df$taxa), "[",-3)
df <- tibble(x = rep(-Inf, length(tbs)), 
             y = rep(Inf, length(tbs)), 
             taxa = levels(as.factor(hist_muMax_df$taxa)), 
             tbl = tbs)


hist_muMax_plot <-  hist.mu.plot + 
  theme(legend.position='none', 
        axis.title.x=element_blank())  +
  geom_table(data = df, aes(x = x, y = y, label = tbl),
             hjust=-.3,vjust = 1.5, 
             table.theme = ttheme_gtlight(base_size = 40, padding = unit(c(1,2),'char')),
             fill='white') 
hist_muMax_plot


#### Fv.Fm ####

##### read in data #####
histFire <- read.csv("../rawData/histData/histFire.csv")


###### remove outliers ######
clean.df <- function(df, org, culture){
  high <- filter(df,File == paste(org,'pFe19', sep='')) %>% 
    bind_cols(Treatment="High Fe", Organism=culture)
  high <- remove.outliers(high,'Fv.Fm')
  low <- filter(df,File == paste(org,'pFe21.9', sep='')) %>% 
    bind_cols(Treatment="Low Fe", Organism=culture)
  low <- remove.outliers(low, 'Fv.Fm')
  df <- bind_rows(high,low) %>% 
    select(c('Organism', 'Treatment', 'Fv.Fm','File')) %>%
    group_by(Treatment) 
  
}
## remove treatments not used in experiment
histFire = histFire[(str_detect(histFire$File, '[[:digit:]]{2}.*-N')==F),]

hist.fv.fm <- ungroup(hist.fv.fm)

##### summarize fv.fm #####
hist.avg.fv.fm <- hist.fv.fm %>% group_by(Treatment, Organism) %>%
  summarize('h.fv.fm'=mean(Fv.Fm), 'sd.fv.fm'=sd(Fv.Fm)) %>%
  mutate(taxa= case_when(Organism %in% c('G. oceanica','G. huxleyi')~'Coccolithophore',
                         Organism %in% c('C. closterium 4','C. closterium 8') ~'Diatom'))

#####test fv.fm and plot#####
hist.fv.fm.t <- function(organism){
  df <- filter(hist.fv.fm,Organism==eval(quote(organism)))
  high.fe <- filter(df,Treatment=="High Fe")
  low.fe <- filter(df,Treatment=="Low Fe")
   #run normality check if there are enough rows
    normal.h <-shapiro_test(high.fe$Fv.Fm)
    normal.l <- shapiro_test(low.fe$Fv.Fm)
    if (normal.h$p.value < 0.05 || normal.l$p.value < 0.05){ #if not normal,...
      print("Sample not normal, run a wilcoxon rank test")
      stat.test <- wilcox_test(df,Fv.Fm~Treatment) %>%
        adjust_pvalue(method='fdr') %>%
        add_significance('p.adj')
      stat.test
    }else{
      print("Sample is normal, run an f test")
      f.test <- var.test(Fv.Fm~Treatment,df)  # to test equal variance between treatments, we run an F test
      print(f.test)
      if(f.test$p.value >0.05){
        # the H0: sample variances are equal, so if p-value > 0.05, run a t test with
        # equal variances. 
        print("run a t.test with equal variance")
        stat.test <- t_test(df, Fv.Fm~Treatment) %>% 
          adjust_pvalue(method='fdr') %>%
          add_significance('p.adj')
        stat.test
      }else{
        #if p-value < 0.05, run a t test with unequal variances 
        print("run a t test with unequal variance")  
        stat.test <- t_test(df, Fv.Fm~Treatment, var.equal = FALSE) %>% 
          adjust_pvalue(method='fdr') %>%
          add_significance('p.adj')
        stat.test
      }}
  stat.test$Organism <- organism
  stat.test 
}


h.fv.fm.4.sig <- hist.fv.fm.t('C. closterium 4')
h.fv.fm.8.sig <-hist.fv.fm.t('C. closterium 8')
h.fv.fm.6.sig <-hist.fv.fm.t('G. oceanica')
h.fv.fm.13.sig <-t_test(ungroup(fv.fm13.h),Fv.Fm~Treatment, var.equal=F) %>%
  adjust_pvalue(method='fdr') %>%
  add_significance('p.adj')

hist.fv.fm.sig <- data.frame(
  'Organism'=c('C. closterium 4','C. closterium 8','G. oceanica','G. huxleyi'),
  'sig'=c('NS','NS','**','NS'),'Treatment'='High Fe')
hist.avg.fv.fm <- full_join(hist.avg.fv.fm,hist.fv.fm.sig)

hist.avg.fv.fm$Organism <- factor(hist.avg.fv.fm$Organism, levels=c('G. oceanica','G. huxleyi',
                                                          'C. closterium 8','C. closterium 4'))
hist.fv.fm.plot <- ggplot(hist.avg.fv.fm, aes(Organism, h.fv.fm, fill=Treatment))+
  geom_bar(size=2, color='black',stat='identity', position = position_dodge2(0.8))+
  geom_errorbar(aes(ymin=h.fv.fm-sd.fv.fm, ymax=h.fv.fm+sd.fv.fm),
                position = position_dodge(width=0.9),
                width=0.2, size=2) +
  geom_text(aes(label=sig, y=h.fv.fm+sd.fv.fm, 
                vjust=-1.6),
            size=c(40/.pt),
            show.legend=FALSE) +
  facet_grid(~taxa, scales="free_x", space="free_x", switch="x")+
  labs(y='Fv.Fm')+
  theme_pubr()+
  theme(text = element_text(size=40),
        axis.title.x=element_blank(), 
        axis.text.x = element_text(face="italic"),
        legend.position = 'none',
        strip.placement = 'outside',
        strip.text = element_text(size=40, face='bold'),
        strip.background = element_rect(fill='white', color='white'))  +
  scale_fill_manual(values=c("black","white")) +
  scale_x_discrete(labels=c(
    "C. closterium 8" = "C. closterium UGA8\n(Inner)",
    "G. oceanica" = "G. oceanica\n(Inner)",
    "C. closterium 4" = "C. closterium UGA4\n(Outer)",
    "G. huxleyi" = "G. huxleyi\n(Outer)")) +
  scale_y_continuous(expand=expansion(mul=c(0,0.2)))

hist.fv.fm.plot

#### Tqa####

##### remove outliers #####
clean.tqa <- function(df, org, culture){
  high <- filter(df,File == paste(org,'pFe19', sep='')) %>% 
    bind_cols(Treatment="High Fe", Organism=culture)
  high <- remove.outliers(high,'TauAv1')
  low <- filter(df,File == paste(org,'pFe21.9', sep='')) %>% 
    bind_cols(Treatment="Low Fe", Organism=culture)
  low <- remove.outliers(low, 'TauAv1')
  df <- bind_rows(high,low) %>% 
    select(c('Organism', 'Treatment', 'TauAv1','File')) %>%
    group_by(Treatment) 
  
}

tq4.h <- clean.tqa(fv.fm.04, '04', 'C. closterium 4')
tq8.h <- clean.tqa(fv.fm.08, '08', 'C. closterium 8')
tq6.h <- clean.tqa(fv.fm.06, '06', 'G. oceanica')
tq13.h <- clean.tqa(fv.fm.13, '13', 'G. huxleyi')

hist.tq <- bind_rows(tq4.h, tq8.h, tq6.h, tq13.h)
hist.tq <- ungroup(hist.tq)

#### summarize tqa #####

hist.avg.tq <- hist.tq %>% group_by(Organism, Treatment) %>%
  summarize('h.tauAv1'=mean(TauAv1), 'sd.tauAv1'=sd(TauAv1)) %>% mutate(taxa='Diatom') %>% 
  mutate(taxa= case_when(Organism %in% c('G. oceanica','G. huxleyi')~'Coccolithophore',
                         Organism %in% c('C. closterium 4','C. closterium 8') ~'Diatom'))


##### test and plot #####
hist.tq.t <- function(organism){
  df <- filter(hist.tq,Organism==eval(quote(organism)))
  high.fe <- filter(df,Treatment=="High Fe")
  low.fe <- filter(df,Treatment=="Low Fe")
  #run normality check if there are enough rows
  normal.h <-shapiro_test(high.fe$TauAv1)
  normal.l <- shapiro_test(low.fe$TauAv1)
  if (normal.h$p.value < 0.05 || normal.l$p.value < 0.05){ #if not normal,...
    print("Sample not normal, run a wilcoxon rank test")
    stat.test <- wilcox_test(df,TauAv1~Treatment) %>%
      adjust_pvalue(method='fdr') %>%
      add_significance('p.adj')
    stat.test
  }else{
    print("Sample is normal, run an f test")
    f.test <- var.test(TauAv1~Treatment,df)  # to test equal variance between treatments, we run an F test
    print(f.test)
    if(f.test$p.value >0.05){
      # the H0: sample variances are equal, so if p-value > 0.05, run a t test with
      # equal variances. 
      print("run a t.test with equal variance")
      stat.test <- t_test(df, TauAv1~Treatment) %>% 
        adjust_pvalue(method='fdr') %>%
        add_significance('p.adj')
      stat.test
    }else{
      #if p-value < 0.05, run a t test with unequal variances 
      print("run a t test with unequal variance")  
      stat.test <- t_test(df, TauAv1~Treatment, var.equal = FALSE) %>% 
        adjust_pvalue(method='fdr') %>%
        add_significance('p.adj')
      stat.test
    }}
  stat.test$Organism <- organism
  stat.test 
}

h.tq.4.sig <- hist.tq.t('C. closterium 4')
h.tq.8.sig <-hist.tq.t('C. closterium 8')
h.tq.6.sig <-hist.tq.t('G. oceanica')
h.tq.13.sig <-t_test(ungroup(tq13.h),TauAv1~Treatment, var.equal=F) %>%
  adjust_pvalue(method='fdr') %>%
  add_significance('p.adj')

hist.tq.sig <- data.frame(
  'Organism'=c('C. closterium 4','C. closterium 8','G. oceanica','G. huxleyi'),
  'sig.tq'=c('NS','*','NS','NS'),'Treatment'='Low Fe')
hist.avg.tq <- full_join(hist.avg.tq,hist.tq.sig)
hist.avg.tq$Organism <- factor(hist.avg.tq$Organism, levels=c('G. oceanica','G. huxleyi',
                                                                    'C. closterium 8','C. closterium 4'))
hist.tq.plot <- ggplot(hist.avg.tq, aes(Organism, h.tauAv1, fill=Treatment))+
  geom_bar(size=2, color='black',stat='identity', position = position_dodge2(0.8))+
  geom_errorbar(aes(ymin=h.tauAv1-sd.tauAv1, ymax=h.tauAv1+sd.tauAv1),
                position = position_dodge(width=0.9),
                width=0.2, size=2) +
  geom_text(aes(label=sig.tq, y=h.tauAv1+sd.tauAv1, 
                vjust=-1.6),
            size=c(40/.pt),
            show.legend=FALSE) +
  facet_grid(~taxa, scales="free_x", space="free_x", switch="x")+
  labs(y='tQa')+
  theme_pubr()+
  theme(text = element_text(size=40),
        axis.title.x=element_blank(), 
        axis.text.x = element_text(face="italic"),
        legend.position = 'right',
        strip.placement = 'outside',
        strip.text = element_text(size=40, face='bold'),
        strip.background = element_rect(fill='white', color='white'))  +
  scale_fill_manual(values=c("black","white")) +
  scale_x_discrete(labels=c(
    "C. closterium 8" = "C. closterium UGA8\n(Inner)",
    "G. oceanica" = "G. oceanica\n(Inner)",
    "C. closterium 4" = "C. closterium UGA4\n(Outer)",
    "G. huxleyi" = "G. huxleyi\n(Outer)")) +
  scale_y_continuous(expand=expansion(mul=c(0,0.2)))

hist.tq.plot

hist.physio.plot <-hist_muMax_plot+hist.fv.fm.plot+hist.tq.plot +
  plot_layout(ncol = 2, nrow=2) &
  theme(legend.direction = 'vertical',legend.text = element_text(size=30),
        legend.title = element_text(size=50), legend.key.size =unit(2,'cm'))

hist.physio.plot <- hist.physio.plot + 
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size=50, face='bold'))


hist.physio.plot

