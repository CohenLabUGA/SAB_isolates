

#### check list of required packages is installed, if not installs and loads them
### function modified from following stackoverflow:
# https://stackoverflow.com/questions/15155814/check-if-r-package-is-installed-then-load-library

checkAndLoadPackages <- function(...,silent=FALSE){
  
  #check names and run 'require' function over if the given package is installed
  requirePkg<- function(pkg){if(length(setdiff(pkg,rownames(installed.packages())))==0)
    suppressPackageStartupMessages(require(pkg, quietly = TRUE,character.only = TRUE))
  }
  
  packages <- as.vector(unlist(list(...)))
  if(!is.character(packages))stop("No numeric allowed! Input must contain package names to install and load")
  
  if (length(setdiff(packages,rownames(installed.packages()))) > 0 )
    install.packages(setdiff(packages,rownames(installed.packages())),
                     repos="https://cloud.r-project.org")
  
  res<- unlist(sapply(packages, requirePkg))
  
  if(silent == FALSE && !is.null(res)) {cat("\nBellow Packages Successfully Loaded:\n\n")
    print(res)
  }
}

anova.fire <- function(organism, var){
  aov.fire <- fire %>% filter(Organism==organism) 
  aov.fire <- data.frame(Treatment=aov.fire$Treatment, 
                         Param = aov.fire[[var]], 
                         Organism= aov.fire$Organism)
  print(aov.fire)
  #1. run anova using aov
  anova.org <- aov(Param~as.factor(Treatment), data=aov.fire)
  #print(summary(anova.org))
  
  #2. create summary table of organism for plotting later
  summary.org <- aov.fire %>% 
    group_by(Treatment) %>%
    summarize('av.Param' = mean(Param),'sd'=sd(Param))%>%
    arrange(desc('av.Param')) ##!!! this is important so letters from cld match right later
  
  #3. run tukey's test using anova output
  tukey.org <- TukeyHSD(anova.org)
  #print(tukey.org)
  
  #4. compact letters asociated with tukey's output
  cld.org <- multcompLetters4(anova.org,tukey.org)
  #print(cld.org)
  cld.org <- as.data.frame.list(cld.org$`as.factor(Treatment)`)
  
  #5. add the results to summary table in a new column
  # this will make row binding later easier
  summary.org$Tukey <- cld.org$Letters
  summary.org <- mutate(summary.org, Organism=organism)
  summary.org
}


barplot.aov <- function(df.sig){
  df.sig <- df.sig %>% group_by(Organism, Treatment)
  plot = ggplot(df.sig, aes(Organism ,av.Param, fill=Treatment)) +
    geom_bar(color="black",stat='identity', position='dodge2') + 
    geom_errorbar(aes(ymin=av.Param-sd, ymax=av.Param+sd),width=0.2,
                  position = position_dodge(width=0.9))+
    geom_text(aes(label=Tukey, y=av.Param+sd,vjust=-1.8,group=Treatment),
              show.legend=FALSE,position=position_dodge(width=0.9)) +
    geom_shadowtext(aes(x=Organism, label='-----', y=h.AvgFv.Fm,group=Treatment),
                    bg.r=.08, bg.colour='white',color='black',
                    position=position_dodge2(0.9)) +
    facet_grid(~taxa, scales="free_x", space="free_x", switch="x") +
    theme_pubr()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(face="italic"),legend.position='none',
          strip.placement = 'outside',
          strip.text = element_text(size=15,face='bold'),
          strip.background = element_rect(fill='white', color='white')) + 
    scale_x_discrete(labels=c(
      "C. closterium UGA8" = "C. closterium UGA8\n(Inner)",
      "G. oceanica" = "G. oceanica\n(Inner)",
      "C. closterium UGA4" = "C. closterium UGA4\n(Outer)",
      "G. huxleyi" = "G. huxleyi\n(Outer)")) +
    scale_fill_manual(labels=c('High Iron','Low Iron','Iron Amendment'),
                      values=c("black","white","grey"))+
    labs(y="Fv/Fm", fill='Treatment') + 
    scale_y_continuous(expand=expansion(mul=c(0,0.2)))
  plot
}