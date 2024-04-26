#Making tables for thesis
source('functions.R')

checkAndLoadPackages('tidyverse','flextable','officer')

set_flextable_defaults(font.size=12, font.family = 'Times New Roman')

####---Physiology Table---####
avgGrowth <- read.csv('../output/AvgGrowth.csv')
avgGrowth$Treatment <- str_replace(avgGrowth$Treatment, 'Fe','Iron')
avgChla <- read.csv('../output/Avg.Chla.csv')
avgChla$Treatment <- str_replace(avgChla$Treatment, 'Fe','Iron')
avgFire <- read.csv('../output/AvgFire.csv')

phys_table <- left_join(avgFire,avgGrowth, by=c('Organism','Treatment')) 
phys_table <- left_join(phys_table,avgChla, by=c('Organism','Treatment')) 
phys_table$Organism <- factor(phys_table$Organism, 
                              levels = c('C. closterium UGA8','G. oceanica',
                                         'C. closterium UGA4','G. huxleyi')) 
phys_table <- arrange(phys_table, Organism)


phys_table <- mutate(phys_table, 'Growth'=(paste(signif(mu,2),'±',round(mu.sd,2))),
       'Chla'=(paste(signif(chla,2),'±',round(chla.sd,2))),
       'Fv/Fm'=(paste(signif(AvgFv.Fm,2),'±',round(SdFv.Fm,2))),
       'tQa'=(paste(signif(AvgtQa,2),'±',round(SdtQa,2))),
       'sigma'=(paste(signif(AvgSigma,2),'±',round(SdSigma,2))))

phys_table <- select(phys_table, c(1,2,13:17))
phys_table <- pivot_longer(phys_table,cols=!c('Organism','Treatment'),names_to = 'Parameter',
                           values_to = 'mean')
phys_table$mean <-  str_replace_na(phys_table$mean, '--')
phys_table$mean <- str_replace(phys_table$mean, '0 ± NA','--')
phys_table$Treatment <- str_replace(phys_table$Treatment, 'Fe','Iron')
phys_table <- pivot_wider(phys_table, names_from = Organism, values_from = mean)
phys_table <- group_by(phys_table, Parameter)
phys_table <- phys_table[,c(2,1,3,5,4,6)]
phys_table <- phys_table%>%arrange(match(Parameter,c('Growth','Chla','Fv/Fm','tQa','sigma')))
phys_table$Parameter <- str_replace_all(phys_table$Parameter, c('Chla'='Chl a (µg/L)',
                                                                'Growth'='Growth Rate (day^-1)',
                                                                'sigma'='Sigma'))


ps <- as_grouped_data(phys_table, groups=c('Parameter'))
exp_tab <- as_flextable(ps, hide_grouplabel = T) %>% 
  style(j=1, i=~!is.na(Parameter), pr_t=fp_text_default(bold=TRUE),
        pr_p=officer::fp_par(text.align='left', padding=5)) %>%
  prepend_chunks(i = ~is.na(Parameter), j = 1, as_chunk("\t")) %>%
  set_header_labels(values=list(Treatment='',
                                'C. closterium 8'='C. closterium UGA8',
                                'C. closterium 4'='C. closterium UGA4',
                                'G. oceanica'='G. oceanica',
                                'G. huxleyi'='G. huxleyi')) %>%
  hline(i=c(4,8,12), part='body') %>% 
  add_header_row(values=c('','Inner Shelf','Outer Shelf'), 
                 colwidth=c(1,2,2))%>%
  italic(part='header', i=2) %>%autofit()
exp_tab


hist_table <- full_join(hist.mu.table, h.fv.fm)
hist_table$Treatment <- str_remove(hist_table$Treatment, ' Fe|Fe ')
hist_table <- pivot_longer(hist_table, cols=!c(Organism,Treatment), names_to = 'Parameter', 
                           values_to='xx')
hist_table <- pivot_wider(hist_table, names_from = Organism, values_from = xx)
hist_table <- group_by(hist_table, Parameter)
hist_table <- hist_table[,c(2,1,4,6,3,5)]
hist_table <- hist_table%>%arrange(match(Parameter,c('GrowthRate','Chla','Fv/Fm','tQa')))
hist_table$Parameter <- str_replace_all(hist_table$Parameter, 
                                        c('Chla'='Chl a (µg/L)',
                                          'GrowthRate'='Growth Rate (day^-1)'))

hs <- as_grouped_data(hist_table, groups=c('Parameter'))
hist_tab <- as_flextable(hs, hide_grouplabel = T) %>% 
  style(j=1, i=~!is.na(Parameter), pr_t=fp_text_default(bold=TRUE),
        pr_p=officer::fp_par(text.align='left', padding=5)) %>%
  prepend_chunks(i = ~is.na(Parameter), j = 1, as_chunk("\t")) %>%
  set_header_labels(values=list(Treatment='',
                                'C. closterium 8'='C. closterium UGA8',
                                'C. closterium 4'='C. closterium UGA4',
                                'G. oceanica'='G. oceanica',
                                'G. huxleyi'='G. huxleyi')) %>%
  hline(i=c(3,6), part='body') %>% 
  add_header_row(values=c('','Inner Shelf','Outer Shelf'), 
                 colwidth=c(1,2,2))%>%
  italic(part='header', i=2) %>% autofit()

hist_tab


#### Sequencing Table ####
illumina <- read.csv('./Illumina/illumina_stat.csv')
illumina$Sample.Name <- str_remove(illumina$Sample.Name, 'R.')
illumina$X..GC <- str_remove(illumina$X..GC, '[:punct:]')
illumina$Length <- str_remove(illumina$Length, 'bp')
illumina$X..GC <- as.numeric(illumina$X..GC)
illumina$Length <- as.numeric(illumina$Length)
illumina <- group_by(illumina, Organism)%>%
  summarise(across(starts_with(c('X..GC','Length')),~mean(.x)),
            'M.seqs'=sum(M.Seqs), 'Samples'=n_distinct(Sample.Name))
colnames(illumina) <- c('Organism','%GC','Mean length (bp)','M Sequences','Samples (n)')
illumina <- pivot_longer(illumina, cols=!Organism, names_to = 'Parameter', values_to = 'stat') %>% mutate('type'='Illumina')
  
  
pac <- read.csv('./PacBio/PacBio_stat.csv')[1:4,]
pac <- pac[,c(1,5,6,8)]
colnames(pac) <- c('Organism', 'Hifi reads', 'Mean Hifi length (bp)','Hifi yield (Mbp)')
pac <- pivot_longer(pac, cols=!Organism, names_to = 'Parameter', values_to = 'stat') %>% mutate('type'='PacBio')

seqs <- full_join(illumina,pac)
seqs$stat <- round(seqs$stat, digits=2)
seqs <- pivot_wider(seqs, names_from = Organism, values_from = stat)


seq_tab <- as_grouped_data(seqs, groups=c('type'))
seq_tab <- as_flextable(seq_tab, hide_grouplabel = T) %>% 
  style(j=1, i=~!is.na(type), pr_t=fp_text_default(bold=TRUE),
        pr_p=officer::fp_par(text.align='left', padding=5)) %>%
  prepend_chunks(i = ~is.na(type), j = 1, as_chunk("\t")) %>%
  italic(part='header', j=c(2:5)) %>%
  hline(i=c(5), part='body') %>% 
  set_header_labels(values=list(Parameter=''))%>%
                                  autofit()
seq_tab

save_as_image(seq_tab, path='../output/Assembly/sequences_table.png')

#### DEGs Ko Table ####
de_ko <- read.csv('./deseq/de_ko_stat.csv')

de_ko <- mutate(de_ko, 'new'=paste(de_ko$count, ' (',de_ko$percent,'%)',sep=''))

de_ko <- de_ko %>% select(!c('count', 'percent')) %>%
  pivot_wider(names_from = Organism, values_from = new)

de_ko_tab <- flextable(de_ko) %>%
  italic(j=c(2:5), part='header') %>% bold(j=1, part='body') %>%
  set_header_labels(values=list(Treatment.comparison='')) %>%
  autofit()
de_ko_tab

save_as_image(de_ko_tab, path='../output/Assembly/de_ko.png')

#### DEGs All Table ####
de <- read.csv('./deseq/de_stat.csv')

de <- mutate(de, 'new'=paste(de$count, ' (',de$percent,'%)',sep=''))

de <- de %>% select(!c('count', 'percent')) %>%
  pivot_wider(names_from = Organism, values_from = new)

de_tab <- flextable(de) %>%
  italic(j=c(2:5), part='header') %>% bold(j=1, part='body') %>%
  set_header_labels(values=list(Treatment.comparison='')) %>%
   autofit()
de_tab
save_as_image(de_tab, path='../output/Assembly/de.png')

#### Isolation location Table ####
small_border <- fp_border(color='lightgrey',width=0.5)
big_border <- fp_border(color='grey', width=1)

isolate <- read.csv('../output/Physio/isolation_tab.csv')[1:4,1:7]
iso_tab <- flextable(isolate) %>%
  set_header_labels(values=list(Cruise.date='Cruise datae',
                                Shelf.location='SAB Shelf zone',
                                Lat.Lot='Lon/Lat',
                                Temperature='Temperature (degC)'))%>%
italic(j=1, part='body') %>% bold(part='header') %>%
border_remove() %>%
hline_top(border=big_border,part='all') %>% hline_bottom(border=small_border) %>%
  autofit()

iso_tab <- merge_v(
  iso_tab, j=c('Cruise.date','Lat.Lot', 'Shelf.location','Depth','Temperature')) %>%
 border_inner(part='body', border=small_border) 
 
iso_tab

save_as_image(iso_tab, '../output/Physio/isolation_tab.png')

#### Assembly stats Table ####

assembly <- read.csv('../output/Assembly/assembly_tab_stats.csv')[,1:5]
aa <- assembly %>% select(!'Shelf') %>%
  pivot_wider(names_from = Organism,values_from = value)

assembly_tab <- as_grouped_data(aa, groups=c('Parameter')) 
at=as_flextable(assembly_tab, hide_grouplabel = T) %>% 
  style(j=1, i=~!is.na(Parameter), pr_t=fp_text_default(bold=TRUE),
        pr_p=officer::fp_par(text.align='left', padding=5)) %>%
  prepend_chunks(i = ~is.na(Parameter), j = 1, as_chunk("\t")) %>%
  set_header_labels(values=list(Stat='')) %>%
   align(align='center',part='header')%>%
  add_header_row(values=c('','Inner Shelf','Outer Shelf'), colwidth=c(1,2,2))%>%
  italic(part='header', i=2) %>%
  merge_v(part='header') %>%
   autofit()
at

save_as_image(at,'../output/Assembly/assembly_stats_table.png')

#### Annotation Table ####

annotation <- read.csv('../output/Assembly/annotation_tab_stats.csv')[,1:4]
annotation <- annotation %>%  pivot_wider(names_from = Organism,values_from = value)
anno_tab <- as_grouped_data(annotation, groups=c('Parameter')) 
anno_tab=as_flextable(anno_tab, hide_grouplabel = T) %>% 
  style(j=1, i=~!is.na(Parameter), pr_t=fp_text_default(bold=TRUE),
        pr_p=officer::fp_par(text.align='left', padding=5)) %>%
  prepend_chunks(i = ~is.na(Parameter), j = 1, as_chunk("\t")) %>%
  set_header_labels(values=list(Stat='')) %>%
  align(align='center',part='header')%>%
  add_header_row(values=c('','Inner Shelf','Outer Shelf'), colwidth=c(1,2,2))%>%
  italic(part='header', i=2) %>%
  merge_v(part='header') %>%
  autofit()
anno_tab

save_as_image(anno_tab, '../output/Assembly/annotation_stats_table.png')
