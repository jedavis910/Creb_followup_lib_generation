library(tidyverse)
library(stringr)
library(viridis)
library(cowplot)

figurefont_theme <- theme(text = element_text(size = 8)) +
  theme(axis.title = element_text(size = 8)) +
  theme(legend.title = element_text(size = 9)) +
  theme(axis.text = element_text(size = 8))

#Import random backgrounds

bgs_noRE_gc <- read_csv('bgs_noRE_gc.csv', 
                       col_names = c('sequence', 'gc_content'),
                       col_types  = 'cn')

#Round gc content calc to 2 digits, take 3 of backgrounds of each gc content
#between 0.35 and 0.65

bgs_top3 <- bgs_noRE_gc %>%
  mutate(gc_content = round(gc_content, digits = 2)) %>%
  filter(gc_content >= 0.35 & gc_content <= 0.65) %>%
  group_by(gc_content) %>%
  slice(1:3)

total_bg_number <- bgs_top3 %>%
  group_by(gc_content) %>%
  summarize(backgrounds = n())

#Export as chosen backgrounds

write_csv(bgs_top3, 'bgs_top3.csv', col_names = FALSE)


#Checking renaming in python with background number and gc content

name_gc_content <- read_csv('bgs_top3_name.csv', 
                            col_names = c('name', 'sequence'))

check_name <- left_join(name_gc_content, bgs_top3, by = 'sequence')




#GC content of new 6-site library

sixsite <- read_csv('cre_test_sixsite.csv',
                    col_names = c('name', 'sequence'))

sixsite_split <- sixsite %>%
  separate(name, into = c("subpool", "site1", "site2", "site3", "site4", 
                          "site5", "site6", "fluff", 'background', 'fluff2',
                          'gc'), sep = "_", convert = TRUE) %>%
  select(-fluff, -fluff2) %>%
  mutate(consensus = str_detect(site1, "consensus") + 
           str_detect(site2, "consensus") + 
           str_detect(site3, "consensus") + 
           str_detect(site4, "consensus") + 
           str_detect(site5, "consensus") + 
           str_detect(site6, "consensus"))

test <- sixsite_split %>%
  filter(consensus == 1) %>%
  arrange(gc)

gc_distr <- ggplot(sixsite_split, aes(gc)) +
  facet_grid(consensus ~ ., scales = 'free') +
  geom_histogram(binwidth = 0.01) +
  scale_x_continuous('GC content', 
                     breaks = seq(from = 0.35, to = 0.65, by = 0.05)) +
  scale_y_continuous('# Sequences') +
  figurefont_theme +
  panel_border(colour = 'black') +
  theme(strip.background = element_rect(colour="black", fill="white"))

ggsave('sixsite_gccontent_distr.pdf', gc_distr, units = 'in', width = 3.5,
       height = 4.5)







