library(readxl)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

ICR96 <- read_excel("~/Dropbox/Master_UOC/TFM/ICR96/ICR96.xlsx", sheet = "CNVs")
ICR96$type <- paste(ICR96$ExonCNVType, ICR96$ExonCNVSize)
ICR96_count <- ICR96 %>% group_by(Gene, type) %>% dplyr::summarize(n=n())
ICR96_count2 <- ICR96 %>% group_by(Gene, ExonCNVType, ExonCNVSize) %>% dplyr::summarize(n=n())

display.brewer.pal(n = 8, name = 'Spectral')
colors <- brewer.pal(n = 8, name = "Spectral")[c(2,3,7:8)]

tiff("~/Dropbox/Master_UOC/TFM/TFM_code/analysis/evaluation/graphs/cnvs_ICR96.tiff", units="in", width=10, height=6, res=300)

ggplot(data=ICR96_count, aes(x=Gene, y = n, fill = type, color = type)) +
  geom_bar(stat = "identity", alpha = 0.2) + 
  theme_classic() +
  ylab("Number of samples") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
  theme(legend.title = element_blank(), text = element_text(size=18)) +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors)

dev.off()
