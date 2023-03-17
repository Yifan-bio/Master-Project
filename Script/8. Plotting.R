# Script for all plotting involved in section 8 Methylation analysis
# Last modified: 17 March 2023

library(tidyverse)
library(gridExtra)

#### Plot8A Scatter Plot compare methylation and accessibility changes ####
ggplot(myDiff %>% arrange(desc(con)),aes(y = meth.diff,x=Acc.Fold)) + 
  geom_point(aes(colour = con)) +
  geom_smooth(method = lm,color = "black") + 
  theme_bw() + 
  scale_x_continuous(limits = c(-6,6),expand = c(0,0)) +
  scale_y_continuous(limits = c(-100,100),expand = c(0,0)) +
  labs(x = expression("Accessibility change (Log"[2]*" Fold Change)"),y = "Methylation changes (Percentage)") + 
  scale_color_manual(values = c("#B9181A","#C0C0C0")) +
  theme(text = element_text(size = 20),legend.position = "none",rect = element_rect(fill = "transparent"))
ggsave("test.png",width = 9.5,height = 7.9,dpi = 1200)