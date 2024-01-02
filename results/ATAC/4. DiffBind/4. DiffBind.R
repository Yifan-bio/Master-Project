# Plot 1A

x = Diff.bind
x$condition = "Unchanged"
x$condition[x$Fold > 1 & x$FDR < 0.05] = "Opened"
x$condition[x$Fold < -1 & x$FDR < 0.05] = "closed"

ggplot(data = x,aes(x=Fold,y=-log10(FDR),color = condition)) +
  geom_point() + 
  scale_color_manual(values=c("red4","forestgreen",'grey')) + 
  theme(panel.grid = element_blank()) +
  labs(title = "Plot 1A) Change in accessibility in peak regions across the samples",
       x = 'Fold change',
       y = '-log10(FDR)',
       colour = 'Accessibility change')
ggsave(dpi = 900,filename = "Volcano Plot",height = 6,width = 10)


#Visualisation
plot(dbaObject)
dba.plotPCA(dbaObject)





boxplot(C1,C2,T1,T2,
        main = "Distribution of accessible peak",
        at = c(1,2,4,5),
        names = c("Control 1", "Control 2", "Treated 1", "Treated 2"),
        ylab = 'Sample',
        xlab = 'Size of peak(bp)',
        col = c("orange","red"),
        border = "brown",
        horizontal = T,
        notch = TRUE
)

test = function(input) {
  x = dbaObject[["peaks"]][[input]]
  x$num = x$End - x$Start
  x$num[x$num > 1500] = 1500
  x = x$num
  #return(x$num)
}

# Plot 1B) DESeq2 vocano plot
x = Diff.bind %>% select(Fold,FDR) 
x$con = "non-significant"
x$con[x$Fold > 1 & x$FDR < 0.05] = "upregulated"
x$con[x$Fold < -1 & x$FDR < 0.05] = "downregulated"
ggplot(x,aes(x=Fold,y=-log10(FDR),color = con)) + geom_point() +
  scale_colour_manual(values = c("green","grey","red")) +                       # Change the colour of dots
  theme_bw() +                                                                  # make it bw rather then grey background
  labs(x = "Fold change in chromatin accessibility",y = "-log10(FDR) value",color = "Gene status") + #naming
  scale_y_continuous(expand = c(0,0)) +                                         # Remove the gap between the results and x margin
  geom_vline(xintercept = 1) + geom_hline(yintercept = -log10(0.05)) + geom_vline(xintercept = -1)
scale_fill