# boxplot for peak size distribution
x = read_delim('./Input/peak/0hr_rep1_hmm_peaks.gappedPeak',col_names = F)
x$X5 = x$X3 - x$X2
x$X5 = x$X3 - x$X2
x = x%>%filter(X4 != "HighCoveragePeak_0")

x = data.frame(name = "Treat_1",width = x$X5)

z = rbind(x,z)
ggplot(z,aes(x = name,y=width)) + geom_boxplot() +
  coord_flip() + theme_bw() +
  labs(x = "Sample",y = "Peak size (bp)")

# Peak distribution between replicates (if same peak is opened in replicates)
x = data.frame(name = c("Unique","Two replicate","Three replicate","All replicate"),number = c(14006,18348,9053,21053))
ggplot(x,aes(x = name,y=number)) + geom_bar(stat = "identity") +
  labs(x="Number of replicates in all conditions containing the peak",y="Number of peaks") +
  theme_bw() +
  scale_y_continuous(expand = c(0,500))

