# ------------------------------ Plot2A ---------------------------------------
# Plotting the alignment

x = read.delim("../DiffBind/test.txt")
colnames(x) = c("Name","Untreat Replicate 1","Untreat Replicate 2","Treat Replicate 1","Treat Replicate 2")
x = x %>% pivot_longer(!Name,names_to = "Condition",values_to = "read")
y=c("Final","Filter","Duplicate","MT","Unmapped")
ggplot(x,aes(x=read,y=Condition,fill=forcats::fct_rev(Name))) + geom_bar(stat = "identity",position="stack") +
  scale_x_continuous(labels = unit_format(unit = "",scale = 1e-6),expand = c(0,0),limits = c(0,1.5e8)) +
  labs(x = "Read Pairs (Million)",y="Sample") +
  theme_bw()
  scale_fill_manual(values=c("red4","forestgreen","lightblue","deeppink4","blue4"))
  



# ------------------------------ Plot2B ---------------------------------------
# Plotting