## ------------------------------------------------------------------------
library(LC50)
head(toxicity)

## ----fig.width=7,fig.height=6--------------------------------------------
library(ggplot2)
ggplot(toxicity,aes(x=conc,y=alive/total,colour=group))+
  geom_point()+
  facet_grid(temperature~salinity)

