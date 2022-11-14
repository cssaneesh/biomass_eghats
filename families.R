library(ggplot2)

data <- data.frame(
  category=c("Poaceae",
             'Fabaceae', 
             "Asteraceae", 
             "Acanthaceae", 
             'Rubiaceae', 
             'Euphorbiaceae', 
             'Others'),
  count=c(18,15,8,4,4,3,17 )
)

data$fraction = data$count / sum(data$count)
data$ymax = cumsum(data$fraction)
data$ymin = c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$category, "\n Species: ", data$count)

families.pie <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect(col='darkgrey', linetype= 'dotted')+
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")+
  geom_text( x=2, aes(y=labelPosition, label=label),
             #color='black',
             color='brown',
             #color= c('#2b8cbe', '#7bccc4', '#ccebc5','#f0f9e8','#08589e', '#a8ddb5', '#4eb3d3'), 
             size=4)
families.pie

ggsave('families.jpg',
       width = 10,
       height = 6,
       dpi = 300)

 


library(plotly)
data <- data.frame(
  category=c("Poaceae",
             'Fabaceae', 
             "Asteraceae", 
             "Acanthaceae", 
             'Rubiaceae', 
             'Euphorbiaceae', 
             'Others'),
  count=c(18,15,8,4,4,3,17 )
)

fig <- data %>% plot_ly(labels= ~ category, values= ~ count)
fig <- fig %>% add_pie(hole= 0.4)
fig <- fig %>% layout(title = "% of families",  showlegend = T)

fig
