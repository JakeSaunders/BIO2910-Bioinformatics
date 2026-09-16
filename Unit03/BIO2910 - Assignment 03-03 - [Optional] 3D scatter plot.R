# Install the package if not already installed
# install.packages("plotly")

# Load the library
library(plotly)

finches <- finches %>% 
  mutate( Sex = as.factor(Sex))


finches %>% 
  plot_ly(
    x = ~Weight,
    y = ~Tarsus,
    z = ~Wing,
    color = ~Sex,
    colors = c('#BF382A', '#0C4B8E','darkgreen'),
    type = "scatter3d",
    mode = "markers"
  ) %>% 
  layout(
    scene = list(
      xaxis = list(title = 'Weight (g)'),
      yaxis = list(title = 'Tarus Length (mm)'),
      zaxis = list(title = 'Wing Length (mm)')
    )
  )
