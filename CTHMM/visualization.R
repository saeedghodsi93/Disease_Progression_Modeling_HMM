
# visualization
visualization.func <- function(H, a, tau.true, tau.obs, z.true, z.obs, y.obs, u.obs) {
  
  fig1 <- plot_ly(x = tau.true[1,], y = z.true[1,], type = 'scatter', mode = 'lines+markers') 
  fig2 <- plot_ly(x = tau.obs[1,], y = y.obs[1,]+1, type = 'scatter', mode = 'lines+markers') 
  fig3 <- plot_ly(x = tau.obs[1,], y = u.obs[1,]+1, type = 'scatter', mode = 'lines+markers') 
  fig <- subplot(fig1, fig2, fig3, nrows = 3) %>% 
    layout(title = list(text = "Stacked Subplots"),
           plot_bgcolor='#e5ecf6', 
           xaxis = list( 
             zerolinecolor = '#ffff', 
             zerolinewidth = 2, 
             gridcolor = 'ffff'), 
           yaxis = list( 
             zerolinecolor = '#ffff', 
             zerolinewidth = 2, 
             gridcolor = 'ffff')) 
  fig
  
}

