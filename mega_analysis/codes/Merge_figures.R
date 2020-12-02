#merging figures
library(multipanelfigure)

fig1 <- multi_panel_figure(width = c(25, 25),
                      height = c(25,25),
                      figure_name = "MDD",
                      # # row_spacing = c(0,2),
                      # column_spacing = c(2,0,3,4),
                      unit = "cm") %>% 
  fill_panel("./results/figures/CrossDisease_scatterplotMat.pdf",column = 1, row = 1) %>% 
  fill_panel("./results/figures/RRHO_allComparisons.pdf",column = 2,row = 1) %>% 
  fill_panel("./results/figures/figure2c.pdf", column = 1, row = 2) %>% 
  fill_panel("./results/figures/figure2d.pdf", column = 2,row = 2)

save_multi_panel_figure(figure =fig1,filename = "./results/figures/fig2.pdf",dpi = 300)
