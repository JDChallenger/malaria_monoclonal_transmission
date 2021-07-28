#Join the plots together
library(gridExtra)
library(grid)
library(ggpubr)

#Show all 6 panels, as per Figure 1. Note: you have to run the 3 regression scripts first
grid.arrange(pl1,pl2,pl3,pl4,pl5,pl6,nrow=3)

#Alternative version:

dm <- text_grob("DMFA", face = "bold")
sm <- text_grob("SMFA", face = "bold")
blnk <- rectGrob(gp=gpar(fill="white",col="white"))
pfs230 <- text_grob("Pfs230", face = "bold")
pfs25 <- text_grob("Pfs25", face = "bold")
pfs48 <- text_grob("    Pfs48/45", face = "bold")

lay <- rbind(c(1,2,2,2,2,3,3,3,3),
              c(4,5,5,5,5,6,6,6,6),
              c(4,5,5,5,5,6,6,6,6),
              c(4,5,5,5,5,6,6,6,6),
              c(4,5,5,5,5,6,6,6,6),
              c(7,8,8,8,8,9,9,9,9),
              c(7,8,8,8,8,9,9,9,9),
              c(7,8,8,8,8,9,9,9,9),
              c(7,8,8,8,8,9,9,9,9),
              c(10,11,11,11,11,12,12,12,12),
              c(10,11,11,11,11,12,12,12,12),
              c(10,11,11,11,11,12,12,12,12),
              c(10,11,11,11,11,12,12,12,12)
)

g <-arrangeGrob(grobs = list(blnk,sm,dm,pfs25,pl1,pl2,pfs48,pl3,pl4,pfs230,pl5,pl6),
                  layout_matrix = lay)
#Note: correctly outputting the greek symbols can be an issue on some computers
ggsave(file="Figure1.pdf", g,height = 14, width = 13)#,  device = cairo_pdf) #saves g




