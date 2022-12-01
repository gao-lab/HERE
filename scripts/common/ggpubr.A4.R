library(ggpubr)

ggsave.A4 <- function(filename, plot, width.r, height.r){
    ggsave(filename=filename, plot=plot, device="png", width=width.r*21, height=height.r*29.7, units="cm", dpi=300)
}
