library(ggplot2)
library(dplyr)

mitopathy_Tc_df <- tibble::tribble(
                        ~Sample, ~Pyr_Uri, ~celltype, ~Viabilty, ~Pct_division, ~Pro_Index,
                      "Healthy",     TRUE,    "CD4+",      64.3,          92.7,       2.71,
                      "Healthy",    FALSE,    "CD4+",        54,          93.4,       2.69,
                          "PT1",     TRUE,    "CD4+",      48.3,          27.1,       1.16,
                          "PT1",    FALSE,    "CD4+",      16.7,          13.2,          1,
                          "PT2",     TRUE,    "CD4+",      27.6,          19.9,       1.08,
                          "PT2",    FALSE,    "CD4+",      0.22,             0,          1,
                          "PT3",     TRUE,    "CD4+",      44.6,          15.4,       1.23,
                          "PT3",    FALSE,    "CD4+",        32,          12.2,          1,
                      "Healthy",     TRUE,    "CD8+",      91.1,          94.8,       2.92,
                      "Healthy",    FALSE,    "CD8+",      86.7,          92.1,       2.85,
                          "PT1",     TRUE,    "CD8+",      66.1,          18.4,       1.23,
                          "PT1",    FALSE,    "CD8+",      21.2,          16.2,          1,
                          "PT2",     TRUE,    "CD8+",      42.1,          25.1,        1.1,
                          "PT2",    FALSE,    "CD8+",      0.97,             0,          1,
                          "PT3",     TRUE,    "CD8+",      75.1,          22.4,       1.18,
                          "PT3",    FALSE,    "CD8+",      62.3,            18,       1.16
                      )
library(BuenColors)
                    
t1 <- mitopathy_Tc_df %>% 
  filter(Sample != "Healthy") %>% 
  ggplot(aes(x = Pyr_Uri, y = Pct_division, group=Sample, color=Sample))+
  geom_line() + geom_point() + scale_color_manual(values=c( "black", "dodgerblue3","firebrick")) +
   pretty_plot(fontsize = 6)+ facet_wrap(~celltype) + xlab("P&U") + ylab("% divided") +
  theme(legend.position = "none") +
  scale_y_continuous()
cowplot::ggsave2(t1, file = "divided.pdf", width = 1.8, height = 1)

t2 <- mitopathy_Tc_df %>% 
  filter(Sample != "Healthy") %>% 
  ggplot(aes(x = Pyr_Uri, y = Viabilty, group=Sample, color=Sample))+
  geom_line() + geom_point() + scale_color_manual(values=c("black", "dodgerblue3","firebrick")) +
  ylim(0, NA) + pretty_plot(fontsize = 6)+ facet_wrap(~celltype) + xlab("P&U") +
  theme(legend.position = "none") + ylab("% viable")+
  scale_y_continuous()
cowplot::ggsave2(t2, file = "viability.pdf", width = 1.8, height = 1)

mitopathy_Tc_df %>% 
  ggplot(aes(x = Pyr_Uri, y = Pro_Index, group=Sample, color=Sample))+
  geom_line() + geom_point() + scale_color_manual(values=c("green", "black", "blue", "red")) +
  ylim(0, NA) + theme_bw() + facet_wrap(~celltype) + xlab("1mM Pyruvate + 200nM Uridine") + ylab("Proliferartion Index")

# CD4 viability
t.test(c(48.3, 27.6, 44.6), c(16.7, 0.22, 32), paired = TRUE)

# CD8 viability
t.test(c(66.1, 42.1, 75.1), c(21.2, 0.97, 62.3), paired = TRUE)

# CD4 pct division
t.test(c(27.1, 19.9, 15.4), c(13.2, 0, 12.2), paired = TRUE)

# CD8 pct division
t.test(c(18.4, 25.1, 22.4), c(1.62, 0, 18), paired = TRUE)

# CD4 proliferation index
t.test(c(1.16, 1.08, 1.12), c(1, 1, 1), paired = TRUE)

# CD8 proliferation index
t.test(c(1.23, 1.1, 1.18), c(1, 1, 1.16), paired = TRUE)


