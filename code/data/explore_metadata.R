## ---- setup_env
library("here")

source(here("code", "data", "create_ubiopred_eset.R"), echo = TRUE)

library("reshape2")
library("ggplot2")
library("corrr")


## ---- explore_gx_pca
es_ubiopred_wb %>% pData %$% lm(gxPC1 ~ cohort) %>% tidy
es_ubiopred_wb %>% pData %$% lm(gxPC2 ~ cohort) %>% tidy

es_ubiopred_wb %>% pData %$% lm(gxPC1 ~ sex) %>% tidy
es_ubiopred_wb %>% pData %$% lm(gxPC2 ~ sex) %>% tidy

es_ubiopred_wb %>% pData %$% lm(gxPC1 ~ race) %>% tidy
es_ubiopred_wb %>% pData %$% lm(gxPC2 ~ race) %>% tidy

es_ubiopred_wb %>% pData %$% lm(gxPC1 ~ site) %>% tidy
es_ubiopred_wb %>% pData %$% lm(gxPC2 ~ site) %>% tidy

# pdf(file = here("reports", "diagnostics", "gxPCA.pdf"), width = 8, height = 8, pointsize = 8)

screeplot(res_pca_ubiopred_wb$pcs, type = "lines", main = "U-BIOPRED WB gxPCs")

pca2DPlot(pcaObj = res_pca_ubiopred_wb,
          plot.dim = c(1, 2),
          labelVariable = "cohort",
          title = "U-BIOPRED WB [Cohort]",
          plotOutPutFlag = FALSE,
          equalRange = FALSE,
          xlab = "PC1",
          ylab = "PC2",
          cex.legend = 0.6,
          cex = 1,
          cex.lab = 1,
          cex.axis = 1,
          legendPosition = "topleft"
)

pca2DPlot(pcaObj = res_pca_ubiopred_wb,
          plot.dim = c(1, 2),
          labelVariable = "sex",
          title = "U-BIOPRED WB [Sex]",
          plotOutPutFlag = FALSE,
          equalRange = FALSE,
          xlab = "PC1",
          ylab = "PC2",
          cex.legend = 0.6,
          cex = 1,
          cex.lab = 1,
          cex.axis = 1,
          legendPosition = "topleft"
)

pca2DPlot(pcaObj = res_pca_ubiopred_wb,
          plot.dim = c(1, 2),
          labelVariable = "race",
          title = "U-BIOPRED WB [Race]",
          plotOutPutFlag = FALSE,
          equalRange = FALSE,
          xlab = "PC1",
          ylab = "PC2",
          cex.legend = 0.6,
          cex = 1,
          cex.lab = 1,
          cex.axis = 1,
          legendPosition = "topleft"
)

pca2DPlot(pcaObj = res_pca_ubiopred_wb,
          plot.dim = c(1, 2),
          labelVariable = "site",
          title = "U-BIOPRED WB [Site]",
          plotOutPutFlag = FALSE,
          equalRange = FALSE,
          xlab = "PC1",
          ylab = "PC2",
          cex.legend = 0.6,
          cex = 1,
          cex.lab = 1,
          cex.axis = 1,
          legendPosition = "topleft"
)

# dev.off()


## ---- explore_cbcs
# pdf(file = here("reports", "diagnostics", "CBC.pdf"), width = 8, height = 8, pointsize = 8)

# plot estimated CBCs by asthma severity
pData(es_ubiopred_wb) %>%
    select(cohort, starts_with("EST")) %>% 
    melt(id.vars = "cohort") %>%
    mutate(variable = str_replace(variable, "EST", "")) %>%
    ggplot(aes(x = variable, y = value)) +
    geom_boxplot(aes(fill = cohort), notch = TRUE) +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "top",
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, title = "")) +
    labs(
        y = "Estimated Proportion",
        x = "Blood Cell Type",
        title = "U-BIOPRED WB") +
    ylim(0, 100)

# plot measured CBCs by asthma severity
pData(es_ubiopred_wb) %>%
    select(cohort, starts_with("PCT")) %>% 
    melt(id.vars = "cohort") %>%
    mutate(variable = str_replace(variable, "PCT", "")) %>%
    ggplot(aes(x = variable, y = value)) +
    geom_boxplot(aes(fill = cohort), notch = TRUE) +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "top",
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, title = "")) +
    labs(
        y = "Measured Proportion",
        x = "Blood Cell Type",
        title = "U-BIOPRED WB") +
    ylim(0, 100)

# plot measured vs. estimated CBCs
pData(es_ubiopred_wb) %>% 
    select(cohort, starts_with("PCT"), starts_with("EST")) %>% 
    ggplot(aes(x = ESTLYMPH, y = PCTLYMPH)) + 
    geom_point() +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "top",
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, title = "")) +
    labs(
        y = "Measured Proportion",
        x = "Estimated Proportion",
        title = "LYMPH") +
    ylim(0, 100) +
    xlim(0, 100) +
    geom_abline(intercept = 0, slope = 1, colour = "gray60", size = 0.25, linetype = "dashed") + 
    geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, size = 0.25) 

pData(es_ubiopred_wb) %>% 
    select(cohort, starts_with("PCT"), starts_with("EST")) %>% 
    ggplot(aes(x = ESTMONO, y = PCTMONO)) + 
    geom_point() +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "top",
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, title = "")) +
    labs(
        y = "Measured Proportion",
        x = "Estimated Proportion",
        title = "MONO") +
    ylim(0, 100) +
    xlim(0, 100) +
    geom_abline(intercept = 0, slope = 1, colour = "gray60", size = 0.25, linetype = "dashed") + 
    geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, size = 0.25) 

pData(es_ubiopred_wb) %>% 
    select(cohort, starts_with("PCT"), starts_with("EST")) %>% 
    ggplot(aes(x = ESTNEUT, y = PCTNEUT)) + 
    geom_point() +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "top",
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, title = "")) +
    labs(
        y = "Measured Proportion",
        x = "Estimated Proportion",
        title = "NEUT") +
    ylim(0, 100) +
    xlim(0, 100) +
    geom_abline(intercept = 0, slope = 1, colour = "gray60", size = 0.25, linetype = "dashed") + 
    geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, size = 0.25) 

# dev.off()

# measure correlations of estimated and measured CBCs
pData(es_ubiopred_wb) %$% cor.test(ESTLYMPH, PCTLYMPH) %>% tidy
pData(es_ubiopred_wb) %$% cor.test(ESTMONO, PCTMONO) %>% tidy
pData(es_ubiopred_wb) %$% cor.test(ESTNEUT, PCTNEUT) %>% tidy

pData(es_ubiopred_wb) %>% use_series(ESTMONO) %>% summary
pData(es_ubiopred_wb) %>% use_series(PCTMONO) %>% summary

# measure all pair-wise CBC correlations
pData(es_ubiopred_wb) %>% select(starts_with("PCT"), starts_with("EST")) %>% correlate %>% fashion
pData(es_ubiopred_wb) %>% select(starts_with("PCT"), starts_with("EST")) %>% correlate %>% stretch %>% fashion %>% filter(x != y, grepl("^PCT", x))

# measure correlations of estimated CBCs
pData(es_ubiopred_wb) %>% select(starts_with("PCT")) %>% correlate %>% stretch %>% fashion %>% filter(x != y) %>% arrange(r) %>% filter(duplicated(r))

pData(es_ubiopred_wb) %$% cor.test(PCTNEUT, PCTEOS) %>% tidy
pData(es_ubiopred_wb) %$% cor.test(PCTNEUT, PCTLYMPH) %>% tidy
pData(es_ubiopred_wb) %$% cor.test(PCTNEUT, PCTMONO) %>% tidy

# check relationships of CBCs with disease severity
pData(es_ubiopred_wb) %$% lm(PCTEOS ~ cohort) %>% tidy

pData(es_ubiopred_wb) %$% lm(PCTLYMPH ~ cohort) %>% tidy
pData(es_ubiopred_wb) %$% lm(ESTLYMPH ~ cohort) %>% tidy

pData(es_ubiopred_wb) %$% lm(PCTMONO ~ cohort) %>% tidy
pData(es_ubiopred_wb) %$% lm(ESTMONO ~ cohort) %>% tidy

pData(es_ubiopred_wb) %$% lm(PCTNEUT ~ cohort) %>% tidy
pData(es_ubiopred_wb) %$% lm(ESTNEUT ~ cohort) %>% tidy
