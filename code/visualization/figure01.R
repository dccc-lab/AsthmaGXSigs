## ---- setup_env
library("here")

source(here("code", "models", "run_gsea.R"), echo = TRUE)


## ---- make_figure01
trem1_pathways <- fgseaRes1 %>% arrange(desc(NES)) %>% filter(padj < 0.05) %>% use_series(pathway)

pdf(file = here("reports", "publication", "Figure01.pdf"), width = 8, height = 7, pointsize = 8)
plotGseaTable(pathways = gsc_c7_control[trem1_pathways],
              stats = dge1_rnks,
              fgseaRes = fgseaRes1,
              colwidths = c(8.0, 3.0, 0.8, 1.0, 1.0),
              gseaParam = 1)
dev.off()


## ---- answer_main_text_queries
# How many cases and controls?
es_ubiopred_wb %>% pData %$% table(asthma)
es_ubiopred_wb %>% pData %$% table(cohort, asthma)
es_ubiopred_wb %>% pData %$% table(severe, asthma)
es_ubiopred_wb %>% pData %$% table(smoker, asthma)
es_ubiopred_wb %>% pData %>% filter(asthma == TRUE) %$% table(severe)
es_ubiopred_wb %>% pData %>% filter(asthma == TRUE) %$% table(smoker, severe)

es_ubiopred_wb %>% pData %$% table(cohort)
es_ubiopred_wb[, !es_ubiopred_wb$smoker] %>% pData %$% table(cohort)
es_ubiopred_wb[, !es_ubiopred_wb$smoker & es_ubiopred_wb$asthma] %>% pData %$% table(cohort)
es_ubiopred_wb[, es_ubiopred_wb$cohort %in% c("Healthy, non-smoking", "Severe asthma, non-smoking")] %>% pData %$% table(cohort)
es_ubiopred_wb[, es_ubiopred_wb$cohort %in% c("Healthy, non-smoking", "Severe asthma, non-smoking")] %>% pData %$% table(cohort)
es_ubiopred_wb[, es_ubiopred_wb$cohort %in% c("Healthy, non-smoking", "Moderate asthma, non-smoking")] %>% pData %$% table(cohort)

# How many genes tested for differential expression?
dge1 %>% use_series(frame) %>% nrow

# How many differentially expressed genes in each comparison?
dge1 %>% use_series(frame) %>% filter(p.adj < 0.05) %>% nrow # cases vs. ctrls. (all subjects)
dge2 %>% use_series(frame) %>% filter(p.adj < 0.05) %>% nrow # cases vs. ctrls. (no smokers)
dge3 %>% use_series(frame) %>% filter(p.adj < 0.05) %>% nrow # severe cases vs. moderate cases (no smokers)
dge4 %>% use_series(frame) %>% filter(p.adj < 0.05) %>% nrow # severe cases vs. ctrls. (no smokers)
dge5 %>% use_series(frame) %>% filter(p.adj < 0.05) %>% nrow # moderate cases vs. ctrls. (no smokers)

# How many TREM1 gene sets were tested?
fgseaRes1 %>% use_series(pathway) %>% length

# How large were the TREM-1 gene sets tested?
gsc_c7_control_trem1 %>% sapply(length) %>% summary

# How many TREM1 gene sets were enriched?
trem1_pathways %>% length

# What are the core enriched TREM1 pathway genes?
core_genes <- fgseaRes1 %>% arrange(NES) %>% filter(padj < 0.05) %>% use_series(leadingEdge) %>% unlist %>% unique %>% sort
core_genes %>% print
core_genes %>% length

# What is the overlap of this core with highlighted TREM1 genes from Table 3 of Croteau-Chonka et al. (2017)?
table3_genes <- c("CCL23", "OLIG1", "OLIG2", "GFOD1", "RHOBTB3", "HSD3B7")
intersect(core_genes, table3_genes)
dge1 %>% use_series(frame) %>% filter(geneSymbols %in% table3_genes)

# Which primary analysis gene sets were enriched in secondary analyses?
fgseaRes1 %>% arrange(desc(NES)) %>% select(-leadingEdge) %>% filter(padj < 0.05) %>% mutate_at(c("pval", "padj", "ES", "NES"), funs(signif(., 3)))
fgseaRes2 %>% filter(pathway %in% trem1_pathways) %>% select(-leadingEdge) %>% filter(padj < 0.05) %>% mutate_at(c("pval", "padj", "ES", "NES"), funs(signif(., 3)))
fgseaRes3 %>% filter(pathway %in% trem1_pathways) %>% select(-leadingEdge) %>% filter(padj < 0.05) %>% mutate_at(c("pval", "padj", "ES", "NES"), funs(signif(., 3)))
fgseaRes4 %>% filter(pathway %in% trem1_pathways) %>% select(-leadingEdge) %>% filter(padj < 0.05) %>% mutate_at(c("pval", "padj", "ES", "NES"), funs(signif(., 3)))
fgseaRes5 %>% filter(pathway %in% trem1_pathways) %>% select(-leadingEdge) %>% filter(padj < 0.05) %>% mutate_at(c("pval", "padj", "ES", "NES"), funs(signif(., 3)))

# Which primary analysis gene sets were NOT enriched in secondary analyses?
fgseaRes3 %>% filter(pathway %in% trem1_pathways) %>% select(-leadingEdge) %>% filter(padj > 0.05) %>% mutate_at(c("pval", "padj", "ES", "NES"), funs(signif(., 3)))
fgseaRes4 %>% filter(pathway %in% trem1_pathways) %>% select(-leadingEdge) %>% filter(padj > 0.05) %>% mutate_at(c("pval", "padj", "ES", "NES"), funs(signif(., 3)))


## ---- clean_up_env
session_info()