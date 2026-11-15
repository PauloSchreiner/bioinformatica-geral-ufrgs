# ============================================================
# ETAPA 04 - ANOTAÇÃO DOS CLUSTERS, COMPOSIÇÃO E FIBROSE
# ============================================================
#
# Fluxograma: (10) anotação dos clusters  + discussão biológica
#
# Entrada : cache/03_objeto_com_scores.rds
#           results/04_anotacao_curada.csv  <- EDITÁVEL À MÃO
# Saída   : figuras finais, tabelas de composição e objeto anotado.
#
# A anotação NÃO é automática: a etapa 03 gera uma PROPOSTA
# (03_anotacao_proposta_automatica.csv) a partir dos scores de painéis;
# esta etapa aplica a versão CURADA, conferida contra os marcadores
# encontrados de novo pelo FindAllMarkers.
# ============================================================

source(file.path(dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[
  grep("^--file=", commandArgs(trailingOnly = FALSE))])), "00_config.R"))

step("ETAPA 04 - Anotação final")
mem_report("início")

combined <- readRDS(file.path(cache_dir, "03_objeto_com_scores.rds"))
Idents(combined) <- "seurat_clusters"

anot_file <- file.path(results_dir, "04_anotacao_curada.csv")
if (!file.exists(anot_file))
  stop("Falta o arquivo de anotação curada: ", anot_file)

anot <- read.csv(anot_file, stringsAsFactors = FALSE)
anot$cluster <- as.character(anot$cluster)
stopifnot(all(levels(combined$seurat_clusters) %in% anot$cluster))

mapa <- setNames(anot$tipo_celular, anot$cluster)
# unname() é essencial: o vetor vindo de mapa[...] carrega os números dos
# clusters como names, e o Seurat tenta casar esses names com os nomes das
# células ("No cell overlap between new meta data and Seurat object").
combined$tipo_celular <- factor(
  unname(mapa[as.character(combined$seurat_clusters)]),
  levels = unique(anot$tipo_celular[order(anot$ordem)])
)

mapa_g <- setNames(anot$compartimento, anot$cluster)
combined$compartimento <- factor(
  unname(mapa_g[as.character(combined$seurat_clusters)]),
  levels = unique(anot$compartimento[order(anot$ordem)])
)

Idents(combined) <- "tipo_celular"
print(table(combined$tipo_celular))

# Paleta estável por tipo celular, agrupada por compartimento (00_config.R).
tipos <- levels(combined$tipo_celular)
pal <- paleta_por_compartimento(anot)
pal <- pal[tipos]

pal_comp <- setNames(
  sapply(familias_compartimento[levels(combined$compartimento)],
         function(a) colorRampPalette(a)(3)[2]),
  levels(combined$compartimento))

# ------------------------------------------------------------
# UMAP ANOTADO
# ------------------------------------------------------------
p_anot <- DimPlot(combined, reduction = "umap", group.by = "tipo_celular",
                  label = TRUE, label.size = 3.8, repel = TRUE,
                  cols = pal, raster = FALSE) +
  NoLegend() +
  labs(title = "UMAP anotado - tipos celulares do pulmão murino",
       subtitle = paste0(ncol(combined), " células | ",
                         length(tipos), " populações"))
save_plot(p_anot, "05_Anotacao", "03_UMAP_anotado.png", width = 10, height = 8.5)

p_anot_leg <- DimPlot(combined, reduction = "umap", group.by = "tipo_celular",
                      cols = pal, raster = FALSE) +
  labs(title = "UMAP anotado (com legenda)")
save_plot(p_anot_leg, "05_Anotacao", "04_UMAP_anotado_legenda.png",
          width = 11, height = 8)

p_comp <- DimPlot(combined, reduction = "umap", group.by = "compartimento",
                  cols = pal_comp, raster = FALSE) +
  labs(title = "UMAP por compartimento tecidual")
save_plot(p_comp, "05_Anotacao", "05_UMAP_compartimentos.png",
          width = 9, height = 7)

p_anot_split <- DimPlot(combined, reduction = "umap", group.by = "tipo_celular",
                        split.by = "sample", cols = pal, ncol = 4,
                        raster = FALSE) +
  labs(title = "Tipos celulares ao longo do tempo")
save_plot(p_anot_split, "05_Anotacao", "06_UMAP_anotado_por_timepoint.png",
          width = 22, height = 6)

# DotPlot final: tipos celulares x marcadores que os definem.
marcadores_finais <- unique(unlist(strsplit(anot$marcadores_chave, "\\s*;\\s*")))
marcadores_finais <- intersect(marcadores_finais, rownames(combined))

p_dot_final <- DotPlot(combined, features = marcadores_finais,
                       group.by = "tipo_celular") +
  RotatedAxis() +
  scale_colour_gradientn(colours = c("#EFEFEF", "#4C78A8", "#1B3A5C")) +
  labs(title = "Marcadores que sustentam cada anotação",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(size = 8))
save_plot(p_dot_final, "05_Anotacao", "07_dotplot_anotacao_final.png",
          width = 18, height = 7)

save_table(anot, "04_anotacao_final_aplicada.csv")

# ------------------------------------------------------------
# COMPOSIÇÃO CELULAR AO LONGO DO TEMPO
# ------------------------------------------------------------
# NOTA METODOLÓGICA (importante para a discussão):
# há UMA biblioteca por timepoint. Portanto não existem réplicas
# biológicas independentes e NÃO cabe teste estatístico de diferença
# de composição. As proporções abaixo são DESCRITIVAS e refletem as
# células recuperadas após dissociação, QC e subamostragem — não a
# abundância real no pulmão.
step("Composição celular por timepoint")

comp <- as.data.frame(table(sample = combined$sample,
                            tipo_celular = combined$tipo_celular)) %>%
  rename(n_celulas = Freq) %>%
  group_by(sample) %>%
  mutate(proporcao = n_celulas / sum(n_celulas),
         pct = round(100 * proporcao, 2)) %>%
  ungroup()
save_table(comp, "04_composicao_por_timepoint.csv")

# Tabela larga, mais fácil de ler num slide.
comp_wide <- comp %>%
  select(tipo_celular, sample, pct) %>%
  tidyr::pivot_wider(names_from = sample, values_from = pct)
save_table(comp_wide, "04_composicao_por_timepoint_tabela.csv")
print(as.data.frame(comp_wide))

p_stack <- ggplot(comp, aes(sample, proporcao, fill = tipo_celular)) +
  geom_col(width = .75) +
  scale_fill_manual(values = pal, name = "Tipo celular") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(labels = sample_labels) +
  labs(title = "Composição celular em cada timepoint",
       subtitle = "5.000 células por biblioteca (subamostragem uniforme)",
       x = NULL, y = "Proporção das células recuperadas") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_stack, "05_Anotacao", "08_composicao_empilhada.png",
          width = 10, height = 7)

p_lines <- ggplot(comp, aes(sample, pct, group = tipo_celular,
                            colour = tipo_celular)) +
  geom_line(linewidth = .9) + geom_point(size = 2) +
  scale_colour_manual(values = pal, guide = "none") +
  facet_wrap(~ tipo_celular, scales = "free_y", ncol = 5) +
  scale_x_discrete(labels = c(Control = "Ctrl", D21 = "D21",
                              D42 = "D42", D63 = "D63")) +
  labs(title = "Dinâmica de cada população ao longo da resolução da fibrose",
       subtitle = "eixo y livre por painel | descritivo, sem réplicas biológicas",
       x = NULL, y = "% das células da biblioteca") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(size = 8, face = "bold"))
save_plot(p_lines, "05_Anotacao", "09_dinamica_por_populacao.png",
          width = 16, height = 9)

# Mesma informação como heatmap de variação relativa ao controle.
comp_fc <- comp %>%
  select(tipo_celular, sample, pct) %>%
  tidyr::pivot_wider(names_from = sample, values_from = pct) %>%
  mutate(across(all_of(samples[-1]),
                ~ log2((.x + 0.05) / (Control + 0.05)))) %>%
  select(-Control) %>%
  tidyr::pivot_longer(-tipo_celular, names_to = "sample", values_to = "log2FC")

p_fc <- ggplot(comp_fc, aes(sample, tipo_celular, fill = log2FC)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%+.1f", log2FC)), size = 3) +
  scale_fill_gradientn(colours = c("#2C7BB6", "white", "#D7191C"),
                       limits = c(-3, 3), oob = scales::squish,
                       name = "log2(vs\ncontrole)") +
  labs(title = "Variação da abundância relativa em relação ao controle",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11)
save_plot(p_fc, "05_Anotacao", "10_variacao_vs_controle.png",
          width = 8, height = 8)

# ------------------------------------------------------------
# GENES DE FIBROSE
# ------------------------------------------------------------
step("Exploração dos genes de fibrose")

fibrose <- list(
  "Colageno_MEC"   = c("Col1a1","Col1a2","Col3a1","Fn1","Eln","Postn","Sparc"),
  "Miofibroblasto" = c("Acta2","Tagln","Cthrc1","Tnc","Thbs2"),
  "Sinalizacao"    = c("Tgfb1","Pdgfrb","Ctgf","Ccn2","Serpine1"),
  "Remodelamento"  = c("Mmp2","Mmp14","Timp1","Lox","Loxl1"),
  "Epitelio_alveolar" = c("Sftpc","Ager","Krt8","Areg")
)
fibrose <- lapply(fibrose, intersect, y = rownames(combined))
fibrose <- fibrose[sapply(fibrose, length) > 0]

p_fib_dot <- DotPlot(combined, features = fibrose, group.by = "tipo_celular") +
  RotatedAxis() +
  scale_colour_gradientn(colours = c("#EFEFEF", "#E45756", "#7F1416")) +
  labs(title = "Genes de fibrose e remodelamento por tipo celular",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p_fib_dot, "06_Fibrose", "01_dotplot_genes_fibrose.png",
          width = 18, height = 7)

genes_destaque <- intersect(c("Col1a1","Acta2","Cthrc1","Postn",
                              "Sftpc","Ager","Krt8","Pecam1"),
                            rownames(combined))

p_feat <- FeaturePlot(combined, features = genes_destaque, ncol = 4,
                      order = TRUE, raster = FALSE,
                      cols = c("grey90", "#D7191C")) &
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank())
save_plot(p_feat, "06_Fibrose", "02_featureplot_genes_chave.png",
          width = 17, height = 8)

# Score de fibrose por célula e sua evolução temporal.
combined <- AddModuleScore(combined, features = fibrose["Colageno_MEC"],
                           name = "score_mec", seed = 123, ctrl = 50)
combined <- AddModuleScore(combined, features = fibrose["Miofibroblasto"],
                           name = "score_mio", seed = 123, ctrl = 50)

p_score_umap <- FeaturePlot(combined, features = c("score_mec1", "score_mio1"),
                            reduction = "umap", ncol = 2, order = TRUE,
                            raster = FALSE,
                            cols = c("grey90", "#4C78A8", "#7F1416")) &
  theme(axis.title = element_blank())
save_plot(p_score_umap, "06_Fibrose", "03_score_fibrose_umap.png",
          width = 13, height = 6)

# Evolução do score dentro de cada tipo celular.
score_df <- combined@meta.data %>%
  select(sample, tipo_celular, score_mec1, score_mio1)

p_score_time <- ggplot(score_df, aes(sample, score_mec1, fill = sample)) +
  geom_violin(scale = "width", colour = NA, alpha = .85) +
  geom_boxplot(width = .13, outlier.shape = NA, fill = "white", alpha = .85) +
  scale_fill_manual(values = sample_colors, guide = "none") +
  facet_wrap(~ tipo_celular, scales = "free_y", ncol = 5) +
  scale_x_discrete(labels = c(Control = "Ctrl", D21 = "D21",
                              D42 = "D42", D63 = "D63")) +
  labs(title = "Score de colágeno/matriz extracelular por tipo celular e timepoint",
       x = NULL, y = "score de módulo (MEC)") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(size = 8, face = "bold"))
save_plot(p_score_time, "06_Fibrose", "04_score_MEC_por_tipo_e_tempo.png",
          width = 16, height = 9)

score_resumo <- score_df %>%
  group_by(tipo_celular, sample) %>%
  summarise(score_MEC_medio = round(mean(score_mec1), 4),
            score_mio_medio = round(mean(score_mio1), 4),
            n = n(), .groups = "drop")
save_table(score_resumo, "04_scores_fibrose_por_tipo_e_timepoint.csv")

# Foco no compartimento mesenquimal (onde a fibrose acontece).
mes <- grep("ibroblast|ericit|esot|sculo liso", tipos, value = TRUE)
if (length(mes) > 0) {
  sub_mes <- subset(combined, subset = tipo_celular %in% mes)
  message("Células mesenquimais: ", ncol(sub_mes))

  p_mes <- VlnPlot(sub_mes,
                   features = intersect(c("Col1a1","Cthrc1","Acta2","Postn"),
                                        rownames(sub_mes)),
                   group.by = "sample", pt.size = 0, ncol = 4,
                   cols = unname(sample_colors)) &
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 9))
  save_plot(p_mes, "06_Fibrose", "05_mesenquima_genes_por_timepoint.png",
            width = 15, height = 5)

  mes_tab <- sub_mes@meta.data %>%
    group_by(sample, tipo_celular) %>%
    summarise(n = n(), .groups = "drop")
  save_table(mes_tab, "04_mesenquima_por_timepoint.csv")
  rm(sub_mes); gc(verbose = FALSE)
}

# ------------------------------------------------------------
step("Salvando objeto final e metadados")
saveRDS(combined, file.path(cache_dir, "04_objeto_anotado.rds"))

meta_out <- combined@meta.data
meta_out$celula <- rownames(meta_out)
save_table(meta_out, "04_metadados_das_celulas.csv")

writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "sessionInfo.txt"))
mem_report("fim da etapa 04")
message("\nETAPA 04 CONCLUÍDA.")
