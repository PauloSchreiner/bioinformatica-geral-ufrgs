# ============================================================
# ETAPA 03 - GENES MARCADORES DOS CLUSTERS
# ============================================================
#
# Fluxograma: (9) identificação de genes marcadores de clusters
#
# Entrada : cache/02_objeto_processado.rds
# Saída   : tabelas de marcadores, scores de painéis canônicos e uma
#           PROPOSTA automática de anotação (revisada na etapa 04).
# ============================================================

source(file.path(dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[
  grep("^--file=", commandArgs(trailingOnly = FALSE))])), "00_config.R"))

step("ETAPA 03 - Marcadores")
mem_report("início")

combined <- readRDS(file.path(cache_dir, "02_objeto_processado.rds"))
Idents(combined) <- "seurat_clusters"
mem_report("objeto carregado")

# ------------------------------------------------------------
# (9) FindAllMarkers
# ------------------------------------------------------------
# Compara cada cluster contra TODAS as outras células (one-vs-rest)
# com o teste de Wilcoxon (não-paramétrico; não assume normalidade,
# adequado a dados esparsos e com zeros).
#
# only.pos = TRUE        -> só genes ENRIQUECIDOS no cluster; são esses
#                           que servem para rotular o cluster.
# min.pct = 0.25         -> o gene tem de ser detectado em >=25% das
#                           células de um dos lados; descarta genes
#                           raríssimos cujo p-valor é instável.
# logfc.threshold = 0.25 -> exige diferença mínima de ~1,19x; evita
#                           milhares de "significâncias" triviais que
#                           o n enorme de células produz.
#
# LEITURA CORRETA DOS p-VALORES: cada célula é tratada como réplica,
# o que infla a significância. Aqui o p-valor serve para ORDENAR genes
# candidatos, não como evidência estatística de diferença entre
# condições experimentais.
step("(9) FindAllMarkers - Wilcoxon, one-vs-rest")
if (requireNamespace("presto", quietly = TRUE))
  message("presto disponível: Wilcoxon acelerado.")

# Etapa cara (~3 min): se a tabela já existe, reaproveita.
markers_file <- file.path(results_dir, "03_marcadores_todos.csv")

if (file.exists(markers_file)) {
  message("Reaproveitando ", basename(markers_file),
          " (apague o arquivo para recalcular).")
  markers <- read.csv(markers_file, stringsAsFactors = FALSE)
  markers$cluster <- factor(markers$cluster,
                            levels = levels(combined$seurat_clusters))
} else {
  t0 <- Sys.time()
  markers <- FindAllMarkers(combined, assay = "RNA", only.pos = TRUE,
                            min.pct = marker_min_pct,
                            logfc.threshold = marker_logfc,
                            verbose = FALSE)
  message("Tempo: ", round(difftime(Sys.time(), t0, units = "min"), 1), " min")
  mem_report("após FindAllMarkers")

  markers <- markers %>%
    mutate(pct_diff = round(pct.1 - pct.2, 3)) %>%
    arrange(cluster, desc(avg_log2FC))

  save_table(markers, "03_marcadores_todos.csv")
}
message("Marcadores encontrados: ", nrow(markers),
        " em ", nlevels(combined$seurat_clusters), " clusters")

# Top 10 por cluster, ordenados por magnitude do efeito (não por p).
top10 <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()
save_table(top10, "03_marcadores_top10_por_cluster.csv")

# Versão "cartão de visita": uma linha por cluster, fácil de ler.
resumo <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(
    n_marcadores = n(),
    top10_genes = paste(head(gene[order(-avg_log2FC)], 10), collapse = ", "),
    .groups = "drop"
  ) %>%
  left_join(as.data.frame(table(cluster = combined$seurat_clusters)),
            by = "cluster") %>%
  rename(n_celulas = Freq) %>%
  relocate(n_celulas, .after = cluster)
save_table(resumo, "03_marcadores_resumo_por_cluster.csv")
print(as.data.frame(resumo[, c("cluster", "n_celulas", "top10_genes")]))

# ------------------------------------------------------------
# HEATMAP DOS MARCADORES
# ------------------------------------------------------------
# Subamostramos 100 células por cluster: o heatmap fica legível e
# evitamos criar uma matriz densa de 20.000 colunas.
step("Heatmap dos marcadores")
top5 <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5, with_ties = FALSE) %>%
  pull(gene) %>% unique()

set.seed(123)
sub <- subset(combined, downsample = heatmap_cells_per_cluster)
sub <- ScaleData(sub, features = top5, verbose = FALSE)

p_heat <- DoHeatmap(sub, features = top5, size = 3.2, raster = TRUE) +
  scale_fill_gradientn(colours = c("#2C7BB6", "white", "#D7191C"),
                       na.value = "white", name = "expressão\nescalonada") +
  theme(axis.text.y = element_text(size = 6.5)) +
  labs(title = "Top 5 marcadores por cluster",
       subtitle = paste0(heatmap_cells_per_cluster,
                         " células por cluster (subamostragem para visualização)"))

save_plot(p_heat, "04_Marcadores", "01_heatmap_top5_marcadores.png",
          width = 13, height = 15)
rm(sub); gc(verbose = FALSE)

# Dotplot com os 3 melhores marcadores de cada cluster.
top3 <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 3, with_ties = FALSE) %>%
  pull(gene) %>% unique()

p_dot_top <- DotPlot(combined, features = top3, cluster.idents = FALSE) +
  RotatedAxis() +
  scale_colour_gradientn(colours = c("#EFEFEF", "#4C78A8", "#1B3A5C")) +
  labs(title = "Top 3 marcadores por cluster",
       x = NULL, y = "Cluster") +
  theme(axis.text.x = element_text(size = 7))

save_plot(p_dot_top, "04_Marcadores", "02_dotplot_top3_marcadores.png",
          width = 20, height = 7)

# ------------------------------------------------------------
# PAINÉIS DE MARCADORES CANÔNICOS DO PULMÃO DE CAMUNDONGO
# ------------------------------------------------------------
# A anotação NÃO pode depender de um gene só. Cada tipo celular é
# representado por um painel de marcadores consagrados na literatura
# de pulmão murino (Tabula Muris, LungMAP, atlas de fibrose por bleomicina).
marker_panels <- list(
  "Fibroblastos"        = c("Col1a1","Col1a2","Dcn","Pdgfra","Lum","Col3a1","Gsn"),
  "Miofibroblastos"     = c("Acta2","Tagln","Postn","Cthrc1","Tnc","Eln"),
  "Macrofagos_alveolares" = c("Siglecf","Ear2","Marco","Atp6v0d2","Krt79","Chil3"),
  "Macrofagos_intersticiais" = c("C1qa","C1qb","C1qc","Mrc1","Csf1r","Adgre1"),
  "Monocitos"           = c("Ccr2","Ly6c2","Plac8","Ms4a4c","Chil3"),
  "Neutrofilos"         = c("S100a8","S100a9","Retnlg","Mmp9","Il1b"),
  "Celulas_dendriticas" = c("Cd209a","H2-Ab1","Cd74","Itgax","Flt3"),
  "Celulas_T"           = c("Cd3d","Cd3e","Cd3g","Trbc2","Lck","Thy1"),
  "Celulas_B"           = c("Cd79a","Cd79b","Ms4a1","Ighm","Ebf1"),
  "Celulas_NK"          = c("Nkg7","Klrd1","Gzma","Ncr1","Prf1"),
  "Epitelio_AT1"        = c("Ager","Pdpn","Hopx","Akap5","Rtkn2"),
  "Epitelio_AT2"        = c("Sftpc","Sftpb","Sftpa1","Lamp3","Abca3"),
  "Celulas_Club"        = c("Scgb1a1","Scgb3a2","Cyp2f2","Cckar"),
  "Celulas_ciliadas"    = c("Foxj1","Tmem212","Dynlrb2","Pifo"),
  "Endotelio_capilar"   = c("Pecam1","Cdh5","Cldn5","Ptprb","Gpihbp1"),
  "Endotelio_linfatico" = c("Mmrn1","Prox1","Lyve1","Ccl21a"),
  "Pericitos_mural"     = c("Rgs5","Pdgfrb","Cspg4","Notch3"),
  "Mesotelio"           = c("Msln","Upk3b","Wt1","Lrrn4"),
  "Eritrocitos"         = c("Hba-a1","Hbb-bs","Alas2"),
  "Proliferativas"      = c("Mki67","Top2a","Birc5","Ccnb1")
)

presentes <- lapply(marker_panels, intersect, y = rownames(combined))
ausentes  <- lapply(marker_panels, setdiff,   y = rownames(combined))

# Para o DotPlot cada gene só pode aparecer uma vez (ver 00_config.R).
# Os scores de módulo, abaixo, continuam usando os painéis COMPLETOS.
presentes_dot <- dedup_paineis(presentes)

save_table(
  data.frame(
    painel = names(marker_panels),
    genes_presentes = sapply(presentes, paste, collapse = ", "),
    genes_ausentes  = sapply(ausentes,  paste, collapse = ", ")
  ),
  "03_paineis_marcadores_canonicos.csv"
)

# DotPlot de todos os painéis: a figura central da anotação.
p_dot_panel <- DotPlot(combined, features = presentes_dot) +
  RotatedAxis() +
  scale_colour_gradientn(colours = c("#EFEFEF", "#4C78A8", "#1B3A5C")) +
  labs(title = "Painéis de marcadores canônicos do pulmão murino",
       subtitle = "tamanho = % de células que expressam | cor = expressão média escalonada",
       x = NULL, y = "Cluster") +
  theme(axis.text.x = element_text(size = 7),
        strip.text.x = element_text(size = 7, angle = 90, hjust = 0))

save_plot(p_dot_panel, "05_Anotacao", "01_dotplot_marcadores_canonicos.png",
          width = 26, height = 8)

# ------------------------------------------------------------
# SCORE DE MÓDULO POR PAINEL  ->  PROPOSTA DE ANOTAÇÃO
# ------------------------------------------------------------
# AddModuleScore calcula, por célula, a expressão média do painel menos
# a de um conjunto-controle de genes com expressão basal semelhante.
# Isso corrige o viés de profundidade: um score alto significa
# "este painel está acima do esperado NESTA célula".
step("Scores de módulo por painel celular")
paineis_ok <- presentes[sapply(presentes, length) >= 2]

combined <- AddModuleScore(combined, features = paineis_ok,
                           name = "PANEL_", seed = 123, ctrl = 50)

score_cols <- paste0("PANEL_", seq_along(paineis_ok))
colnames(combined@meta.data)[match(score_cols, colnames(combined@meta.data))] <-
  paste0("score_", names(paineis_ok))

score_mat <- combined@meta.data %>%
  select(cluster = seurat_clusters, starts_with("score_")) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean), .groups = "drop")

names(score_mat) <- sub("^score_", "", names(score_mat))
save_table(score_mat, "03_scores_medios_por_cluster.csv")

# Proposta automática: para cada cluster, os 3 painéis de maior score.
mat <- as.matrix(score_mat[, -1])
rownames(mat) <- score_mat$cluster

proposta <- do.call(rbind, lapply(seq_len(nrow(mat)), function(i) {
  o <- order(mat[i, ], decreasing = TRUE)
  data.frame(
    cluster = rownames(mat)[i],
    proposta_1 = colnames(mat)[o[1]], score_1 = round(mat[i, o[1]], 3),
    proposta_2 = colnames(mat)[o[2]], score_2 = round(mat[i, o[2]], 3),
    proposta_3 = colnames(mat)[o[3]], score_3 = round(mat[i, o[3]], 3)
  )
}))

proposta <- proposta %>%
  left_join(resumo[, c("cluster", "n_celulas", "top10_genes")], by = "cluster") %>%
  relocate(n_celulas, .after = cluster)

save_table(proposta, "03_anotacao_proposta_automatica.csv")
print(proposta[, c("cluster", "n_celulas", "proposta_1", "score_1", "proposta_2")])

# Heatmap dos scores: cluster x painel. Mostra de relance quem é quem.
if (requireNamespace("pheatmap", quietly = TRUE)) {
  z <- scale(mat)
  png(file.path(figures_dir, "05_Anotacao", "02_heatmap_scores_paineis.png"),
      width = 2200, height = 1500, res = 200)
  pheatmap::pheatmap(
    t(z), cluster_cols = TRUE, cluster_rows = TRUE,
    color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
    main = "Score médio (z-score) de cada painel celular por cluster",
    fontsize = 9, border_color = NA)
  dev.off()
  message("  [fig] 05_Anotacao/02_heatmap_scores_paineis.png")
}

saveRDS(combined, file.path(cache_dir, "03_objeto_com_scores.rds"))
mem_report("fim da etapa 03")
message("\nETAPA 03 CONCLUÍDA.")
