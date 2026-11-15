# ============================================================
# ETAPA 02 - NORMALIZAÇÃO, HVG, SCALING, PCA, CLUSTERS, UMAP
# ============================================================
#
# Fluxograma: (3) normalização -> (4) genes altamente variáveis ->
#             (5) scaling -> (6) PCA -> (7) clusterização -> (8) UMAP
#
# Entrada : cache/01_objetos_filtrados.rds (4 objetos pequenos)
# Saída   : cache/02_objeto_processado.rds + figuras + tabelas
# ============================================================

source(file.path(dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[
  grep("^--file=", commandArgs(trailingOnly = FALSE))])), "00_config.R"))

step("ETAPA 02 - Clusterização")
mem_report("início")

objs <- readRDS(file.path(cache_dir, "01_objetos_filtrados.rds"))

# ------------------------------------------------------------
# MERGE DAS 4 BIBLIOTECAS
# ------------------------------------------------------------
# add.cell.ids prefixa o nome de cada célula com a biblioteca de
# origem, evitando colisão de barcodes (o mesmo barcode 10x pode
# aparecer em bibliotecas diferentes).
#
# ATENÇÃO CONCEITUAL: isto é um MERGE, não uma INTEGRAÇÃO.
# Não estamos corrigindo batch effect — e a etapa 02b mostra por quê
# isso é uma decisão deliberada neste desenho experimental.
combined <- merge(objs[[1]], y = objs[-1],
                  add.cell.ids = names(objs), project = "GSE328925")
rm(objs); gc(verbose = FALSE)

# No Seurat v5+, o merge mantém uma "layer" de contagens por objeto.
# JoinLayers funde tudo numa única matriz, necessária para normalizar
# e testar marcadores no conjunto todo.
combined <- JoinLayers(combined, assay = "RNA")
DefaultAssay(combined) <- "RNA"
combined$sample <- factor(combined$sample, levels = samples)

message("Objeto combinado: ", nrow(combined), " genes x ",
        ncol(combined), " células")
print(table(combined$sample))
mem_report("após merge")

# ------------------------------------------------------------
# (3) NORMALIZAÇÃO
# ------------------------------------------------------------
# Problema: células diferentes foram sequenciadas a profundidades
# diferentes. Uma célula com 20.000 UMIs tem contagens maiores que uma
# com 2.000 simplesmente por isso, não por biologia.
#
# LogNormalize: conta_do_gene / total_de_UMIs_da_célula * 10.000, e
# depois log1p(). O fator 10.000 ("counts per 10k") é a convenção do
# Seurat; o log comprime a cauda de genes muito expressos, impedindo
# que uns poucos genes dominem a PCA.
step("(3) Normalização - LogNormalize, scale.factor = 10000")
combined <- NormalizeData(combined, normalization.method = "LogNormalize",
                          scale.factor = 10000, verbose = FALSE)
mem_report("após normalização")

# ------------------------------------------------------------
# (4) GENES ALTAMENTE VARIÁVEIS (HVG)
# ------------------------------------------------------------
# A maior parte dos 24.709 genes é ruído ou é expressa igualmente em
# todas as células: não ajuda a separar tipos celulares e só adiciona
# dimensões. O método "vst" modela a relação média-variância esperada
# e escolhe os genes cuja variância EXCEDE o esperado para a sua média.
# 2000 é o padrão do Seurat e funciona bem para tecidos heterogêneos.
step("(4) Genes altamente variáveis - vst, 2000 genes")
combined <- FindVariableFeatures(combined, selection.method = "vst",
                                 nfeatures = n_variable_features,
                                 verbose = FALSE)

top20 <- head(VariableFeatures(combined), 20)
message("Top 20 HVG: ", paste(top20, collapse = ", "))

p_hvg <- LabelPoints(VariableFeaturePlot(combined), points = top20,
                     repel = TRUE, max.overlaps = 30) +
  labs(title = "Genes altamente variáveis (vst)",
       subtitle = paste(n_variable_features,
                        "genes selecionados; 20 mais variáveis rotulados")) +
  theme(legend.position = "bottom")

save_plot(p_hvg, "02_HVG_PCA", "01_genes_altamente_variaveis.png",
          width = 10, height = 7)

hvg_df <- HVFInfo(combined)
hvg_df$gene <- rownames(hvg_df)
hvg_df$hvg_selecionado <- hvg_df$gene %in% VariableFeatures(combined)
hvg_df$rank_hvg <- match(hvg_df$gene, VariableFeatures(combined))
save_table(hvg_df[order(hvg_df$rank_hvg), ], "02_genes_variaveis.csv")

# ------------------------------------------------------------
# (5) SCALING
# ------------------------------------------------------------
# A PCA é sensível à escala: sem padronizar, genes muito expressos
# (ex.: colágenos, actina) teriam variância absoluta enorme e
# capturariam as primeiras PCs sozinhos. ScaleData centraliza cada gene
# em média 0 e desvio-padrão 1 — assim cada gene contribui pelo seu
# PADRÃO, não pela sua magnitude.
#
# MEMÓRIA: counts e data são ESPARSAS (93% de zeros, que não ocupam
# espaço). scale.data NÃO pode ser esparsa: ao centrar na média, cada zero
# vira -media/desvio, um número explícito. A matriz fica densa.
# Escalar os 24.709 genes custaria 24709 x 20000 x 8 bytes = 3,95 GB.
# Restringindo aos 2000 HVG — que é o que o RunPCA lê de qualquer forma —
# o custo cai para 320 MB.
step("(5) Scaling - apenas os 2000 HVG (economia de memória)")
combined <- ScaleData(combined, features = VariableFeatures(combined),
                      verbose = FALSE)
mem_report("após scaling")

# ------------------------------------------------------------
# (6) PCA
# ------------------------------------------------------------
# Reduz 2000 dimensões (genes) a algumas dezenas de componentes, cada
# um sendo uma combinação linear de genes que covariam — na prática,
# "programas" de expressão. É linear, determinística e serve de base
# tanto para o grafo de vizinhos quanto para o UMAP.
step("(6) PCA - 50 componentes")
combined <- RunPCA(combined, features = VariableFeatures(combined),
                   npcs = n_pcs, seed.use = 123, verbose = FALSE)

save_plot(VizDimLoadings(combined, dims = 1:4, reduction = "pca", ncol = 4),
          "02_HVG_PCA", "02_PCA_loadings.png", width = 14, height = 7)

p_pca <- DimPlot(combined, reduction = "pca", group.by = "sample",
                 cols = sample_colors, shuffle = TRUE, raster = FALSE) +
  labs(title = "PCA (PC1 x PC2) por biblioteca")
save_plot(p_pca, "02_HVG_PCA", "03_PCA_por_biblioteca.png",
          width = 8, height = 6.5)

# ElbowPlot: desvio-padrão explicado por PC. O "cotovelo" indica a
# partir de onde as PCs capturam sobretudo ruído técnico.
sdev <- Stdev(combined, reduction = "pca")
elbow_df <- data.frame(PC = seq_along(sdev), sd = sdev)

p_elbow <- ggplot(elbow_df, aes(PC, sd)) +
  geom_point(size = 1.8, colour = "#4C78A8") +
  geom_vline(xintercept = max(dims_use), linetype = "dashed",
             colour = "#E45756") +
  annotate("text", x = max(dims_use) + 1, y = max(sdev) * .85,
           label = paste0("PCs usadas = ", max(dims_use)),
           hjust = 0, colour = "#E45756", size = 4) +
  labs(title = "ElbowPlot - variância capturada por componente principal",
       subtitle = "a curva achata depois de ~20 PCs",
       x = "Componente principal", y = "Desvio-padrão") +
  theme_classic(base_size = 12)

save_plot(p_elbow, "02_HVG_PCA", "04_PCA_elbow.png", width = 8, height = 6)
save_table(elbow_df, "02_PCA_desvio_padrao.csv")

# Genes que mais carregam nas primeiras PCs (útil para interpretar).
loadings <- Loadings(combined, "pca")
pc_top <- do.call(rbind, lapply(1:6, function(i) {
  l <- sort(loadings[, i], decreasing = TRUE)
  data.frame(PC = paste0("PC", i),
             positivos = paste(names(head(l, 10)), collapse = ", "),
             negativos = paste(names(head(rev(l), 10)), collapse = ", "))
}))
save_table(pc_top, "02_PCA_genes_por_componente.csv")

# ------------------------------------------------------------
# (7) CLUSTERIZAÇÃO
# ------------------------------------------------------------
# FindNeighbors constrói um grafo KNN no espaço das 20 PCs e o refina
# para SNN (compartilhamento de vizinhos), que é mais robusto a ruído.
# FindClusters aplica Louvain: procura comunidades densamente conectadas.
# A "resolution" é o botão de granularidade — maior = mais clusters.
step("(7) Clusterização - grafo SNN + Louvain")
combined <- FindNeighbors(combined, dims = dims_use, verbose = FALSE)

# Resolução principal + resoluções extras (para mostrar estabilidade).
for (res in sort(unique(c(cluster_resolution, extra_resolutions)))) {
  combined <- FindClusters(combined, resolution = res,
                           random.seed = 123, verbose = FALSE)
}

res_col <- paste0("RNA_snn_res.", cluster_resolution)
combined$seurat_clusters <- combined[[res_col]][, 1]
Idents(combined) <- "seurat_clusters"

n_clusters_df <- data.frame(
  resolucao = sort(unique(c(cluster_resolution, extra_resolutions))),
  n_clusters = sapply(sort(unique(c(cluster_resolution, extra_resolutions))),
                      function(r) nlevels(combined[[paste0("RNA_snn_res.", r)]][, 1]))
)
print(n_clusters_df)
save_table(n_clusters_df, "02_clusters_por_resolucao.csv")

message("Resolução usada: ", cluster_resolution, " -> ",
        nlevels(combined$seurat_clusters), " clusters")
print(table(combined$seurat_clusters))
mem_report("após clusterização")

# ------------------------------------------------------------
# (8) UMAP
# ------------------------------------------------------------
# Projeção NÃO-LINEAR das mesmas 20 PCs em 2D, apenas para visualizar.
# Preserva bem a vizinhança local (quem é parecido fica perto) e mal a
# distância global — a distância entre duas ilhas distantes NÃO é
# proporcional à diferença biológica entre elas.
step("(8) UMAP")
combined <- RunUMAP(combined, dims = dims_use, seed.use = 123, verbose = FALSE)

p_umap_cl <- DimPlot(combined, reduction = "umap", label = TRUE,
                     label.size = 4.5, repel = TRUE, raster = FALSE) +
  NoLegend() +
  labs(title = paste0("UMAP - ", nlevels(combined$seurat_clusters),
                      " clusters (resolução ", cluster_resolution, ")"))

p_umap_s <- DimPlot(combined, reduction = "umap", group.by = "sample",
                    cols = sample_colors, shuffle = TRUE, raster = FALSE) +
  labs(title = "UMAP por biblioteca / timepoint")

save_plot(p_umap_cl, "03_UMAP", "01_UMAP_clusters.png", width = 8.5, height = 7)
save_plot(p_umap_s,  "03_UMAP", "02_UMAP_por_biblioteca.png", width = 8.5, height = 7)
save_plot(p_umap_cl | p_umap_s, "03_UMAP", "03_UMAP_lado_a_lado.png",
          width = 15, height = 6.5)

p_split <- DimPlot(combined, reduction = "umap", split.by = "sample",
                   label = TRUE, label.size = 3.5, repel = TRUE,
                   ncol = 4, raster = FALSE) +
  NoLegend() +
  labs(title = "UMAP separado por timepoint",
       subtitle = "mesma projeção; cada painel mostra só as células daquela biblioteca")
save_plot(p_split, "03_UMAP", "04_UMAP_split_por_biblioteca.png",
          width = 18, height = 5.5)

# Densidade por amostra: destaca onde cada timepoint concentra células.
umap_df <- data.frame(Embeddings(combined, "umap"),
                      sample = combined$sample,
                      cluster = combined$seurat_clusters)
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")

p_dens <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
  geom_point(data = transform(umap_df, sample = NULL),
             colour = "grey88", size = .25) +
  geom_point(aes(colour = sample), size = .3) +
  scale_colour_manual(values = sample_colors, guide = "none") +
  facet_wrap(~ sample, nrow = 1, labeller = as_labeller(sample_labels)) +
  theme_void(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 12)) +
  labs(title = "Ocupação do UMAP por timepoint",
       subtitle = "cinza = todas as células; cor = células daquela biblioteca")

save_plot(p_dens, "03_UMAP", "05_UMAP_densidade_por_biblioteca.png",
          width = 16, height = 4.8)

# ------------------------------------------------------------
# (opcional) AVALIAÇÃO DE BATCH EFFECT / INTEGRAÇÃO
# ------------------------------------------------------------
# O enunciado pede: "se necessário, realizar o passo de integração".
# Aqui rodamos Harmony APENAS como diagnóstico, para mostrar o efeito
# e justificar a decisão de NÃO usá-lo na análise principal:
# neste desenho há exatamente UMA biblioteca por condição, então
# "batch" e "condição biológica" são a mesma variável. Corrigir o batch
# apagaria justamente o sinal que queremos observar.
if (requireNamespace("harmony", quietly = TRUE)) {
  step("(extra) Harmony - diagnóstico de batch effect")
  combined <- harmony::RunHarmony(combined, group.by.vars = "sample",
                                  reduction.use = "pca",
                                  reduction.save = "harmony", verbose = FALSE)
  combined <- RunUMAP(combined, reduction = "harmony", dims = dims_use,
                      reduction.name = "umap_harmony", seed.use = 123,
                      verbose = FALSE)

  p_h1 <- DimPlot(combined, reduction = "umap", group.by = "sample",
                  cols = sample_colors, shuffle = TRUE, raster = FALSE) +
    labs(title = "Sem integração (usado na análise)")
  p_h2 <- DimPlot(combined, reduction = "umap_harmony", group.by = "sample",
                  cols = sample_colors, shuffle = TRUE, raster = FALSE) +
    labs(title = "Com Harmony (apenas diagnóstico)")

  save_plot(p_h1 | p_h2, "03_UMAP", "06_UMAP_sem_vs_com_harmony.png",
            width = 15, height = 6.5)
  mem_report("após harmony")
}

# ------------------------------------------------------------
step("Salvando objeto processado")
saveRDS(combined, file.path(cache_dir, "02_objeto_processado.rds"))
mem_report("fim da etapa 02")
message("\nETAPA 02 CONCLUÍDA.")
