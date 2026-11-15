# ============================================================
# ETAPA 01 - MATRIZ DE CONTAGENS, QC E FILTRAGEM
# ============================================================
#
# Fluxograma: (1) matriz de contagens -> (2) QC e filtragem
#
# ESTRATÉGIA DE MEMÓRIA (máquina com 8 GB de RAM):
# as 4 bibliotecas somam 160 milhões de valores não-zero. Carregar as
# quatro simultaneamente custaria ~2 GB só nas contagens brutas, e o
# pico durante a leitura do .mtx é ~2x isso. Por isso cada biblioteca é
# lida, filtrada, subamostrada e descartada da memória ANTES da próxima.
# Só os metadados (tabela pequena) de TODAS as células são acumulados,
# de modo que as estatísticas de QC descrevem o dataset completo.
# ============================================================

source(file.path(dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[
  grep("^--file=", commandArgs(trailingOnly = FALSE))])), "00_config.R"))

step("ETAPA 01 - QC e filtragem")
mem_report("início")

# ------------------------------------------------------------
# DOWNLOAD (só roda se os arquivos ainda não estiverem no disco)
# ------------------------------------------------------------
# O GEO entrega os arquivos como "GSM9694319_control_barcodes.tsv.gz".
# O Read10X() exige exatamente barcodes/features/matrix dentro de uma
# pasta por amostra — por isso renomeamos no próprio download.
gsm_ids  <- c(Control = "GSM9694319", D21 = "GSM9694320",
              D42     = "GSM9694321", D63 = "GSM9694322")
gsm_sufx <- c(Control = "control", D21 = "d21", D42 = "d42", D63 = "d63")

# Um .gz truncado só daria erro muito depois, no meio da importação.
gzip_ok <- function(path) {
  file.exists(path) &&
    suppressWarnings(system2("gzip", c("-t", shQuote(path)),
                             stdout = FALSE, stderr = FALSE)) == 0
}

for (s in samples) {
  dir.create(file.path(data_dir, s), recursive = TRUE, showWarnings = FALSE)

  arquivos <- c(barcodes = "barcodes.tsv.gz",
                features = "features.tsv.gz",
                matrix   = "matrix.mtx.gz")

  for (tipo in names(arquivos)) {
    destino <- file.path(data_dir, s, arquivos[[tipo]])
    if (gzip_ok(destino)) next

    origem <- paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9694nnn/", gsm_ids[[s]],
      "/suppl/", gsm_ids[[s]], "_", gsm_sufx[[s]], "_", tipo,
      if (tipo == "matrix") ".mtx.gz" else ".tsv.gz")

    message("  baixando ", s, "/", arquivos[[tipo]])
    download.file(origem, destino, mode = "wb", method = "libcurl", quiet = TRUE)

    if (!gzip_ok(destino))
      stop("Download falhou ou o arquivo está corrompido: ", destino)
  }
}
message("Dados disponíveis em: ", data_dir)

meta_all   <- list()   # metadados de TODAS as células (antes do filtro)
objs_small <- list()   # objetos já filtrados e subamostrados
import_log <- list()

for (s in samples) {

  step(paste("Biblioteca:", s, "-", sample_labels[[s]]))

  # --- (1) MATRIZ DE CONTAGENS ------------------------------
  # Read10X lê o trio barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz
  # e devolve uma matriz ESPARSA (genes x células). Esparsa = os zeros,
  # que são ~93% da matriz, não ocupam memória.
  counts <- Read10X(data.dir = file.path(data_dir, s))

  # min.cells = 3 descarta genes vistos em menos de 3 células: são quase
  # sempre ruído de ambiente/PCR e encolhem bastante a matriz.
  # min.features = 0 porque o filtro por célula é decidido DEPOIS de
  # olhar as distribuições de QC, não antes.
  seu <- CreateSeuratObject(counts = counts, project = s,
                            min.cells = min_cells_per_gene,
                            min.features = 0)
  rm(counts); gc(verbose = FALSE)

  seu$sample <- s

  # --- (2) MÉTRICAS DE QC -----------------------------------
  # percent.mt  : % de UMIs de genes mitocondriais (prefixo "mt-" no mm10).
  #               Alta => membrana rompida, citoplasma vazou, célula morrendo.
  # percent.ribo: % de proteínas ribossomais (Rps*/Rpl*). Muito alta pode
  #               indicar célula com RNA degradado ou estresse de dissociação.
  # percent.hb  : % de hemoglobina. Pulmão é muito vascularizado; hemácias
  #               e contaminação por sangue aparecem aqui.
  seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern = "^mt-")
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^Rp[sl]")
  seu[["percent.hb"]]   <- PercentageFeatureSet(seu, pattern = "^Hb[ab]-")

  n_before <- ncol(seu)
  n_genes  <- nrow(seu)

  meta_all[[s]] <- seu@meta.data[, c("sample", "nCount_RNA", "nFeature_RNA",
                                     "percent.mt", "percent.ribo", "percent.hb")]

  # --- FILTRAGEM --------------------------------------------
  keep <- seu$nFeature_RNA >= min_features &
          seu$nFeature_RNA <= max_features &
          seu$nCount_RNA   >= min_counts &
          seu$percent.mt   <  max_mito

  seu <- subset(seu, cells = colnames(seu)[keep])
  n_after <- ncol(seu)

  # --- SUBAMOSTRAGEM (restrição de memória) -----------------
  # Amostragem aleatória simples, com semente fixa => reproduzível.
  set.seed(123)
  if (ncol(seu) > cells_per_sample) {
    seu <- subset(seu, cells = sample(colnames(seu), cells_per_sample))
  }
  n_kept <- ncol(seu)

  import_log[[s]] <- data.frame(
    sample          = s,
    label           = sample_labels[[s]],
    genes_detectados = n_genes,
    celulas_brutas   = n_before,
    celulas_pos_QC   = n_after,
    pct_aprovadas    = round(100 * n_after / n_before, 1),
    celulas_usadas   = n_kept
  )

  message(sprintf("  %s: %d células brutas -> %d pós-QC (%.1f%%) -> %d usadas",
                  s, n_before, n_after, 100 * n_after / n_before, n_kept))

  objs_small[[s]] <- seu
  rm(seu); gc(verbose = FALSE)
  mem_report(paste("após", s))
}

# ------------------------------------------------------------
# TABELAS DE QC
# ------------------------------------------------------------
meta_all <- bind_rows(meta_all)
meta_all$sample <- factor(meta_all$sample, levels = samples)

qc_stats <- meta_all %>%
  group_by(sample) %>%
  summarise(
    n_celulas    = n(),
    genes_med    = median(nFeature_RNA),
    genes_q05    = quantile(nFeature_RNA, 0.05),
    genes_q95    = quantile(nFeature_RNA, 0.95),
    umis_med     = median(nCount_RNA),
    umis_q05     = quantile(nCount_RNA, 0.05),
    umis_q95     = quantile(nCount_RNA, 0.95),
    mito_med     = round(median(percent.mt), 2),
    mito_q95     = round(quantile(percent.mt, 0.95), 2),
    ribo_med     = round(median(percent.ribo), 2),
    hb_med       = round(median(percent.hb), 2),
    .groups = "drop"
  )

save_table(qc_stats, "01_QC_metricas_por_biblioteca.csv")
save_table(bind_rows(import_log), "01_QC_resumo_filtragem.csv")
print(as.data.frame(qc_stats))
print(bind_rows(import_log))

# ------------------------------------------------------------
# FIGURAS DE QC  (feitas a partir da tabela de metadados, não do
# objeto Seurat -> não exige ter as 4 matrizes na memória)
# ------------------------------------------------------------
vln <- function(var, titulo, ylab, hlines = NULL, log = FALSE) {
  p <- ggplot(meta_all, aes(x = sample, y = .data[[var]], fill = sample)) +
    geom_violin(scale = "width", trim = TRUE, colour = NA, alpha = .85) +
    geom_boxplot(width = .12, outlier.shape = NA, fill = "white", alpha = .8) +
    scale_fill_manual(values = sample_colors, guide = "none") +
    scale_x_discrete(labels = sample_labels) +
    labs(title = titulo, x = NULL, y = ylab) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
  if (!is.null(hlines))
    p <- p + geom_hline(yintercept = hlines, linetype = "dashed",
                        colour = "grey25", linewidth = .4)
  if (log) p <- p + scale_y_log10()
  p
}

p1 <- vln("nFeature_RNA", "Genes detectados por célula",
          "nFeature_RNA", c(min_features, max_features), log = TRUE)
p2 <- vln("nCount_RNA", "UMIs (contagens) por célula",
          "nCount_RNA", min_counts, log = TRUE)
p3 <- vln("percent.mt", "% de contagens mitocondriais",
          "percent.mt", max_mito)
p4 <- vln("percent.ribo", "% de contagens ribossomais", "percent.ribo")

save_plot((p1 | p2) / (p3 | p4), "01_QC",
          "01_QC_antes_da_filtragem.png", width = 12, height = 9)

# Relação UMIs x genes, colorida pela % mitocondrial: o "L" clássico.
# Células boas ficam na diagonal; gotículas vazias no canto inferior
# esquerdo; células moribundas aparecem em vermelho.
set.seed(123)
sub_meta <- meta_all[sample(nrow(meta_all), min(25000, nrow(meta_all))), ]

p_scatter <- ggplot(sub_meta, aes(nCount_RNA, nFeature_RNA, colour = percent.mt)) +
  geom_point(size = .3, alpha = .5) +
  scale_colour_gradientn(colours = c("#2C7BB6", "#FFFFBF", "#D7191C"),
                         limits = c(0, 25), oob = scales::squish,
                         name = "% mito") +
  scale_x_log10() + scale_y_log10() +
  geom_hline(yintercept = c(min_features, max_features),
             linetype = "dashed", colour = "grey30") +
  geom_vline(xintercept = min_counts, linetype = "dashed", colour = "grey30") +
  facet_wrap(~ sample, nrow = 1, labeller = as_labeller(sample_labels)) +
  labs(title = "UMIs x genes detectados por célula",
       subtitle = "linhas tracejadas = limites de filtragem aplicados",
       x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") +
  theme_bw(base_size = 11)

save_plot(p_scatter, "01_QC", "02_QC_umis_vs_genes.png",
          width = 14, height = 4.5)

# Barras antes/depois da filtragem.
il <- bind_rows(import_log)
bar_df <- data.frame(
  sample = factor(rep(il$sample, 2), levels = samples),
  etapa  = rep(c("Antes do QC", "Após o QC"), each = nrow(il)),
  n      = c(il$celulas_brutas, il$celulas_pos_QC)
)

p_bar <- ggplot(bar_df, aes(sample, n, fill = etapa)) +
  geom_col(position = position_dodge(.75), width = .7) +
  geom_text(aes(label = n), position = position_dodge(.75),
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("Antes do QC" = "grey70",
                               "Após o QC" = "#4C78A8")) +
  scale_x_discrete(labels = sample_labels) +
  labs(title = "Células antes e depois do controle de qualidade",
       subtitle = sprintf(
         "filtros: %d <= nFeature <= %d | nCount >= %d | %%mito < %d",
         min_features, max_features, min_counts, max_mito),
       x = NULL, y = "Número de células", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

save_plot(p_bar, "01_QC", "03_QC_celulas_antes_depois.png",
          width = 9, height = 6)

# ------------------------------------------------------------
# SAÍDA: objeto pequeno, já filtrado e subamostrado
# ------------------------------------------------------------
step("Salvando objetos filtrados")
saveRDS(objs_small, file.path(cache_dir, "01_objetos_filtrados.rds"))

message("\nCélulas totais que seguem para a etapa 02: ",
        sum(sapply(objs_small, ncol)))
mem_report("fim da etapa 01")

writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "sessionInfo.txt"))
message("\nETAPA 01 CONCLUÍDA.")
