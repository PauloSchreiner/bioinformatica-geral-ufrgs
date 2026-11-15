# ============================================================
# 00. CONFIGURAÇÃO COMUM A TODAS AS ETAPAS
# ============================================================
#
# GSE328925 - Longitudinal scRNA-seq do pulmão de camundongo
# durante a RESOLUÇÃO espontânea da fibrose induzida por bleomicina.
#
# Desenho: 4 bibliotecas 10x Genomics 3' (CellRanger 7.2.0, mm10)
#   Control = salina (dia 0)      GSM9694319
#   D21     = bleomicina, dia 21  GSM9694320  (fibrose estabelecida)
#   D42     = bleomicina, dia 42  GSM9694321  (resolução em curso)
#   D63     = bleomicina, dia 63  GSM9694322  (resolução avançada)
#
# Este arquivo é carregado (source) pelas etapas 01, 02 e 03.
# Nenhuma análise acontece aqui: só caminhos, parâmetros e utilitários.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

set.seed(123)

# ------------------------------------------------------------
# Caminhos
# ------------------------------------------------------------
# Os arquivos 10x (barcodes/features/matrix) ocupam ~580 MB e ficam
# FORA do repositório git. Os resultados e figuras ficam DENTRO.
data_dir    <- path.expand("~/Bioinformatica/GSE328925/data")
# Diretório do projeto: resolvido a partir do caminho do próprio script
# (funciona tanto com Rscript quanto com source() no RStudio).
.this_file <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) normalizePath(f) else NULL
}, error = function(e) NULL)

if (!is.null(.this_file)) {
  project_dir <- normalizePath(file.path(dirname(.this_file), ".."))
} else {
  project_dir <- path.expand(
    "~/Cadeiras/6o-sem/bioinfo-geral/bioinformatica-geral-ufrgs/tarefa-3/GSE328925"
  )
}

results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")
cache_dir   <- path.expand("~/Bioinformatica/GSE328925/cache")  # RDS grandes

# Subpastas de figuras, numeradas na ordem da apresentação.
figure_subdirs <- c(
  "01_QC",
  "02_HVG_PCA",
  "03_UMAP",
  "04_Marcadores",
  "05_Anotacao",
  "06_Fibrose"
)

for (d in c(results_dir, cache_dir,
            file.path(figures_dir, figure_subdirs))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------------------------------------
# Amostras
# ------------------------------------------------------------
samples <- c("Control", "D21", "D42", "D63")

sample_labels <- c(
  Control = "Salina (controle)",
  D21     = "Bleomicina D21",
  D42     = "Bleomicina D42",
  D63     = "Bleomicina D63"
)

# Paleta fixa por amostra (usada em todas as figuras).
sample_colors <- c(
  Control = "#4C78A8",
  D21     = "#E45756",
  D42     = "#F58518",
  D63     = "#54A24B"
)

# ------------------------------------------------------------
# PARÂMETROS DA ANÁLISE  (todos justificados no RELATORIO.md)
# ------------------------------------------------------------

# -- Etapa 01: QC -------------------------------------------
min_cells_per_gene <- 3     # gene precisa aparecer em >=3 células
min_features       <- 500   # células com <500 genes: debris / gotícula vazia
max_features       <- 6000  # células com >6000 genes: suspeita de doublet
max_mito           <- 10    # % mitocondrial máxima (célula estressada/morta)
min_counts         <- 1000  # UMIs mínimos por célula

# Subamostragem: nº máximo de células mantidas POR BIBLIOTECA.
# Motivo principal: a máquina tem 8 GB de RAM e o dataset completo
# (72.757 células x 32.286 genes, 160 milhões de valores não-zero)
# exigiria ~5 GB só para as matrizes de contagem + normalizada.
# Efeito colateral desejável: iguala o "peso" das 4 bibliotecas,
# evitando que a maior domine a clusterização.
cells_per_sample <- 5000

# -- Etapa 02: redução de dimensionalidade e clusters --------
n_variable_features <- 2000  # genes altamente variáveis (HVG)
n_pcs               <- 50    # PCs calculadas (para o ElbowPlot)
dims_use            <- 1:20  # PCs efetivamente usadas (ver ElbowPlot)
cluster_resolution  <- 0.5   # resolução principal do algoritmo de Louvain
# Resoluções extras só para mostrar a estabilidade do agrupamento.
extra_resolutions   <- c(0.2, 0.3, 0.8, 1.0)

# -- Etapa 03: marcadores ------------------------------------
marker_min_pct   <- 0.25
marker_logfc     <- 0.25
heatmap_cells_per_cluster <- 100  # subamostragem só para desenhar o heatmap

# ------------------------------------------------------------
# Utilitários
# ------------------------------------------------------------

save_plot <- function(plot, subdir, filename, width = 8, height = 6) {
  path <- file.path(figures_dir, subdir, filename)
  ggsave(path, plot = plot, width = width, height = height,
         dpi = 200, bg = "white", limitsize = FALSE)
  message("  [fig] ", file.path(subdir, filename))
  invisible(path)
}

save_table <- function(df, filename) {
  path <- file.path(results_dir, filename)
  write.csv(df, path, row.names = FALSE)
  message("  [tab] ", filename)
  invisible(path)
}

# DotPlot do Seurat exige genes ÚNICOS quando recebe uma lista de
# painéis: se o mesmo gene aparece em dois painéis, ele cria níveis de
# fator duplicados e falha. Esta função mantém a primeira ocorrência.
dedup_paineis <- function(lst) {
  vistos <- character(0)
  out <- lapply(lst, function(g) {
    g <- setdiff(g, vistos)
    vistos <<- c(vistos, g)
    g
  })
  out[sapply(out, length) > 0]
}

# Paleta por compartimento: cada compartimento tecidual recebe uma FAMÍLIA
# de cor (epitélio = roxo, endotélio = vermelho/laranja, mesênquima = verde,
# mieloide = azul, linfoide = ciano, proliferativas = oliva, técnico = cinza)
# e os tipos dentro dela variam só em luminosidade. Assim, numa figura com
# 25 categorias, dá para ler o compartimento de relance.
familias_compartimento <- list(
  "Epitelio"       = c("#4A1486", "#C994D6"),
  "Endotelio"      = c("#99000D", "#FDAE6B"),
  "Mesenquima"     = c("#00441B", "#A1D99B"),
  "Mieloide"       = c("#08306B", "#9ECAE1"),
  "Proliferativas" = c("#8C8C00", "#D9D96B"),
  "Linfoide"       = c("#00494D", "#7FD4D4"),
  "Tecnico"        = c("#525252", "#BDBDBD")
)

paleta_por_compartimento <- function(anot) {
  anot <- anot[order(anot$ordem), ]
  tipos <- unique(anot$tipo_celular)
  comp  <- anot$compartimento[match(tipos, anot$tipo_celular)]
  cores <- character(0)
  for (cp in unique(comp)) {
    idx <- which(comp == cp)
    anchors <- familias_compartimento[[cp]]
    if (is.null(anchors)) anchors <- c("#333333", "#CCCCCC")
    cores <- c(cores, colorRampPalette(anchors)(max(length(idx), 2))[seq_along(idx)])
  }
  setNames(cores, tipos)
}

# Imprime o uso de memória do processo R (ajuda a monitorar os 8 GB).
mem_report <- function(tag = "") {
  g <- gc(verbose = FALSE, full = TRUE)
  used_mb <- sum(g[, "used"] * c(56, 8)) / 1e6  # aproximação
  rss <- tryCatch({
    as.numeric(
      sub("[^0-9]*([0-9]+).*", "\\1",
          grep("^VmRSS", readLines("/proc/self/status"), value = TRUE))
    ) / 1024
  }, error = function(e) NA_real_)
  message(sprintf("  [mem] %-28s R=%6.0f MB | RSS=%6.0f MB",
                  tag, used_mb, rss))
  invisible(rss)
}

# Marca o início/fim de uma etapa no log.
step <- function(msg) {
  message("\n", strrep("=", 62))
  message("  ", msg, "   (", format(Sys.time(), "%H:%M:%S"), ")")
  message(strrep("=", 62))
}

# Evita que o pacote future tente paralelizar e duplique a memória.
options(future.globals.maxSize = 2 * 1024^3)
