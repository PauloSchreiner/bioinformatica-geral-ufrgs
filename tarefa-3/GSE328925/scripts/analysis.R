# GSE328925 — Análise de scRNA-seq
# Pipeline exploratório com Seurat v5
#
# Organismo: Mus musculus
# Tecido: pulmão
# Condições: controle (salina), bleomicina D21, D42 e D63
#
# Objetivos:
# 1. Importar as matrizes de expressão do GEO.
# 2. Avaliar e filtrar células de baixa qualidade.
# 3. Normalizar os dados e identificar genes variáveis.
# 4. Reduzir a dimensionalidade e identificar clusters.
# 5. Identificar genes marcadores dos clusters.
# 6. Explorar a identidade celular e padrões relacionados à fibrose.
#
# IMPORTANTE:
# - Este é um estudo exploratório.
# - Há uma biblioteca GEO por condição/timepoint.
# - O número de animais que contribuiu para cada biblioteca não está
#   esclarecido nos metadados disponíveis.
# - Portanto, comparações entre bibliotecas serão descritivas.
# - Não trataremos células individuais como réplicas biológicas.
#
# Recomendação: executar seção por seção e examinar os gráficos de QC
# antes de decidir os filtros definitivos.


# ============================================================
# 00. CONFIGURAÇÕES INICIAIS
# ============================================================

# Instalar os pacotes necessários, caso ainda não estejam instalados.
# Executar esta instalação apenas uma vez.

# install.packages("Seurat")
# install.packages("ggplot2")
# install.packages("patchwork")
# install.packages("dplyr")

# Carregar bibliotecas.
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Verificar a versão do Seurat.
packageVersion("Seurat")

# Este pipeline foi escrito para Seurat v5.
if (packageVersion("Seurat") < "5.0.0") {
  stop("Este script requer Seurat v5 ou superior.")
}

# Definir semente para tornar operações aleatórias reproduzíveis.
set.seed(123)

# Diretório principal do projeto.
# Ajustar este caminho caso tenhas criado a pasta em outro local.
project_dir <- path.expand("~/Bioinformatica/GSE328925")

# Criar diretórios de trabalho.
dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(project_dir, "data"),
           recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(project_dir, "results"),
           recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(project_dir, "figures"),
           recursive = TRUE, showWarnings = FALSE)

# Subdiretórios para organizar as figuras.
figure_dirs <- c(
  "QC",
  "PCA",
  "UMAP",
  "Markers",
  "Fibrosis"
)

for (d in figure_dirs) {
  dir.create(
    file.path(project_dir, "figures", d),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# Função auxiliar para salvar gráficos.
save_plot <- function(plot, filename, width = 8, height = 6) {
  
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}


# ============================================================
# 01. DOWNLOAD DOS DADOS DO GEO
# ============================================================

# O dataset GSE328925 contém quatro bibliotecas.
#
# Cada biblioteca possui três arquivos:
# - barcodes.tsv.gz: identificadores das células.
# - features.tsv.gz: identificadores dos genes/features.
# - matrix.mtx.gz: matriz esparsa de contagens.
#
# O formato é compatível com a função Read10X() do Seurat.

samples <- c(
  Control = "GSM9694319",
  D21     = "GSM9694320",
  D42     = "GSM9694321",
  D63     = "GSM9694322"
)

# Sufixos usados nos nomes originais dos arquivos.
suffixes <- c(
  Control = "control",
  D21     = "d21",
  D42     = "d42",
  D63     = "d63"
)

# Função para testar a integridade de um arquivo gzip.
check_gzip <- function(path) {
  
  if (!file.exists(path)) {
    return(FALSE)
  }
  
  status <- suppressWarnings(
    system2("gzip", args = c("-t", shQuote(path)))
  )
  
  return(status == 0)
}

# Baixar os arquivos de cada amostra.
for (sample_name in names(samples)) {
  
  gsm <- samples[[sample_name]]
  suffix <- suffixes[[sample_name]]
  
  sample_dir <- file.path(
    project_dir,
    "data",
    sample_name
  )
  
  dir.create(
    sample_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  base_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/",
    "GSM9694nnn/",
    gsm,
    "/suppl/"
  )
  
  original_files <- c(
    barcodes = paste0(gsm, "_", suffix, "_barcodes.tsv.gz"),
    features = paste0(gsm, "_", suffix, "_features.tsv.gz"),
    matrix   = paste0(gsm, "_", suffix, "_matrix.mtx.gz")
  )
  
  for (file_type in names(original_files)) {
    
    file_url <- paste0(
      base_url,
      original_files[[file_type]]
    )
    
    # Read10X espera estes nomes de arquivos.
    destination <- file.path(
      sample_dir,
      paste0(
        file_type,
        ifelse(file_type == "matrix", ".mtx.gz", ".tsv.gz")
      )
    )
    
    # Não baixar novamente um arquivo que já está íntegro.
    if (check_gzip(destination)) {
      message("Arquivo já disponível: ", destination)
      next
    }
    
    message("Baixando: ", original_files[[file_type]])
    
    download.file(
      url = file_url,
      destfile = destination,
      mode = "wb",
      method = "libcurl"
    )
    
    # Verificar se o download foi concluído corretamente.
    if (!check_gzip(destination)) {
      stop(
        "Falha ao baixar ou verificar o arquivo: ",
        destination
      )
    }
  }
}

# Conferir os arquivos baixados.
for (sample_name in names(samples)) {
  
  sample_dir <- file.path(
    project_dir,
    "data",
    sample_name
  )
  
  cat("\nAmostra:", sample_name, "\n")
  print(file.info(
    list.files(sample_dir, full.names = TRUE)
  )[, "size", drop = FALSE])
}


# ============================================================
# 02. IMPORTAÇÃO DAS MATRIZES
# ============================================================

# Read10X() importa os três arquivos de cada biblioteca.
#
# A matriz contém:
# - linhas: genes;
# - colunas: células;
# - valores: contagens de RNA/UMIs.
#
# Como a matriz é esparsa, o Seurat não precisa armazenar
# explicitamente todos os valores iguais a zero.
#
# Criamos um objeto separado por biblioteca para preservar
# a identificação da origem de cada célula.

seurat_list <- list()

for (sample_name in names(samples)) {
  
  message("Importando amostra: ", sample_name)
  
  sample_dir <- file.path(
    project_dir,
    "data",
    sample_name
  )
  
  # Importar a matriz.
  counts <- Read10X(data.dir = sample_dir)
  
  # Criar objeto Seurat.
  #
  # min.cells = 3:
  # descarta genes detectados em menos de três células.
  #
  # min.features = 0:
  # não aplicar ainda um filtro de genes por célula.
  # Queremos avaliar a distribuição de QC antes de decidir
  # os limites definitivos.
  
  seu <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 3,
    min.features = 0
  )
  
  # Adicionar metadado identificando a biblioteca.
  seu$sample <- sample_name
  
  # Guardar o objeto.
  seurat_list[[sample_name]] <- seu
  
  # Liberar a matriz temporária da memória.
  rm(counts)
  gc()
}

# Conferir dimensões de cada objeto.
for (sample_name in names(seurat_list)) {
  
  cat("\n", sample_name, "\n")
  
  print(dim(seurat_list[[sample_name]]))
}

# As dimensões são:
# número de genes x número de células.
#
# Se alguma amostra tiver dimensão inesperada ou zero células,
# interromper a análise e verificar os arquivos importados.


# ============================================================
# 03. MÉTRICAS DE CONTROLE DE QUALIDADE (QC)
# ============================================================

# Vamos avaliar três métricas principais:
#
# nFeature_RNA:
# número de genes detectados por célula.
#
# nCount_RNA:
# número total de contagens/UMIs por célula.
#
# percent.mt:
# porcentagem das contagens atribuídas a genes mitocondriais.
#
# Interpretação geral:
#
# - Poucos genes: pode indicar gotículas vazias ou células
#   de baixa qualidade.
# - Muitas contagens/genes: pode indicar doublets,
#   embora também possa representar células grandes ou
#   biologicamente muito ativas.
# - Alta porcentagem mitocondrial: pode indicar células
#   danificadas ou em processo de morte.
#
# Em camundongos, os símbolos de genes mitocondriais
# geralmente começam com "mt-".

for (sample_name in names(seurat_list)) {
  
  seu <- seurat_list[[sample_name]]
  
  seu[["percent.mt"]] <- PercentageFeatureSet(
    seu,
    pattern = "^mt-"
  )
  
  seurat_list[[sample_name]] <- seu
}

# Combinar os objetos em um único objeto.
#
# add.cell.ids adiciona um prefixo aos nomes das células,
# evitando colisões entre bibliotecas diferentes.
#
# Não estamos realizando integração entre amostras.
# Apenas reunimos as células para uma análise conjunta.

combined <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "GSE328925"
)

# Definir o ensaio principal.
DefaultAssay(combined) <- "RNA"

# Conferir número de células por biblioteca.
table(combined$sample)

# Salvar os números antes da filtragem.
cells_before <- table(combined$sample)

# Visualizar distribuições de QC por amostra.

p_qc_features <- VlnPlot(
  combined,
  features = "nFeature_RNA",
  group.by = "sample",
  pt.size = 0
) +
  ggtitle("Genes detectados por célula")

p_qc_counts <- VlnPlot(
  combined,
  features = "nCount_RNA",
  group.by = "sample",
  pt.size = 0
) +
  ggtitle("Contagens por célula")

p_qc_mito <- VlnPlot(
  combined,
  features = "percent.mt",
  group.by = "sample",
  pt.size = 0
) +
  ggtitle("Porcentagem de contagens mitocondriais")

p_qc <- p_qc_features / p_qc_counts / p_qc_mito

print(p_qc)

save_plot(
  p_qc,
  file.path(
    project_dir,
    "figures",
    "QC",
    "QC_before_filtering.png"
  ),
  width = 10,
  height = 12
)

# Gráficos de dispersão para investigar relações entre métricas.

p_scatter1 <- FeatureScatter(
  combined,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
) +
  ggtitle("Contagens versus genes detectados")

p_scatter2 <- FeatureScatter(
  combined,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
) +
  ggtitle("Contagens versus porcentagem mitocondrial")

print(p_scatter1)
print(p_scatter2)

save_plot(
  p_scatter1,
  file.path(
    project_dir,
    "figures",
    "QC",
    "QC_counts_vs_features.png"
  )
)

save_plot(
  p_scatter2,
  file.path(
    project_dir,
    "figures",
    "QC",
    "QC_counts_vs_mito.png"
  )
)

# Resumo numérico das métricas de QC por biblioteca.
#
# Estes valores ajudam a escolher limites de filtragem
# de maneira fundamentada, em vez de usar cortes arbitrários.

qc_by_sample <- combined@meta.data %>%
  group_by(sample) %>%
  summarise(
    n_cells = n(),
    
    features_median = median(nFeature_RNA),
    features_q05 = quantile(nFeature_RNA, 0.05),
    features_q95 = quantile(nFeature_RNA, 0.95),
    
    counts_median = median(nCount_RNA),
    counts_q05 = quantile(nCount_RNA, 0.05),
    counts_q95 = quantile(nCount_RNA, 0.95),
    
    mito_median = median(percent.mt),
    mito_q95 = quantile(percent.mt, 0.95),
    
    .groups = "drop"
  )

print(qc_by_sample)

write.csv(
  qc_by_sample,
  file.path(
    project_dir,
    "results",
    "QC_metrics_by_sample.csv"
  ),
  row.names = FALSE
)

# IMPORTANTE:
# Antes de prosseguir, examinar os gráficos e a tabela.
#
# Os limites abaixo são provisórios.
# Se uma biblioteca apresentar uma distribuição muito diferente,
# investigar antes de aplicar os mesmos limites a todas.


# ============================================================
# 04. FILTRAGEM DE CÉLULAS
# ============================================================

# Definir os limites de filtragem.
#
# Estes valores são pontos de partida, não regras universais.
#
# min_features:
# limite inferior de genes detectados.
#
# max_features:
# limite superior de genes detectados.
# Valores excepcionalmente elevados podem sugerir doublets.
#
# max_mito:
# limite máximo de porcentagem mitocondrial.
#
# Rever estes parâmetros após examinar os gráficos de QC.

min_features <- 200
max_features <- 6000
max_mito <- 15

# Aplicar filtros.
combined <- subset(
  combined,
  subset =
    nFeature_RNA >= min_features &
    nFeature_RNA <= max_features &
    percent.mt < max_mito
)

# Conferir quantas células restaram.
cells_after <- table(
  factor(
    combined$sample,
    levels = names(cells_before)
  )
)

# Comparar número de células antes e depois.
qc_summary <- data.frame(
  sample = names(cells_before),
  before = as.integer(cells_before),
  after = as.integer(cells_after)
)

qc_summary$removed <- (
  qc_summary$before - qc_summary$after
)

qc_summary$percent_retained <- round(
  100 * qc_summary$after / qc_summary$before,
  1
)

print(qc_summary)

write.csv(
  qc_summary,
  file.path(
    project_dir,
    "results",
    "QC_summary.csv"
  ),
  row.names = FALSE
)

# Visualizar novamente as métricas após filtragem.
p_qc_after <- VlnPlot(
  combined,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "sample",
  pt.size = 0,
  ncol = 3
)

print(p_qc_after)

save_plot(
  p_qc_after,
  file.path(
    project_dir,
    "figures",
    "QC",
    "QC_after_filtering.png"
  ),
  width = 14,
  height = 5
)

# Salvar checkpoint após QC.
saveRDS(
  combined,
  file.path(
    project_dir,
    "results",
    "GSE328925_QC_filtered.rds"
  )
)

# OBSERVAÇÃO:
# Não estamos realizando uma identificação formal de doublets.
# Os limites de nFeature_RNA ajudam a remover algumas células
# suspeitas, mas não substituem ferramentas específicas
# de detecção de doublets.
#
# Também não devemos interpretar as proporções de células
# recuperadas como abundâncias exatas no tecido original.


# ============================================================
# 05. PREPARAÇÃO DAS MATRIZES PARA ANÁLISE CONJUNTA
# ============================================================

# No Seurat v5, a combinação de objetos pode manter camadas
# separadas para cada biblioteca.
#
# Como não faremos integração entre amostras, vamos unir
# as camadas do ensaio RNA antes da normalização conjunta.

combined <- JoinLayers(
  combined,
  assay = "RNA"
)

DefaultAssay(combined) <- "RNA"

# Conferir o objeto.
combined

# Salvar checkpoint.
saveRDS(
  combined,
  file.path(
    project_dir,
    "results",
    "GSE328925_QC_joined.rds"
  )
)


# ============================================================
# 06. NORMALIZAÇÃO
# ============================================================

# A matriz original contém contagens brutas de RNA.
#
# Células diferentes podem ter profundidades de sequenciamento
# distintas, isto é, diferentes números totais de UMIs.
#
# LogNormalize:
# 1. Divide as contagens de cada gene pelo total de contagens
#    da respectiva célula.
# 2. Multiplica pelo fator de escala (10.000).
# 3. Aplica transformação logarítmica.
#
# Isso permite comparar perfis de expressão entre células
# com diferentes profundidades de sequenciamento.

combined <- NormalizeData(
  combined,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# Conferir se a normalização foi realizada.
DefaultAssay(combined)


# ============================================================
# 07. IDENTIFICAÇÃO DE GENES ALTAMENTE VARIÁVEIS
# ============================================================

# Nem todos os genes são igualmente informativos para distinguir
# células.
#
# FindVariableFeatures identifica genes com variação de expressão
# relativamente elevada entre as células.
#
# Esses genes serão utilizados nas etapas de PCA e clustering.
#
# A seleção "vst" é uma abordagem padrão do Seurat.

combined <- FindVariableFeatures(
  combined,
  selection.method = "vst",
  nfeatures = 2000
)

# Visualizar genes variáveis.
p_variable <- VariableFeaturePlot(combined)

top10 <- head(
  VariableFeatures(combined),
  10
)

p_variable_labeled <- LabelPoints(
  plot = p_variable,
  points = top10,
  repel = TRUE
)

print(p_variable_labeled)

save_plot(
  p_variable_labeled,
  file.path(
    project_dir,
    "figures",
    "PCA",
    "Variable_features.png"
  ),
  width = 10,
  height = 7
)

# Registrar os genes variáveis.
write.csv(
  data.frame(
    gene = VariableFeatures(combined)
  ),
  file.path(
    project_dir,
    "results",
    "Variable_features.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 08. ESCALONAMENTO DOS DADOS
# ============================================================

# ScaleData centraliza e escala a expressão gênica.
#
# Para cada gene:
# - centraliza os valores em torno da média;
# - divide pela variabilidade (desvio-padrão).
#
# Isso evita que genes com escalas de expressão muito diferentes
# dominem a PCA.
#
# Para economizar memória, escalamos somente os 2.000 genes
# variáveis, e não todos os genes do transcriptoma.

combined <- ScaleData(
  combined,
  features = VariableFeatures(combined)
)


# ============================================================
# 09. ANÁLISE DE COMPONENTES PRINCIPAIS (PCA)
# ============================================================

# PCA (Principal Component Analysis) é uma técnica de redução
# de dimensionalidade.
#
# Em vez de representar cada célula por milhares de genes,
# a PCA cria componentes principais que resumem padrões
# de variação da expressão gênica.
#
# Cada componente representa uma combinação de genes.
#
# Vamos calcular 30 componentes principais.

combined <- RunPCA(
  combined,
  features = VariableFeatures(combined),
  npcs = 30,
  seed.use = 123
)

# Visualizar genes que mais contribuem para os componentes 1 e 2.
p_pca_loadings <- VizDimLoadings(
  combined,
  dims = 1:2,
  reduction = "pca"
)

print(p_pca_loadings)

save_plot(
  p_pca_loadings,
  file.path(
    project_dir,
    "figures",
    "PCA",
    "PCA_loadings.png"
  ),
  width = 10,
  height = 7
)

# Visualizar células no espaço das duas primeiras PCs.
p_pca <- DimPlot(
  combined,
  reduction = "pca",
  group.by = "sample"
) +
  ggtitle("PCA colorida por biblioteca")

print(p_pca)

save_plot(
  p_pca,
  file.path(
    project_dir,
    "figures",
    "PCA",
    "PCA_by_sample.png"
  )
)

# ElbowPlot:
# ajuda a avaliar quantos componentes principais são úteis.
#
# Não existe um número universal de PCs adequado para todos
# os datasets. A escolha deve considerar a queda da variância
# explicada e a estrutura observada nos dados.

p_elbow <- ElbowPlot(
  combined,
  ndims = 30
)

print(p_elbow)

save_plot(
  p_elbow,
  file.path(
    project_dir,
    "figures",
    "PCA",
    "PCA_elbow.png"
  )
)

# Escolha inicial de 20 componentes.
#
# Rever após observar o ElbowPlot.
dims_use <- 1:20


# ============================================================
# 10. CONSTRUÇÃO DA REDE DE VIZINHOS E CLUSTERING
# ============================================================

# FindNeighbors:
# identifica células com perfis semelhantes no espaço
# dos componentes principais.
#
# FindClusters:
# utiliza a rede de vizinhos para identificar comunidades
# de células com perfis de expressão semelhantes.
#
# A resolução controla a granularidade do agrupamento.
# Valores maiores geralmente produzem mais clusters.
#
# Começaremos com resolução 0.5.
# Posteriormente poderemos avaliar outras resoluções.

combined <- FindNeighbors(
  combined,
  dims = dims_use
)

combined <- FindClusters(
  combined,
  resolution = 0.5,
  random.seed = 123
)

# Conferir quantidade de células em cada cluster.
table(Idents(combined))

# Visualizar clusters no espaço da PCA.
p_pca_clusters <- DimPlot(
  combined,
  reduction = "pca",
  group.by = "seurat_clusters",
  label = TRUE
) +
  ggtitle("Clusters celulares na PCA")

print(p_pca_clusters)

save_plot(
  p_pca_clusters,
  file.path(
    project_dir,
    "figures",
    "PCA",
    "PCA_clusters.png"
  )
)


# ============================================================
# 11. UMAP
# ============================================================

# UMAP (Uniform Manifold Approximation and Projection) é uma
# técnica de redução de dimensionalidade usada para representar
# células em duas dimensões.
#
# Células com perfis semelhantes tendem a aparecer próximas.
#
# IMPORTANTE:
# - Distâncias no UMAP não equivalem diretamente a distâncias
#   biológicas ou diferenças quantitativas de expressão.
# - A separação visual de grupos não demonstra causalidade.
# - UMAP é uma visualização, não um teste estatístico.

set.seed(123)

combined <- RunUMAP(
  combined,
  dims = dims_use,
  seed.use = 123
)

# UMAP colorido por cluster.
p_umap_clusters <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Clusters celulares")

# UMAP colorido por biblioteca.
p_umap_samples <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "sample"
) +
  ggtitle("Distribuição das células por biblioteca")

# Visualizar lado a lado.
p_umap <- p_umap_clusters | p_umap_samples

print(p_umap)

save_plot(
  p_umap_clusters,
  file.path(
    project_dir,
    "figures",
    "UMAP",
    "UMAP_clusters.png"
  ),
  width = 9,
  height = 7
)

save_plot(
  p_umap_samples,
  file.path(
    project_dir,
    "figures",
    "UMAP",
    "UMAP_samples.png"
  ),
  width = 9,
  height = 7
)

# Visualização adicional:
# dividir o UMAP por biblioteca.
#
# Isso permite observar a distribuição dos clusters
# em cada condição/timepoint.

p_umap_split <- DimPlot(
  combined,
  reduction = "umap",
  split.by = "sample",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
)

print(p_umap_split)

save_plot(
  p_umap_split,
  file.path(
    project_dir,
    "figures",
    "UMAP",
    "UMAP_split_by_sample.png"
  ),
  width = 16,
  height = 7
)


# ============================================================
# 12. IDENTIFICAÇÃO DE GENES MARCADORES DOS CLUSTERS
# ============================================================

# Agora queremos descobrir quais genes caracterizam cada cluster.
#
# FindAllMarkers compara a expressão de cada cluster com
# as demais células.
#
# only.pos = TRUE:
# mantém genes com expressão relativamente maior no cluster.
#
# min.pct = 0.25:
# exige que o gene seja detectado em pelo menos 25% das células
# de pelo menos um dos grupos comparados.
#
# logfc.threshold = 0.25:
# aplica um limite mínimo de diferença de expressão.
#
# ATENÇÃO:
# Os p-valores dessa análise não representam evidência
# de diferenças entre condições experimentais.
# São resultados exploratórios de expressão diferencial
# entre células agrupadas por cluster.
#
# Também não substituem uma análise com réplicas biológicas.

# Garantir que as identidades atuais sejam os clusters.
Idents(combined) <- "seurat_clusters"

# Encontrar marcadores.
markers <- FindAllMarkers(
  combined,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Ordenar marcadores por cluster e significância.
markers <- markers %>%
  arrange(cluster, p_val_adj, desc(avg_log2FC))

# Conferir os resultados.
head(markers)

# Salvar tabela completa.
write.csv(
  markers,
  file.path(
    project_dir,
    "results",
    "cluster_markers.csv"
  ),
  row.names = FALSE
)

# Salvar os 10 principais marcadores de cada cluster.
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_min(
    order_by = p_val_adj,
    n = 10,
    with_ties = FALSE
  )

write.csv(
  top_markers,
  file.path(
    project_dir,
    "results",
    "top10_markers_per_cluster.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. VISUALIZAÇÃO DOS MARCADORES
# ============================================================

# Heatmap:
# permite observar padrões de expressão dos genes marcadores
# entre os clusters.
#
# Usaremos os cinco principais marcadores de cada cluster.
#
# O heatmap utiliza os valores escalonados do objeto.

top5_markers <- markers %>%
  group_by(cluster) %>%
  slice_min(
    order_by = p_val_adj,
    n = 5,
    with_ties = FALSE
  ) %>%
  pull(gene) %>%
  unique()

# Conferir quais genes estão disponíveis.
top5_markers <- intersect(
  top5_markers,
  rownames(combined)
)

if (length(top5_markers) > 0) {
  
  p_heatmap <- DoHeatmap(
    combined,
    features = top5_markers,
    group.by = "seurat_clusters"
  ) +
    NoLegend()
  
  print(p_heatmap)
  
  save_plot(
    p_heatmap,
    file.path(
      project_dir,
      "figures",
      "Markers",
      "Cluster_marker_heatmap.png"
    ),
    width = 12,
    height = 10
  )
}


# ============================================================
# 14. PAINEL DE MARCADORES PARA ANOTAÇÃO CELULAR
# ============================================================

# A identificação de tipos celulares depende da expressão
# de conjuntos de genes marcadores conhecidos.
#
# Não devemos atribuir identidade a um cluster com base
# em apenas um gene.
#
# Os genes abaixo são candidatos para investigar populações
# celulares comuns no pulmão de camundongos.
#
# A presença de um marcador não garante, sozinha,
# a identidade de um cluster.
#
# A anotação deverá considerar:
# - vários marcadores concordantes;
# - marcadores de populações alternativas;
# - genes marcadores encontrados em FindAllMarkers;
# - coerência biológica do conjunto.

marker_panels <- list(
  
  Fibroblastos = c(
    "Col1a1", "Col1a2", "Dcn", "Pdgfra",
    "Lum", "Col3a1"
  ),
  
  Miofibroblastos = c(
    "Acta2", "Tagln", "Postn",
    "Cthrc1", "Tnc"
  ),
  
  Macrofagos = c(
    "Adgre1", "C1qa", "C1qb",
    "Lyz2", "Csf1r"
  ),
  
  Monocitos = c(
    "Ccr2", "Ly6c2", "S100a8",
    "Lyz2"
  ),
  
  Neutrofilos = c(
    "S100a8", "S100a9", "Ly6g",
    "Mpo", "Retnlg"
  ),
  
  Celulas_T = c(
    "Cd3d", "Cd3e", "Trac",
    "Cd4", "Cd8a"
  ),
  
  Celulas_B = c(
    "Cd79a", "Ms4a1", "Cd74",
    "Cd37"
  ),
  
  Celulas_NK = c(
    "Nkg7", "Klrd1", "Klrk1",
    "Prf1"
  ),
  
  Epiteliais_AT1 = c(
    "Ager", "Pdpn", "Aqp5",
    "Hopx"
  ),
  
  Epiteliais_AT2 = c(
    "Sftpc", "Sftpb", "Abca3",
    "Sftpa1"
  ),
  
  Celulas_Club = c(
    "Scgb1a1", "Cyp2f2",
    "Scgb3a1"
  ),
  
  Celulas_ciliadas = c(
    "Foxj1", "Tubb4b",
    "Pifo"
  ),
  
  Epiteliais_basal = c(
    "Krt5", "Krt14",
    "Trp63"
  ),
  
  Endoteliais = c(
    "Pecam1", "Cdh5", "Kdr",
    "Emcn"
  ),
  
  Pericitos = c(
    "Rgs5", "Pdgfrb",
    "Cspg4"
  ),
  
  Mesoteliais = c(
    "Msln", "Upk3b",
    "Wt1"
  )
)

# Manter somente genes presentes na matriz.
all_panel_genes <- unique(
  unlist(marker_panels)
)

present_panel_genes <- intersect(
  all_panel_genes,
  rownames(combined)
)

missing_panel_genes <- setdiff(
  all_panel_genes,
  present_panel_genes
)

cat("\nGenes do painel ausentes na matriz:\n")
print(missing_panel_genes)

# Salvar a lista de genes ausentes para documentação.
write.csv(
  data.frame(gene = missing_panel_genes),
  file.path(
    project_dir,
    "results",
    "marker_genes_not_found.csv"
  ),
  row.names = FALSE
)

# DotPlot:
# - tamanho do ponto: proporção de células que expressam o gene;
# - cor: expressão média do gene no grupo.
#
# Esse gráfico é especialmente útil para comparar
# vários marcadores entre clusters.

if (length(present_panel_genes) > 0) {
  
  p_dotplot <- DotPlot(
    combined,
    features = present_panel_genes,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis() +
    ggtitle("Painel de marcadores celulares")
  
  print(p_dotplot)
  
  save_plot(
    p_dotplot,
    file.path(
      project_dir,
      "figures",
      "Markers",
      "Cell_type_marker_dotplot.png"
    ),
    width = 16,
    height = 9
  )
}


# ============================================================
# 15. EXPLORAÇÃO DE GENES RELACIONADOS À FIBROSE
# ============================================================

# A fibrose pulmonar envolve remodelamento da matriz
# extracelular, deposição de colágeno e alterações
# na atividade de fibroblastos e outras populações.
#
# Vamos explorar alguns genes relacionados a esses processos.
#
# ATENÇÃO:
# A expressão de um gene isolado não comprova a presença
# de fibrose nem permite concluir que uma condição causou
# determinada alteração.
#
# A interpretação deve considerar o tipo celular,
# o contexto temporal e as limitações experimentais.

fibrosis_genes <- c(
  "Col1a1",
  "Col1a2",
  "Col3a1",
  "Acta2",
  "Tagln",
  "Fn1",
  "Tgfb1",
  "Postn",
  "Cthrc1",
  "Sftpc",
  "Ager"
)

# Verificar quais genes estão presentes.
fibrosis_genes_present <- intersect(
  fibrosis_genes,
  rownames(combined)
)

fibrosis_genes_missing <- setdiff(
  fibrosis_genes,
  fibrosis_genes_present
)

cat("\nGenes de interesse ausentes:\n")
print(fibrosis_genes_missing)

# Visualizar expressão dos genes no UMAP,
# separando as células por biblioteca.
#
# A escala de cores é mantida entre bibliotecas para facilitar
# a comparação visual dos níveis de expressão.

for (gene in fibrosis_genes_present) {
  
  p <- FeaturePlot(
    combined,
    features = gene,
    reduction = "umap",
    split.by = "sample",
    keep.scale = "all"
  ) +
    plot_annotation(
      title = paste("Expressão de", gene)
    )
  
  print(p)
  
  save_plot(
    p,
    file.path(
      project_dir,
      "figures",
      "Fibrosis",
      paste0("FeaturePlot_", gene, ".png")
    ),
    width = 14,
    height = 7
  )
}

# DotPlot dos genes relacionados à fibrose por cluster.
if (length(fibrosis_genes_present) > 0) {
  
  p_fibrosis_dot <- DotPlot(
    combined,
    features = fibrosis_genes_present,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis() +
    ggtitle("Genes relacionados à fibrose por cluster")
  
  print(p_fibrosis_dot)
  
  save_plot(
    p_fibrosis_dot,
    file.path(
      project_dir,
      "figures",
      "Fibrosis",
      "Fibrosis_genes_by_cluster.png"
    ),
    width = 12,
    height = 7
  )
}


# ============================================================
# 16. DISTRIBUIÇÃO DOS CLUSTERS ENTRE AS BIBLIOTECAS
# ============================================================

# Vamos descrever a proporção de células de cada cluster
# em cada biblioteca.
#
# Isso permite observar padrões como:
# - clusters encontrados em todas as bibliotecas;
# - clusters mais representados em determinado timepoint;
# - populações aparentemente restritas a uma biblioteca.
#
# IMPORTANTE:
# Essas proporções representam as células recuperadas e
# analisadas, não necessariamente a abundância real no tecido.
#
# Além disso, sem réplicas biológicas independentes conhecidas,
# não devemos aplicar testes inferenciais de diferenças
# de composição celular entre condições.

cluster_counts <- as.data.frame(
  table(
    sample = combined$sample,
    cluster = combined$seurat_clusters
  )
)

colnames(cluster_counts) <- c(
  "sample",
  "cluster",
  "n_cells"
)

# Calcular proporção de cada cluster dentro de cada biblioteca.
cluster_proportions <- cluster_counts %>%
  group_by(sample) %>%
  mutate(
    proportion = n_cells / sum(n_cells)
  ) %>%
  ungroup()

# Salvar tabela.
write.csv(
  cluster_proportions,
  file.path(
    project_dir,
    "results",
    "cluster_proportions_by_sample.csv"
  ),
  row.names = FALSE
)

# Gráfico de barras empilhadas.
p_cluster_proportions <- ggplot(
  cluster_proportions,
  aes(
    x = sample,
    y = proportion,
    fill = cluster
  )
) +
  geom_col() +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    title = "Distribuição dos clusters por biblioteca",
    x = "Biblioteca",
    y = "Proporção de células",
    fill = "Cluster"
  ) +
  theme_classic()

print(p_cluster_proportions)

save_plot(
  p_cluster_proportions,
  file.path(
    project_dir,
    "figures",
    "UMAP",
    "Cluster_proportions_by_sample.png"
  ),
  width = 9,
  height = 7
)


# ============================================================
# 17. PREPARAR TABELA PARA ANOTAÇÃO MANUAL DOS CLUSTERS
# ============================================================

# A anotação celular não deve ser automatizada sem inspeção.
#
# Vamos criar uma tabela para preencher após examinar:
# - marcadores identificados por FindAllMarkers;
# - DotPlot de marcadores conhecidos;
# - heatmap;
# - distribuição dos genes no UMAP.
#
# A coluna cell_type deverá ser preenchida com a identidade
# celular mais plausível.
#
# A coluna evidence deverá registrar os genes que sustentam
# a anotação.
#
# Exemplos de anotações possíveis:
# Fibroblasto, Macrófago, Célula AT2, Endotélio etc.

cluster_ids <- levels(
  factor(combined$seurat_clusters)
)

annotation_template <- data.frame(
  cluster = cluster_ids,
  cell_type = rep(NA_character_, length(cluster_ids)),
  evidence = rep(NA_character_, length(cluster_ids)),
  confidence = rep(NA_character_, length(cluster_ids))
)

write.csv(
  annotation_template,
  file.path(
    project_dir,
    "results",
    "cluster_annotation_template.csv"
  ),
  row.names = FALSE
)

# Não preencher automaticamente os tipos celulares.
# A anotação será feita após interpretar os marcadores.


# ============================================================
# 18. SALVAR OBJETO FINAL E METADADOS
# ============================================================

# Salvar objeto Seurat completo.
#
# Ele contém:
# - contagens;
# - dados normalizados;
# - genes variáveis;
# - metadados;
# - PCA;
# - clusters;
# - UMAP.

saveRDS(
  combined,
  file.path(
    project_dir,
    "results",
    "GSE328925_processed.rds"
  )
)

# Exportar metadados das células.
write.csv(
  combined@meta.data,
  file.path(
    project_dir,
    "results",
    "cell_metadata.csv"
  )
)

# Registrar informações da sessão:
# versões do R e dos pacotes utilizados.
capture.output(
  sessionInfo(),
  file = file.path(
    project_dir,
    "results",
    "sessionInfo.txt"
  )
)

# Conferir estrutura final do objeto.
combined

# Conferir células por biblioteca.
table(combined$sample)

# Conferir células por cluster.
table(combined$seurat_clusters)

message("Pipeline concluído.")
