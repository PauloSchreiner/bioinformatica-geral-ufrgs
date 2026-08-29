# Tarefa 1
## Do sequenciamento ao alinhamento

### 1. Obter genoma de referência 
Escolhi trabalhar com o genoma de Saccharomyces cerevisiae. Acessei o NCBI Genomes via https://www.ncbi.nlm.nih.gov/datasets/genome/ e pesquisei Saccharomyces cerevisiae. Escolhi a montagem 	R64 da cepa S288C por ser a montagem de referência. Identifiquei o código do RefSeq dessa montagem: GCF_000146045

Usei a ferramenta ncbi-datasets-cli para poder baixar o genoma apenas com o código RefSeq:

```bash
# preparei um ambiente no mamba
mamba create -n bioinfo-geral
mamba activate bioinfo-geral
mamba install -c bioconda ncbi-datasets-cli -y

# usei a ferramenta para baixar o genoma
datasets download genome accession GCF_000146045.2 --include genome --filename sc_ref_genome.zip
```

### 2. Obter reads do SRA 
O genoma que escolhi pegar do SRA foi o genoma de "Wyeast 3068, a Weihenstephan Weizen tetraploid S. cerevisiae ale strain frequently used to brew Hefeweizen style beers", obtido fazendo a busca por "Saccharomyces cerevisiae" AND ("serial repitching" OR "repitching") AND "Illumina" no SRA. Com isso, obti a run ID SRR25850990.

Decidi baixar um subconjunto de reads usando fastq-dump, que descobri ser uma ferramenta da suíte sra-tools. Usei a flag -X para limitar o número de fragmentos para 50000, a fim de obter a cobertura média de 1.2X (considerando que o genoma de S. cerevisiae tem 12MB e que as reads terão ~150pb). 

```bash
mamba install -c bioconda -c conda-forge sra-tools

fastq-dump -X 50000 -F --outdir reads --split-files SRR25850990
# --split-files faz com que as reads F e R venham em arquivos separados, facilitando o alinhamento
```

### 3. Indexar o genoma de referência

**Observação: Li na documentação do bioconda que é necessário rodar os comandos abaixo para sempre priorizar o canal conda-forge e depois o canal bioconda ao instalar ferramentas:**

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# assim, nem é mais necessário explicitar os canais -c bioconda nas instalações

mamba install bwa samtools -y
```

**Observação: Vi que é necessário descompactar o genoma de referência. Quando baixei, veio em formato .zip. Descompactei com unzip, movi o .fna para a pasta certa, /ref, e renomeei para ref_genome.fna**

```bash
# indexei o genoma de referencia:

bwa index sc_ref_genome.fna 

# foram gerados os arquivos:
ls ref/
ref_genome.fna      ref_genome.fna.ann  ref_genome.fna.pac
ref_genome.fna.amb  ref_genome.fna.bwt  ref_genome.fna.sa

```



### 4. Alinhar as reads contra o genoma

```bash
bwa mem -t 4 ref/ref_genome.fna reads/SRR25850990_1.fastq  reads/SRR25850990_2.fastq > aligned.sam
```

### 5. Gerar e indexar o BAM

```bash
cat aligned.sam | samtools sort > aligned.sorted.bam

samtools index aligned.sorted.bam

ls
# aligned.sam         aligned.sorted.bam.bai  reads
# aligned.sorted.bam  README.md               ref

```

## 6. Visualização no IGV

Constatei que a cobertura estava baixa. Logo, resolvi rodar o fastq-dump de novo com -X 200000 e repeti os comandos acima

```bash
fastq-dump -X 200000 -F --outdir reads --split-files SRR25850990

```
