# Crie um diretório para os FASTQ subamostrados
mkdir -p subamostras

# Loop para processar uma amostra por vez
for srr in SRR37459106 SRR37459107 SRR37459108 SRR37459109 SRR37459110 SRR37459111; do
    echo "========================================"
    echo "Processando $srr..."
    echo "========================================"
    
    # Converte apenas 7.5M spots (= 15M reads paired-end)
    fastq-dump --split-files -X 5000000 -O subamostras "$srr"
    
    # Comprime os FASTQ gerados
    gzip subamostras/${srr}_1.fastq subamostras/${srr}_2.fastq
    
    # Remove o arquivo .sra para liberar espaço
    rm -rf "$srr"
    
    # Mostra o espaço restante
    echo "Espaço após $srr:"
    df -h | tail -1
    echo ""
done
