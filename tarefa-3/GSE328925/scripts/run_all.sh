#!/usr/bin/env bash
# ============================================================
# Executa a análise completa do GSE328925, etapa por etapa.
#
# Cada etapa roda num PROCESSO SEPARADO de R. Isso é intencional:
# ao terminar, o processo morre e devolve TODA a memória ao sistema.
# Numa máquina de 8 GB isso é a diferença entre rodar e travar.
#
# Uso:
#   bash scripts/run_all.sh          # roda tudo
#   bash scripts/run_all.sh 03 04    # roda só as etapas 03 e 04
#
# Pico de memória observado: ~4,5 GB (etapa 02).
# Tempo total: ~15 min num notebook com 8 GB e 8 threads.
# ============================================================
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p results/logs

ETAPAS=("01:01_qc_filtragem.R" "02:02_clusterizacao.R" \
        "03:03_marcadores.R"   "04:04_anotacao_final.R")

ALVOS=("$@")
roda() { [ ${#ALVOS[@]} -eq 0 ] && return 0
         for a in "${ALVOS[@]}"; do [ "$a" = "$1" ] && return 0; done; return 1; }

for e in "${ETAPAS[@]}"; do
  num="${e%%:*}"; script="${e##*:}"
  roda "$num" || continue

  log="results/logs/${num}_$(basename "$script" .R).log"
  echo "▶  Etapa $num  ($script)  →  $log"

  # Roda o R em segundo plano e acompanha o pico de RSS.
  Rscript "scripts/$script" > "$log" 2>&1 &
  pid=$!; pico=0
  while kill -0 $pid 2>/dev/null; do
    rss=$(awk '/VmRSS/{print $2}' "/proc/$pid/status" 2>/dev/null || echo 0)
    [ -n "$rss" ] && [ "$rss" -gt "$pico" ] && pico=$rss
    sleep 5
  done
  wait $pid || { echo "✗  Etapa $num falhou. Veja $log"; exit 1; }

  echo "   PICO_RSS_MB=$((pico/1024))" | tee -a "$log"
  echo "✓  Etapa $num concluída"
done

echo
echo "Pronto. Figuras em figures/ e tabelas em results/."
