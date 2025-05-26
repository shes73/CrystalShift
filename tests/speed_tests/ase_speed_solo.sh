#!/bin/bash

CIF_DIR="./ase_speed_test_solo"
SCRIPT_FILE="ase_speed_solo.py"
LOG_FILE="./logfile.log"

cd "$CIF_DIR" || { echo "Ошибка: не удалось перейти в директорию $CIF_DIR"; exit 1; }

shopt -s nullglob
cif_files=(*.cif)

if [ ${#cif_files[@]} -eq 0 ]; then
    echo "В директории $CIF_DIR нет CIF-файлов."
    exit 1
fi

> "$LOG_FILE"

for CIF_FILE in "${cif_files[@]}"; do
    START_TIME=$(date +%s%N)

    python ../"$SCRIPT_FILE" "$CIF_FILE"

    END_TIME=$(date +%s%N)
    ELAPSED_TIME=$((END_TIME - START_TIME))

    echo "$ELAPSED_TIME, $CIF_FILE" >> "$LOG_FILE"
done

echo "Операция завершена. Результаты записаны в $LOG_FILE"
