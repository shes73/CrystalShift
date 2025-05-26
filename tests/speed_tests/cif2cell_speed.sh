#!/bin/bash

CIF_DIR="./cif2cell_speed_test"
LOG_FILE="./logfile.log"

# Проверяем, что команда cif2cell доступна
if ! command -v cif2cell &> /dev/null; then
    echo "Ошибка: cif2cell не найден. Убедитесь, что он установлен и доступен в PATH."
    exit 1
fi

# Переходим в директорию с CIF-файлами
cd "$CIF_DIR" || { echo "Ошибка: не удалось перейти в директорию $CIF_DIR"; exit 1; }

# Получаем список CIF-файлов
CIF_FILES=$(ls *.cif 2>/dev/null)
if [ -z "$CIF_FILES" ]; then
    echo "В директории $CIF_DIR нет CIF-файлов."
    exit 1
fi

> "$LOG_FILE"

for CIF_FILE in $CIF_FILES; do
    BASE_NAME=$(basename "$CIF_FILE" .cif)
    POSCAR_FILE="${BASE_NAME}_POSCAR"

    START_TIME=$(date +%s%N)

    # Выполняем команду cif2cell с нужными параметрами
    cif2cell -f "$CIF_FILE" -p VASP -o "$POSCAR_FILE"

    END_TIME=$(date +%s%N)
    ELAPSED_TIME=$((END_TIME - START_TIME))

    echo "$ELAPSED_TIME, $CIF_FILE" >> "$LOG_FILE"
done

echo "Операция завершена. Результаты записаны в $LOG_FILE"
