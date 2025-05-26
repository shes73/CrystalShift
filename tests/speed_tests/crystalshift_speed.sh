#!/bin/bash

CIF_DIR="./cif_speed_test"
LOG_FILE="./logfile.log"
SPEED_TEST_EXE="./speed_test.exe"

# Проверяем, существует ли исполняемый файл
if [ ! -f "$SPEED_TEST_EXE" ]; then
    echo "Ошибка: Исполняемый файл $SPEED_TEST_EXE не найден."
    exit 1
fi

# Сохраняем абсолютный путь к исполняемому файлу
SPEED_TEST_EXE="$(realpath "$SPEED_TEST_EXE")"

# Переходим в директорию с CIF-файлами
cd "$CIF_DIR" || exit

CIF_FILES=$(ls *.cif 2>/dev/null)
> "$LOG_FILE"

for CIF_FILE in $CIF_FILES; do
    BASE_NAME=$(basename "$CIF_FILE" .cif)
    POSCAR_FILE="${BASE_NAME}_POSCAR"

    START_TIME=$(date +%s%N)
    {
        echo "$CIF_FILE"
        echo "$POSCAR_FILE"
    } | "$SPEED_TEST_EXE"
    END_TIME=$(date +%s%N)

    ELAPSED_TIME=$((END_TIME - START_TIME))
    echo "$ELAPSED_TIME, $CIF_FILE" >> "$LOG_FILE"
done

echo "Операция завершена. Результаты записаны в $LOG_FILE"
