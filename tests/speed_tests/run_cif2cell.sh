#!/bin/bash

paths=(
  "$HOME/AppData/Local/Programs/Python/Python39-32/Scripts/cif2cell.exe"
  "$HOME/AppData/Local/Programs/Python/Python311/Scripts/cif2cell.exe"
  "$HOME/AppData/Roaming/Python/Python39/Scripts/cif2cell.exe"
)

for path in "${paths[@]}"; do
    if [[ -f "$path" ]]; then
        cif2cell="$path"
        break
    fi
done

if [[ -z "$cif2cell" ]]; then
    echo "Ошибка: cif2cell не найден. Убедитесь, что он установлен через pip."
    exit 1
fi

"$cif2cell" "$@"
