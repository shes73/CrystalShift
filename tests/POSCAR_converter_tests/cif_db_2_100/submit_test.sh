#!/bin/bash

program1="./format_converter_test.exe"
program2="./statistics.exe"
logfile="test2.log"
file_directory="."

echo "FILENAME       TEST   RMSE(x)   RMSE(y)   RMSE(z)" > "$logfile"

run_program1() {
    local filename=$1
    local filename_POSCAR="${filename%.cif}_POSCAR"
    local filename_CIF_from_POSCAR="${filename%.cif}_from_POSCAR.cif"

    /usr/bin/expect <<EOF
    spawn $program1
    expect "Enter the path to the CIF file: "
    send "$filename\r"
    expect "Choose the order for writing atoms to POSCAR:\n"
    expect "1. From lightest to heaviest;\n"
    expect "2. From heaviest to lightest;\n"
    expect "3. Custom order (input symbols separated by spaces).\n"
    expect "Enter your choice: "
    send "1\r"
    expect "Enter the path to the POSCAR file: "
    send "$filename_POSCAR\r"
    expect "Enter the path to the new CIF file: "
    send "$filename_CIF_from_POSCAR\r"
    expect eof
    set exit_code [lindex [wait] 3]
    set log [open "$logfile" "a"]
    if {\$exit_code != 0} {
        puts \$log "$filename FAIL"
    }
    close \$log
EOF
}

run_program2() {
    local filename=$1
    local filename_CIF_from_POSCAR="${filename%.cif}_from_POSCAR.cif"

    /usr/bin/expect <<EOF
    spawn $program2
    expect "Enter the path to the first CIF file: "
    send "$filename\r"
    expect "Enter the path to the second CIF file: "
    send "$filename_CIF_from_POSCAR\r"
    expect eof
    set exit_code [lindex [wait] 3]
    set log [open "$logfile" "a"]
    if {\$exit_code == 0} {
        puts \$log "$filename SUCCESS"
    } else {
        puts \$log "$filename FAIL"
    }
    close \$log
EOF

    output=$($program2 <<EOF
$filename
$filename_CIF_from_POSCAR
EOF
)
    last_line=$(echo "$output" | tail -n 1)
    rmse_x=$(echo "$last_line" | awk '{print $(NF-2)}')
    rmse_y=$(echo "$last_line" | awk '{print $(NF-1)}')
    rmse_z=$(echo "$last_line" | awk '{print $NF}')

    # Записываем RMSE в log-файл
    sed -i "\$s/$/ $rmse_x $rmse_y $rmse_z/" "$logfile"
}

for file in $file_directory/*.cif; do
    if [ -f "$file" ]; then
        run_program1 "$file"
        run_program2 "$file"
    fi
done
