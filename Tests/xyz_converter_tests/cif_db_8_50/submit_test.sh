#!/bin/bash

program1="./format_converter_xyz_test.exe"
program2="./statistics.exe"
logfile="test8.log"
file_directory="."

echo "FILENAME       TEST   RMSE(x)   RMSE(y)   RMSE(z)" > "$logfile"

run_program1() {
    local filename=$1
    local filename_xyz="${filename%.cif}.xyz"
    local filename_CIF_from_xyz="${filename%.cif}_from_xyz.cif"

    /usr/bin/expect <<EOF
    spawn $program1
    expect "Enter the path to the CIF file: "
    send "$filename\r"
    expect "Choose the format for writing xyz file:\n"
    expect "1. Standard xyz format;\n"
    expect "2. Extended xyz format (with lattice parameters in the comment line);\n"
    expect "Enter your choice: "
    send "2\r"
    expect "Enter the path to the xyz file: "
    send "$filename_xyz\r"
    expect "Enter the path to the new CIF file: "
    send "$filename_CIF_from_xyz\r"
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
    local filename_CIF_from_xyz="${filename%.cif}_from_xyz.cif"

    /usr/bin/expect <<EOF
    spawn $program2
    expect "Enter the path to the first CIF file: "
    send "$filename\r"
    expect "Enter the path to the second CIF file: "
    send "$filename_CIF_from_xyz\r"
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
$filename_CIF_from_xyz
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
