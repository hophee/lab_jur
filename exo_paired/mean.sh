#!/bin/bash

if [ -z "$1" ]; then
    echo "Ошибка: Укажите файл с параметрами."
    echo "Использование: $0 <файл_с_параметрами.txt> <vcf_file.vcf.gz>"
    exit 1
fi

if [ -z "$2" ]; then
    echo "Ошибка: Укажите VCF файл."
    echo "Использование: $0 <файл_с_параметрами.txt> <vcf_file.vcf.gz>"
    exit 1
fi

PARAM_FILE="$1"
VCF_FILE="$2"

PARAMS=$(cat "$PARAM_FILE")

for PARAM in $PARAMS; do
    VALUES=$(bcftools query -f "%${PARAM}\n" "$VCF_FILE")

    MEAN=$(echo "$VALUES" | awk '
    {
        if ($1 != ".") {
            sum += $1;
            count++;
        }
    }
    END {
        if (count > 0) {
            printf "%.2f", sum / count;
        } else {
            print "N/A";  # Если значений нет
        }
    }')

    echo "PARAM: \"$PARAM\""
    echo "MEAN: \"$MEAN\""
    echo ""
done