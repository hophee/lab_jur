#!/bin/bash

# 1. Конвертируем каждый EPS в PDF
for eps_file in ECZV_*.eps; do
    if [ -f "$eps_file" ]; then
        echo "Конвертируем $eps_file в PDF..."
        gs -q -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dEPSCrop -sOutputFile="${eps_file%.*}.pdf" "$eps_file"
    else
        echo "Файлов ECZV_*.eps не найдено!"
        exit 1
    fi
done

# 2. Объединяем все PDF в один
if ls ECZV_*.pdf 1> /dev/null 2>&1; then
    echo "Объединяем PDF-файлы в merged.pdf..."
    # Проверяем, установлен ли pdfunite (быстрее)
    if command -v pdfunite &> /dev/null; then
        pdfunite ECZV_*.pdf merged.pdf
    else
        gs -q -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -sOutputFile=merged.pdf ECZV_*.pdf
    fi
    echo "Готово! Итоговый файл: merged.pdf"
else
    echo "Не удалось создать PDF-файлы для объединения."
    exit 1
fi
