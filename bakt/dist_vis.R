library(ComplexHeatmap)
library(tidyverse)
library(reshape2)

# Чтение файла
text <- readLines("/home/iz-user/Documents/proteo_lab/bakta_zvl/folding/dist.mat")[-169]

# Извлечение названий
names <- sub("^>", "", grep("^>", text, value = TRUE))[-85]

# Фильтруем только строки с числами
data_lines <- text[!grepl("^>", text) & text != ""]

# Преобразуем каждую строку в числовой вектор
split_lines <- lapply(data_lines, function(line) as.numeric(strsplit(trimws(line), "\\s+")[[1]]))

# Проверка правильности структуры треугольника
expected_lengths <- seq_along(split_lines)
actual_lengths <- sapply(split_lines, length)

if (!all(actual_lengths == expected_lengths)) {
  stop("Формат нижнего треугольника нарушен. Строки не имеют ожидаемой длины.")
}

# Объединяем в один вектор
dist_values <- unlist(split_lines)

# Получаем размер матрицы
n <- length(split_lines) + 1  # т.к. в нижнем треугольнике n-1 строк

# Создаем и заполняем матрицу
dist_matrix <- matrix(0, nrow = n, ncol = n, byrow = TRUE)
n <- nrow(dist_matrix)

# Получаем индексы нижнего треугольника в матрице (строка и столбец)
idx <- which(lower.tri(dist_matrix), arr.ind = TRUE)

# Сортируем индексы по строкам, затем по столбцам (то есть построчно)
idx_ordered <- idx[order(idx[,1], idx[,2]), ]

# Теперь запишем dist_values в матрицу по этим индексам построчно
for (i in seq_len(nrow(idx_ordered))) {
  dist_matrix[idx_ordered[i, "row"], idx_ordered[i, "col"]] <- dist_values[i]
}

# Делаем симметричную матрицу
dist_matrix <- dist_matrix + t(dist_matrix)


# Присваиваем имена
rownames(dist_matrix) <- colnames(dist_matrix) <- names


# Тепловая карта
Heatmap(dist_matrix, row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5), heatmap_legend_param = list(title = 'base pairs \n distance'))
#heatmap(dist_matrix, symm = TRUE, main = "Distance Matrix Heatmap", cexRow = 0.5, cexCol = 0.5)

#Дерево
hc <- hclust(as.dist(dist_matrix))
plot(hc, main = "Hierarchical Clustering", xlab = '', sub = '')