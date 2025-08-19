library(circlize)
library(tidyverse)

eczv <- read_tsv('bakt_files/ECZV.tsv', skip=2) %>% janitor::clean_names()
narna <- eczv %>%
  filter(gene == 'naRNA4') %>%
  select(number_sequence_id, start, stop, strand, gene)


circos.initializeWithIdeogram(plotType = NULL)
circos.genomicInitialize(eczv %>% select(number_sequence_id, start, stop, strand, gene))


circos.genomicTrackPlotRegion(
  narna, 
  ylim = c(0, 1), 
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = "red", border = "darkred", lwd = 2, ...)
  }
)


# Подготовка данных с координатами по Y
narna_plot <- narna %>%
  mutate(
    ytop = ifelse(strand == "+", 1.5, 0.5),
    ybottom = ifelse(strand == "+", 1.0, 0.0),
    color = ifelse(strand == "+", "red", "blue"),
    border_color = ifelse(strand == "+", "darkred", "darkblue")
  )

circos.initializeWithIdeogram(plotType = NULL)
circos.genomicInitialize(eczv %>% select(number_sequence_id, start, stop, strand, gene))

circos.genomicTrackPlotRegion(
  narna_plot, 
  ylim = c(0, 2), 
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, 
                       col = value$color, 
                       border = value$border_color, 
                       lwd = 2, 
                       ytop = value$ytop, 
                       ybottom = value$ybottom, 
                       ...)
  }
)

# Добавляем подписи
circos.genomicTrackPlotRegion(
  narna_plot,
  ylim = c(0, 2),
  panel.fun = function(region, value, ...) {
    y_position <- ifelse(value$strand == "+", 1.25, 0.25)
    circos.genomicText(region, value, labels = value$gene, 
                       facing = "clockwise", adj = c(0, 0.5), 
                       cex = 0.7, y = y_position, ...)
  }
)

circos.clear()


