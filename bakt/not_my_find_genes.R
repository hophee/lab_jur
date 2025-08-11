# AI redacted script
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(org.EcK12.eg.db))

# [1]-*.tsv path, [2] - *.pdf output
arg <- commandArgs(trailingOnly = TRUE)
tsv_path <- arg[1]
print(paste0("!!!!!!!", tsv_path, ("!!!!!!!")))
lines <- readLines(tsv_path)
hash_lines <- grep("^#", lines)
skip_lines <- if (length(hash_lines) > 0) max(hash_lines) - 1 else 0
eczv <- read_tsv(tsv_path, skip = skip_lines) %>% janitor::clean_names()

if ("gene" %in% colnames(eczv)) {
    print("Gene col is here!")
} else {
    stop("No Gene col - no plots!")
}
if (!("naRNA4" %in% eczv$gene)) {
    stop("naRNA4 not found in gene list")
}
eczv <- eczv %>%
    mutate(db_xrefs = strsplit(db_xrefs, ",\\s*")) %>%
    unnest(db_xrefs) %>%
    separate(db_xrefs, into = c("db", "id"), sep = ":", extra = "merge") %>%
    distinct(number_sequence_id, type, start, stop, strand, locus_tag, gene, product, db, .keep_all = TRUE) %>%
    pivot_wider(
        names_from = db,
        values_from = id,
        values_fill = NA,
        values_fn = ~ paste(unique(.x), collapse = "; ")
    )
entr <- bitr(
    eczv$gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.EcK12.eg.db
) %>% dplyr::rename("gene" = "SYMBOL")
eczv <- left_join(eczv, entr, by = "gene")
keggs <- bitr_kegg(na.omit(eczv$ENTREZID), fromType = "ncbi-geneid", toType = "kegg", organism = "eco", drop = TRUE) %>% dplyr::rename("ENTREZID" = "ncbi-geneid")
eczv <- left_join(eczv, keggs, by = "ENTREZID")

get_neighbor <- function(nd) {
    idxes <- which(eczv$gene == "naRNA4")
    right_neighbors <- eczv %>%
        dplyr::slice(outer((1:nd), idxes, "+") %>% as.vector()) %>%
        dplyr::filter(gene != "naRNA4") %>%
        mutate(side = "right")
    left_neighbors <- eczv %>%
        dplyr::slice(outer((-nd:-1), idxes, "+") %>% as.vector()) %>%
        dplyr::filter(gene != "naRNA4") %>%
        mutate(side = "left")
    neighbors <- rbind(right_neighbors, left_neighbors) %>%
        dplyr::arrange(start) %>%
        distinct(start, .keep_all = TRUE)
    return(neighbors)
}

get_goplot <- function(genelist, ont_type, side_num, nd, pcut_gp = 0.05, qcut_gp = 1, sc = 15) {
    # Фильтруем NA заранее
    genelist <- genelist[!is.na(genelist)]

    if (length(genelist) == 0) {
        return(ggplot() +
            ggtitle("No genes for GO analysis"))
    }

    types_int <- c("BP", "MF", "CC")
    onts <- c("Biological Process", "Molecular Function", "Cellular Component")
    side <- c("Right", "Left")

    tryCatch(
        {
            ego_result <- enrichGO(
                gene = genelist,
                OrgDb = org.EcK12.eg.db,
                keyType = "ENTREZID",
                ont = types_int[ont_type],
                pvalueCutoff = pcut_gp,
                pAdjustMethod = "BH",
                qvalueCutoff = qcut_gp,
                readable = TRUE
            )

            if (!is.null(ego_result) && nrow(ego_result@result) > 0) {
                plot <- goplot(ego_result, showCategory = sc) +
                    ggtitle(paste0(side[side_num], " genes, ", onts[ont_type], " ontology, nd=", as.character(nd))) +
                    theme(plot.title = element_text(hjust = 0.5))
                return(plot)
            } else {
                return(ggplot() +
                    ggtitle(paste0("No GO terms found for ", side[side_num], " genes, nd=", nd)))
            }
        },
        error = function(e) {
            return(ggplot() +
                ggtitle(paste0("GO plot error: ", e$message)))
        }
    )
}

# nd - neighbor distance
get_plots <- function(nd, pcut = 1, qcut = 1) {
    neighbors <- get_neighbor(nd)

    # Проверка наличия соседей
    if (nrow(neighbors) == 0) {
        empty_plot <- ggplot() +
            ggtitle(paste0("No neighbors found for nd=", nd))
        return(replicate(10, empty_plot, simplify = FALSE))
    }

    go_right <- neighbors %>%
        dplyr::filter(side == "right") %>%
        pull(ENTREZID) %>%
        na.omit() %>%
        as.vector()
    go_right <- go_right[!is.na(go_right)]

    go_left <- neighbors %>%
        dplyr::filter(side == "left") %>%
        pull(ENTREZID) %>%
        na.omit() %>%
        as.vector()
    go_left <- go_left[!is.na(go_left)]

    # Проверка наличия генов
    if (length(go_right) == 0 && length(go_left) == 0) {
        empty_plot <- ggplot() +
            ggtitle(paste0("No genes found for neighbors, nd=", nd))
        return(replicate(10, empty_plot, simplify = FALSE))
    }

    # GO analysis для правых генов
    right_dots <- tryCatch(
        {
            if (length(go_right) > 0) {
                ego_right <- enrichGO(
                    gene = go_right,
                    OrgDb = org.EcK12.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pvalueCutoff = pcut,
                    pAdjustMethod = "BH",
                    qvalueCutoff = qcut,
                    readable = TRUE
                )

                if (!is.null(ego_right) && nrow(ego_right@result) > 0) {
                    dotplot(ego_right,
                        showCategory = 10,
                        font.size = 10,
                        color = "p.adjust",
                        split = "ONTOLOGY",
                        title = paste0("Right genes, nd=", as.character(nd))
                    ) + facet_grid(ONTOLOGY ~ ., scales = "free")
                } else {
                    ggplot() +
                        ggtitle(paste0("No GO terms for right genes, nd=", nd))
                }
            } else {
                ggplot() +
                    ggtitle(paste0("No right genes found, nd=", nd))
            }
        },
        error = function(e) {
            ggplot() +
                ggtitle(paste0("GO right error: ", e$message))
        }
    )

    # GO plots для правых генов
    right_go_bp <- get_goplot(go_right, 1, 1, nd, pcut_gp = pcut)
    right_go_mf <- get_goplot(go_right, 2, 1, nd, pcut_gp = pcut)
    right_go_cc <- get_goplot(go_right, 3, 1, nd, pcut_gp = pcut)

    # GO analysis для левых генов
    left_dots <- tryCatch(
        {
            if (length(go_left) > 0) {
                ego_left <- enrichGO(
                    gene = go_left,
                    OrgDb = org.EcK12.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pvalueCutoff = pcut,
                    pAdjustMethod = "BH",
                    qvalueCutoff = qcut,
                    readable = TRUE
                )

                if (!is.null(ego_left) && nrow(ego_left@result) > 0) {
                    dotplot(ego_left,
                        showCategory = 10,
                        font.size = 10,
                        color = "p.adjust",
                        split = "ONTOLOGY",
                        title = paste0("Left genes, nd=", as.character(nd))
                    ) + facet_grid(ONTOLOGY ~ ., scales = "free")
                } else {
                    ggplot() +
                        ggtitle(paste0("No GO terms for left genes, nd=", nd))
                }
            } else {
                ggplot() +
                    ggtitle(paste0("No left genes found, nd=", nd))
            }
        },
        error = function(e) {
            ggplot() +
                ggtitle(paste0("GO left error: ", e$message))
        }
    )

    # GO plots для левых генов
    left_go_bp <- get_goplot(go_left, 1, 2, nd, pcut_gp = pcut)
    left_go_mf <- get_goplot(go_left, 2, 2, nd, pcut_gp = pcut)
    left_go_cc <- get_goplot(go_left, 3, 2, nd, pcut_gp = pcut)

    # KEGG analysis
    right_kegg <- neighbors %>%
        dplyr::filter(side == "right") %>%
        pull(kegg) %>%
        na.omit() %>%
        as.vector()
    right_kegg <- right_kegg[!is.na(right_kegg)]

    left_kegg <- neighbors %>%
        dplyr::filter(side == "left") %>%
        pull(kegg) %>%
        na.omit() %>%
        as.vector()
    left_kegg <- left_kegg[!is.na(left_kegg)]

    dot_kegg_right <- tryCatch(
        {
            if (length(right_kegg) > 0) {
                ekegg_right <- enrichKEGG(
                    gene = right_kegg,
                    organism = "eco",
                    pvalueCutoff = 1
                )

                if (nrow(ekegg_right@result[ekegg_right@result$p.adjust <= 0.05, ]) > 0) {
                    dotplot(ekegg_right, title = paste0("Right genes, KEGG, nd=", nd))
                } else {
                    ggplot() +
                        ggtitle(paste0("No KEGG enrichment found for right genes, nd=", nd))
                }
            } else {
                ggplot() +
                    ggtitle(paste0("No KEGG IDs for right genes, nd=", nd))
            }
        },
        error = function(e) {
            ggplot() +
                ggtitle(paste0("KEGG right error: ", e$message))
        }
    )

    dot_kegg_left <- tryCatch(
        {
            if (length(left_kegg) > 0) {
                ekegg_left <- enrichKEGG(
                    gene = left_kegg,
                    organism = "eco",
                    pvalueCutoff = 1
                )

                if (nrow(ekegg_left@result[ekegg_left@result$p.adjust <= 0.05, ]) > 0) {
                    dotplot(ekegg_left, title = paste0("Left genes, KEGG, nd=", nd))
                } else {
                    ggplot() +
                        ggtitle(paste0("No KEGG enrichment found for left genes, nd=", nd))
                }
            } else {
                ggplot() +
                    ggtitle(paste0("No KEGG IDs for left genes, nd=", nd))
            }
        },
        error = function(e) {
            ggplot() +
                ggtitle(paste0("KEGG left error: ", e$message))
        }
    )

    return(list(right_dots, right_go_mf, right_go_cc, right_go_bp, left_dots, left_go_mf, left_go_cc, left_go_bp, dot_kegg_right, dot_kegg_left))
}

# Функция для вывода всех графиков в PDF
print_plots <- function(plots_list) {
    for (plot in plots_list) {
        tryCatch(
            {
                print(plot)
            },
            error = function(e) {
                print(ggplot() +
                    ggtitle(paste("Plot print error:", e$message)))
            }
        )
    }
}

# merged pdf of nds-s
pdf(arg[2], width = 15, height = 12)
tryCatch(
    {
        result1 <- get_plots(nd = 1, pcut = 0.5)
        print_plots(result1)

        result3 <- get_plots(nd = 3, pcut = 0.3)
        print_plots(result3)
    },
    error = function(e) {
        print(ggplot() +
            ggtitle(paste("Critical error:", e$message)))
    }
)
dev.off()
