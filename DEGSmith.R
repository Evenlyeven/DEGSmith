suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(writexl)))
suppressMessages(suppressWarnings(library(rmarkdown)))
suppressMessages(suppressWarnings(library(plotly)))

## ===== define a function of the analysis ===== ##
write_degsmith_report_rmd <- function(rmd_path,
                                      rdata_path,
                                      report_title = "DEGSmith Volcano Plots",
                                      report_mode = c("per_celltype", "combined", "single"),
                                      volcano_logfc = 0.25,
                                      volcano_padj = 0.05) {
  report_mode <- match.arg(report_mode)
  
  lines <- c(
    "---",
    paste0("title: \"", report_title, "\""),
    "output:",
    "  html_document:",
    "    theme: flatly",
    "    toc: false",
    "    df_print: paged",
    "params:",
    "  cell_type: ~",
    paste0("  report_mode: \"", report_mode, "\""),
    "---",
    "",
    "<style>",
    ".main-container { max-width: 100% !important; }",
    ".plotly.html-widget { height: 100vh; width: 100%; }",
    "</style>",
    "",
    "```{r echo=FALSE, message=FALSE, warning=FALSE, fig.width=14, fig.height=8, out.width='100%'}",
    "library(ggplot2)",
    "library(plotly)",
    "library(magrittr)",
    "library(tidyverse)",
    "",
    paste0("load('", rdata_path, "')"),
    "stopifnot(exists('results'))",
    "plot_dir <- file.path(results$output_dir, 'plots')",
    "dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)",
    "",
    "`%||%` <- function(x, y) if (is.null(x)) y else x",
    "",
    "cond1 <- results$params$condition_1",
    "cond2 <- results$params$condition_2",
    "table_choice <- results$volcano_table %||% 'specific'",
    "",
    paste0("VOLCANO_LOGFC <- ", volcano_logfc),
    paste0("VOLCANO_PADJ  <- ", volcano_padj),
    "",
    "volcano_plot <- function(degs, title) {",
    "  if (is.null(degs) || nrow(degs) == 0) return(NULL)",
    "  degs <- degs %>% as_tibble()",
    "  if (!all(c('gene','avg_log2FC','p_val_adj') %in% colnames(degs))) {",
    "    # if gene is rownames",
    "    if (!'gene' %in% colnames(degs) ) degs <- tibble::rownames_to_column(degs, 'gene')",
    "  }",
    "  if (!all(c('gene','avg_log2FC','p_val_adj') %in% colnames(degs))) return(NULL)",
    "",
    "  degs <- degs %>% mutate(p_val_adj = ifelse(p_val_adj == 0, .Machine$double.xmin, p_val_adj))",
    "",
    "  up_lab   <- paste0('Upregulated in ', cond1)",
    "  down_lab <- paste0('Downregulated in ', cond1)",
    "",
    "  degs_vp <- degs %>%",
    "    mutate(significance = case_when(",
    "      p_val_adj <= VOLCANO_PADJ & avg_log2FC >=  VOLCANO_LOGFC ~ up_lab,",
    "      p_val_adj <= VOLCANO_PADJ & avg_log2FC <= -VOLCANO_LOGFC ~ down_lab,",
    "      TRUE ~ 'Not significant'",
    "    ))",
    "",
    "  p <- ggplot(",
    "    degs_vp,",
    "    aes(",
    "      x = avg_log2FC,",
    "      y = -log10(p_val_adj),",
    "      color = significance,",
    "      text = paste0(",
    "        'Gene: ', gene,",
    "        '<br>avg_log2FC: ', signif(avg_log2FC, 4),",
    "        '<br>-log10(p_val_adj): ', signif(-log10(p_val_adj), 4),",
    "        '<br>p_val_adj: ', signif(p_val_adj, 4),",
    "        '<br>Significance: ', significance",
    "      )",
    "    )",
    "  ) +",
    "    geom_point(alpha = 0.6, size = 0.8) +",
    "    scale_color_manual(values = setNames(",
    "      c('red', 'forestgreen', 'grey'),",
    "      c(up_lab, down_lab, 'Not significant')",
    "    )) +",
    "    theme_minimal() +",
    "    labs(title = title, x = 'Average Log2 Fold Change', y = '-Log10 Adjusted P-value') +",
    "    theme(text = element_text(size = 16), legend.title = element_blank(), legend.position = 'right') +",
    "    geom_vline(xintercept = c(-VOLCANO_LOGFC, VOLCANO_LOGFC), linetype = 'dashed', color = 'black') +",
    "    geom_hline(yintercept = -log10(VOLCANO_PADJ), linetype = 'dashed', color = 'black')",
    "",
    "  ggplotly(p, tooltip = 'text') |>",
    "    layout(",
    "      autosize = TRUE,",
    "      margin = list(l = 60, r = 20, b = 60, t = 80),",
    "      legend = list(orientation = 'v')",
    "    ) |>",
    "    config(displayModeBar = TRUE, responsive = TRUE) |>",
    "    htmlwidgets::onRender(\"",
    "      function(el) {",
    "        var t = document.createElement('div');",
    "        t.style.cssText = 'position:fixed;top:20px;left:50%;transform:translateX(-50%);background:rgba(40,40,40,.95);color:#fff;padding:12px 28px;border-radius:10px;font:600 15px/1 monospace;z-index:9999;opacity:0;transition:opacity .3s;pointer-events:none;box-shadow:0 4px 12px rgba(0,0,0,.3)';",
    "        document.body.appendChild(t);",
    "        el.on('plotly_click', function(d) {",
    "          if (!d||!d.points||!d.points[0]) return;",
    "          var p = d.points[0];",
    "          var txt = Array.isArray(p.data.text) ? p.data.text[p.pointNumber] : p.data.text;",
    "          if (!txt) return;",
    "          var m = txt.match(/Gene:\\\\s*([^<]+)/);",
    "          if (!m) return;",
    "          var gene = m[1].trim();",
    "          function show() {",
    "            t.textContent = '\\u2705 Copied: ' + gene;",
    "            t.style.opacity = '1';",
    "            setTimeout(function(){ t.style.opacity='0'; }, 1800);",
    "          }",
    "          if (navigator.clipboard && navigator.clipboard.writeText) {",
    "            navigator.clipboard.writeText(gene).then(show).catch(function(){",
    "              var a=document.createElement('textarea');a.value=gene;",
    "              a.style.cssText='position:fixed;left:-9999px';",
    "              document.body.appendChild(a);a.select();",
    "              try{document.execCommand('copy')}catch(e){}",
    "              document.body.removeChild(a);show();",
    "            });",
    "          }",
    "        });",
    "      }",
    "    \")",
    "    htmlwidgets::onRender(\"",
    "      function(el) {",
    "        var t = document.createElement('div');",
    "        t.style.cssText = 'position:fixed;top:20px;left:50%;transform:translateX(-50%);background:rgba(40,40,40,.95);color:#fff;padding:12px 28px;border-radius:10px;font:600 15px/1 monospace;z-index:9999;opacity:0;transition:opacity .3s;pointer-events:none;box-shadow:0 4px 12px rgba(0,0,0,.3)';",
    "        document.body.appendChild(t);",
    "        el.on('plotly_click', function(d) {",
    "          if (!d||!d.points||!d.points[0]) return;",
    "          var p = d.points[0];",
    "          var txt = Array.isArray(p.data.text) ? p.data.text[p.pointNumber] : p.data.text;",
    "          if (!txt) return;",
    "          var m = txt.match(/Gene:\\\\s*([^<]+)/);",
    "          if (!m) return;",
    "          var gene = m[1].trim();",
    "          function show() {",
    "            t.textContent = '\\u2705 Copied: ' + gene;",
    "            t.style.opacity = '1';",
    "            setTimeout(function(){ t.style.opacity='0'; }, 1800);",
    "          }",
    "          if (navigator.clipboard && navigator.clipboard.writeText) {",
    "            navigator.clipboard.writeText(gene).then(show).catch(function(){",
    "              var a=document.createElement('textarea');a.value=gene;",
    "              a.style.cssText='position:fixed;left:-9999px';",
    "              document.body.appendChild(a);a.select();",
    "              try{document.execCommand('copy')}catch(e){}",
    "              document.body.removeChild(a);show();",
    "            });",
    "          }",
    "        });",
    "      }",
    "    \")",
    "}",
    "",
    "plot_list <- if (table_choice == 'specific') results$deg_specific_list else results$deg_list",
    "mode <- params$report_mode %||% 'per_celltype'",
    "ct_param <- params$cell_type",
    "",
    "if (!is.null(ct_param) && nzchar(ct_param)) {",
    "  degs <- plot_list[[ct_param]]",
    "  w <- volcano_plot(degs, paste0('DEGSmith: ', ct_param, ' | ', cond1, ' vs ', cond2))",
    "",
    "  cat('## ', ct_param, '\\n\\n', sep='')",
    "",
    "  if (!is.null(w)) {",
    "    htmlwidgets::saveWidget(",
    "      widget = w,",
    "      file = file.path(plot_dir, paste0(gsub('[^A-Za-z0-9._-]+', '_', ct_param), '_volcano.html')),",
    "      selfcontained = FALSE",
    "    )",
    "    print(w)",
    "  } else {",
    "    cat('No DEGs to plot.\\n')",
    "  }",
    "",
    "} else if (mode == 'combined') {",
    "  degs <- results$deg_all_tbl",
    "  w <- volcano_plot(degs, paste0('DEGSmith: all cell types | ', cond1, ' vs ', cond2))",
    "  if (!is.null(w)) print(w) else cat('No DEGs to plot.\\n')",
    "",
    "} else {",
    "  for (ct in names(plot_list)) {",
    "    degs <- plot_list[[ct]]",
    "    if (is.null(degs) || nrow(degs) == 0) next",
    "    w <- volcano_plot(degs, paste0('DEGSmith: ', ct, ' | ', cond1, ' vs ', cond2))",
    "    cat('## ', ct, '\\n\\n', sep='')",
    "    if (!is.null(w)) print(w) else cat('No DEGs to plot.\\n')",
    "    cat('\\n\\n')",
    "  }",
    "}",
    "```"
  )
  
  writeLines(lines, con = rmd_path)
}

render_degsmith_per_celltype <- function(out_dir,
                                         prefix,
                                         rmd_template,
                                         cell_types,
                                         safe_tag) {
  # Volcano widgets will still be written to results$output_dir/plots by the Rmd itself.
  # We render the *wrapper* HTML files into a temp directory so they don't clutter out_dir.
  
  tmp_dir <- file.path(out_dir, "__tmp_rmd_render__")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (ct in cell_types) {
    wrapper_html <- paste0(prefix, safe_tag(ct), ".html")  # this is the unwanted file
    rmarkdown::render(
      input = rmd_template,
      output_file = wrapper_html,
      output_dir  = tmp_dir,
      params = list(cell_type = ct, report_mode = "single"),
      quiet = TRUE
    )
  }
  
  # delete wrapper outputs in one shot
  unlink(tmp_dir, recursive = TRUE, force = TRUE)
}

cellTypeDE <- function(seurat_obj,
                       seurat_format = c("auto", "rds", "rdata", "qs2"),
                       output_dir,
                       cell_types,
                       condition_1,
                       condition_2,
                       cond_col = "orig.ident",
                       celltype_col = "annotation",
                       de_assay = NA,
                       de_slot = NA,
                       marker_assay = "SCT",
                       marker_slot = "data",
                       test_use = "wilcox",
                       de_min_pct = 0.01,
                       de_logfc_threshold = 0.1,
                       marker_min_pct = 0.01,
                       marker_logfc_threshold = 0.1,
                       marker_only_pos = TRUE,
                       marker_logfc_cutoff = 0.5,
                       marker_padj_cutoff = 0.05,
                       celltype_deg_logfc_cutoff = 0.25,
                       celltype_deg_padj_cutoff  = 0.05,
                       min_cells_per_group = 10,
                       saveRData = TRUE,
                       write_combined_workbook = TRUE,
                       volcano_logfc = 0.25,
                       volcano_padj = 0.05,
                       volcano_table = c("specific", "deg"),
                       prefix = "",
                       make_report = FALSE,
                       render_report = FALSE,
                       report_title = "DEGSmith Volcano Plots",
                       report_mode = c("per_celltype", "combined")) {
  seurat_format <- match.arg(seurat_format)
  
  ## ---- read seurat object (rds / RData / qs2) ---- ##
  in_path <- seurat_obj
  if (!file.exists(in_path))
    stop("Input file does not exist: ", in_path)
  
  ext <- tolower(tools::file_ext(in_path))
  
  # Decide format
  fmt <- seurat_format
  if (fmt == "auto") {
    fmt <- dplyr::case_when(
      ext == "rds"                    ~ "rds",
      ext %in% c("RData", "rda")      ~ "RData",
      ext == "qs2"                     ~ "qs2",
      TRUE                            ~ NA_character_
    )
    if (is.na(fmt)) {
      stop(
        "Could not infer file type from extension: .",
        ext,
        "\nUse --seurat_format rds|RData|qs2"
      )
    }
  }
  
  # Read based on format
  if (fmt == "rds") {
    seurat_obj <- readRDS(in_path)
    
  } else if (fmt == "RData") {
    temp_env <- new.env()
    load(in_path, envir = temp_env)
    
    seurat_vars <- ls(temp_env)
    seurat_objs <- Filter(function(x)
      inherits(temp_env[[x]], "Seurat"), seurat_vars)
    
    if (length(seurat_objs) == 0) {
      stop("No Seurat object found in the RData file.")
    } else if (length(seurat_objs) > 1) {
      warning("Multiple Seurat objects found. Using the first one: ",
              seurat_objs[1])
    }
    
    seurat_obj <- temp_env[[seurat_objs[1]]]
    message("Using Seurat object from .RData: ", seurat_objs[1])
    
  } else if (fmt == "qs2") {
    if (!requireNamespace("qs2", quietly = TRUE)) {
      stop(
        "Input is qs2 but package 'qs2' is not installed. Install it or use a different format."
      )
    }
    seurat_obj <- qs2::qs_read(in_path)
    
  } else {
    stop("Unsupported --seurat_format: ", fmt)
  }
  
  volcano_table <- match.arg(volcano_table) #accepts exactly one of those choices
  report_mode <- match.arg(report_mode)
  
  ## ---- checks ---- ##
  if (!inherits(seurat_obj, "Seurat"))
    stop("The input must be a Seurat object.")
  
  if (!cond_col %in% colnames(seurat_obj@meta.data)) {
    stop(
      "cond_col='",
      cond_col,
      "' not found in meta.data. Available columns:\n",
      paste(colnames(seurat_obj@meta.data), collapse = ", ")
    )
  }
  if (!celltype_col %in% colnames(seurat_obj@meta.data)) {
    stop(
      "celltype_col='",
      celltype_col,
      "' not found in meta.data. Available columns:\n",
      paste(colnames(seurat_obj@meta.data), collapse = ", ")
    )
  }
  
  if (!marker_assay %in% Assays(seurat_obj)) {
    stop(
      "marker_assay='",
      marker_assay,
      "' not found. Available assays: ",
      paste(Assays(seurat_obj), collapse = ", ")
    )
  }
  
  # Parse cell types input:
  # - comma-separated list, OR
  # - the literal 'all' (case-insensitive) to use all cell types present in celltype_col
  cell_types_raw <- cell_types
  
  if (length(cell_types_raw) == 1 &&
      tolower(trimws(cell_types_raw)) == "all") {
    cell_types <- seurat_obj@meta.data[[celltype_col]] %>%
      as.character() %>%
      trimws() %>%
      discard( ~ is.na(.x) || .x == "")
    
    cell_types <- sort(unique(cell_types))
    
    message(
      "cell_types='all' detected. Using all cell types (n=",
      length(cell_types),
      "): ",
      paste(cell_types, collapse = ", ")
    )
    
  } else {
    # normal behavior: allow comma-separated string or vector
    if (length(cell_types_raw) == 1 && grepl(",", cell_types_raw)) {
      cell_types <- strsplit(cell_types_raw, ",")[[1]] %>% trimws()
    } else {
      cell_types <- as.character(cell_types_raw) %>% trimws()
    }
    
    cell_types <- cell_types %>% discard( ~ is.na(.x) || .x == "")
  }
  
  # pre-filter cell types that don't have enough cells in BOTH conditions
  counts_by_ct_cond <- seurat_obj@meta.data %>%
    as_tibble() %>%
    mutate(ct = as.character(.data[[celltype_col]]), cond = as.character(.data[[cond_col]])) %>%
    filter(!is.na(ct), ct != "", !is.na(cond), cond != "") %>%
    count(ct, cond, name = "n")
  
  eligible_ct <- counts_by_ct_cond %>%
    filter(cond %in% c(condition_1, condition_2)) %>%
    tidyr::pivot_wider(
      names_from = cond,
      values_from = n,
      values_fill = 0
    ) %>%
    filter(.data[[condition_1]] >= min_cells_per_group, .data[[condition_2]] >= min_cells_per_group) %>%
    pull(ct)
  
  dropped <- setdiff(cell_types, eligible_ct)
  if (length(dropped) > 0) {
    message(
      "Skipping ",
      length(dropped),
      " cell types due to min_cells_per_group in BOTH conditions: ",
      paste(dropped, collapse = ", ")
    )
  }
  cell_types <- eligible_ct
  
  # File-name safe tag
  safe_tag <- function(x) {
    x <- as.character(x) %>% trimws()
    x <- gsub("&", "and", x, fixed = TRUE)
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    ifelse(nzchar(x), x, "NA")
  }
  
  ## ---- output folder (timestamped) ---- ##
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_dir <- file.path(output_dir, paste0(prefix, "DEGSmith_", timestamp))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("Output folder: ", out_dir)
  message("Cell types: ", paste(cell_types, collapse = ", "))
  message("Compare within cell type: ",
          condition_1,
          " vs ",
          condition_2,
          " (",
          cond_col,
          ")")
  message("Cell type column: ", celltype_col)
  
  ## ---- interpret assay/slot NA -> NULL (so Seurat defaults if not set) ---- ##
  de_assay_use <- if (is.na(de_assay) ||
                      de_assay == "NA" || de_assay == "")
    NULL
  else
    de_assay
  de_slot_use  <- if (is.na(de_slot)  ||
                      de_slot  == "NA" || de_slot  == "")
    NULL
  else
    de_slot
  
  ## ---- main loops ---- ##
  deg_list <- list()
  marker_list <- list()
  marker_keep_list <- list()
  deg_specific_list <- list()
  deg_specific_filtered_list <- list()
  
  meta_tbl <- seurat_obj@meta.data %>% as_tibble(rownames = "cell_id")
  
  for (ct in cell_types) {
    message("----- ", ct, " -----")
    
    # Cells for within-cell-type DE
    cells_1 <- meta_tbl %>%
      filter(.data[[celltype_col]] == ct, .data[[cond_col]] == condition_1) %>%
      pull(cell_id)
    cells_2 <- meta_tbl %>%
      filter(.data[[celltype_col]] == ct, .data[[cond_col]] == condition_2) %>%
      pull(cell_id)
    
    if (length(cells_1) < min_cells_per_group ||
        length(cells_2) < min_cells_per_group) {
      warning(
        sprintf(
          "[%s] Too few cells: %s=%d, %s=%d. DEG will be empty.",
          ct,
          condition_1,
          length(cells_1),
          condition_2,
          length(cells_2)
        )
      )
      deg_tbl <- tibble()
    } else {
      deg_tbl <- FindMarkers(
        object = seurat_obj,
        ident.1 = cells_1,
        ident.2 = cells_2,
        assay = de_assay_use,
        slot = de_slot_use,
        test.use = test_use,
        min.pct = de_min_pct,
        logfc.threshold = de_logfc_threshold
      ) %>%
        rownames_to_column("gene") %>%
        arrange(desc(avg_log2FC))
    }
    
    deg_list[[ct]] <- deg_tbl
    
    # Markers: cell type vs rest
    Idents(seurat_obj) <- seurat_obj@meta.data[[celltype_col]]
    
    marker_tbl <- FindMarkers(
      object = seurat_obj,
      ident.1 = ct,
      ident.2 = NULL,
      assay = marker_assay,
      slot = marker_slot,
      test.use = test_use,
      only.pos = marker_only_pos,
      min.pct = marker_min_pct,
      logfc.threshold = marker_logfc_threshold
    ) %>%
      rownames_to_column("gene") %>%
      arrange(desc(avg_log2FC))
    
    marker_list[[ct]] <- marker_tbl
    
    marker_keep <- marker_tbl %>%
      filter(avg_log2FC >= marker_logfc_cutoff,
             p_val_adj <= marker_padj_cutoff)
    marker_keep_list[[ct]] <- marker_keep
    
    # Specific DEGs = overlap(DEGs, marker_keep)
    deg_specific <- deg_tbl %>%
      filter(gene %in% marker_keep$gene)
    deg_specific_list[[ct]] <- deg_specific
    
    deg_specific_filtered <- deg_specific %>%
      filter(
        abs(avg_log2FC) >= celltype_deg_logfc_cutoff,
        p_val_adj <= celltype_deg_padj_cutoff
      )
    deg_specific_filtered_list[[ct]] <- deg_specific_filtered
    
    # Write per-celltype files
    ct_tag <- safe_tag(ct)
    
    writexl::write_xlsx(deg_tbl, path = file.path(
      out_dir,
      paste0(
        prefix,
        "DEG_",
        ct_tag,
        "_",
        condition_1,
        "_vs_",
        condition_2,
        ".xlsx"
      )
    ))
    
    writexl::write_xlsx(
      list(markers_all = marker_tbl, markers_keep = marker_keep),
      path = file.path(
        out_dir,
        paste0(prefix, "Markers_", ct_tag, "_vs_rest.xlsx")
      )
    )
    
    writexl::write_xlsx(deg_specific, path = file.path(
      out_dir,
      paste0(
        prefix,
        "DEG_specific_",
        ct_tag,
        "_",
        condition_1,
        "_vs_",
        condition_2,
        ".xlsx"
      )
    ))
    
    writexl::write_xlsx(deg_specific_filtered, path = file.path(
      out_dir,
      paste0(
        prefix,
        "DEG_specific_filtered_",
        ct_tag,
        "_",
        condition_1,
        "_vs_",
        condition_2,
        "_logfc",
        celltype_deg_logfc_cutoff,
        "_padj",
        celltype_deg_padj_cutoff,
        ".xlsx"
      )
    ))
  }
  
  # Choose which table is used for volcano/report
  volcano_source_list <- if (volcano_table == "specific")
    deg_specific_list
  else
    deg_list
  
  # Combined table (used only if report_mode == "combined")
  deg_all_tbl <- purrr::imap_dfr(volcano_source_list, function(df, ct) {
    if (is.null(df) || nrow(df) == 0)
      return(tibble())
    df %>% mutate(cell_type = ct)
  })
  
  ## ---- save RData / combined workbook ---- ##
  results <- list(
    deg_list = deg_list,
    marker_list = marker_list,
    marker_keep_list = marker_keep_list,
    deg_specific_list = deg_specific_list,
    deg_specific_filtered_list = deg_specific_filtered_list,
    deg_all_tbl = deg_all_tbl,
    volcano_table = volcano_table,
    params = list(
      cell_types = cell_types,
      condition_1 = condition_1,
      condition_2 = condition_2,
      cond_col = cond_col,
      celltype_col = celltype_col,
      test_use = test_use,
      de_assay = de_assay_use,
      de_slot = de_slot_use,
      marker_assay = marker_assay,
      marker_slot = marker_slot,
      de_min_pct = de_min_pct,
      de_logfc_threshold = de_logfc_threshold,
      marker_min_pct = marker_min_pct,
      marker_logfc_threshold = marker_logfc_threshold,
      marker_only_pos = marker_only_pos,
      marker_logfc_cutoff = marker_logfc_cutoff,
      marker_padj_cutoff = marker_padj_cutoff,
      min_cells_per_group = min_cells_per_group
    ),
    output_dir = out_dir
  )
  
  if (saveRData) {
    save(results, file = file.path(out_dir, paste0(prefix, "DEGSmith_results.RData")))
  }
  
  ## ---- generate / render HTML volcano report ---- ##
  if (isTRUE(make_report) || isTRUE(render_report)) {
    if (!saveRData) {
      stop(
        "Report generation needs saveRData=TRUE because the Rmd loads DEGSmith_results.RData."
      )
    }
    
    rdata_path <- file.path(out_dir, paste0(prefix, "DEGSmith_results.RData"))
    rmd_path   <- file.path(out_dir, paste0(prefix, "DEGSmith_volcano_report.Rmd"))
    
    # write template
    write_degsmith_report_rmd(
      rmd_path = rmd_path,
      rdata_path = rdata_path,
      report_title = report_title,
      report_mode = report_mode,
      # per_celltype / combined
      volcano_logfc = volcano_logfc,
      volcano_padj = volcano_padj
    )
    message("Wrote report template: ", rmd_path)
    
    if (isTRUE(render_report)) {
      if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        stop(
          "Package 'rmarkdown' is required to render HTML. Install it or run with --make_report only."
        )
      }
      
      if (report_mode == "per_celltype") {
        # one HTML per cell type (filename = cell type)
        render_degsmith_per_celltype(
          out_dir = out_dir,
          prefix = prefix,
          rmd_template = rmd_path,
          cell_types = names(volcano_source_list),
          # plots are based on volcano_source_list
          safe_tag = safe_tag
        )
        message("Rendered per-cell-type HTML files in: ", out_dir)
        
      } else {
        # single combined HTML
        plot_dir <- file.path(out_dir, "plots")
        dir.create(plot_dir,
                   recursive = TRUE,
                   showWarnings = FALSE)
        
        rmarkdown::render(
          input = rmd_path,
          output_dir = plot_dir,
          params = list(report_mode = "combined"),
          quiet = TRUE
        )
        message("Rendered combined HTML report in: ", out_dir)
      }
    }
  }
  
  if (write_combined_workbook) {
    sheets <- list()
    for (ct in names(deg_list)) {
      ct_tag <- safe_tag(ct)
      sheets[[paste0(ct_tag, "_DEG")]] <- deg_list[[ct]]
      sheets[[paste0(ct_tag, "_MarkersKeep")]] <- marker_keep_list[[ct]]
      sheets[[paste0(ct_tag, "_SpecificDEG")]] <- deg_specific_list[[ct]]
      sheets[[paste0(ct_tag, "_SpecificDEG_Filtered")]] <- deg_specific_filtered_list[[ct]]
    }
    writexl::write_xlsx(sheets, path = file.path(out_dir, paste0(prefix, "combined_cellTypeDE.xlsx")))
  }
  
  message("Done.")
  invisible(results)
}

## ===== define options for the script ===== ##
description_text <- "This script performs cell-type-specific differential expression (DE) using Seurat FindMarkers.

For each requested cell type:
1) Within-cell-type DE: (cell type, condition_1) vs (cell type, condition_2)
2) Cell-type markers: (cell type) vs (all other cell types)
3) Specific DEGs: overlap between (1) and filtered markers from (2)

Outputs:
- DEG_<celltype>_<cond1>_vs_<cond2>.xlsx
- Markers_<celltype>_vs_rest.xlsx (sheets: markers_all, markers_keep)
- DEG_specific_<celltype>_<cond1>_vs_<cond2>.xlsx
- optional: DEGSmith_results.RData
- optional: combined_cellTypeDE.xlsx

Usage:
Rscript cellTypeDE.R --seurat_obj <rds/RData/qs2> --output_dir <dir> --cell_types <comma-separated> --condition_1 <A> --condition_2 <B> [other options...]"

option_list <- list(
  make_option(
    c("--seurat_obj"),
    type = "character",
    default = NULL,
    help = "Path to the Seurat object file (rds, RData or qs2)"
  ),
  make_option(
    c("--seurat_format"),
    type = "character",
    default = "auto",
    help = "Input format: 'auto' (default), 'rds', 'RData', or 'qs2'"
  ),
  make_option(
    c("--output_dir"),
    type = "character",
    default = "./",
    help = "Directory to save the output files"
  ),
  make_option(
    c("--cell_types"),
    type = "character",
    default = NULL,
    help = "Comma-separated cell types (e.g., Cycling_Hi_Tumor1,Cycling_Hi_Tumor4) OR 'all' to use all cell types in celltype_col"
  ),
  
  make_option(
    c("--condition_1"),
    type = "character",
    default = NULL,
    help = "Condition A label (e.g., Pd1_C5aR)"
  ),
  make_option(
    c("--condition_2"),
    type = "character",
    default = NULL,
    help = "Condition B label (e.g., Ctrl_IgG)"
  ),
  
  make_option(
    c("--cond_col"),
    type = "character",
    default = "orig.ident",
    help = "Metadata column for condition (default: orig.ident)"
  ),
  make_option(
    c("--celltype_col"),
    type = "character",
    default = "annotation",
    help = "Metadata column for cell type (default: annotation)"
  ),
  
  make_option(
    c("--test_use"),
    type = "character",
    default = "wilcox",
    help = "DE test for both steps (default: wilcox)"
  ),
  
  # within-cell-type DE
  make_option(
    c("--de_assay"),
    type = "character",
    default = "NA",
    help = "Assay for within-cell-type DE (default: NA = Seurat default)"
  ),
  make_option(
    c("--de_slot"),
    type = "character",
    default = "NA",
    help = "Slot for within-cell-type DE (default: NA = Seurat default)"
  ),
  make_option(
    c("--de_min_pct"),
    type = "double",
    default = 0.01,
    help = "min.pct for DE (default: 0.01)"
  ),
  make_option(
    c("--de_logfc_threshold"),
    type = "double",
    default = 0.1,
    help = "logfc.threshold for DE (default: 0.1)"
  ),
  
  # markers
  make_option(
    c("--marker_assay"),
    type = "character",
    default = "SCT",
    help = "Assay for marker calling (default: SCT)"
  ),
  make_option(
    c("--marker_slot"),
    type = "character",
    default = "data",
    help = "Slot for marker calling (default: data)"
  ),
  make_option(
    c("--marker_only_pos"),
    type = "logical",
    default = TRUE,
    help = "only.pos for marker calling (default: TRUE)"
  ),
  make_option(
    c("--marker_min_pct"),
    type = "double",
    default = 0.01,
    help = "min.pct for markers (default: 0.01)"
  ),
  make_option(
    c("--marker_logfc_threshold"),
    type = "double",
    default = 0.1,
    help = "logfc.threshold for markers (default: 0.1)"
  ),
  make_option(
    c("--marker_logfc_cutoff"),
    type = "double",
    default = 0.5,
    help = "Filter markers: avg_log2FC >= cutoff (default: 0.5)"
  ),
  make_option(
    c("--marker_padj_cutoff"),
    type = "double",
    default = 0.05,
    help = "Filter markers: p_val_adj <= cutoff (default: 0.05)"
  ),
  make_option(
    c("--celltype_deg_logfc_cutoff"),
    type = "double",
    default = 0.25,
    help = "Filter cell-type-specific DEGs: abs(avg_log2FC) >= cutoff (default: 0.25)"
  ),
  make_option(
    c("--celltype_deg_padj_cutoff"),
    type = "double",
    default = 0.05,
    help = "Filter cell-type-specific DEGs: p_val_adj <= cutoff (default: 0.05)"
  ),
  
  # render interactive volcano plot
  make_option(
    c("--make_report"),
    action = "store_true",
    default = FALSE,
    help = "If set, generate an Rmd report template in output folder"
  ),
  make_option(
    c("--render_report"),
    action = "store_true",
    default = FALSE,
    help = "If set, render the report to HTML (implies --make_report)"
  ),
  make_option(
    c("--report_title"),
    type = "character",
    default = "DEGSmith Volcano Plots",
    help = "Title for the HTML report"
  ),
  make_option(
    c("--report_mode"),
    type = "character",
    default = "per_celltype",
    help = "What to plot: 'combined' (all cell types) or 'per_celltype'"
  ),
  make_option(
    c("--volcano_logfc"),
    type = "double",
    default = 0.25,
    help = "Log2FC threshold used for volcano significance"
  ),
  make_option(
    c("--volcano_padj"),
    type = "double",
    default = 0.05,
    help = "Adjusted p-value threshold used for volcano significance"
  ),
  make_option(
    c("--volcano_table"),
    type = "character",
    default = "specific",
    help = "Which DE table to plot in volcano: 'deg' or 'specific' (default: specific)"
  ),
  
  # misc
  make_option(
    c("--min_cells_per_group"),
    type = "integer",
    default = 10,
    help = "Minimum cells per condition within cell type (default: 10)"
  ),
  make_option(
    c("--saveRData"),
    type = "logical",
    default = TRUE,
    help = "Whether to save results as RData (default: TRUE)"
  ),
  make_option(
    c("--write_combined_workbook"),
    type = "logical",
    default = TRUE,
    help = "Whether to write a combined workbook (default: TRUE)"
  ),
  make_option(
    c("--prefix"),
    type = "character",
    default = "",
    help = "Filename prefix (default: empty)"
  )
)

opt_parser <- OptionParser(option_list = option_list, description = description_text)
opt <- parse_args(opt_parser)

## ===== check the input parameters ===== ##
required <- c("seurat_obj",
              "cell_types",
              "condition_1",
              "condition_2",
              "output_dir")
missing <- required[map_lgl(required, ~ is.null(opt[[.x]]) ||
                              !nzchar(opt[[.x]]))]
if (any(missing)) {
  print_help(opt_parser)
  stop("Missing required arguments: ",
       paste(required[missing], collapse = ", "))
}

## ===== run ===== ##
cellTypeDE(
  seurat_obj = opt$seurat_obj,
  seurat_format = opt$seurat_format,
  output_dir = opt$output_dir,
  cell_types = opt$cell_types,
  condition_1 = opt$condition_1,
  condition_2 = opt$condition_2,
  cond_col = opt$cond_col,
  celltype_col = opt$celltype_col,
  de_assay = opt$de_assay,
  de_slot = opt$de_slot,
  marker_assay = opt$marker_assay,
  marker_slot = opt$marker_slot,
  test_use = opt$test_use,
  de_min_pct = opt$de_min_pct,
  de_logfc_threshold = opt$de_logfc_threshold,
  marker_min_pct = opt$marker_min_pct,
  marker_logfc_threshold = opt$marker_logfc_threshold,
  marker_only_pos = opt$marker_only_pos,
  marker_logfc_cutoff = opt$marker_logfc_cutoff,
  marker_padj_cutoff = opt$marker_padj_cutoff,
  celltype_deg_logfc_cutoff = opt$celltype_deg_logfc_cutoff,
  celltype_deg_padj_cutoff  = opt$celltype_deg_padj_cutoff,
  min_cells_per_group = opt$min_cells_per_group,
  saveRData = opt$saveRData,
  write_combined_workbook = opt$write_combined_workbook,
  prefix = opt$prefix,
  make_report = opt$make_report,
  render_report = opt$render_report,
  report_title = opt$report_title,
  report_mode = opt$report_mode,
  volcano_logfc = opt$volcano_logfc,
  volcano_padj = opt$volcano_padj,
  volcano_table = opt$volcano_table
)
