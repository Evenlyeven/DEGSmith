<p align="center">
  <img src="hex.png" width="180"/>
</p>

# DEGSmith: Cell-Type–Specific Differential Expression Pipeline for Seurat

**DEGSmith** is a command-line R pipeline for performing cell-type–resolved differential expression (DE) analysis from Seurat objects. It supports two comparison modes:

- **`condition` mode** — compares two experimental conditions within each cell type
- **`population` mode** — directly compares two cell populations against each other

For each selected cell type (in `condition` mode), it computes:  
	1.	**Within-cell-type DE**  
(cell type, condition_1) vs (cell type, condition_2)  
	2.	**Cell-type markers**  
(cell type) vs (all other cell types)  
	3.	**Cell-type–specific DEGs**  
Overlap between DEGs and high-confidence markers  

In `population` mode, it performs a single direct DE comparison between two named populations, along with markers for each population vs rest.

Results are exported as Excel files and optionally rendered into interactive Plotly volcano HTML reports.

## 📦 Features

	•	Accepts Seurat objects in `.rds`, `.RData`, or `.qs2` format
	•	Two comparison modes: `condition` (within cell type) and `population` (cross-population)
	•	Performs within-cell-type DE using `FindMarkers()`
	•	Identifies cell-type markers (vs rest)
	•	Derives cell-type–specific DEGs via marker overlap
	•	Automatically filters low-cell-count groups
	•	Exports per-cell-type Excel files
	•	Generates a combined multi-sheet workbook
	•	Optional gene annotation table with automatic ZFIN hyperlink generation
	•	Optionally produces interactive volcano HTML reports
	•	Click any point on the volcano plot to copy the gene name to clipboard
	•	Fully command-line configurable
	•	Timestamped output directory (prevents overwriting)

## 🚀 Quick Start

### 🔧 Requirements

Install required R packages:

```r
install.packages(c(
  "tidyverse",
  "Seurat",
  "optparse",
  "magrittr",
  "openxlsx",
  "rmarkdown",
  "plotly"
))
```

If using `.qs2` input:
```r
install.packages("qs2")
```

If using `--annotation_table`:
```r
install.packages("readxl")
```

### 🖥️ Basic Usage

**Condition mode** — DE between two conditions within each cell type:
```bash
Rscript DEGSmith.R \
  --seurat_obj path/to/seurat_object.rds \
  --output_dir results/ \
  --cell_types T_cell,B_cell \
  --condition_1 Treated \
  --condition_2 Ctrl
```

Use all cell types:
```bash
Rscript DEGSmith.R \
  --seurat_obj path/to/seurat_object.qs2 \
  --output_dir results/ \
  --cell_types all \
  --condition_1 Treated \
  --condition_2 Control
```

**Population mode** — direct DE between two cell populations:
```bash
Rscript DEGSmith.R \
  --seurat_obj path/to/seurat_object.rds \
  --output_dir results/ \
  --comparison_mode population \
  --population_1 T_cell \
  --population_2 B_cell
```

Generate interactive volcano report:
```bash
Rscript DEGSmith.R \
  --seurat_obj path/to/seurat_object.qs2 \
  --output_dir results/ \
  --cell_types all \
  --condition_1 Treated \
  --condition_2 Control \
  --make_report \
  --render_report
```

### 📝 Parameters

**Required (condition mode)**
| Parameter | Description |
|-----------|-------------|
| `--seurat_obj` | Path to Seurat object file (`.rds`, `.RData`, `.qs2`) |
| `--cell_types` | Comma-separated list of cell types or `"all"` |
| `--condition_1` | Condition A label used for DE comparison |
| `--condition_2` | Condition B label used for DE comparison |
| `--output_dir` | Directory where output files will be saved |

**Required (population mode)**
| Parameter | Description |
|-----------|-------------|
| `--seurat_obj` | Path to Seurat object file (`.rds`, `.RData`, `.qs2`) |
| `--output_dir` | Directory where output files will be saved |
| `--comparison_mode` | Set to `population` |
| `--population_1` | Name of cell population 1 (value in `--celltype_col`) |
| `--population_2` | Name of cell population 2 (value in `--celltype_col`) |

**Comparison Mode**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--comparison_mode` | `condition` | Comparison mode: `condition` or `population` |
| `--population_1` | `NULL` | Cell population 1 for population-vs-population DE |
| `--population_2` | `NULL` | Cell population 2 for population-vs-population DE |

**Input Format**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--seurat_format` | `auto` | Input format: `auto`, `rds`, `RData`, or `qs2` |

**Metadata Columns**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--cond_col` | `orig.ident` | Metadata column containing condition labels |
| `--celltype_col` | `annotation` | Metadata column containing cell type annotations |

**Within-Cell-Type DE**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--test_use` | `wilcox` | Statistical test used in `FindMarkers()` |
| `--de_assay` | Seurat default | Assay used for within-cell-type DE |
| `--de_slot` | Seurat default | Slot used for within-cell-type DE |
| `--de_min_pct` | `0.01` | Minimum percent expression required for testing |
| `--de_logfc_threshold` | `0.1` | Log2 fold-change threshold for testing |

**Marker Detection (Cell Type vs Rest)**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--marker_assay` | `SCT` | Assay used for marker detection |
| `--marker_slot` | `data` | Slot used for marker detection |
| `--marker_only_pos` | `TRUE` | If TRUE, only positive markers are returned |
| `--marker_min_pct` | `0.01` | Minimum percent expression for marker detection |
| `--marker_logfc_threshold` | `0.1` | Log2FC threshold for marker testing |
| `--marker_logfc_cutoff` | `0.5` | Log2FC cutoff used to retain high-confidence markers |
| `--marker_padj_cutoff` | `0.05` | Adjusted p-value cutoff for marker filtering |

**Specific DEG Filtering**
| Parameter                        | Default | Description |
|-----------------------------------|---------|-------------|
| `--celltype_deg_logfc_cutoff`     | `0.25`  | Minimum absolute log2FC for cell-type–specific DEGs |
| `--celltype_deg_padj_cutoff`      | `0.05`  | Adjusted p-value cutoff for cell-type–specific DEGs |
| `--min_cells_per_group`           | `10`    | Minimum cells per condition within a cell type |

**Gene Annotation**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--annotation_table` | `NULL` | Path to `.xlsx` annotation file with a required `gene` column; any additional columns are appended to all DEG/marker tables. If a `zfin_id` column is present, clickable ZFIN hyperlinks are generated automatically. |

**Volcano Report**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--make_report` | `FALSE` | Generate volcano Rmd template |
| `--render_report` | `FALSE` | Render interactive HTML volcano report |
| `--report_title` | `DEGSmith Volcano Plots` | Title used in the HTML report |
| `--report_mode` | `per_celltype` | Report mode: `per_celltype` or `combined` |
| `--volcano_logfc` | `0.25` | Log2FC threshold used for volcano classification |
| `--volcano_padj` | `0.05` | Adjusted p-value threshold used for volcano classification |
| `--volcano_table` | `specific` | Which table to plot: `deg` or `specific` |

**Output Controls**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--saveRData` | `TRUE` | Save results as `DEGSmith_results.RData` |
| `--write_combined_workbook` | `TRUE` | Write combined multi-sheet Excel workbook |
| `--prefix` | `""` | Prefix added to all output filenames |

## 📂 Output Structure

Each run creates a timestamped directory: `DEGSmith_20260212_153045/`  
<br>
**Per Cell Type (condition mode)**  
	•	DEG_<celltype>_<cond1>_vs_<cond2>.xlsx  
	•	Markers_<celltype>_vs_rest.xlsx  
	•	DEG_specific_<celltype>_<cond1>_vs_<cond2>.xlsx  
	•	DEG_specific_filtered_<celltype>_...xlsx  
  <br>
**Population mode**  
	•	DEG_<pop1>_vs_<pop2>.xlsx  
	•	Markers_<pop1>_and_<pop2>_vs_rest.xlsx  
	•	DEG_specific_<pop1>_vs_<pop2>.xlsx  
	•	DEG_specific_filtered_<pop1>_vs_<pop2>_...xlsx  
  <br>
**Combined workbook sheets**  
	•	<celltype>_DEG  
	•	<celltype>_MarkersKeep_1  
	•	<celltype>_MarkersKeep_2  
	•	<celltype>_SpecificDEG  
	•	<celltype>_SpecificDEG_Filtered  
  <br>
**Optional**  
	•	DEGSmith_results.RData  
	•	Interactive volcano HTML report(s)  

## 📌 Notes

 - The Seurat object must contain the specified metadata columns.  
 - Only `SCT` assay is supported for marker detection by default.  
 - Cell types with insufficient cells per condition are automatically skipped.  
 - If multiple Seurat objects are present in an `.RData` file, the first one is used.  
 - Volcano plots are interactive (Plotly) with hover tooltips. **Click any point to copy the gene name to clipboard** — a confirmation toast will appear briefly.  
 - In `population` mode, `--cell_types`, `--condition_1`, and `--condition_2` are not required.  
 - Each run creates a new timestamped output folder to prevent overwriting.  
