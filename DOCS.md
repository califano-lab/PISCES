# Table of Contents

- [Table of Contents](#table-of-contents)
- [`Utils` Class](#utils-class)
  - [Overview](#overview)
  - [Public Methods](#public-methods)
    - [1. `init()`](#1-init)
    - [2. `compute_distance_matrix()`](#2-compute_distance_matrix)
    - [3. `do_paths_exist()`](#3-do_paths_exist)
    - [4. `create_directories()`](#4-create_directories)
  - [Example Usage](#example-usage)
- [`Converter` Class](#converter-class)
  - [Overview](#overview-1)
  - [Public Methods](#public-methods-1)
    - [1. `initialize()`](#1-initialize)
    - [2. `convert_ensembl_to_gene_names()`](#2-convert_ensembl_to_gene_names)
  - [Example Usage](#example-usage-1)
- [`Loader` Class](#loader-class)
  - [Overview](#overview-2)
  - [Public Methods](#public-methods-2)
    - [1. `initialize()`](#1-initialize-1)
    - [2. `load_data()`](#2-load_data)
  - [Example Usage](#example-usage-2)
- [`Preprocessor` Class](#preprocessor-class)
  - [Overview](#overview-3)
  - [Public Methods](#public-methods-3)
    - [1. `initialize()`](#1-initialize-2)
    - [2. `preprocess_data()`](#2-preprocess_data)
  - [Example Usage](#example-usage-3)
- [`Integrator` Class](#integrator-class)
  - [Overview](#overview-4)
  - [Public Methods](#public-methods-4)
    - [1. `initialize()`](#1-initialize-3)
    - [2. `integrate_data()`](#2-integrate_data)
  - [Example Usage](#example-usage-4)
- [`Clusterer` Class](#clusterer-class)
  - [Overview](#overview-5)
  - [Public Methods](#public-methods-5)
    - [1. `initialize()`](#1-initialize-4)
    - [2. `run_clustering()`](#2-run_clustering)
    - [3. `calc_silhouette_scores()`](#3-calc_silhouette_scores)
    - [4. `set_best_clusters()`](#4-set_best_clusters)
    - [5. `find_top_genes()`](#5-find_top_genes)
  - [Example Usage](#example-usage-5)
- [`Plotter` Class](#plotter-class)
  - [Overview](#overview-6)
  - [Public Methods](#public-methods-6)
    - [1. `initialize()`](#1-initialize-5)
    - [2. `plot_violin()`](#2-plot_violin)
    - [3. `plot_silhouette_scores()`](#3-plot_silhouette_scores)
    - [4. `plot_umap_clusters()`](#4-plot_umap_clusters)
    - [5. `plot_umap_with_labels()`](#5-plot_umap_with_labels)
    - [6. `plot_gene_heatmap()`](#6-plot_gene_heatmap)
    - [7. `plot_cluster_freq_by()`](#7-plot_cluster_freq_by)
  - [Example Usage](#example-usage-6)
- [`MetacellGenerator` Class](#metacellgenerator-class)
  - [Overview](#overview-7)
  - [Public Methods](#public-methods-7)
    - [1. `initialize()`](#1-initialize-6)
    - [2. `generate_metacell_matrices()`](#2-generate_metacell_matrices)
  - [Example Usage](#example-usage-7)
- [`Runner` Class](#runner-class)
  - [Overview](#overview-8)
  - [Public Methods](#public-methods-8)
    - [1. `initialize()`](#1-initialize-7)
    - [2. `run_aracne()`](#2-run_aracne)
    - [3. `run_viper()`](#3-run_viper)
  - [Example Usage](#example-usage-8)
- [`TableGenerator` Class](#tablegenerator-class)
  - [Overview](#overview-9)
  - [Public Methods](#public-methods-9)
    - [1. `initialize()`](#1-initialize-8)
    - [2. `generate_cell_count_summary()`](#2-generate_cell_count_summary)
  - [Example Usage](#example-usage-9)

For any additional questions or contributions, please refer to the repository’s [README.md](./README.md) or open an issue/pull request on GitHub.

# `Utils` Class

The **`Utils`** class provides a set of utility methods that perform various common tasks within the analysis pipeline. Most of these tasks revolve around directory setup, path checking, and data manipulation (like generating distance matrices). 

---

## Overview

- **Author**: *Luka Jovanović / Obradović Lab / Columbia University*
- **Created**: *January, 2025*
- **Dependencies**: `R6`, `stats`, `utils`

The `Utils` class is built on the R6 OOP system. An instance of this class is typically created only once at the beginning of the pipeline to initialize paths for subsequent modules and perform sanity checks on user-defined paths. However, you are free to instantiate it any time you need its functionality.

---

## Public Methods

### 1. `init()`

**Description**  
Initializes key directories for plotting, ARACNe, and VIPER output, ensuring they exist. Also checks if certain critical user-defined paths (like the base data path and ARACNe binary) exist.

**Usage**  
```r
utils$init(
  base_data_path,
  base_output_path,
  aracne_binary_path,
  regulator_dir_path
)
```

**Arguments**  
- `base_data_path`  
  A **character** string specifying the path where the raw data files or patient directories are located.  
  *Example:* `"/home/user/scRNAseq_data"`

- `base_output_path`  
  A **character** string specifying where you want to write all output files and folders.  
  *Example:* `"/home/user/output"`

- `aracne_binary_path`  
  A **character** string specifying the path to the ARACNe3 binary executable.  
  *Example:* `"/home/user/bin/ARACNe3_app_release"`

- `regulator_dir_path`  
  A **character** string specifying the path to the directory containing regulator .txt files.  
  *Example:* `"/home/user/regulators/human_hugo"`

**Details**  
1. Constructs subdirectories under `base_output_path`:
   - `plots`
   - `plots/reclustered_viper_plots`
   - `aracne_results`
   - `viper_results`
2. Calls [`create_directories()`](#4-create_directories) to ensure these directories exist.
3. Calls [`do_paths_exist()`](#3-do_paths_exist) to verify that `base_data_path`, `base_output_path`, `aracne_binary_path`, and `regulator_dir_path` are valid paths on the file system.

**Value / Return**  
Returns a **named list** with the following paths:
```r
list(
  plot_output_path = "path/to/plots",
  plot_viper_output_path = "path/to/plots/reclustered_viper_plots",
  aracne_output_path = "path/to/aracne_results",
  viper_output_path = "path/to/viper_results"
)
```

**Side Effects**  
- Creates directories if they do not exist.  
- Throws an error if any required path does not exist.

---

### 2. `compute_distance_matrix()`

**Description**  
Computes a **distance matrix** from a given data matrix using **1 - Pearson correlation**. This is useful in hierarchical clustering or other methods that require a distance measure rather than a similarity measure.

**Usage**  
```r
utils$compute_distance_matrix(dat_mat)
```

**Arguments**  
- `dat_mat`  
  A **matrix** (or a data frame that can be coerced to a matrix) containing gene expression or other numerical data. Rows typically represent genes (features), and columns represent samples (cells or metacells).  
  *Example:* A 10000 (genes) x 500 (samples) matrix.

**Details**  
1. If `dat_mat` is not already a matrix, it is converted via `as.matrix(dat_mat)`.  
2. The function computes the Pearson correlation matrix using `cor(dat_mat, method = "pearson")`.  
3. It then transforms this correlation matrix into a distance matrix by doing `1 - correlation`.  
4. Finally, it wraps the resulting distance values in a `dist` object (as provided by base R).

**Value / Return**  
An object of class **`dist`** suitable for clustering functions like `hclust()` or for use in other dimensionality reduction or clustering workflows.

---

### 3. `do_paths_exist()`

**Description**  
Verifies that a set of paths exist on the local file system. If any path does not exist, the method prints a message and **stops** execution with an error.

**Usage**  
```r
utils$do_paths_exist(paths)
```

**Arguments**  
- `paths`  
  A **character vector** of file/directory paths to be checked.  
  *Example:* `c("/path/to/data", "/path/to/ARACNe3_app_release", "/some/other/path")`

**Details**  
1. Each path in `paths` is normalized using `normalizePath()`, handling different OS path separators.  
2. Checks each path using `file.exists()`.  
3. Collects and prints any non-existent paths.  
4. If at least one path does not exist, the function calls `stop()`—halting execution with an error message.

**Value / Return**  
A **character vector** of non-existent paths. In most usage scenarios, this vector will be empty (and the function will continue silently). However, if non-existent paths are found, the function will raise an error instead of returning.

**Side Effects**  
- Throws an error if any path fails to exist, thus halting the pipeline until the user corrects the input paths.

---

### 4. `create_directories()`

**Description**  
Ensures that a list of directories all exist. If a directory does not exist, it is created (recursively if necessary).

**Usage**  
```r
utils$create_directories(dir_paths)
```

**Arguments**  
- `dir_paths`  
  A **character vector** specifying one or more directories to create if they do not already exist.  
  *Example:* `c("/home/user/output/plots", "/home/user/output/aracne_results")`

**Details**  
1. Iterates over each directory path in `dir_paths`.  
2. Uses `dir.exists()` to check if the directory is present.  
3. If it is missing, calls `dir.create(dir_path, recursive = TRUE)` to create it (and any missing parent directories) automatically.

**Value / Return**  
None (invisibly returns `NULL`). This is a utility method and is primarily invoked for its side effects.

**Side Effects**  
- Creates new directories on the file system.

---

## Example Usage

```r
# ---------------------------------------------------------------
# Example: Initializing the pipeline using the Utils class
# ---------------------------------------------------------------

# 1. Instantiate the Utils class
utils <- Utils$new()

# 2. Define necessary paths
base_data_path     <- "/path/to/data"
base_output_path   <- "/path/to/output"
aracne_binary_path <- "/path/to/ARACNe3_app_release"
regulator_dir_path <- "/path/to/regulators"

# 3. Initialize paths, which creates the required output subdirectories
paths <- utils$init(
  base_data_path,
  base_output_path,
  aracne_binary_path,
  regulator_dir_path
)
# 'paths' is a list containing subdirectory paths for plots, ARACNe, etc.

# 4. Compute a distance matrix from a sample matrix
dummy_matrix <- matrix(rnorm(1000), nrow = 100)  # 100 genes x 10 samples
distance_mat <- utils$compute_distance_matrix(dummy_matrix)
# 'distance_mat' is now a dist object

# 5. Check if certain paths exist
tryCatch({
  utils$do_paths_exist(c(base_data_path, "/some/random/nonexistent/path"))
}, error = function(e) {
  message("Caught error (as expected): ", e$message)
})

# 6. (Manually) create any additional directories you might need
utils$create_directories(c("/home/user/output/additional_plots"))

# 7. Optionally, clean up
rm(utils)
```
[Back to Table of Contents](#table-of-contents)


# `Converter` Class

The **`Converter`** class enables users to transform non-Read10X single-cell RNA-seq data into a format compatible with the Read10X workflow. Its main function is converting ENSEMBL IDs into gene symbols (HUGO), allowing downstream analysis tools (e.g., Seurat) to recognize and work with consistent feature names.

## Overview

- **Author**: *Luka Jovanović / Obradović Lab / Columbia University*
- **Created**: *January, 2025* 
- **Dependencies**: `R6`, `AnnotationDbi`, `org.Hs.eg.db`, `tools`  

An instance of `Converter` is typically needed only if you have non-Read10X data. It scans `.rds` files within a specified directory (`base_data_path`), attempts to convert ENSEMBL IDs to gene symbols, and then saves the processed data back as `.rds`.

---

## Public Methods

### 1. `initialize()`

**Description**  
Constructor method that stores the `base_data_path` for your data directory. This path is used when scanning for `.rds` files to convert.

**Usage**  
```r
converter <- Converter$new(base_data_path)
```

### 2. `convert_ensembl_to_gene_names()`

**Description**  
Converts ENSEMBL IDs to gene symbols (HUGO) for **all `.rds` files** located in `base_data_path`. This is particularly useful if you have single-cell data saved with ENSEMBL IDs as row names, and you need to switch to readable gene symbols for downstream tools (e.g., Seurat).

**Usage**  
```r
converter$convert_ensembl_to_gene_names()
```

**Arguments**  
_None._  

**Details**  
1. **Scan for `.rds` files** in `base_data_path`:  
   - Gathers all file paths ending with the `.rds` extension.  
2. **Attempt to read** each file using `readRDS()`.  
3. **Extract gene expression** data from the object—assumed to be in `@assayData$exprs`.  
4. **Convert ENSEMBL IDs** to gene symbols using `org.Hs.eg.db`. Rows without a corresponding SYMBOL entry are excluded.  
5. **Save the modified data** back to an `.rds` file (overwriting the original, or creating a new file with the same name).  
6. If any file fails to process (due to a malformed `.rds` or unexpected structure), the method prints a message and raises an error via `stop()`.

**Value / Return**  
- **None**. This method primarily operates by writing updated `.rds` files to disk, returning invisibly.

**Side Effects**  
- Replaces or overwrites `.rds` files with new versions containing gene symbols as row names.  
- Terminates execution (`stop()`) if any file fails to process, preventing silent errors in the pipeline.

---

## Example Usage

```r
# ---------------------------------------------------------------
# Example: Converting ENSEMBL IDs to Gene Symbols
# ---------------------------------------------------------------

# 1. Instantiate the Converter with your base data path
converter <- Converter$new("/path/to/non_read10x_data")

# 2. Run the conversion process
#    This scans .rds files in /path/to/non_read10x_data, 
#    converts ENSEMBL IDs to gene symbols, and saves
#    updated files back to the same directory.
converter$convert_ensembl_to_gene_names()

# 3. Optionally, remove or reassign the converter object
rm(converter)
```
[Back to Table of Contents](#table-of-contents)


# `Loader` Class

The **`Loader`** class is responsible for loading patient data from one or more directories, constructing Seurat objects, and attaching relevant metadata. It automates the process of reading in multiple patients’ data, whether you have it as `.rds` files or from 10X feature-barcode matrices.

## Overview

- **Author**: *Luka Jovanović / Obradović Lab / Columbia University*  
- **Created**: *January, 2025*  
- **Dependencies**: `R6`, `Seurat`  

The `Loader` class can handle two main data-loading scenarios:  
1. **`.rds` files** (which might already be partially processed).  
2. **Feature-barcode matrices** in the standard 10X format.

Each loaded dataset becomes a Seurat object. Patient metadata is automatically attached as object metadata.

---

## Public Methods

### 1. `initialize()`

**Description**  
Constructs a new `Loader` object. Stores a reference to your patient list and the base/relative paths for data so it can seamlessly locate patient files or directories.

**Usage**  
```r
loader <- Loader$new(
  patients,
  base_data_path,
  patient_data_path
)
```

**Arguments**  
- `patients`  
  A **list** of patient metadata. Each element is typically a named list or a small dictionary containing fields like `id`, `type`, etc.  
  *Example:*  
  ```r
  list(
    list(id = "Patient1", type = "Early"),
    list(id = "Patient2", type = "Late")
  )
  ```

- `base_data_path`  
  A **character** string specifying the path that contains patient directories (or `.rds` files).  
  *Example:* `"/home/user/scRNAseq_data"`

- `patient_data_path`  
  A **character** string specifying the **relative** path, within each patient’s folder, that leads to the feature-barcode matrix directory if you’re using 10X data.  
  *Example:* `"analysis/cellranger-count-default/cellranger_count_outs/filtered_feature_bc_matrix"`

**Details**  
1. **Stores** the `patients`, `base_data_path`, and `patient_data_path` internally for use by other methods (particularly `load_data()`).  
2. **Does not** immediately load any data. That step happens when you call `load_data()`.

**Value / Return**  
- Returns a **`Loader`** object (invisibly) that can then be used to load data in subsequent method calls.

**Side Effects**  
- None at this stage. No data is actually read until `load_data()` is invoked.

---

### 2. `load_data()`

**Description**  
Iterates over each patient entry in `patients`, checks for an `.rds` file named `"<patient_id>.rds"` in `base_data_path`, and if found, loads that as a Seurat object. Otherwise, it tries to load the standard 10X feature-barcode matrix under `[base_data_path]/[patient_id]/[patient_data_path]`. In either case, the method then attaches all patient metadata to the Seurat object.

**Usage**  
```r
patient_seurat_list <- loader$load_data()
```

**Arguments**  
_None._ (This method uses arguments from the class fields set during `initialize()`.)

**Details**  
1. **Scans** the `base_data_path` for each patient’s ID to see if a corresponding `<ID>.rds` file exists.  
2. If the `.rds` file is found, it loads the data via a private helper `load_rds_file()`. Otherwise, it uses `Read10X()` (by default) from the directory `[patient_id]/[patient_data_path]`.  
3. **Creates** a Seurat object with a minimal filter (e.g., `min.features = 200`, `min.cells = 50`) to remove very sparse cells.  
4. **Renames** the cells with `add.cell.id = patient$id`.  
5. **Attaches** each key-value pair in the patient metadata list (e.g., `type`, `id`, etc.) to the Seurat object’s metadata.  
6. **Returns** a list containing one Seurat object per patient.

**Value / Return**  
- A **list** of Seurat objects, each corresponding to a patient.  
- Each object has the raw counts or expression data loaded and the relevant metadata fields attached.

**Side Effects**  
- Prints status messages (e.g., “Loading patient: XYZ”) unless suppressed.  
- If any patient data fails to load, the method **stops** with an error.  
- On successful completion, you have a ready-to-use list of Seurat objects for immediate downstream analysis.

---

## Example Usage

```r
# ---------------------------------------------------------------
# Example: Loading multiple patients' data into Seurat objects
# ---------------------------------------------------------------

# 1. Define your patient metadata
patients <- list(
  list(id = "Patient1", type = "Early"),
  list(id = "Patient2", type = "Late")
)

# 2. Instantiate the Loader
loader <- Loader$new(
  patients = patients,
  base_data_path = "/path/to/data",
  patient_data_path = "analysis/cellranger-count-default/cellranger_count_outs/filtered_feature_bc_matrix"
)

# 3. Load data for each patient
patient_seurat_list <- loader$load_data()

# Now you can explore or combine these Seurat objects.
# For instance:
# patient_seurat_list[[1]]  # Seurat object for Patient1
# patient_seurat_list[[2]]  # Seurat object for Patient2

# 4. Clean up or keep the loader in your workspace
rm(loader)
```
[Back to Table of Contents](#table-of-contents)


# `Preprocessor` Class

The **`Preprocessor`** class handles a series of common data-cleaning and normalization steps on a list of Seurat objects. Its workflow includes:

1. Calculating mitochondrial gene percentages per cell.
2. Filtering cells by mitochondrial gene content and total RNA count.
3. Normalizing data with **SCTransform**.
4. Annotating cells using the **SingleR** approach with a blueprint-encode reference.

## Overview

- **Author**: *Luka Jovanović / Obradović Lab / Columbia University*  
- **Created**: *January, 2025*  
- **Dependencies**: `R6`, `Seurat`, `SingleR`, `celldex`

A typical usage scenario involves providing one or more Seurat objects (e.g., from multiple patients) with varying quality or coverage. The `Preprocessor` methods then unify and clean the data, returning improved Seurat objects ready for downstream analysis.

---

## Public Methods

### 1. `initialize()`

**Description**  
Constructs a new `Preprocessor` object, storing references to your Seurat objects as well as user-defined thresholds for filtering.

**Usage**  
```r
preprocessor <- Preprocessor$new(
  seurat_list,
  mt_threshold = 25,
  min_rna = 1000,
  max_rna = 15000,
  verbose = FALSE
)
```

**Arguments**  
- `seurat_list`  
  A **list of Seurat objects** to preprocess. Typically, these objects come from an earlier data-loading step.  
  *Example:*  
  ```r
  list(
    SeuratObject_patient1,
    SeuratObject_patient2
  )
  ```
- `mt_threshold`  
  A **numeric** threshold (as a percentage) for mitochondrial gene content. Cells with `percent.mt` above this value are discarded. Defaults to `25`.
- `min_rna`  
  A **numeric** specifying the minimum `nCount_RNA` for retaining cells. Defaults to `1000`.
- `max_rna`  
  A **numeric** specifying the maximum `nCount_RNA` for retaining cells. Defaults to `15000`.
- `verbose`  
  A **logical** indicating whether to print detailed log messages. Defaults to `FALSE`.

**Details**  
1. **Instantiates** a `Preprocessor` object with all user-defined thresholds.  
2. **Retrieves** blueprint-encode reference data from `celldex::BlueprintEncodeData()` for SingleR-based annotation.

**Value / Return**  
- A **`Preprocessor`** object used to call `preprocess_data()`.

**Side Effects**  
- None immediately. The actual preprocessing occurs in `preprocess_data()`.

---

### 2. `preprocess_data()`

**Description**  
Performs the entire preprocessing workflow on each Seurat object in `seurat_list`, including:

1. Calculating `percent.mt` (mitochondrial gene percentage).
2. Filtering cells by `percent.mt`, `nCount_RNA` minimum/maximum.
3. Running **SCTransform** for normalization.
4. Annotating cells using **SingleR** with the blueprint-encode reference.

**Usage**  
```r
processed_seurat_list <- preprocessor$preprocess_data()
```

**Arguments**  
_None._ (Uses the fields set in `initialize()`.)

**Details**  
1. **Iterates** over each Seurat object in `seurat_list`.  
2. **Calculates** mitochondrial gene percentages (`percent.mt`) via a private method.  
3. **Filters** cells that exceed `mt_threshold` or fall outside `min_rna` to `max_rna`.  
4. **Normalizes** each object with `SCTransform()`, controlling for total RNA count and mitochondrial content.  
5. **Annotates** each cell with SingleR’s predicted “blueprint_labels” and “blueprint_pvals”.  
6. Returns the updated list of Seurat objects, each now containing additional metadata.

**Value / Return**  
- A **list** of **preprocessed** Seurat objects, each having new metadata fields (`percent.mt`, `blueprint_labels`, `blueprint_pvals`, etc.).

**Side Effects**  
- Prints diagnostic messages if `verbose` is `TRUE`.  
- Raises an error (via `stop()`) if any single Seurat object fails in the process.

---

## Example Usage

```r
# ---------------------------------------------------------------
# Example: Preprocessing a list of Seurat objects
# ---------------------------------------------------------------

# Suppose we already have a list of Seurat objects (patient_seurat_list):
# patient_seurat_list <- list(seurat_obj_1, seurat_obj_2, ...)

# 1. Instantiate the Preprocessor with desired thresholds
preprocessor <- Preprocessor$new(
  seurat_list = patient_seurat_list,
  mt_threshold = 25,   # Drop cells with > 25% mt genes
  min_rna = 1000,      # Retain cells with >= 1000 RNA counts
  max_rna = 15000,     # Retain cells with <= 15000 RNA counts
  verbose = TRUE       # Print detailed logs
)

# 2. Run the preprocessing pipeline
processed_list <- preprocessor$preprocess_data()

# 3. Each element in 'processed_list' is a fully preprocessed Seurat object,
#    ready for downstream steps (integration, clustering, etc.)
# processed_list[[1]]
# processed_list[[2]]

# 4. Cleanup or keep the preprocessor for reference
rm(preprocessor)
```
[Back to Table of Contents](#table-of-contents)


# `Integrator` Class

The **`Integrator`** class handles the important task of combining multiple Seurat objects into a single, integrated dataset. It supports two primary methods for integration:

1. **Seurat-based Integration** (SCT-based workflow).
2. **FastMNN** (via **batchelor** and **SeuratWrappers**).

## Overview

- **Author**: *Luka Jovanović / Obradović Lab / Columbia University*  
- **Created**: *January, 2025*  
- **Dependencies**: `R6`, `Seurat`, `SeuratWrappers`, `batchelor`

A typical usage scenario is to bring together multiple batches or patient datasets, correct for batch effects, and create a more cohesive representation of the combined data. The choice between “seurat” and “fastMNN” integration methods can depend on user preference or performance considerations.

---

## Public Methods

### 1. `initialize()`

**Description**  
Constructs a new `Integrator` object. Users can control the integration approach—either **Seurat** or **fastMNN**—and specify key parameters like the number of features to use or the reference dataset index (for Seurat-based integration).

**Usage**  
```r
integrator <- Integrator$new(
  seurat_list,
  nfeatures = 4000,
  reference = 1,
  integration_type = "seurat",
  verbose = FALSE
)
```

**Arguments**  
- `seurat_list`  
  A **list** of Seurat objects that will be integrated.  
  *Example:*  
  ```r
  list(
    SeuratObject_patient1,
    SeuratObject_patient2,
    ...
  )
  ```
- `nfeatures`  
  A **numeric** specifying how many variable features to select for integration. Defaults to `4000`.
- `reference`  
  A **numeric** index (e.g., `1`, `2`) indicating which Seurat object should serve as the reference dataset when running Seurat-based integration. Defaults to `1`.
- `integration_type`  
  A **character** string specifying the integration method—either `"seurat"` or `"fastMNN"`. Defaults to `"seurat"`.
- `verbose`  
  A **logical** indicating whether to print detailed messages during integration. Defaults to `FALSE`.

**Details**  
- **Initialization** merely captures these parameters; the actual integration is triggered when you call `integrate_data()`.
- The choice of `reference` is relevant **only** if `integration_type` is `"seurat"`.
- If `integration_type` is `"fastMNN"`, the `reference` argument is ignored.

**Value / Return**  
A configured **`Integrator`** object.

**Side Effects**  
_None._ No integration is performed during object instantiation.

---

### 2. `integrate_data()`

**Description**  
Executes the integration process on the provided list of Seurat objects using one of the two available methods: **Seurat SCT-based** or **fastMNN**. Returns a single integrated Seurat object.

**Usage**  
```r
integrated_seurat <- integrator$integrate_data()
```

**Arguments**  
_None._ (This method relies on class fields specified at initialization.)

**Details**  
1. **Checks** the selected `integration_type`:
   - If `"seurat"`, runs a **Seurat-based** integration workflow:
     1. **SelectIntegrationFeatures** to find top variable features.  
     2. **PrepSCTIntegration** and `RunPCA` on each object.  
     3. **FindIntegrationAnchors** using `rpca` reduction.  
     4. **IntegrateData** to produce a merged Seurat object (SCT assay).  
   - If `"fastMNN"`, runs **fastMNN** using `RunFastMNN()` on normalized data from each object.  
   - If neither, raises an error.  
2. **Returns** the integrated Seurat object, which you can then further analyze (e.g., run UMAP, clustering, differential expression, etc.).

**Value / Return**  
A **Seurat object** containing the integrated assay and metadata.

**Side Effects**  
- Prints messages about the progress if `verbose = TRUE`.  
- Raises an error (`stop()`) if an unsupported integration method is selected.

---

## Example Usage

```r
# ---------------------------------------------------------------
# Example: Integrating multiple Seurat objects
# ---------------------------------------------------------------

# 1. Suppose we already have a list of preprocessed Seurat objects
#    from different patients/batches
preprocessed_seurat_list <- list(
  seurat_obj_batch1,
  seurat_obj_batch2,
  ...
)

# 2. Instantiate the Integrator for a Seurat-based integration
integrator <- Integrator$new(
  seurat_list = preprocessed_seurat_list,
  nfeatures = 3000,            # e.g., select 3000 features
  reference = 1,               # Use the first object as reference
  integration_type = "seurat", # alternatively, fastMNN by setting
  verbose = TRUE
)

# 3. Integrate the data
integrated_seurat <- integrator$integrate_data()

rm(integrator)
```
[Back to Table of Contents](#table-of-contents)


# `Clusterer` Class

The **`Clusterer`** class facilitates clustering operations on integrated Seurat objects. It allows you to:

1. Run dimensionality reductions (PCA, UMAP).
2. Perform neighbor finding and Louvain clustering across multiple resolutions.
3. Identify an optimal clustering resolution via silhouette scores.
4. Assign “best” clusters to the Seurat object.
5. Detect top marker genes for each cluster.

## Overview

- **Author**: *Luka Jovanović / Obradović Lab / Columbia University*  
- **Created**: *January, 2025*  
- **Dependencies**: `R6`, `Seurat`, `dplyr`, `cluster`, `umap`, `pheatmap`, `ggplot2`

The `Clusterer` class leverages some functionality from the [`Utils` class](#utils-class) for computing distance matrices (used in silhouette-score evaluation). It’s typically used after data integration to examine different granularity levels of clusters and to identify biologically relevant groupings.

---

## Public Methods

### 1. `initialize()`

**Description**  
Sets up a `Clusterer` object with a target Seurat object, verbosity settings, a range of resolutions to test, and the PCA dimensions to use in clustering.

**Usage**  
```r
clusterer <- Clusterer$new(
  seurat_obj,
  verbose = FALSE,
  resolutions = seq(0.1, 1, by = 0.1),
  dims = 1:50
)
```

**Arguments**  
- `seurat_obj`  
  A **Seurat** object containing integrated single-cell data.  
  *Example:*  
  ```r
  integrated_seurat
  ```
- `verbose`  
  A **logical** indicating whether to print detailed progress messages. Defaults to `FALSE`.
- `resolutions`  
  A **numeric vector** of possible resolution values for **Louvain clustering**. Defaults to `seq(0.1, 1, by = 0.1)`.
- `dims`  
  A **numeric vector** specifying which PCA dimensions to use (for example, `1:30` to use the first 30 PCs). Defaults to `1:50`.

**Details**  
- Storing multiple resolutions lets you evaluate different cluster granularities via **silhouette scores** and pick the best resolution automatically.

**Value / Return**  
- A **`Clusterer`** object ready to run the clustering pipeline on your `seurat_obj`.

**Side Effects**  
_None at this stage. Actual clustering occurs in `run_clustering()`._

---

### 2. `run_clustering()`

**Description**  
Executes a standard clustering workflow on the Seurat object. It runs **PCA**, **UMAP**, identifies neighbors, and then runs Louvain clustering at each specified resolution.

**Usage**  
```r
clustered_seurat <- clusterer$run_clustering()
```

**Arguments**  
_None._ (Uses the fields in `initialize()`.)

**Details**  
1. **RunPCA** using the variable features.  
2. **RunUMAP** to embed the data into a 2D space for visualization.  
3. **FindNeighbors** on the PCA embeddings (dimensions set in `dims`).  
4. **FindClusters** across all `resolutions`, assigning each cell to multiple cluster columns in the Seurat object’s metadata (e.g., `integrated_snn_res.0.1`, `integrated_snn_res.0.2`, etc.).

**Value / Return**  
- The **Seurat object** (same as `self$seurat_obj`), updated with:
  - PCA embeddings (`seurat_obj@reductions$pca`).
  - UMAP embeddings (`seurat_obj@reductions$umap`).
  - Neighbor graph.
  - Clustering results for each tested resolution.

**Side Effects**  
- Mutates `self$seurat_obj` in place with new dimensional reductions and cluster assignments.

---

### 3. `calc_silhouette_scores()`

**Description**  
Evaluates how well each resolution separates cells into clusters by computing **silhouette scores** across the tested resolutions (e.g., `0.1`, `0.2`, ... , `1`). Repeated subsampling of cells ensures a robust score estimate.

**Usage**  
```r
sil_results <- clusterer$calc_silhouette_scores()
best_res <- sil_results$best_resolution
mean_scores <- sil_results$mean_scores
sd_scores <- sil_results$sd_scores
```

**Arguments**  
_None._ (Uses the resolutions from initialization and the PCA dimensions via `self$dims`.)

**Details**  
1. Extracts the cluster assignments for each resolution.  
2. Subsamples cells (max 1000) multiple times (100 by default).  
3. Computes **1 - Pearson correlation** distance on the chosen PCA embeddings.  
4. Calculates silhouette widths for each resolution.  
5. Aggregates into mean ± SD silhouette scores, returning:
   - **`best_resolution`**: the resolution with the highest mean silhouette.  
   - **`mean_scores`** and **`sd_scores`**: named vectors keyed by resolution.

**Value / Return**  
A **list** with elements:
- `best_resolution`  
- `mean_scores`  
- `sd_scores`

**Side Effects**  
- None, aside from printing progress if `verbose` is on.

---

### 4. `set_best_clusters()`

**Description**  
Assigns a single “best” resolution as the **active identity** in the Seurat object, making it easier to refer to or plot a single cluster solution.

**Usage**  
```r
seurat_obj <- clusterer$set_best_clusters(best_resolution)
```

**Arguments**  
- `best_resolution`  
  A **numeric** indicating which resolution to select. Typically, you supply the `best_resolution` from `calc_silhouette_scores()`.

**Details**  
1. Looks up the appropriate metadata column (e.g., `integrated_snn_res.0.3`) based on `best_resolution`.  
2. Moves that column to `seurat_clusters` for standard usage in Seurat.  
3. Sets **`Idents`** to the same cluster IDs.

**Value / Return**  
- The updated **Seurat object**, with the chosen resolution set as the main clustering identity.

**Side Effects**  
- Overrides the `seurat_clusters` column in your Seurat object metadata.

---

### 5. `find_top_genes()`

**Description**  
Identifies the top marker genes for each cluster by performing a simple log fold-change analysis between each cluster and all other cells. Useful for quick inspection of distinctive genes in each cluster.

**Usage**  
```r
top_genes_df <- clusterer$find_top_genes(
  assay_name = "SCT",
  n_top_genes = 5,
  logfc_threshold = 0.25
)
```

**Arguments**  
- `assay_name`  
  A **character** string specifying which assay to use (e.g., `"SCT"`, `"RNA"`). Defaults to `"SCT"`.
- `n_top_genes`  
  A **numeric** specifying how many top genes to return for each cluster. Defaults to `5`.
- `logfc_threshold`  
  A **numeric** threshold for the minimum average log2 fold change. Genes must exceed this difference to be considered markers. Defaults to `0.25`.

**Details**  
1. **Extracts** scaled data from the specified assay.  
2. **Separates** cells into cluster vs. non-cluster groups for each cluster.  
3. **Computes** average logFC for each gene (mean cluster expression minus mean “all other” expression).  
4. **Filters** by `logfc_threshold` and returns the top `n_top_genes` based on the largest logFC.

**Value / Return**  
A **data frame** with columns:
- `gene`
- `cluster`
- `avg_log2FC`

**Side Effects**  
- None. It merely returns a table of candidate marker genes.

---

## Example Usage

```r
# ---------------------------------------------------------------
# Example: Clustering an integrated Seurat object
# ---------------------------------------------------------------

# 1. Instantiate the Clusterer
clusterer <- Clusterer$new(
  seurat_obj = integrated_seurat,
  verbose = TRUE,           # Print progress
  resolutions = seq(0.1, 1, 0.1),
  dims = 1:30               # Use first 30 PCs
)

# 2. Run the clustering pipeline
clustered_seurat <- clusterer$run_clustering()

# 3. Calculate silhouette scores and pick best resolution
sil_results <- clusterer$calc_silhouette_scores()
best_res <- sil_results$best_resolution

# 4. Assign best resolution clusters as active identity
integrated_seurat <- clusterer$set_best_clusters(best_res)

# 5. Find top marker genes for each cluster
top_markers <- clusterer$find_top_genes(
  assay_name = "SCT",
  n_top_genes = 5,
  logfc_threshold = 0.25
)

# 6. top_markers now has a row for each cluster’s top genes
head(top_markers)

# 7. Clean up if desired
rm(clusterer)
```
[Back to Table of Contents](#table-of-contents)


# `Plotter` Class

The **`Plotter`** class centralizes various plotting routines for visualizing clustering results from a Seurat object. It can produce:

- Violin plots of gene expression by group,
- Silhouette score diagnostic plots,
- UMAP cluster maps (with refined labeling options),
- Gene heatmaps,
- Cluster frequency comparisons across metadata categories.

## Overview

- **Author**: *Luka Jovanović / Obradović Lab / Columbia University*  
- **Created**: *January, 2025*  
- **Dependencies**: `R6`, `Seurat`, `ggplot2`, `plyr`, `pheatmap`

By delegating plotting logic to this class, you keep your main workflow scripts more concise. `Plotter` relies on a properly clustered Seurat object and, optionally, a user-defined set of cluster labels to make final figures publication-ready.

---

## Public Methods

### 1. `initialize()`

**Description**  
Constructs a `Plotter` object tied to a specific Seurat object, an output directory for saving plots, and optional cluster labels.

**Usage**  
```r
plotter <- Plotter$new(
  seurat_obj,
  plot_output_path,
  cluster_labels = NULL,
  resolutions = seq(0.1, 1, by = 0.1)
)
```

**Arguments**  
- `seurat_obj`  
  A **Seurat** object containing clustering and dimension-reduction (UMAP) results.  
- `plot_output_path`  
  A **character** specifying where to save generated plots (e.g., `"/path/to/plots"`).  
- `cluster_labels`  
  A **character vector** of labels to assign to clusters (e.g., `c("CD8 T-cell", "B-cell", ...)`).  
  If `NULL`, numeric cluster IDs are used.  
- `resolutions`  
  A **numeric vector** of resolution values used for clustering. Defaults to `seq(0.1, 1, by = 0.1)`.

**Details**  
- This simply stores references. Actual plots are generated by calling the appropriate methods.
- If `cluster_labels` is shorter (or longer) than the number of clusters, warnings may arise and any extra clusters remain unlabeled (or trimmed).

**Value / Return**  
A **`Plotter`** object for generating various cluster-related figures.

**Side Effects**  
_None._ No plots are created until you call the methods below.

---

### 2. `plot_violin()`

**Description**  
Creates a **violin plot** of one or more features across a grouping variable. Saves the result to **`plot_output_path`**.

**Usage**  
```r
plotter$plot_violin(
  features = "geneX",
  group_by = "seurat_clusters",
  pt_size = 0
)
```

**Arguments**  
- `features`  
  A **character vector** of features (e.g., gene names) to plot.
- `group_by`  
  A **metadata column** in `seurat_obj` to group cells by (e.g., `"seurat_clusters"`, `"type"`, etc.).
- `pt_size`  
  A **numeric** controlling the overlaid point size. Defaults to `0` (no points).

**Details**  
- Internally uses **`VlnPlot`** from Seurat.
- Automatically saves a PNG file named `"<feature>_violin_plot.png"` in the `plot_output_path`.

**Value / Return**  
_None._ Plot is saved to disk. In R, this is an invisible side effect.

**Side Effects**  
- A PNG file is written to **`plot_output_path`**.

---

### 3. `plot_silhouette_scores()`

**Description**  
Visualizes **silhouette scores** (mean ± standard deviation) against multiple resolutions. Highlights the “best” resolution with the highest mean silhouette.

**Usage**  
```r
plotter$plot_silhouette_scores(
  mean_scores,
  sd_scores
)
```

**Arguments**  
- `mean_scores`  
  A **numeric vector** of mean silhouette scores for each resolution.
- `sd_scores`  
  A **numeric vector** of standard deviations corresponding to `mean_scores`.

**Details**  
- Creates an **error-bar** style plot using base graphics (`errbar`) plus a line connecting the mean scores.
- Displays a **legend** with “Best Resolution = X”.
- Saves the result as `silhouette_scores.png`.

**Value / Return**  
_None._ The silhouette scores plot is saved to disk.

**Side Effects**  
- Writes a `silhouette_scores.png` file to `plot_output_path`.

---

### 4. `plot_umap_clusters()`

**Description**  
Draws a **UMAP** layout colored by cluster IDs. Optionally updates cluster labels if you provided `cluster_labels`.

**Usage**  
```r
updated_seurat <- plotter$plot_umap_clusters()
```

**Arguments**  
_None._

**Details**  
1. Checks how many clusters exist vs. how many labels were provided.  
2. Maps cluster IDs to user-specified labels if possible, or appends generic labels (e.g., `"Cluster1"`) if the user labels are insufficient.  
3. **Saves** a file named `umap_clustering_results.png` in `plot_output_path`.

**Value / Return**  
- Returns the **Seurat object** with updated `seurat_clusters` metadata if label remapping occurs.

**Side Effects**  
- Creates a UMAP cluster plot (PNG) on disk.

---

### 5. `plot_umap_with_labels()`

**Description**  
Generates a UMAP plot colored by **refined labels**—an additional metadata column named `"refined_labels"`. This can be helpful after using SingleR or other label-refinement methods.

**Usage**  
```r
plotter$plot_umap_with_labels()
```

**Arguments**  
_None._

**Details**  
1. Calls a private function to **filter** blueprint labels by p-value and frequency.  
2. Plots a labeled UMAP using `"refined_labels"`.  
3. Saves as `umap_refined_labels.png` in `plot_output_path`.

**Value / Return**  
_None._ The plot is saved to disk.

**Side Effects**  
- **Stops** with an error if `blueprint_labels` or `blueprint_pvals` are missing.  
- Writes a `umap_refined_labels.png` file.

---

### 6. `plot_gene_heatmap()`

**Description**  
Creates a **heatmap** of selected genes (optionally grouped by cluster). It can also incorporate top genes per cluster, row-wise scaling, and color annotations.

**Usage**  
```r
plotter$plot_gene_heatmap(
  genes,
  genes_by_cluster = TRUE,
  n_top_genes_per_cluster = 5,
  color_palette = NULL,
  scaled = TRUE
)
```

**Arguments**  
- `genes`  
  A **character vector** of genes to display in the heatmap.  
- `genes_by_cluster`  
  A **logical** indicating whether to group or annotate genes by cluster membership. Defaults to `TRUE`.  
- `n_top_genes_per_cluster`  
  A **numeric** controlling how many top genes to select for each cluster if `genes_by_cluster` is `TRUE`.  
- `color_palette`  
  A **character vector** of custom colors for cluster annotation. By default, uses `hue_pal` from **`scales`**.  
- `scaled`  
  A **logical** indicating if data is already scaled. If `FALSE`, row-wise z-score scaling is applied within the function.

**Details**  
1. Subsets the **SCT**-based scaled data (or raw data, then scales if `scaled=FALSE`).  
2. Orders columns by cluster, coloring columns by cluster identity.  
3. Optionally uses a row annotation if `genes_by_cluster` is `TRUE`.  
4. Generates a **pheatmap** with user-defined or auto-generated color palette.  
5. Saves the heatmap to **`gene_heatmap.png`**.

**Value / Return**  
- Returns the **pheatmap object** (invisibly).  
- Also writes the heatmap to disk.

**Side Effects**  
- May issue warnings if some requested genes are not found in the expression matrix.

---

### 7. `plot_cluster_freq_by()`

**Description**  
Visualizes how often clusters appear across different **group_by** categories (e.g., early vs. late patients). Outputs either a **dot plot** or a **box plot**.

**Usage**  
```r
plotter$plot_cluster_freq_by(
  col_names,
  plot_title,
  group_by,
  plot_type = "dot",
  binwidth = 0.01
)
```

**Arguments**  
- `col_names`  
  A **character vector** specifying column names for early/late data merges. E.g.:
  ```r
  c("Early_p1", "Early_p2", ..., "Late_p1", "Late_p2", ...)
  ```
- `plot_title`  
  A **character** string for the plot’s title.
- `group_by`  
  A **metadata column** used to define grouping (“type”, “condition”, etc.).
- `plot_type`  
  A **character** either `"dot"` or `"box"` specifying the plot style. Defaults to `"dot"`.
- `binwidth`  
  A **numeric** controlling the bin width for the dot plot. Defaults to `0.01`.

**Details**  
1. Validates required metadata columns (`id`, `seurat_clusters`, and `group_by`).  
2. Constructs a **frequency table** of cluster occurrences under each group.  
3. Merges early/late (or other categories) into a single data frame.  
4. Draws either a **dot plot** or **box plot**, saving to `plot_output_path`.

**Value / Return**  
_None._ A figure (`cluster_frequencies_by_<group_by>_<plot_type>.png`) is saved.

**Side Effects**  
- Halts with an error if the relevant metadata columns are missing.  
- Writes a PNG to disk in `plot_output_path`.

---

## Example Usage

```r
# ---------------------------------------------------------------
# Example: Using Plotter on an integrated, clustered Seurat object
# ---------------------------------------------------------------

# 1. Instantiate the Plotter
plotter <- Plotter$new(
  seurat_obj = integrated_seurat,
  plot_output_path = "/path/to/plots",
  cluster_labels = c("Cluster1", "Cluster2", "Tumor", "B-cells", "T-cells")
)

# 2. Plot a violin of geneX, grouped by seurat_clusters
plotter$plot_violin(features = "geneX", group_by = "seurat_clusters", pt_size = 0)

# 3. Suppose we have silhouette scores from the Clusterer
mean_scores <- c(0.35, 0.40, 0.42, 0.38, 0.30, 0.29, 0.25, 0.23, 0.20)
sd_scores   <- c(0.05, 0.04, 0.03, 0.06, 0.05, 0.05, 0.03, 0.02, 0.01)
plotter$plot_silhouette_scores(mean_scores, sd_scores)

# 4. Plot UMAP with cluster labels
plotter$plot_umap_clusters()

# 5. If SingleR-based refined labels exist, visualize them
plotter$plot_umap_with_labels()

# 6. Plot a heatmap of some top marker genes
top_markers <- c("CD3D", "CD3E", "MS4A1", "CD79A", "EPCAM")
plotter$plot_gene_heatmap(genes = top_markers, scaled = TRUE)

# 7. Display cluster frequencies across "type" (e.g., Early vs. Late)
col_names <- c("Early_p1","Early_p2","Late_p1","Late_p2")
plotter$plot_cluster_freq_by(col_names, plot_title = "Cluster Frequency",
                             group_by = "type", plot_type = "dot")

# 8. Cleanup
rm(plotter)
```
[Back to Table of Contents](#table-of-contents)


# `MetacellGenerator` Class

The **`MetacellGenerator`** class produces **metacell matrices** from single-cell data stored in a Seurat object. These metacell matrices are especially useful for ARACNe analysis, where they can help reduce noise from individual cells and highlight gene–gene regulatory interactions more robustly.

## Overview

- **Author**: *Luka Jovanović / Obradović Lab / Columbia University*  
- **Created**: *January, 2025*  
- **Dependencies**: `R6`, `Seurat`, `ggplot2`, `reshape2`  
- **Utilities**: Relies on an internal instance of **`Utils`** for computing distance matrices (used in the metacell construction process).

A “metacell” here is an aggregated cell profile formed by combining the expression of a cell’s nearest neighbors in order to **smooth** and **denoise** single-cell data. You can specify parameters such as how many neighbors to consider, and how large each metacell subset should be for downstream network inference.

---

## Public Methods

### 1. `initialize()`

**Description**  
Instantiates a `MetacellGenerator` object, storing references to a Seurat object and a base output path for saving results.

**Usage**  
```r
metacell_gen <- MetacellGenerator$new(
  seurat_obj,
  base_output_path
)
```

**Arguments**  
- `seurat_obj`  
  A **Seurat** object containing single-cell data (ideally with SCT counts in `seurat_obj[["SCT"]]@counts`).  
  *Example:*  
  ```r
  integrated_seurat
  ```
- `base_output_path`  
  A **character** string specifying the path where all outputs (e.g., metacell files) will be saved.  
  *Example:*  
  ```r
  "/path/to/output"
  ```

**Details**  
- Keeps a reference to the `seurat_obj` for extracting cluster assignments and expression data.  
- Uses an internal `Utils$new()` object to compute distance matrices for metacell creation.

**Value / Return**  
A **`MetacellGenerator`** object.

**Side Effects**  
_None._ No metacell generation occurs until `generate_metacell_matrices()` is called.

---

### 2. `generate_metacell_matrices()`

**Description**  
Produces one metacell matrix **per cluster** from the Seurat object. Can optionally filter out small clusters and specify how many neighbors to use for each metacell. This results in two versions of the matrix per cluster:  
1. **All cells** (metacells built from the entire cluster).  
2. **Subset** (metacells with a maximum of `sub_size` columns, if the cluster is large).

**Usage**  
```r
metacell_matrices <- metacell_gen$generate_metacell_matrices(
  out_name = "metacell",
  size_thresh = 50,
  num_neighbors = 5,
  sub_size = 200
)
```

**Arguments**  
- `out_name`  
  A **character** string used as the file prefix when saving metacell matrices. Defaults to `"metacell"`.
- `size_thresh`  
  A **numeric** specifying the minimum cluster size to process. Clusters smaller than this are skipped. Defaults to `50`.
- `num_neighbors`  
  A **numeric** controlling how many nearest neighbors to aggregate for each cell in building a metacell. Defaults to `5`.
- `sub_size`  
  A **numeric** specifying the maximum number of columns (cells) in the final “subset” matrix. If the cluster’s metacell matrix has more columns than `sub_size`, it is randomly downsampled. Defaults to `200`.

**Details**  
1. **Extracts** expression data (`["SCT"]@counts`) and cluster labels (`seurat_clusters`).  
2. **Groups** cells by cluster, skipping those below `size_thresh`.  
3. **Builds** metacells by summing the expression of each cell with its `num_neighbors` nearest neighbors (based on **1 - Pearson correlation**).  
4. **Generates** and saves two versions for each cluster:
   - “all” metacells (no column downsampling).  
   - “sub” version (downsampled if more than `sub_size` columns).  
   - Each is saved both as `.rds` and `.tsv` (ARACNe-friendly).  
5. Returns a **list** of the final subset metacell matrices (one per cluster).

**Value / Return**  
- A **list** of metacell matrices (the “subset” versions) for each cluster.

**Side Effects**  
- Multiple files are written to **`base_output_path`**:
  - `"{out_name}_clust-{i}-metaCells_all.rds"`, `"{out_name}_clust-{i}-metaCells_all.tsv"`, …  
  - `"{out_name}_clust-{i}-metaCells_sub.rds"`, `"{out_name}_clust-{i}-metaCells_sub.tsv"`, …
- Prints warnings or skips clusters that fail the size threshold or other checks.

---

## Example Usage

```r
# ---------------------------------------------------------------
# Example: Generating metacell matrices for ARACNe
# ---------------------------------------------------------------

# 1. Suppose we have an integrated, clustered Seurat object
#    with default 'seurat_clusters' and SCT counts
my_seurat <- integrated_seurat  # or any other Seurat object

# 2. Instantiate the MetacellGenerator
metacell_gen <- MetacellGenerator$new(
  seurat_obj = my_seurat,
  base_output_path = "/path/to/output"
)

# 3. Create metacell matrices:
#    - Only consider clusters >= 50 cells
#    - Each cell is aggregated with 5 neighbors
#    - Subset final matrix columns to at most 200
metacell_matrices <- metacell_gen$generate_metacell_matrices(
  out_name = "my_metacells",
  size_thresh = 50,
  num_neighbors = 5,
  sub_size = 200
)

# 4. Check the returned list of subsetted metacell matrices
names(metacell_matrices)  # Should list cluster indices or IDs
head(metacell_matrices[[1]])  # Inspect the first cluster's matrix

# 5. Cleanup if desired
rm(metacell_gen)
```
[Back to Table of Contents](#table-of-contents)


# `Runner` Class

The **`Runner`** class orchestrates the **ARACNe** and **VIPER** analyses on metacell expression data. It encapsulates the logic for:

1. Preparing expression files and regulator lists for ARACNe.
2. Executing ARACNe3 in a systematic way (optionally in parallel threads).
3. Processing ARACNe outputs into regulon objects compatible with VIPER.
4. Running VIPER to compute activity scores for regulators.

## Overview

- **Author**: *Luka Jovanović / Obradović Lab / Columbia University*  
- **Created**: *January, 2025*  
- **Dependencies**: `R6`, `Seurat`, `viper` (for creating regulons and running VIPER)

A typical usage scenario is that you’ve already generated **metacell matrices** (one per cluster) with the [`MetacellGenerator`](#metacellgenerator-class), and now you want to:  
1. Run **ARACNe** on those matrices to reconstruct gene regulatory networks.  
2. Convert the results to **regulon** objects.  
3. Run **VIPER** to get regulon activity scores, which can be integrated into your Seurat object for further analysis.

---

## Public Methods

### 1. `initialize()`

**Description**  
Constructs a new `Runner` object with all the parameters required for ARACNe–VIPER analysis, including:

- Paths to ARACNe binary, regulator files, and output directories.  
- The metacell matrix (or list of them) to be analyzed.  
- The Seurat object from which you’ll extract expression data for VIPER.

**Usage**  
```r
runner <- Runner$new(
  metacell_mat,
  base_output_path,
  aracne_binary_path,
  aracne_output_path,
  regulator_files,
  seurat_obj,
  viper_output_path,
  threads = 4,
  seed = 42
)
```

**Arguments**  
- `metacell_mat`  
  A **list** (or single object) representing the metacell matrices. Each element typically corresponds to a specific cluster’s metacell data, if you have multiple.  
- `base_output_path`  
  A **character** string for the overall output folder.  
- `aracne_binary_path`  
  A **character** string for the location of the ARACNe3 executable.  
- `aracne_output_path`  
  A **character** specifying where ARACNe outputs (networks, consolidated results) should be written.  
- `regulator_files`  
  A **named list** of regulator file paths. Example:  
  ```r
  list(
    tfs = "path/to/tfs-hugo.txt",
    cotfs = "path/to/cotfs-hugo.txt",
    ...
  )
  ```  
- `seurat_obj`  
  A **Seurat object** from which VIPER will pull expression data (in the `"integrated"` assay by default, layer = `"scale.data"`).  
- `viper_output_path`  
  A **character** specifying where VIPER outputs should be written.  
- `threads`  
  A **numeric** indicating how many parallel threads ARACNe can use. Defaults to `4`.  
- `seed`  
  A **numeric** seed for reproducibility. Defaults to `42`.

**Details**  
- Stores all parameters in fields for subsequent steps.  
- If you have multiple regulator files, ARACNe will be run for each regulator set.

**Value / Return**  
A **`Runner`** object configured to run ARACNe and VIPER.

**Side Effects**  
_None._ Actual computation happens in `run_aracne()` and `run_viper()`.

---

### 2. `run_aracne()`

**Description**  
Executes the ARACNe pipeline on your metacell matrices. For each expression file and each regulator file, ARACNe is invoked. Results (network files, `.tsv`) are saved to **`aracne_output_path`** under subfolders named after the regulator set + matrix name.

**Usage**  
```r
runner$run_aracne()
```

**Arguments**  
_None._ (Uses fields set during `initialize()`.)

**Details**  
1. **Verifies** that `metacell_mat` is not empty.  
2. **Locates** or processes expression files (`_all_all.txt.tsv` by default).  
3. **Runs** ARACNe using a system call to `aracne_binary_path`.  
4. Creates an output directory for each combination of expression file x regulator set.  
5. The user should check the console/logs for ARACNe’s progress and results.

**Value / Return**  
_None._ Results are written to disk in subfolders of `aracne_output_path`.

**Side Effects**  
- If `metacell_mat` is empty, raises an error and stops.

---

### 3. `run_viper()`

**Description**  
Conducts VIPER analysis on the ARACNe output to infer regulator activity. It generates **regulon objects** from ARACNe `.tsv` files and then calls `viper()` on them. Results are saved under `viper_output_path`.

**Usage**  
```r
runner$run_viper()
```

**Arguments**  
_None._ (Uses fields set during `initialize()`.)

**Details**  
1. **Extracts** the integrated expression matrix from `seurat_obj` (via `GetAssayData`).  
2. **Scans** `aracne_output_path` for ARACNe result files (by pattern).  
3. **Converts** each ARACNe file into a regulon object (`aracne2regulon`).  
4. **Prunes** the regulon (optional step for refinement).  
5. **Runs** `viper()` to compute regulatory activity scores.  
6. Saves the final viper scores as **`viper_results.rds`** in `viper_output_path`.

**Value / Return**  
_None._ The main result is `viper_results.rds` on disk.

**Side Effects**  
- If no ARACNe output is found, it halts with an error.  
- The user can read the results into R by:  
  ```r
  viper_res <- readRDS(file.path(viper_output_path, "viper_results.rds"))
  ```
- After re-clustering by VIPER results, you can integrate them back into your Seurat workflows.

---

## Example Usage

```r
# ---------------------------------------------------------------
# Example: Running ARACNe and VIPER on metacell data
# ---------------------------------------------------------------

# 1. Suppose we have a list of metacell matrices, each from different clusters
#    or different conditions
metacell_matrices <- list(
  cluster1_matrix,
  cluster2_matrix,
  ...
)

# 2. And we have a Seurat object for which we want to run VIPER
seurat_obj <- integrated_seurat

# 3. Regulator files to be used by ARACNe
regulator_files <- list(
  tfs = "/path/to/tfs-hugo.txt",
  cotfs = "/path/to/cotfs-hugo.txt",
  ...
)

# 4. Instantiate the Runner
runner <- Runner$new(
  metacell_mat = metacell_matrices,
  base_output_path = "/path/to/output",
  aracne_binary_path = "/path/to/ARACNe3_app_release",
  aracne_output_path = "/path/to/output/aracne_results",
  regulator_files = regulator_files,
  seurat_obj = seurat_obj,
  viper_output_path = "/path/to/output/viper_results",
  threads = 4,
  seed = 42
)

# 5. Run ARACNe
runner$run_aracne()
# ARACNe output is placed in /path/to/output/aracne_results/...

# 6. Run VIPER
runner$run_viper()
# A file 'viper_results.rds' is saved under /path/to/output/viper_results

# 7. Clean up if desired
rm(runner)
```
[Back to Table of Contents](#table-of-contents)


# `TableGenerator` Class

The **`TableGenerator`** class creates various summary tables derived from a **Seurat** object’s metadata. It currently supports generating a **cell count summary**, comparing the number of cells before and after quality control.

## Overview

- **Author**: *Luka Jovanović / Obradović Lab / Columbia University*  
- **Created**: *January, 2025*  
- **Dependencies**: `R6`, `Seurat`  

This class is useful when you want a quick tabular overview of how many cells passed certain QC criteria (e.g., minimum features, maximum counts) per patient. The tables are then saved in **CSV** format to `base_output_path`.

---

## Public Methods

### 1. `initialize()`

**Description**  
Constructs a `TableGenerator` object, holding references to your Seurat object and the output path for generated CSV files.

**Usage**  
```r
table_generator <- TableGenerator$new(
  seurat_obj,
  base_output_path
)
```

**Arguments**  
- `seurat_obj`  
  A **Seurat** object containing single-cell data. This object is expected to have a metadata column `id` indicating the patient or sample ID.  
- `base_output_path`  
  A **character** string specifying where output files (e.g., CSV summaries) will be saved.  

**Details**  
- No table generation is performed here. Call [`generate_cell_count_summary()`](#2-generate_cell_count_summary) to create the actual summary.

**Value / Return**  
A **`TableGenerator`** object.

**Side Effects**  
_None._ Initialization is strictly for storing references.

---

### 2. `generate_cell_count_summary()`

**Description**  
Builds a **cell count summary table** for each patient (`id` in the Seurat metadata). It calculates:  
1. **PreQC**: Number of cells before QC filtering.  
2. **PostQC**: Number of cells that pass two specific thresholds:  
   - `nFeature_RNA > 200`  
   - `nCount_RNA < 2500`  

The resulting table is written to a CSV file named **`cell_counts_summary.csv`** in `base_output_path`.

**Usage**  
```r
summary_df <- table_generator$generate_cell_count_summary()
```

**Arguments**  
_None._ (Uses the Seurat object and output path set in `initialize()`.)

**Details**  
1. Gathers all **unique** patient IDs from `seurat_obj@meta.data$id`.  
2. Subsets the Seurat object for each patient, counting the columns (cells) pre-QC.  
3. Applies two conditions on `nFeature_RNA` and `nCount_RNA` to determine post-QC.  
4. Assembles a data frame with columns:  
   - `Patient`  
   - `PreQC`  
   - `PostQC`  
5. Writes the data frame to `cell_counts_summary.csv` in `base_output_path`.  
6. Returns the data frame to R.

**Value / Return**  
A **data frame** with one row per `Patient`, containing columns: `Patient`, `PreQC`, and `PostQC`.

**Side Effects**  
- If the `id` column is missing from the Seurat metadata, this method may fail or produce unexpected results.  
- A CSV file, `cell_counts_summary.csv`, is created in the specified output path.

---

## Example Usage

```r
# ---------------------------------------------------------------
# Example: Generating a Cell Count Summary
# ---------------------------------------------------------------

# 1. Suppose you have a Seurat object with multiple patients
head(my_seurat@meta.data)

# 2. Instantiate the TableGenerator
table_generator <- TableGenerator$new(
  seurat_obj = my_seurat,
  base_output_path = "/path/to/output"
)

# 3. Create the cell count summary
summary_df <- table_generator$generate_cell_count_summary()

# 4. Check the resulting data frame
print(summary_df)

# 5. A 'cell_counts_summary.csv' file is now saved in /path/to/output
#    with patient-wise cell counts before/after QC.
rm(table_generator)
```
[Back to Table of Contents](#table-of-contents)