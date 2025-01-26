# Table of Contents
- [Table of Contents](#table-of-contents)
  - [Overview](#overview)
  - [Prerequisites](#prerequisites)
    - [R](#r)
      - [Windows](#windows)
      - [macOS](#macos)
      - [Linux](#linux)
      - [GCC Requirement](#gcc-requirement)
    - [R Packages](#r-packages)
    - [ARACNe3](#aracne3)
  - [Setup `main_workflow.r`](#setup-main_workflowr)
    - [Define Your Local Paths](#define-your-local-paths)
    - [Define Metadata](#define-metadata)
    - [Define Other Preferences](#define-other-preferences)
    - [Further Documentation](#further-documentation)
  - [Contributing](#contributing)
    - [How to Contribute](#how-to-contribute)
    - [Code of Conduct](#code-of-conduct)
    - [Reporting Issues](#reporting-issues)
    - [License](#license)
  - [License](#license-1)

## Overview

The **`main_workflow.r`** script provides a **comprehensive pipeline** for analyzing single-cell RNA-seq data—from raw reads all the way through to regulator-level inferences. Specifically, it:

1. **Loads and Converts Data**  
   - Handles data from either 10X feature-barcode matrices or `.rds` files (e.g., if you’ve already converted ensemble IDs to gene symbols).
   - Assigns metadata (e.g., patient ID, sample type) to each cell.

2. **Preprocesses and Normalizes**  
   - Calculates common QC metrics (e.g., mitochondrial gene content, total RNA).
   - Filters cells according to specified thresholds.
   - Uses **SCTransform** for normalization and **SingleR** for preliminary cell-type annotations.

3. **Integrates Multiple Datasets**  
   - Supports either Seurat’s SCT-based integration or a **fastMNN** approach, enabling batch correction across multiple patients or conditions.

4. **Clusters Cells**  
   - Performs PCA/UMAP and identifies clusters at multiple resolutions.
   - Finds an optimal clustering solution via **silhouette scores**.
   - Labels clusters and pinpoints key marker genes.

5. **Generates Metacells**  
   - Constructs aggregated cell profiles (“metacells”) for each cluster, reducing noise and enhancing signal for downstream regulatory network analysis.

6. **Runs ARACNe and VIPER**  
   - Utilizes **ARACNe3** to infer gene regulatory networks from the metacell expression.
   - Converts ARACNe outputs into regulon objects and applies **VIPER** to compute regulator activity scores.

Each major step is encapsulated in a dedicated R6 class—e.g., `Loader` for data loading, `Preprocessor` for QC and normalization, `Clusterer` for clustering, `Plotter` for visualization, `MetacellGenerator` for creating metacells, and `Runner` for ARACNe/VIPER. By editing the **`main_workflow.r`** script’s configuration sections (e.g., paths, metadata), you can tailor every stage of the pipeline to your specific data and experimental needs.

[Back to Table of Contents](#table-of-contents)

## Prerequisites

Before running the script, ensure you have the following installed:

### R

You need R (version 4.0 or later). Follow the instructions below to install R on your operating system.

#### Windows

1. Go to the [CRAN R project page](https://cran.r-project.org/).
2. Click on the "Download R for Windows" link.
3. Click on the "base" link to download the R installer.
4. Run the downloaded installer and follow the on-screen instructions to complete the installation.

#### macOS

1. Go to the [CRAN R project page](https://cran.r-project.org/).
2. Click on the "Download R for macOS" link.
3. Choose the appropriate package file for your macOS version and download it.
4. Open the downloaded package file and follow the on-screen instructions to complete the installation.

#### Linux

For Debian-based distributions (like Ubuntu):

1. Open a terminal.
2. Add the CRAN repository to your sources list:
   ```sh
   sudo sh -c 'echo "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" > /etc/apt/sources.list.d/R.list'
   ```
3. Add the key for the R repository:
   ```sh
   sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
   ```
4. Update your package index and install R:
   ```sh
   sudo apt update
   sudo apt install r-base
   ```

For Red Hat-based distributions (like Fedora):

1. Open a terminal.
2. Add the CRAN repository to your sources list:
   ```sh
   sudo dnf config-manager --add-repo https://cloud.r-project.org/bin/linux/fedora/R.repo
   ```
3. Install R:
   ```sh
   sudo dnf install R
   ```

#### GCC Requirement

Note: Ensure you have GCC version 9.3.1 or later installed, as some R packages may require it. To check your GCC version, run:
```sh
gcc --version
```
If you need to update GCC, you can do so by following the instructions specific to your distribution.

### R Packages

After installing R, you need several R packages. Install the required packages by running the following commands in R:

```R
# Install CRAN packages
install.packages(c("R6", "tools", "dplyr", "cluster", "umap", "pheatmap", "Hmisc", "plyr", "ggplot2", "scales", "reshape2"))

# Install Bioconductor and its packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db", "SingleR", "Seurat", "SeuratWrappers", "batchelor"))

# Install celldex for SingleR annotation reference data
BiocManager::install("celldex")

# Install viper package
BiocManager::install("viper")
```

If you encounter issues installing `SingleR` from `install.packages()`, consider using `devtools` or a related method:
```R
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("dviraran/SingleR")
```

### ARACNe3

The ARACNe3 tool is required for network inference analysis. Follow the steps below to download and set up ARACNe3:

1. Open a terminal.
2. Clone the ARACNe3 repository from GitHub:
   ```sh
   git clone https://github.com/califano-lab/ARACNe3.git
   ```
3. Navigate to the ARACNe3 directory:
   ```sh
   cd ARACNe3
   ```
4. Build the ARACNe3 executable:
   ```sh
   cd build
   cmake ..
   make
   ```
5. After building, the ARACNe3 executable will be located in `build/src/app/`. Ensure you note this path, as you will need it to configure the `aracne_binary_path` in the `main_workflow.r` script.

[Back to Table of Contents](#table-of-contents)

## Setup `main_workflow.r`

To set up the `main_workflow.r` script, you need to define several paths and preferences. Follow the instructions below to configure your script correctly.

**Note**: The `main_workflow.r` script includes detailed comments explaining each step of the process. Pay special attention to the keywords `@todo` and `@note`. The keyword `@todo` indicates that an action is required from you, while the keyword `@note` highlights important information that you should read carefully.

### Define Your Local Paths

Edit the `main_workflow.r` script to define your local paths:

1. **Base Data Path**:
   Define the path of the directory where each of the patient directories is located. For example, if you have patient directories `P1`, `P2`, ..., `Pn`, they would be located at `base_data_path/Pi`.
   ```r
   base_data_path <- "/path/to/your/data"
   ```

2. **Patient Data Path**:
   Define the path of the data directory within each patient directory. For example, if you `cd` into `P1` and find the data in `P1/data`, then `patient_data_path` would be `"data"`.
   ```r
   patient_data_path <- "relative/path/to/patient/data"
   ```

3. **Base Output Path**:
   Define the path of the directory where you want all your output to go.
   ```r
   base_output_path <- "/path/to/your/output"
   ```

4. **ARACNe3 Binary Path**:
   Define the path of your ARACNe3 binary executable on the machine you are running this script on. It must be ARACNe3, not ARACNe2.
   ```r
   aracne_binary_path <- "/path/to/ARACNe3/build/src/app/ARACNe3_app_release"
   ```

5. **Regulator Files**:
   Define paths to regulator files, currently supported only in `.txt` format.
   ```r
   regulator_dir_path <- "/path/to/regulator/files"
   regulator_files <- list(
     cotfs = file.path(regulator_dir_path, "cotfs-hugo.txt"),
     surface = file.path(regulator_dir_path, "surface-hugo.txt"),
     sig = file.path(regulator_dir_path, "sig-hugo.txt"),
     tfs = file.path(regulator_dir_path, "tfs-hugo.txt")
   )
   ```

### Define Metadata

Edit the `main_workflow.r` script to define any metadata you want:

1. **Patient Metadata**:
   Define any metadata for your patients. Ensure you include an `id` for each patient, which should be the name of the topmost directory containing that patient's data.
   ```r
   patients <- list(
     list(id = "P1", type = "Early"),
     list(id = "P2", type = "Early"),
     list(id = "P3", type = "Early"),
     list(id = "P4", type = "Late"),
     list(id = "P5", type = "Late"),
     list(id = "P6", type = "Late")
   )
   ```

### Define Other Preferences

Edit the `main_workflow.r` script to set other preferences:

1. **Verbose Messages**:
   Set to `TRUE` to display verbose messages during the analysis.
   ```r
   my_verbose <- FALSE
   ```

By following these instructions, you should be able to set up and run the `main_workflow.r` script for any dataset. If you encounter any issues or have questions, please refer to the comments in the script or seek assistance from the project maintainers.

### Further Documentation

If you want a deeper dive into the classes and methods used throughout the pipeline (e.g., `Utils`, `Loader`, `Preprocessor`, `Integrator`, etc.), please refer to our [DOCS.md](./DOCS.md). That document provides a class-by-class breakdown of public methods, usage examples, and side effects.

[Back to Table of Contents](#table-of-contents)

## Contributing

We welcome contributions to improve and enhance the `standard-workflow`. If you have suggestions, bug reports, or would like to contribute code, please follow the guidelines below.

### How to Contribute

1. **Fork the Repository**:
   - Navigate to the [GitHub repository](https://github.com/lukagolf/standard-workflow) and fork it to your own GitHub account by clicking the "Fork" button.

2. **Clone Your Fork**:
   - Clone the forked repository to your local machine:
     ```sh
     git clone git@github.com:your-username/standard-workflow.git
     cd standard-workflow
     ```

3. **Create a Branch**:
   - Create a new branch for your feature or bug fix:
     ```sh
     git checkout -b my-feature-branch
     ```

4. **Make Changes**:
   - Implement your changes in the new branch. Make sure your code follows the project's coding standards and includes appropriate comments.

5. **Commit Changes**:
   - Commit your changes with a descriptive commit message:
     ```sh
     git add .
     git commit -m "Description of the changes made"
     ```

6. **Push to GitHub**:
   - Push your changes to your fork on GitHub:
     ```sh
     git push origin my-feature-branch
     ```

7. **Create a Pull Request**:
   - Go to the original repository on GitHub and create a pull request. Provide a clear description of your changes and the problem they solve. If applicable, link any related issues.

### Code of Conduct

To maintain a positive and inclusive community, we require all contributors to adhere to the [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org/). Please read it carefully and ensure that your interactions with the community remain respectful and constructive.

### Reporting Issues

If you encounter any problems or have suggestions for improvements, please open an issue on the GitHub repository. Provide as much detail as possible, including steps to reproduce the issue, expected behavior, and any relevant screenshots or logs.

### License

By contributing to this project, you agree that your contributions will be licensed under the MIT License.

---

Thank you for your interest in contributing to the `standard-workflow`! Your support is greatly appreciated.

[Back to Table of Contents](#table-of-contents)

## License

This project is licensed under the MIT License. See the [LICENSE](https://github.com/your-repo/main_workflow.r/blob/main/LICENSE) file for more information.

[Back to Table of Contents](#table-of-contents)