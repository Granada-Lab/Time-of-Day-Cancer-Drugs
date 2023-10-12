# Circ-TNBC
Collab with Charite about circadian properties of TNBC cell lines

When installing the required packages, be careful that you only have 'umap-learn' installed, not the single-named 'umap'
package. If so please remove the umap package as it will conflict with the import call reserved for umap-learn, namely
'import umap'.

Proposed work plan:

- Download all data relate to the project: 
  - Carolin's data table :white_check_mark:
  - Gene sets definitions
    - Extend genes sets (STRING)
  - CCLE mutation data (check which) :white_check_mark:
  - CCLE transcriptomics data (check which) :white_check_mark:

- Import of data in proper form:
- Preprocess omics data :white_check_mark:
- Extract data for our cell lines :white_check_mark:

- Clustering (and some parameter analysis):
  - PCA expression + color by rythmicity :white_check_mark:
  - PCA mutations + color by rythmicity
  - UMAP expression + color by rythmicity :white_check_mark:
  - UMAP mutations + color by rythmicity
  - tSNE expression + color by rythmicity :white_check_mark:
  - tSNE mutations + color by rythmicity
  - LDA expression + color by subtypes only :white_check_mark:
  - LDA mutations + color by subtypes only
  
- Relogio clustering:
  - hierarchical with pearson as distance:
    - try the different gene sets
  - gene selection as in paper
    - differential expression
    - clustering

- Sort features by correlation (pearson):
  - transcriptomics
  - genomics

- Build predictors:
  - leave-one-out best decision tree -> N trees
