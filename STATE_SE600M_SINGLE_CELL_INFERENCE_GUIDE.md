# STATE SE-600M Single-Cell Inference Guide

## 1. Purpose

This document records how we set up, ran, and interpreted inference with the **STATE State Embedding (SE-600M)** model from Arc Institute.

The evaluation framework focuses on three representations:

1. **Cell embedding**
2. **Gene embedding**
3. **Contextual gene embedding**

For STATE, cell embedding and base gene embedding are confirmed in the source code. A directly exposed contextual gene embedding has **not yet been confirmed** and requires tracing the model's internal Transformer outputs.

## 2. Model and environment

- Repository: `https://github.com/ArcInstitute/state`
- Model repository: `arcinstitute/SE-600M`
- Local repository: `E:\state_fork`
- Package: `arc-state 0.11.3`
- Python: `3.11`

Model files:

```text
models/SE-600M/
├── config.yaml
├── model.safetensors
├── protein_embeddings.pt
└── se600m_epoch16.ckpt
```

The successful inference used `se600m_epoch16.ckpt` plus an explicit `protein_embeddings.pt` override.

## 3. High-level workflow

```text
Single-cell AnnData (.h5ad)
            │
            ├── X: expression matrix
            ├── var: gene information
            └── obs: cell metadata
            │
            ▼
Detect gene identifiers
            │
            ▼
Match genes to STATE protein-embedding vocabulary
            │
            ▼
Detect raw counts vs log1p expression
            │
            ▼
Construct fixed-length cell sentence
            │
            ▼
SE-600M
            │
            ├── cell representation
            └── dataset representation
            │
            ▼
Final cell embedding
```

## 4. Required input

STATE SE inference uses an **AnnData `.h5ad` file**.

| Component | Requirement | Example |
|---|---|---|
| `adata.X` | Cells × genes expression matrix | raw UMI counts |
| `adata.var.index` or `var` column | Valid gene identifiers | `TP53`, `CD3D`, `NKG7` |
| `adata.obs` | Cell metadata | cell IDs, labels, batches |
| Gene identifiers | Must overlap STATE protein embeddings | human gene symbols |
| Matrix representation | Dense or supported sparse matrix | CSR works |

Example:

```text
             TP53   CD3D   NKG7   MS4A1
Cell_001       1      8      0       0
Cell_002       0     12      2       0
Cell_003       1      0     18       0
```

STATE searches `var.index` and columns in `adata.var`, selecting the gene representation with the greatest overlap with its protein embeddings.

## 5. Raw counts vs processed expression

STATE SE has no CLI switch equivalent to scFoundation's `--pre_normalized`.

The loader attempts to detect whether input contains raw UMI counts or already-`log1p` expression. The current source uses `RAW_COUNT_HEURISTIC_THRESHOLD = 35`. If a value exceeds 35, it treats the input as raw counts. For ambiguous inputs it uses an `expm1`-based check.

For reproducible testing, raw UMI counts are a clean input choice when available.

## 6. Real dataset used

We used Scanpy PBMC3K raw data:

`data/pbmc3k/pbmc3k_state.h5ad`

| Property | Value |
|---|---:|
| Cells | 2,700 |
| Genes | 32,738 |
| Matrix | CSR sparse |
| dtype | float32 |
| Values | integer-like raw counts |
| Nonzero entries | 2,286,884 |
| Maximum stored count | 419 |

Gene symbols were available in `var.index`.

## 7. Gene mapping

The downloaded `protein_embeddings.pt` contained **19,790 STATE genes**.

PBMC3K had 32,738 genes, of which 17,855 matched:

- PBMC gene overlap: `17,855 / 32,738 = 54.54%`
- STATE vocabulary coverage: `17,855 / 19,790 = 90.22%`

Inference confirmed:

```text
Auto-detected gene column: var.index
17855 genes mapped to embedding file (out of 32738)
```

## 8. Expression processing

For raw counts STATE derives both expression proportions and log counts.

Example cell:

| Gene | Raw UMI |
|---|---:|
| CD3D | 100 |
| IL7R | 50 |
| NKG7 | 40 |
| LST1 | 10 |
| MS4A1 | 0 |

Total = 200 UMIs.

### Expression proportion

`count_expr_dist = raw_count / total_count`

| Gene | Raw | Proportion |
|---|---:|---:|
| CD3D | 100 | 0.50 |
| IL7R | 50 | 0.25 |
| NKG7 | 40 | 0.20 |
| LST1 | 10 | 0.05 |
| MS4A1 | 0 | 0 |

### log1p

STATE separately applies `log(1 + raw_count)`.

| Gene | Raw | Approx. log1p |
|---|---:|---:|
| CD3D | 100 | 4.615 |
| IL7R | 50 | 3.932 |
| NKG7 | 40 | 3.714 |
| LST1 | 10 | 2.398 |
| MS4A1 | 0 | 0 |

This path does not use scFoundation's normalize-to-10,000 preprocessing.

## 9. Cell sentence construction

STATE does not send every mapped gene into the Transformer. The checkpoint uses:

`pad_length = 2048`

Conceptually:

```text
17,855 mapped genes
        ↓
rank by expression
        ↓
construct fixed sequence
        ↓
[CLS] + up to 2,047 gene positions
        ↓
2,048-token cell sentence
```

The implementation shuffles before sorting to break expression ties. If fewer genes are expressed than available positions, unexpressed genes can be sampled to fill remaining positions. This creates a possible source of stochasticity worth testing.

For selected genes, STATE also constructs an expression feature based on:

`100 × expression_weight`

Thus the model receives gene identity/protein information together with expression information.

## 10. Protein embeddings

STATE uses pretrained protein embeddings as the starting representation for genes:

```text
Gene name
   ↓
protein_embeddings.pt
   ↓
pretrained protein representation
   ↓
STATE gene_embedding_layer
   ↓
STATE gene representation
```

## 11. Three embedding types

### 11.1 Cell embedding

A cell embedding represents the entire transcriptional state of one cell.

```text
Cell expression profile
        ↓
STATE SE-600M
        ↓
Cell vector
```

Generated with:

`state emb transform`

Successful command:

```powershell
uv run state emb transform `
  --checkpoint .\models\SE-600M\se600m_epoch16.ckpt `
  --protein-embeddings .\models\SE-600M\protein_embeddings.pt `
  --input .\data\pbmc3k\pbmc3k_state.h5ad `
  --output .\data\pbmc3k\pbmc3k_se.npy `
  --batch-size 16
```

Output:

`2700 × 2058`

The final representation contains:

```text
2048-D main cell representation
          +
10-D dataset/context representation
          =
2058-D saved representation
```

For evaluation, distinguish:

```python
cell_embedding = X[:, :2048]
dataset_embedding = X[:, -10:]
```

Potential uses include cell-type classification, clustering, cell-state separation, nearest-neighbor retrieval, PCA/UMAP, batch robustness, perturbation-state separation, and linear probes.

### 11.2 Gene embedding

A gene embedding represents the gene itself, independent of a particular cell:

```text
TP53 → one vector
EGFR → one vector
CD3D → one vector
```

STATE exposes:

`Inference.get_gene_embedding(genes)`

which passes the corresponding protein representation through:

`self.model.gene_embedding_layer(protein_embeds)`

Input required: gene names present in STATE's protein-embedding vocabulary.

Potential evaluation includes pathway membership, functional similarity, gene families, interaction networks, gene ontology, and perturbation-target relationships.

### 11.3 Contextual gene embedding

A contextual gene embedding represents a gene conditioned on a particular cell state.

```text
TP53 in T cell     → vector A
TP53 in B cell     → vector B
TP53 in tumor cell → vector C
```

This differs from the base gene embedding because the representation can vary with the surrounding transcriptional context.

A direct contextual-gene-embedding output has **not yet been confirmed** in the STATE inference interface inspected so far. We need to trace the internal model forward path to establish which per-token Transformer tensor, if any, should be treated as the contextual gene representation and how its positions map back to genes.

Current status:

| Embedding | STATE status |
|---|---|
| Cell embedding | Confirmed and generated |
| Gene embedding | Confirmed in source |
| Contextual gene embedding | Requires model-level investigation |

## 12. Embeddings side by side

| Representation | Unit | Input | STATE availability | Example use |
|---|---|---|---|---|
| Cell embedding | vector/cell | whole expression profile | Confirmed | cell identity/state |
| Gene embedding | vector/gene | gene/protein identity | Confirmed | gene relationships |
| Contextual gene embedding | vector/gene/cell | gene + cellular context | Under investigation | context-specific gene behavior |

## 13. `state emb transform` parameters

| Parameter | Purpose | Used |
|---|---|---|
| `--model-folder` | Model checkpoint directory | No in successful run |
| `--checkpoint` | Exact checkpoint | Yes |
| `--config` | Override embedded config | No |
| `--input` | Input `.h5ad` | Yes |
| `--output` | `.npy` or `.h5ad` | Yes |
| `--embed-key` | `.h5ad` embedding key | Default |
| `--protein-embeddings` | Override protein embeddings | Yes |
| `--batch-size` | Cells per forward batch | 16 |
| `--lancedb` | Optional vector database | No |
| `--lancedb-update` | Update LanceDB entries | No |
| `--lancedb-batch-size` | LanceDB write batch size | No |

## 14. Important internal parameters and behavior

| Parameter / behavior | Value / behavior |
|---|---|
| Protein vocabulary | 19,790 genes |
| `pad_length` | 2,048 |
| Main cell representation | 2,048 dimensions |
| Dataset representation | 10 dimensions |
| Final saved representation | 2,058 dimensions |
| Raw/log detection | Automatic heuristic |
| Raw-count threshold | 35 |
| Gene selection | Expression-ranked |
| Tie handling | Random shuffle before sorting |
| Unmapped genes | Excluded |
| Protein embeddings | Used for gene representations |
| Batch size | CLI-overridable |

## 15. Output formats

With `--output output.npy`, STATE writes only the embedding matrix.

With `--output output.h5ad`, STATE preserves the AnnData object and writes embeddings into `obsm` using the embedding key.

## 16. Successful PBMC3K run

```text
Input
2700 cells × 32738 genes
        ↓
17855 genes matched
        ↓
raw counts detected
        ↓
cell sentence construction
        ↓
SE-600M epoch 16
        ↓
2048-D cell representation
+ 10-D dataset representation
        ↓
2700 × 2058 float32 output
```

Runtime:

`1:16:05`

Output statistics:

```text
Shape: (2700, 2058)
dtype: float32
min:  -0.6076075
max:   0.4564924
mean: -0.00018490653
std:   0.026105504
NaN:   0
Inf:   0
```

## 17. Validation plan

### Cell embedding
1. Test repeatability across identical inference runs.
2. Separate the 2048-D cell representation from the 10-D dataset representation.
3. Run PCA/UMAP on the intended representation.
4. Evaluate known cell-type structure.
5. Test k-NN or a linear classifier.
6. Evaluate clustering metrics where labels are available.

### Gene embedding
1. Extract embeddings for known genes.
2. Compare related versus unrelated genes.
3. Test pathway and functional-neighborhood structure.
4. Evaluate known gene sets and interaction relationships.

### Contextual gene embedding
1. Inspect STATE's internal model forward path.
2. Identify per-token hidden states.
3. Verify token-to-gene mapping.
4. Determine whether the representation is before or after contextual Transformer processing.
5. Extract the same gene across multiple cell types/states.
6. Test whether context-dependent changes are biologically meaningful.

## 18. Key comparison with scFoundation

| Feature | scFoundation | STATE SE |
|---|---|---|
| Input | CSV/NPY cells × genes | AnnData `.h5ad` |
| Gene vocabulary | fixed 19,264 genes | protein-embedding vocabulary |
| Our matched vocabulary | fixed alignment | 17,855 PBMC genes matched |
| Raw/log handling | explicit `--pre_normalized` | automatic heuristic |
| Raw preprocessing | normalize to 10k + log1p | proportions + log1p |
| Sequence construction | model-specific 19,266 representation | 2048-token cell sentence |
| Cell output in our test | 3072-D | 2058-D saved output |
| Base gene embedding | supported | supported |
| Contextual gene representation | gene-level model output available | still to verify precisely |


Embedding	Status	Output
Cell embedding	Done	2048-D
Static gene embedding	Done	2048-D
Contextual gene/token embedding	Done	(1, 2047, 2048)