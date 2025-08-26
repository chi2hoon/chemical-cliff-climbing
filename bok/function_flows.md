# Activity Cliff Analyzer - Major Function Flows

## 🔄 Main Application Flow

```mermaid
graph TD
    A[User Uploads Data] --> B[File Validation]
    B --> C{Valid SMILES?}
    C -->|No| D[Show Error & Stop]
    C -->|Yes| E[Convert IC50 to pIC50 if needed]
    E --> F[User Configures Parameters]
    F --> G[Select Detection Method]
    G --> H{Method Type}
    H -->|Basic| I[detect_activity_cliffs]
    H -->|SALI| J[detect_activity_cliffs_sali]
    H -->|k-NN| K[detect_activity_cliffs_knn]
    I --> L[Display Results]
    J --> L
    K --> L
    L --> M[User Clicks LLM Analysis]
    M --> N[Generate Prompt]
    N --> O[Call LLM]
    O --> P[Translate if needed]
    P --> Q[Display AI Interpretation]
```

## 🧬 Core Detection Algorithm Flow

### Basic Threshold Detection (`detect_activity_cliffs`)

```mermaid
graph TD
    A[Input DataFrame] --> B[Convert SMILES to Molecules]
    B --> C[Generate Morgan Fingerprints]
    C --> D[Generate Canonical SMILES]
    D --> E[Create All Pair Combinations]
    E --> F[For Each Pair i,j]
    F --> G{Identical Canonical SMILES?}
    G -->|Yes| H[Skip Pair]
    G -->|No| I[Calculate Tanimoto Similarity]
    I --> J{Similarity >= Threshold?}
    J -->|No| H
    J -->|Yes| K[Calculate Activity Difference]
    K --> L{Activity Diff >= Threshold?}
    L -->|No| H
    L -->|Yes| M[Add to Results]
    M --> N{More Pairs?}
    N -->|Yes| F
    N -->|No| O[Return Results DataFrame]
```

### SALI Ranking Detection (`detect_activity_cliffs_sali`)

```mermaid
graph TD
    A[Input DataFrame] --> B[Prepare Molecules & Fingerprints]
    B --> C[Generate Canonical SMILES]
    C --> D[Optional: Generate Scaffold SMILES]
    D --> E[For Each Pair i,j]
    E --> F{Identical Canonical?}
    F -->|Yes| G[Skip]
    F -->|No| H{Scaffold Constrained?}
    H -->|Yes| I{Same Scaffold?}
    I -->|No| G
    I -->|Yes| J[Calculate Similarity]
    H -->|No| J
    J --> K{Similarity >= Floor?}
    K -->|No| G
    K -->|Yes| L[Calculate Activity Difference]
    L --> M[Calculate SALI = ΔpIC50 / (1-similarity)]
    M --> N[Add to Candidate List]
    N --> O{More Pairs?}
    O -->|Yes| E
    O -->|No| P[Sort by SALI Descending]
    P --> Q[Take Top N Results]
    Q --> R[Return Results DataFrame]
```

### k-NN Based Detection (`detect_activity_cliffs_knn`)

```mermaid
graph TD
    A[Input DataFrame] --> B[Prepare Molecules & Fingerprints]
    B --> C[Generate Canonical SMILES]
    C --> D[Optional: Generate Scaffold SMILES]
    D --> E[Precompute Similarity Matrix]
    E --> F[For Each Molecule i]
    F --> G[Find Top-k Most Similar Neighbors]
    G --> H[For Each Neighbor j]
    H --> I{Scaffold Constrained?}
    I -->|Yes| J{Same Scaffold?}
    J -->|No| K[Skip]
    J -->|Yes| L{Identical Canonical?}
    I -->|No| L
    L -->|Yes| K
    L -->|No| M[Calculate Activity Difference]
    M --> N{Activity Diff >= Threshold?}
    N -->|No| K
    N -->|Yes| O[Add to Pairs Dictionary]
    O --> P{More Neighbors?}
    P -->|Yes| H
    P -->|No| Q{More Molecules?}
    Q -->|Yes| F
    Q -->|No| R[Convert to DataFrame]
    R --> S[Sort by Activity Difference]
    S --> T[Return Results]
```

## 🤖 AI Interpretation Flow

### LLM Integration Flow (`prompts.py`)

```mermaid
graph TD
    A[User Clicks LLM Analysis] --> B[Select Prompt Strategy]
    B --> C{Strategy Type}
    C -->|Basic| D[generate_prompt]
    C -->|Few-shot| E[generate_fewshot_prompt]
    C -->|RAG| F[generate_rag_prompt]
    D --> G[Create System + User Prompt]
    E --> H[Add Few-shot Examples]
    H --> G
    F --> I[detect_diff_type]
    I --> J[retrieve_examples]
    J --> K[Create RAG Context]
    K --> G
    G --> L[call_llm]
    L --> M{Output Language}
    M -->|Korean| N[translate_text]
    M -->|English| O[Format Output]
    N --> O
    O --> P[Display Results]
```

### RAG Retrieval Flow (`rag/retriever.py`)

```mermaid
graph TD
    A[generate_rag_prompt] --> B[detect_diff_type]
    B --> C[Create Query String]
    C --> D[Encode Query with Sentence Transformer]
    D --> E[Search Qdrant Vector Database]
    E --> F[Retrieve Top-k Similar Examples]
    F --> G[Extract Example Text]
    G --> H[Add to Prompt Context]
    H --> I[Send to LLM]
```

## 🎨 Visualization Flow

### Molecular Visualization (`utils.py`)

```mermaid
graph TD
    A[smiles_diff_to_images] --> B[Convert SMILES to RDKit Molecules]
    B --> C[Find Maximum Common Substructure]
    C --> D[Identify Differing Atoms]
    D --> E[Generate Highlighted Images]
    E --> F[Return Image Pair]
    
    G[highlight_canonical_smiles_diff] --> H[Use difflib.SequenceMatcher]
    H --> I[Identify Equal/Replace/Delete/Insert Segments]
    I --> J[Apply HTML Highlighting]
    J --> K[Return Highlighted HTML Strings]
```

## 🔧 Data Processing Flow

### File Processing (`app.py`)

```mermaid
graph TD
    A[File Upload] --> B{File Type}
    B -->|CSV| C[pd.read_csv]
    B -->|Excel| D[pd.read_excel]
    C --> E[Validate SMILES Column]
    D --> E
    E --> F[Check Each SMILES with RDKit]
    F --> G{All Valid?}
    G -->|No| H[Show Invalid SMILES Error]
    G -->|Yes| I{Has pIC50?}
    I -->|Yes| J[Continue]
    I -->|No| K{Has IC50?}
    K -->|Yes| L[compute_pIC50]
    K -->|No| M[Show Error]
    L --> J
    J --> N[Display Data Preview]
```

### Activity Cliff Detection Pipeline

```mermaid
graph TD
    A[Validated DataFrame] --> B[User Selects Method]
    B --> C[User Sets Parameters]
    C --> D[Call Detection Function]
    D --> E[Process Results]
    E --> F[Display Count]
    F --> G[Show Download Button]
    G --> H[For Each Result Pair]
    H --> I[Generate Molecular Images]
    I --> J[Display Side-by-Side]
    J --> K[Show Metrics]
    K --> L[Show Structural Difference Type]
    L --> M[Add LLM Analysis Button]
    M --> N{More Pairs?}
    N -->|Yes| H
    N -->|No| O[End]
```

## 🔍 Key Function Signatures

### Core Detection Functions

```python
# Basic threshold detection
def detect_activity_cliffs(df, sim_thres=0.85, act_thres=1.0):
    """
    Input: DataFrame with SMILES, pIC50 columns
    Output: DataFrame with activity cliff pairs
    """

# SALI ranking detection  
def detect_activity_cliffs_sali(df, top_n=200, radius=2, fp_bits=2048, 
                               sim_floor=0.5, scaffold_constrained=False):
    """
    Input: DataFrame + SALI parameters
    Output: DataFrame ranked by SALI score
    """

# k-NN detection
def detect_activity_cliffs_knn(df, k=8, act_thres=1.0, radius=2, 
                              fp_bits=2048, min_sim=0.0, scaffold_constrained=False):
    """
    Input: DataFrame + k-NN parameters  
    Output: DataFrame with nearest neighbor cliffs
    """
```

### AI Integration Functions

```python
# Prompt generation
def generate_prompt(row):
    """Input: Activity cliff pair row, Output: (system_prompt, user_prompt)"""

def generate_fewshot_prompt(row):
    """Input: Activity cliff pair row, Output: Prompt with examples"""

def generate_rag_prompt(row, client):
    """Input: Activity cliff pair + Qdrant client, Output: RAG-enhanced prompt"""

# LLM calling
def call_llm(system_prompt, user_prompt):
    """Input: Prompts, Output: LLM response"""

def translate_text(text, target_language="Korean"):
    """Input: English text, Output: Translated text"""
```

### Utility Functions

```python
# Molecular processing
def compute_pIC50(ic50_nM):
    """Input: IC50 in nM, Output: pIC50 value"""

def smiles_diff_to_images(smiles1, smiles2, size=(300, 300)):
    """Input: Two SMILES, Output: Highlighted molecular images"""

def highlight_canonical_smiles_diff(a, b):
    """Input: Two canonical SMILES, Output: HTML-highlighted differences"""

# Structural analysis
def detect_diff_type(mol1, mol2):
    """Input: Two RDKit molecules, Output: Difference type classification"""
```

## 📊 Data Flow Summary

1. **Input**: CSV/Excel with SMILES + IC50/pIC50
2. **Validation**: SMILES validation, data conversion
3. **Detection**: One of three algorithms processes all pairs
4. **Results**: Structured DataFrame with cliff pairs
5. **Visualization**: Molecular images and highlighted differences
6. **AI Analysis**: On-demand LLM interpretation
7. **Output**: Downloadable results and AI insights

This flow ensures robust data processing, flexible detection strategies, and comprehensive analysis capabilities for activity cliff identification in drug discovery.

