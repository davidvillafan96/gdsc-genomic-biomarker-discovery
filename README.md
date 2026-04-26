# 🧬 GDSC Project 4 — Genomic Integration & Biomarker Discovery

## 📌 Overview

This project extends a machine learning pipeline for drug response prediction by integrating genomic mutation data to uncover biological drivers of sensitivity and resistance.

The workflow combines pharmacogenomics (GDSC) with genomic mutation data (COSMIC) to move beyond prediction and enable **biomarker discovery and biological interpretation**.

---

## 🎯 Objectives

- Integrate genomic mutation data into drug response models
- Identify biomarkers associated with drug resistance and sensitivity
- Improve model interpretability using statistical and ML approaches
- Explore mutation-driven biological mechanisms

---

## 📊 Data Sources

- **GDSC (Genomics of Drug Sensitivity in Cancer)**
  - Drug response data (LN_IC50)
  - Cancer types and pathways

- **COSMIC (Catalogue Of Somatic Mutations In Cancer)**
  - Genomic mutation profiles
  - Binary mutation encoding (0 = WT, 1 = Mutated)

---

## ⚙️ Environment

- Language: **R**
- Execution: **Google Colab (R runtime)**
- Libraries:
  - tidyverse
  - ggplot2
  - glmnet
  - caret

---

## 🧠 Methodology

### 1. Data Integration

- Merge GDSC drug response with COSMIC mutation data
- Align samples using shared identifiers
- Generate final feature matrix including:
  - Genomic mutations
  - Cancer type
  - Target pathways

---

### 2. Exploratory Analysis

- Gene-level effect analysis
- Example: KRAS mutation impact on LN_IC50
- Visualization:
  - Boxplots (WT vs Mutated)
  - Statistical validation (p-values)

---

### 3. Biomarker Discovery

- Logistic regression to compute:
  - Odds Ratios (OR)
  - Adjusted p-values

- Volcano plots:
  - X-axis: Log Odds Ratio
  - Y-axis: -log10(p-value)

- Identification of key biomarkers

---

### 4. Mutation Signatures

- Detection of:
  - Co-resistance patterns
  - Co-sensitivity relationships

- Network-based interpretation of mutation interactions

---

### 5. Feature Selection (LASSO)

- Cross-validated LASSO regression
- Identification of optimal lambda
- Selection of ~54 predictive features

---

### 6. Model Improvement

- Comparison of models:
  - Without genomic data
  - With genomic integration

- Evaluation metric:
  - RMSE (regression performance)

---

## 📊 Key Results

### 🔬 Top Biomarkers (Resistance)

- KRAS
- TPR
- GLI1
- RBM10
- TP63

---

### 📈 Key Insights

- KRAS mutation significantly alters drug response distribution
- Multiple genes show strong statistical significance
- Odds Ratio > 1 → increased probability of resistance
- Genomic integration improves model performance

---

## 🧬 Biological Interpretation

- Drug response is driven by:
  - Genomic mutations
  - Pathway interactions
  - Cancer-specific context

- Evidence of mutation-driven resistance mechanisms

---

## 🏥 Real-World Impact

- Biomarker-based therapy selection
- Precision oncology applications
- Reduced experimental screening cost
- Foundation for clinical decision-support systems

---

## ⚠️ Limitations

- No transcriptomic or epigenomic data
- Biological noise remains high
- Requires validation in clinical cohorts

---

## 🚀 Future Work

- Multi-omics integration
- Deep learning models for interaction effects
- Clinical validation pipelines
- Deployment as a decision-support tool

---

## 👨‍🔬 Author

David Villafañe  
PhD in Biological Sciences  
Bioinformatics & Machine Learning
