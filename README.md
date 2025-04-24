# NYU-Langone
# Predicting Variant Effects Using Pretrained Nucleotide and Protein Language Models

This project explores the use of pretrained large language models (LLMs) to predict the functional impact of genetic variants, with a focus on non-coding mutations. Conducted as part of a semester-long collaboration with NYU Langone, the work emphasizes **zero-shot evaluation** methods without model fine-tuning.

## Project Overview

We utilize three state-of-the-art models to evaluate pathogenicity scores for nucleotide-level variants:

- **DNABERT** – A BERT-style model trained on genomic sequences
- **Nucleotide Transformer** – A transformer model for DNA sequence representation learning
- **EvoFormer / Evo2** – A protein-language model capturing evolutionary context

These models were assessed using a benchmark dataset of labeled variants (pathogenic vs. benign) under a **zero-shot setting**.

---

## Methods & Metrics

We applied the following methodology:

1. **Zero-shot scoring** of sequences using LLMs
2. **Variant scoring** using log-likelihood ratio (LLR) or delta-log-likelihood
3. **Evaluation metrics**:
   - Area Under the ROC Curve (AUC)
   - Delta / LLR-based ranking
   - Model agreement and comparative performance
