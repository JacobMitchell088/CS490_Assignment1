# DNase Accessible Region Preprocessing
> Author: Jacob Mitchell   
> Date: 3/11/26    

## Overview

**This project preprocesses DNase‑seq accessibility data for machine learning models.**   
The script reads a BED or narrowPeak file containing accessible   
chromosomal regions and automatically determines an optimal sequence   
length **X** that balances data retention and truncation.   

Using this optimal value, the pipeline generates:   

-   **positive.txt** -- DNA sequences from accessible regions
-   **negative.txt** -- DNA sequences sampled from inaccessible genomic gaps

These files can then be used as training inputs for an AI model that
predicts chromatin accessibility.   

------------------------------------------------------------------------   

# How It Works

## 1. Input Data

The program expects a BED‑style file (such as ENCODE **narrowPeak**)   
containing genomic coordinates.   

Example format:   
   
    chr1    10000   10120   peak1   500   
    chr1    10500   10720   peak2   450   
    chr1    11000   11100   peak3   600   

Only the **first three columns** are used:   

-   Chromosome
-   Start position
-   End position

These define accessible genomic intervals.   

------------------------------------------------------------------------   

# Optimal Sequence Length Selection

Accessible regions vary in length. Machine learning models require   
sequences of equal length, so the program determines an optimal length   
**X**   

Instead of simply using the median length, the program evaluates   
multiple candidate values and scores them based on:   

-   **Number of regions retained**
-   **Amount of truncation required**
-   **Overall information preserved**

For each candidate X:   

- kept = number of intervals with length ≥ X
- trim_loss = total bases removed when truncating longer intervals
- score = (kept * X) − penalty * trim_loss

The **X with the highest score** is selected.    

This ensures:   

-   Most accessible regions are preserved
-   Excessive truncation is avoided
-   Training data size remains large

------------------------------------------------------------------------   

# Positive Dataset Creation

Accessible intervals are processed as follows:   

1.  Regions shorter than X are discarded.
2.  Regions longer than X are truncated to `[start, start + X]`.
3.  The resulting BED file is converted to sequences using:

```
bedtools getfasta
```


Example `positive.txt` output:   

    ATGCGTACGTTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC   
    GCTAGCTAGCTAGCTAACGTTAGCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA   
    TTGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA   

Properties:   

-   One DNA sequence per line
-   All sequences are length **X**
-   Only characters **A T G C**
-   All sequences converted to **uppercase**

------------------------------------------------------------------------   

# Negative Dataset Creation

Negative samples are drawn from **gaps between accessible regions**   

Steps:   

1.  Identify genomic gaps between merged accessible intervals.
2.  From each gap, sample windows of length **X**.
3.  Sampling occurs near the **middle of the gap** to avoid edges close
    to accessible regions.
4.  Negative samples are generated so that:

```
|negative| ≈ |positive|
```

(within \~5--10% difference, or a given tolerance value)    

Example `negative.txt`:   

    CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA    
    TTTGGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG    
    AACCGGTTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG    

------------------------------------------------------------------------   

# Requirements

-   Python 3
-   bedtools
-   Reference genome FASTA (e.g., GRCh38)

Example:   

    GRCh38.fa   

------------------------------------------------------------------------   

# Running the Script

Example command:   
```
python preprocess_dnase.py \
--bed experiment1.bed \
--genome /path/to/GRCh38.fa \
--outdir results
```

#### Parameters:

  Parameter       Description   
  `--bed`         Input BED or narrowPeak file   
  `--genome`      Reference genome FASTA   
  `--outdir`      Output directory   
  `--min_x`       Minimum allowed X length   
  `--neg_ratio`   Ratio of negative to positive sequences   

------------------------------------------------------------------------   

# Output Files

The program produces:    
    results/   
     ├── positive.bed   
     ├── negative.bed   
     ├── positive.fa   
     ├── negative.fa   
     ├── positive.txt   
     └── negative.txt   

#### Descriptions:   

  File           Purpose   
  `positive.bed `   Accessible regions trimmed to X   
  `negative.bed  `  Sampled inaccessible windows   
  `positive.fa  `    FASTA sequences of positives   
  `negative.fa  `    FASTA sequences of negatives   
  `positive.txt  `   Training sequences (accessible)   
  `negative.txt  `  Training sequences (inaccessible)   

------------------------------------------------------------------------   

# Data Validation

After execution, the following checks should hold:   

    wc -l positive.txt negative.txt   

Counts should be approximately equal.   

Check sequence length consistency:   

    awk '{print length($0)}' positive.txt | sort -u   

Should return a **single value = X**.   

------------------------------------------------------------------------   

# Summary

This preprocessing pipeline converts DNase accessibility data into   
balanced, fixed‑length DNA sequence datasets suitable for machine learning.  

Key features:   

-   Automatic optimal sequence length selection
-   Balanced positive and negative datasets
-   Clean DNA sequence output
-   Compatible with standard genomic tools such as bedtools
