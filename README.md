# ✍️ Author
- This program was written for a science fair project by **Sanjay Senthilnathan of Adlai E. Stevenson High School** in 10th grade from 2025 to 2026
- Paper title: **Subgroup identification of Kawasaki Disease through a metaanalysis of transcriptome data**
# 📋 Overview
- Uses public transcriptome datasets from GEO (GSE73461, GSE73462, GSE73463, and GSE63881)
- Preprocessing, QC, log2 + Quantile Normalization on data
- Generates
  1. PCA of uncorrected merged data vs PCA of corrected merged data
  2. Silhouette and PAC scores for clusters using hierarchal clustering
  3. Delta area curve for clusters from hierarchal clustering
  4. PCA of KD + Control
  5. PCA of KD cases only colored by phenotype
  6. Heatmap of top 50 genes of cluster vs cluster
  7. Heatmap of top 10 genes of each cluster compared to every other cluster
  8. Volcano plot of cluster vs every other cluster
  9. Heatmap of top pathways for each cluster vs every other cluster
  10. Violin plot of leading edge pathway scores
  11. Dot plots of top 10 up/down regulated pathways of cluster vs healthy control
  12. Phenotype enrichment in each cluster
- Gemini AI used to help with writing, troubleshooting and commenting code
# 💻 System specs
- MacBook Pro, 14 inch, 2021 with 16 GB memory and M1 Pro chip
- Rstudio Version 2025.09.0+387 (2025.09.0+387)
# ❓ Questions
- Please direct all questions to the email address splayfern3@gmail.com 
  
