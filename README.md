# Network of Cannabis Use Characteristics & Psychopathology
R-code to reproduce analyses described in "A network approach to relationships between cannabis use characteristics and psychopathology in the general population" by Linda T. Betz, Nora Penzel, Joseph Kambeitz.  

Code by L. Betz (linda.betz@uk-koeln.de)

There are three files:

* Code_Main_Analysis.R provides code to reproduce results and plots reported in the main manuscript.

* Code_Supplementary_Materials.R provides code to reproduce results and plots reported in the supplementary materials.

* edge_weights_network.csv provides the weighted adjacency matrix containing the coefficients for the Mixed Graphical Model reported in the paper as a .csv-file.

Data required for the analysis (National Comorbidity Survey Data; NCS-1) are available for public use at https://www.icpsr.umich.edu/web/ICPSR/studies/6693.

For optimal reproducibility, we use the R-package "checkpoint", with snapshot date August 1, 2021 (R version 4.1.0). Information on how this works can be found at https://mran.microsoft.com/documents/rro/reproducibility.
