# Network Cannabis Use Characteristics
R-code to reproduce analyses described in "A network approach to relationships between cannabis use characteristics and psychopathology in the general population" by Linda T. Betz, Nora Penzel, Joseph Kambeitz.  

Code by L. Betz (linda.betz@uk-koeln.de)

There are three files:

* Code_Main_Analysis.R provides code to reproduce findings and plots reported in the main manuscript.

* Code_Supplementary_Analysis.R provides code to reproduce findings and plots reported in the supplementary file.

* edge_weights_network.csv provides the weighted adjacency matrix containing the coefficients for the Mixel Graphical Model reported in the paper as a .csv-file.

Data required for the analysis is available for public use at https://www.icpsr.umich.edu/web/ICPSR/studies/6693.

For optimal reproducibility, we use the R-package "checkpoint", with snapshot date April 1, 2020 (R version 3.6.3). Information on how this works can be found at https://mran.microsoft.com/documents/rro/reproducibility.
