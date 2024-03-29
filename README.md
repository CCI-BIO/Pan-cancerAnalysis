<b>A gene regulatory network identifies a cuproplasia-related gene signature pan-cancer to predict patient outcome</b>

This repository includes the scripts and data of the proposed method.

The scripts are in the Script folder and include:</br>
    &nbsp;&nbsp;&nbsp;&nbsp;(1) ProposedMethod_Functions.R - Developed functions used for our method</br>
    &nbsp;&nbsp;&nbsp;&nbsp;(2) Smain.R - Main script for the method, including differential gene expression analysis and gene regulatory network construction</br>
    &nbsp;&nbsp;&nbsp;&nbsp;(3) main02.R - Main script for the method, including critical node identification and Univariate Cox analysis</br>
    &nbsp;&nbsp;&nbsp;&nbsp;(4) main03.R - Risk score model for LGG</br>
    &nbsp;&nbsp;&nbsp;&nbsp;(5) sub.R - Other analyses</br>
    &nbsp;&nbsp;&nbsp;&nbsp;(6) job folder - jobs for HPC
    
The data are in the Data folder and include:</br>
    &nbsp;&nbsp;&nbsp;&nbsp;(1) Input files</br>
    &nbsp;&nbsp;&nbsp;&nbsp;(2) Output folder - Output files for general</br>
    &nbsp;&nbsp;&nbsp;&nbsp;(3) Cancer type folders - Output files for cancer types
    
Note:</br>
    &nbsp;&nbsp;&nbsp;&nbsp;To develop our framework, we used the network control method, which is a library from the paper Liu, Y.-Y., et al. (2011). "Controllability of complex networks." Nature 473: 167. The library can be downloaded from https://scholar.harvard.edu/yyl/code. You may need to contact the authors to get the approval for the usage of the library.

Contact information:</br>
    &nbsp;&nbsp;&nbsp;&nbsp;If there are any issues with the method we encourage you to post on the issues section of this github page. If you are unable to do so please get in touch with Vu Viet Hoang Pham (vu.v.pham@unsw.edu.au).
