# RASDF_GUI
The program of a reliable and adaptive spatiotemporal data fusion method (RASDF), developed by Dizhou Guo.

This program according to the paper Shi W, Guo D*, Zhang H. A reliable and adaptive spatiotemporal data fusion method for blending multi-spatiotemporal-resolution satellite images[J]. Remote Sensing of Environment, 2022,268:112770. https://doi.org/10.1016/j.rse.2021.112770.

The authors greatly thank the great work of FSDAF, SFSDAF and OBSTFM from Xiaolin Zhu, Eileen H. Helmer, Feng Gao, Desheng Liu, Jin Chen, Michael A. Lefsky, Xiaodong Li, Giles M. Foody, Doreen S. Boyd, Yong Ge, Yihang Zhang, Yun Du, Feng Ling and Yue Sun.

Please see the Readme.docx for instructions.


------------------------------------------------------------------------------------------
Update history：

10/16/2021 First version of RASDF GUI developed by MATLAB 2020b GUI

01/10/2022 Modify and redesign the interface on the MATLAB 2020b APP designer for ease of use

01/18/2022 Add TTIF format support and correct a bug in change detection

01/20/2022 Correct a bug caused by too few candidate pixels for unmixing (If the number of pixels available for unmixing is not much larger than the number of classifications, the unmixing step is not performed).
