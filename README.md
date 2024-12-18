# Partial-to-Partial Shape Matching with Geometric Consistency

This repository contains the code of the paper 
"Partial-to-Partial Shape Matching with Geometric Consistency", V.Ehm, M.Gao, P. Roetzer, M. Eisenberger, D. Cremers, F. Bernard. CVPR 2024.
It provides code to solve the partial-to-partial shape matching problem utilizing geometric consistency. 

It consists of two main folders: 
- `processing` (in Python): Preprocessing and Postprocessing of shapes. 
- `matlab`: Matlab code to solve the partial-to-partial shape matching approach, which is based on [1,2,3]. 


## Prerequisites for Processing
Add a conda environment 
````
conda create -n gc-ppsm python=3.8
conda activate gc-ppsm
````
Add additional requirements
````
pip install -r processing/requirements.txt
````

## Prerequisites for Partial-to-Partial Shape Matching
- YALMIP 
- Matlab
- Gurobi (for Matlab) 
- Run the setup_matlab() file (Adjust the path to the repo beforehand)

## Run the Code 
- Run `matlab/call_main.m` to solve the partial-to-partial shape matching problem on low resolution
- Run `processsing/gen_high_res_mat.py` to generate the high resolution solution based on the given low resolution solution

## CUTS24 Dataset
You can find the indices for the CUTS24 dataset [here](https://github.com/vikiehm/geometrically-consistent-partial-partial-shape-matching).

## References 
[1] Windheuser, T., Schlickwei, U., Schimdt, F. R., & Cremers, D. (2011, August). Large‐scale integer linear programming for orientation preserving 3d shape matching. In Computer Graphics Forum (Vol. 30, No. 5, pp. 1471-1480). Oxford, UK: Blackwell Publishing Ltd.

[2] Windheuser, T., Schlickewei, U., Schmidt, F. R., & Cremers, D. (2011, November). Geometrically consistent elastic matching of 3d shapes: A linear programming solution. In 2011 International Conference on Computer Vision (pp. 2134-2141). IEEE.

[3] P. Roetzer, P. Swoboda, D. Cremers, F. Bernard (2021). A Scalable Combinatorial Solver for Elastic Geometrically Consistent 3D Shape Matching. In IEEE Conference on Computer Vision and Pattern Recognition (CVPR). 2022
