% TODO adjust correct path here 
cd ~/Development/gc-ppsm/
addpath(genpath('./matlab/src/constraints/'))
addpath(genpath('./matlab/src/energyComputation/'))
addpath(genpath('./matlab/src/solveShapeMatch/helper/'))
addpath('./matlab/src/import/')
addpath('./matlab/src/solveShapeMatch/')
addpath('./matlab/src/')
addpath('./matlab/')
addpath('./matlab/src/findMeshHoles/')
addpath('./')

addpath(genpath('~/YALMIP')) 
addpath(strjoin(strcat(['~/Repos/gptoolbox/'],{'external','imageprocessing', 'images', 'matrix', 'mesh', 'mex', 'quat','utility','wrappers'}),':'))
addpath ~/gurobi1001/linux64/matlab/         