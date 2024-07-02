function [] = main(curr_animal, num_faces)
close all

double_close = 0;
area_pres = 0;
partial_to_partial = false;
test_vis = false;
coarse_to_fine = true;
closedBound = false;
closed = 0;
warm_start = false;
use_radius = true;
degenerate = true;
relaxed = true;
dummy_elements = false;
curr_type = "holes";
manifold = false;

lambda_mem = 0;
disp(strcat("Starting main with ",curr_animal, " and lambda", string(lambda_mem)));
yalmip('clear')
stk = dbstack;
filepath = which(stk(1).file);
[curr_path,name,~] = fileparts(filepath);
basepath = "/usr/wiss/ehm/Development"
fileinfo = dir(strcat(basepath, "/Unsupervised-Learning-of-Robust-Spectral-Shape-Matching/data/SHREC16_closed/", curr_type, "/off/*",curr_animal ,"*.off"));

fnames = {fileinfo.name};
not_used_shapes = readcell("not_manifold_shapes.csv");
for i=1:size(not_used_shapes,2)
    fnames = strrep(fnames,not_used_shapes(i),'');
end
fnames = fnames(~cellfun('isempty',fnames));
size(fnames)

for i=1:size(fnames,2)
    shape_1 = curr_animal;
    shape_2 = fnames{i};
    splited_shape = split(shape_2,'.');
    shape_2 = splited_shape{1};
    disp(strcat("Start matching ",curr_animal, " with ", shape_2));
    if double_close
        double_close_path = strcat(basepath, "/Unsupervised-Learning-of-Robust-Spectral-Shape-Matching/data/SHREC16_double_closed/",string(num_faces),"_faces_", curr_type, "/");
        feature_folder = "features/";
        file_1 = strcat(double_close_path,curr_animal, '.off');
        file_2 = strcat(double_close_path,shape_2, '.off');

    else
        file_1 = strcat(basepath, "/Unsupervised-Learning-of-Robust-Spectral-Shape-Matching/data/SHREC16_closed/null/off/",curr_animal, '.off');
        file_2 = strcat(basepath, "/Unsupervised-Learning-of-Robust-Spectral-Shape-Matching/data/SHREC16_closed/", curr_type, "/off/",shape_2, '.off');
    end

    [VX,FX]=get_mesh(file_1);
    [VY,FY]=get_mesh(file_2);

    if double_close
        Feats = load(strcat(double_close_path, feature_folder, curr_animal,"_", shape_2, '_feat.mat'));
        Feata = Feats.Feata;
        Featb = Feats.Featb;
    else
        featInfo = load(strcat(basepath, "/Unsupervised-Learning-of-Robust-Spectral-Shape-Matching/results/shrec16_",curr_type, "_fixed/visualization/", curr_animal, "-", shape_2, ".mat"));
    
        FeatAHr = featInfo.featX;
        FeatBHr = featInfo.featY;

        % if not manifold we need to load extra info and reduce only part of the mesh
        if ~manifold
            [VX,~,VY,~, FX,FY, FeatBHr, idx_VX, idx_VY, idx_FX,idx_FY, added_verts] = get_non_manifold_mesh(shape_1, shape_2, num_faces, FeatBHr, basepath);
        else
            [VY, FY,idx_FY,idx_VY] = decimate_libigl(VY,FY,area_partial * num_faces, 'Method', 'naive');
            [VX, FX,idx_FX,idx_VX] = decimate_libigl(VX,FX,num_faces, 'Method', 'naive');

        end
        area_partial = sum(getTriangleAreas(VY, FY));
        mesh = surfaceMesh(VY,FY);
        assert(isEdgeManifold(mesh,true))
        assert(isVertexManifold(mesh))
        assert(isOrientable(mesh))
        if use_radius
           [Feata,Featb] = get_low_res_features(VX,FX,VY,FY,FeatAHr, FeatBHr, idx_VX, idx_VY, curr_animal, shape_2, curr_type, basepath, manifold, added_verts);
            if test_vis    
                vis_features(VX,FX,VY,FY,Feata,Featb);
            end
        else
            Feata = FeatAHr(idx_VX,:);
            Featb = FeatBHr(idx_VY,:);
        end
    end
    
    FY_old = FY;

    non_matched_FY = []; 
    VY_old = VY;


    VX = VX - min(VX);
    VY = VY  - min(VX);
    translation= 0 * [0 1 0];
    VY = VY + translation;
    VY_old = VY_old + translation;

    %lambda_mem = 0;
    % solver options
    %warm_start = load("G_new_pos.mat");
    % might be a good idea to set a time limit sometimes :)
    %timeLimit = 60 * 10; % seconds
    % options.cplex.timelimit = timeLimit;plo
    % options.mosek.MSK_DPAR_OPTIMIZER_MAX_TIME = timeLimit;
    %options.gurobi.TimeLimit = timeLimit;

    shape1 = dictionary('V', {VX}, 'F', {FX}, 'features', {Feata});
    shape2 = dictionary('V', {VY}, 'F', {FY}, 'features', {Featb});

    %filename = 'test_400'
    %keys = ["lambda", "degenerate", "relaxed", "partial-partial", "closedBound", "coarseToFine"];
    %values = [lambda_mem, degenerate, relaxed, partial_to_partial, closedBound,coarse_to_fine];
    keys = ["lambda", "degenerate", "relaxed", "partial-partial", "closedBound", "warmStart", "doubleClosed", "areaPreserve"];
    values = [lambda_mem, degenerate, relaxed, partial_to_partial, closedBound, warm_start, double_close, area_pres];
    disp("Using Gurobi as Solver")
    options = sdpsettings('solver', 'gurobi', 'cachesolver',1, 'verbose', 1, 'debug',1, 'gurobi.Method',3);


    configMatching = dictionary(keys,values);
    curr_filename = strcat(shape_1,'_', shape_2, '_', string(num_faces), '_faces');
    for j = 1:size(configMatching.keys,1)
        curr_key = configMatching.keys{j};
        curr_value = configMatching(curr_key);
        curr_filename = strcat(curr_filename, '_', curr_key, '_', string(curr_value));

    end

    filename = convertStringsToChars(strcat(datestr(now,'mm-dd-yyyy-HH-MM'),curr_filename));
    [G, times, solution] = shapeMatchSolver(shape1,shape2, non_matched_FY, VY_old, closed, dummy_elements, options,filename, configMatching);
end
end


function vis_verts(V,F,show_color, nr, min_c, max_c)
    subplot(2,2,nr)
    patch('Vertices', V, 'Faces', F, 'FaceVertexCData',show_color, 'FaceColor', 'interp', 'EdgeColor', 'interp')
    
    view([90,20]);
    
    colormap(jet(256));
    clim([min_c max_c])
    colorbar;
end

