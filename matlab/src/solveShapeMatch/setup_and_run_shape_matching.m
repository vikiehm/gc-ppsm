function [G, times, solution, combAll, usedComb, vis] = setup_and_run_shape_matching(VX, FX, VY, FY, FeatX, FeatY, lambda_mem, shape1_name, shape2_name, num_faces, oldSol, idxFX, idxFY, curr_nr, output_folder, curr_lambda)
    close all
    partial_to_partial = true;
    closed = 0;
    degenerate = true;
    relaxed = false;
    FY_old = FY;
    VY_old = VY;

    shape1 = dictionary('V', {VX}, 'F', {FX}, ...
        'features', {FeatX}, 'idxFX', {idxFX}, "name", {shape1_name});
    shape2 = dictionary('V', {VY}, 'F', {FY}, ...
        'old_F', {FY_old}, 'old_V', {VY_old}, 'features', {FeatY}, 'idxFY', {idxFY}, "curr_nr", {curr_nr}, "name", {shape2_name}, "output_folder", {output_folder}, "curr_lambda", {curr_lambda});

    keys = ["lambda", "degenerate", "relaxed", "partial-partial", "curr_nr"];
    values = [lambda_mem, degenerate, relaxed, partial_to_partial, curr_nr];

    disp("Using Gurobi as Solver")
    options = sdpsettings('solver', 'gurobi', 'verbose', 1, 'debug', 1, 'gurobi.TimeLimit', 60);

    filepath = fileparts(which('shapeMatchSolver.m'));

    output_path = strcat(filepath, '/..', '/output_new/', shape1_name, '_', shape2_name, '/');

    configMatching = dictionary(keys, values);
    curr_filename = strcat(datestr(now, 'mm-dd-yyyy-HH-MM'), shape1_name, '_', shape2_name, '_', string(num_faces), '_faces');
    add_vars = "";

    keys = ["curr_nr", "relaxed", "partial-partial", "closedBound"];

    for j = 1:size(keys, 1)
        curr_key = keys{j};
        curr_value = configMatching(curr_key);
        add_vars = strcat(add_vars, '_', curr_key, '_', string(curr_value));
    end

    add_vars = strcat(add_vars, "_lp_0");
    filename = convertStringsToChars(strcat(curr_filename, add_vars));

    fileinfo = dir(strcat(output_path, "*", shape1_name, '_', shape2_name, '_', string(num_faces), '_faces', add_vars, "_G.mat"));
    fnamesG = {fileinfo.name}
    if size(fnamesG, 1) > 1
        disp("Loading G from file")
        [FaComboNDeg, FbComboNDeg] = getCombinations(FX, FY, closed);
        [FaComboDegA, FbComboDegA] = getDegenerateCombinations(FX, FY, closed);
        [FbComboDegB, FaComboDegB] = getDegenerateCombinations(FY, FX, closed);
        combAll = [FaComboNDeg, FbComboNDeg; FaComboDegA, FbComboDegA; FaComboDegB, FbComboDegB];
        curr_G_name = fnamesG{end};
        G = load(curr_G_name).G;
        %replace G with used comb
        used_comb_name = strrep(curr_G_name, 'G.mat', 'usedComb.mat');
        vis_name = strrep(curr_G_name, 'G.mat', 'Vis.mat');
        usedComb = load(used_comb_name).usedComb;
        vis = load(vis_name);
        % check if field visDict exists in vis
        checker = isfield(vis, 'visDict');

        if checker
            vis.VY = vis.visDict("VY");
            vis.VX = vis.visDict("VX");
            vis.FY = vis.visDict("FY");
            vis.FX = vis.visDict("FX");
        end
        times = 0;
        solution = 0;
    else
        [G, times, solution, combAll] = shapeMatchSolver(shape1, shape2, closed, options, filename, configMatching, oldSol);
        usedComb = combAll(G, :);
        vis = [];
    end
end
