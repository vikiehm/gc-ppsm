function [G, times, solution, combAll] = shapeMatchSolver(shape1, shape2, closed, varargin)
    %SHAPEMATCHSOLVER Script to solve the shape match problem
    %
    % Syntax:
    % [G, times, solution] = shapeMatchSolver(VX, VY, FX, FY)
    % [G, times, solution] = shapeMatchSolver(VX, VY, FX, FY, solverOps)
    % [G, times, solution] = shapeMatchSolver(VX, VY, FX, FY, solverOps, filename)
    % [G, times, solution] = shapeMatchSolver(VX, VY, FX, FY, solverOps, filename, degenerate)
    %
    % Input Arguments:
    %       - VX:        point coordinates from triangle X
    %       - FX:        triangulation from triangle X
    %       - VY:        point coordinates from triangle Y
    %       - FY:        triangulation from triangle Y
    %       - solverOps: field of created with sdpsettings
    %                    standard solver is cplex
    %       - filename   if filename is given the plot as well as the console
    %                    output is saved to output/filename.txt and
    %                    output/filename.png respectively
    %       - degenerate boolean defining degenerate
    %                    [{true}, false]
    %
    % Output Arguments:
    %       - G:         the matching matrix
    %       - times:     times of the various parts of the shapeMatchSolver
    %                    namely: [tOpt tEnergy tConstraints]
    %       - Gamma:     the optimization variable
    %       - solution:  the solution from the yalmip solver
    %
    % Note: there is an iternal options with which an iterative scheme for
    % solving the shape matching problem can be applied. Simply set:
    % relaxed = true;

    %% Parse Inputs

    VX = cell2mat(shape1('V'));
    VY = cell2mat(shape2('V'));
    FX = cell2mat(shape1('F'));
    FY = cell2mat(shape2('F'));
    corr_infos = varargin{end};

    shapeXName = shape1('name');
    shapeYName = shape2('name');

    curr_nr = cell2mat(shape2('curr_nr'));
    curr_lambda = cell2mat(shape2('curr_lambda'));
    output_folder = shape2('output_folder');

    options = sdpsettings('solver', 'cplex');
    degenerate = false;

    varargin = {varargin{1:end - 1}};

    if nargin > 6
        configMatching = varargin{end};
        partial_to_partial = configMatching("partial-partial");
        degenerate = configMatching("degenerate");

        varargin = {varargin{1:end - 1}};

        for i = 1:2
            if ~iscell(varargin)
                varargin = {varargin};
            end
            if ischar(varargin{end})
                filename = varargin{end};
                varargin = {varargin{1:end - 1}};
            elseif isstruct(varargin{end})
                options = varargin{end};
                varargin = {varargin{1:end - 1}};
            else
                error('Inputs not supported')
            end
            if nargin < 6
                break;
            end
        end

    end

    % switch logging on
    if exist('filename')
        filepath = fileparts(which('shapeMatchSolver.m'));
        output_path = strcat(filepath, '/..', '/', output_folder, '/', shapeXName, '_', shapeYName, '/');
        output_path = output_path{1};
        % check if path exists
        if ~exist(output_path, 'dir')
            mkdir(output_path);
        end
        path = output_path;
        filenameConsole = strcat(path, filename, '.txt');
        diary(filenameConsole);
    end
    disp("+++++++++++++ SHAPE MATCH SOLVER ++++++++++++++ ")

    % Get face combinations
    [FaComboNDeg, FbComboNDeg] = getCombinations(FX, FY);
    [FaComboDegA, FbComboDegA] = getDegenerateCombinations(FX, FY, closed);
    [FbComboDegB, FaComboDegB] = getDegenerateCombinations(FY, FX, closed);

    combAll = [FaComboNDeg, FbComboNDeg; FaComboDegA, FbComboDegA; FaComboDegB, FbComboDegB];

    if size(corr_infos, 1) > 0
        all_matchings = corr_infos('all_matchings');
        all_matchings = all_matchings{1};
        comb1 = combAll(:, [1 4]);
        comb2 = combAll(:, [2 5]);
        comb3 = combAll(:, [3 6]);
        comb1_shift1 = combAll(:, [1 5]);
        comb1_shift2 = combAll(:, [1 6]);
        comb2_shift1 = combAll(:, [2 4]);
        comb2_shift2 = combAll(:, [2 6]);
        comb3_shift1 = combAll(:, [3 4]);
        comb3_shift2 = combAll(:, [3 5]);

        pruneVec = zeros(size(combAll, 1), 1);
        pruneVec = pruneVec | ismember(comb1, all_matchings, 'rows');
        pruneVec = pruneVec | ismember(comb2, all_matchings, 'rows');
        pruneVec = pruneVec | ismember(comb3, all_matchings, 'rows');
        pruneVec = pruneVec | ismember(comb1_shift1, all_matchings, 'rows');
        pruneVec = pruneVec | ismember(comb1_shift2, all_matchings, 'rows');
        pruneVec = pruneVec | ismember(comb2_shift1, all_matchings, 'rows');
        pruneVec = pruneVec | ismember(comb2_shift2, all_matchings, 'rows');
        pruneVec = pruneVec | ismember(comb3_shift1, all_matchings, 'rows');
        pruneVec = pruneVec | ismember(comb3_shift2, all_matchings, 'rows');

        FaComboDegB = FaComboDegB(pruneVec((size(FaComboNDeg, 1) + size(FaComboDegA, 1) + 1):end) == 1, :);
        FbComboDegB = FbComboDegB(pruneVec((size(FbComboNDeg, 1) + size(FbComboDegA, 1) + 1):end) == 1, :);
        FaComboDegA = FaComboDegA(pruneVec((size(FaComboNDeg, 1) + 1):(size(FaComboNDeg, 1) + size(FaComboDegA, 1))) == 1, :);
        FbComboDegA = FbComboDegA(pruneVec((size(FbComboNDeg, 1) + 1):(size(FbComboNDeg, 1) + size(FbComboDegA, 1))) == 1, :);
        FaComboNDeg = FaComboNDeg(pruneVec(1:size(FaComboNDeg, 1)) == 1, :);
        FbComboNDeg = FbComboNDeg(pruneVec(1:size(FbComboNDeg, 1)) == 1, :);
        combAll = combAll(pruneVec, :);
        assert(size(combAll, 1) == (size(FaComboNDeg, 1) + size(FaComboDegA, 1) + size(FaComboDegB, 1)));
        assert(size(combAll, 1) == size(FbComboNDeg, 1) + size(FbComboDegA, 1) + size(FbComboDegB, 1));
    end
    sizeGamma = size(combAll, 1);

    tic

    numFacesX = length(FX);
    numFacesY = length(FY);

    [constraintMatrix, constraintVector, productE] = getConstraints(VX, FX, VY, FY, degenerate, closed, corr_infos);

    % remove unused rows
    if size(corr_infos, 1) > 0
        constraintMatrix = constraintMatrix(:, pruneVec);
        del = constraintMatrix(1:(end - size(FX, 1) - size(FY, 1)), :);
        proj = constraintMatrix((end - size(FX, 1) - size(FY, 1) + 1):end, :);
        % Find the row and column indices of non-zero elements
        [rowIndices, ~] = find(del);

        % Find the unique column indices
        nonRowColumnIndices = unique(rowIndices);
        productE = productE(nonRowColumnIndices, :);
        del = del(nonRowColumnIndices, :);
        constraintMatrix = [del; proj];
    end

    tConstraints = toc;

    %% Construct energy vector
    tic
    E = deformationEnergy(shape1, shape2, combAll);
    tEnergy = toc;
    assert(sum(E < 0) == 0)

    %% Optimize

    % optimization variable (yalmip syntax)
    Gamma = binvar(sizeGamma, 1);

    % find bounding edges for partial shape Y
    [~, ~, ~, bounding_edgesY] = findTriMeshHoles(FY, VY);
    bound_edges = zeros(size(productE, 1), 1);
    size(bounding_edgesY)

    % set product edges to 1 if it contains bounding edge from partial shape
    for i = 1:size(bounding_edgesY, 1)
        curr_b = bounding_edgesY(i, :);
        bound_edges = bound_edges | (sum(abs(double(productE(:, 3:4)) - double(curr_b)), 2) == 0);
        bound_edges = bound_edges | (sum(abs(double(productE(:, 3:4)) - double([curr_b(2), curr_b(1)])), 2) == 0);
    end

    [~, ~, ~, bounding_edgesX] = findTriMeshHoles(FX, VX);

    if partial_to_partial
        for i = 1:size(bounding_edgesX, 1)
            curr_b = bounding_edgesX(i, :);
            bound_edges = bound_edges | (sum(abs(double(productE(:, 1:2)) - curr_b), 2) == 0);
            bound_edges = bound_edges | (sum(abs(double(productE(:, 1:2)) - [curr_b(2), curr_b(1)]), 2) == 0);
        end
    end

    new_del = constraintMatrix(1:end - numFacesX - numFacesY, :);

    new_constr = constraintMatrix(end - numFacesY + 1:end, :);
    constr_X = constraintMatrix((end - numFacesX - numFacesY + 1):(end - numFacesY), :);

    new_compare = constraintVector((end - numFacesY + 1):end, :);
    new_compare_X = constraintVector((end - numFacesX - numFacesY + 1):(end - numFacesY), :);

    constraint = new_constr * Gamma <= new_compare;

    constraint = [constraint, constr_X * Gamma <= new_compare_X];
    constraint = [constraint, (~bound_edges .* new_del) * Gamma == zeros(size(new_del, 1), 1)];
    constraint = [constraint, sum(Gamma) <= (size(FX, 1) + size(FY, 1))];
    if curr_nr > 0
        "Curr nr"
        curr_nr
        constraint = [constraint, sum(Gamma) == curr_nr];
    end
    non_bond_triangles = zeros(3 * size(FX, 1) * size(FY, 1), 1);
    % get non boundary triangles
    for i = 1:size(non_bond_triangles, 1)
        curr_comb = combAll(i, :);
        comb1 = curr_comb(1:3);
        is_bound_X = ismember(comb1(1), bounding_edgesX) || ismember(comb1(2), bounding_edgesX) || ismember(comb1(3), bounding_edgesX);
        comb1 = curr_comb(4:6);
        is_bound_Y = ismember(comb1(1), bounding_edgesY) || ismember(comb1(2), bounding_edgesY) || ismember(comb1(3), bounding_edgesY);
        if ~is_bound_X && ~is_bound_Y
            non_bond_triangles(i) = 1;
        end
    end

    if sum(non_bond_triangles) == 0 %
        disp("No non boundary triangles found")
    else
        disp("ADDING constraint >= 1")
        constraint = [constraint, sum(non_bond_triangles .* Gamma(1:(3 * size(FX, 1) * size(FY, 1)))) >= 1];
    end

    E = E - min(E);
    E = E / max(E);
    disp(curr_lambda)
    size(E')
    size(Gamma)
    objective = (E' * Gamma) - curr_lambda * sum(Gamma);

    tic
    disp("++++++++++ Try to solve for global optimality +++++++++++")

    tic
    disp("++++++++++ Solving for global optimality +++++++++++")
    solution = optimize(constraint, objective, options);
    tOpt = toc;

    tOpt = toc;
    %% Plot Solution
    G = value(Gamma);
    usedE = E' * G;
    %load("test_G.mat")
    % make sure solution really is binary
    G = G > 0.5;
    usedComb = combAll(G, :);
    not_matched_FX = ~(constr_X * G);
    not_matched_FY = ~(new_constr * G);
    fig = plotSolution(G, VX, FX, VY, FY, not_matched_FX, not_matched_FY, shape2, []);

    %% Custom Output
    disp("++++++++++++ END OF SOLVER OUTPUT +++++++++++++")
    disp("Solution contains >>> " + string(nnz(G(1:size(FaComboNDeg, 1), :))) + " <<< NON DEGENERATE matchings")
    disp("Solution contains >>> " + string(nnz(G((size(FaComboNDeg, 1) + 1):end, :))) + " <<< DEGENERATE matchings")
    disp("Solution contains >>> " + string(nnz(G)) + " <<< matchings")
    disp("Solution has wks energy of >>> " + string(E' * G) + " <<< matchings")

    disp("Total time consumption approx: " + string(tOpt + tEnergy + tConstraints))
    disp("      => optimization    : " + string(tOpt) + " s")
    disp("      => energy comp     : " + string(tEnergy) + " s")
    disp("      => constraints comp: " + string(tConstraints) + " s")

    times = [tOpt tEnergy tConstraints];

    visDict = containers.Map;
    visDict("VX") = VX;
    visDict("VY") = VY;
    visDict("FX") = FX;
    visDict("FY") = FY;
    visDict("usedComb") = usedComb;
    %% Save files
    if exist('filename')
        if ~exist(path, 'dir')
            mkdir(path);
        end
        strcat(path, filename, '_G.mat')

        save(strcat(path, filename, '_G.mat'), 'G');
        % save time for speedtest
        save(strcat(path, filename, '_time.mat'), 'tOpt');
        save(strcat(path, filename, '_usedComb.mat'), 'usedComb');
        % save to file
        save(strcat(path, filename, '_Vis.mat'), 'visDict');
        save(strcat(path, filename, '_Vis_py.mat'), "usedComb", "VX", "VY", "FX", "FY");
        save(strcat(path, filename, '_E.mat'), 'usedE');
        save(strcat(path, filename, '_E_vec.mat'), 'E');

        saveas(fig, strcat(path, filename, '.png'));
        saveas(fig, strcat(path, filename, '.fig'));
        diary off;
        disp('Files successfully saved')
    end

    disp("++++++++++ END OF SHAPE MATCH SOLVER ++++++++++")
end
