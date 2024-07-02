function [] = interval_search(curr_animal, shape_1, shape_2, curr_lambda, data_folder)
    curr_animal = char(curr_animal);
    shape_1 = char(shape_1);
    shape_2 = char(shape_2);
    curr_lambda = str2num(curr_lambda);
    disp("Lambda: " + curr_lambda)
    calc_mean = true;
    max_iter = 300;
    max_energy = 200;
    use_lp = 0;
    curr_dataset = "cp2p";
    output_folder = 'output/'

    used_folder = ['./matlab/src/solveShapeMatch/' output_folder '/' shape_1 '_' shape_2 '/'];
    used_energies = max_energy * ones(1, max_iter);
    not_solvable_arr = zeros(1, max_iter);
    L = 1;
    start_step_size = 25;
    min_step = 1;
    all_steps = [min_step:start_step_size:max_iter - 1, max_iter - 1];

    % Set energies in intervals of st
    for i = all_steps
        [curr_G, curr_t, curr_E, not_solvable] = get_G(i + 1, used_folder, curr_animal, shape_1, shape_2, use_lp, calc_mean, output_folder, curr_lambda, curr_dataset, data_folder);
        used_energies(i + 1) = curr_E;
        not_solvable_arr(i + 1) = not_solvable;
        fprintf('%d & %d & %f & %f \n', i, curr_G, curr_t, curr_E);
    end

    checked_all = false;
    new_energies = used_energies;

    while ~checked_all
        checked_all = true;

        % Divide used_energy by i+1
        curr_min = nanmin(used_energies);

        for i = 1:max_iter
            curr_e = used_energies(i);

            % Check next value that is not max_energy
            for j = (i + 1):max_iter
                if used_energies(j) ~= max_energy
                    next_e = used_energies(j);
                    break;
                end
            end
            if size(j, 1) == 0
                continue;
            end
            if not_solvable_arr(j) == 1 || curr_e == max_energy || next_e == max_energy
                continue;
            end

            mean_e = curr_e;
            mean_next_e = next_e;

            % Check if it is possible that there is a value between mean_next_e and mean_e
            % Get values between i+1 and j
            between_vals = i + 1:j - 1;
            if calc_mean
                max_increase_arr = L ./ between_vals;

                max_increase = sum(max_increase_arr);
                curr_increase = mean_next_e - mean_e;
                diff = max_increase - curr_increase;
                min_poss_value = mean_e - diff;
                check_nan = isnan(min_poss_value);

                if (check_nan || min_poss_value < curr_min) && ~isempty(max_increase_arr)
                    checked_all = false;

                    % There is a value between mean_nextx_e and mean_e
                    % Get the value
                    new_i = (j + i) / 2;
                    new_i = round(new_i);

                    % Get the new energy
                    [~, ~, curr_E, not_solvable] = get_G(new_i, used_folder, curr_animal, shape_1, shape_2, use_lp, calc_mean, output_folder, curr_lambda, curr_dataset, data_folder);
                    new_energies(new_i) = curr_E;
                    not_solvable_arr(new_i) = not_solvable;
                end
            end

            used_energies = new_energies;
        end
    end

    function [num_faces, curr_time, usedEnergy, not_solvable] = get_G(i, used_folder, curr_animal, shape_1, shape_2, use_lp, calc_mean, output_folder, curr_lambda, curr_dataset, data_folder)
        fullfile(used_folder, ['*' shape_1 '_' shape_2 '_*_faces_curr_nr_' num2str(i + 1) '_lp_' num2str(use_lp) '_E.mat'])
        files = dir(fullfile(used_folder, ['*' shape_1 '_' shape_2 '_*_faces_curr_nr_' num2str(i + 1) '_lp_' num2str(use_lp) '_E.mat']));
        % check if files are empty
        if isempty(files)
            lambda_mem = 0;
            num_faces = 100;
            curr_nr = i + 1;
            curr_type = "cuts";
            solveOneShape(shape_1, shape_2, curr_animal, lambda_mem, num_faces, curr_type, curr_nr, output_folder, curr_lambda, curr_dataset, data_folder);
            files = dir(fullfile(used_folder, ['*' shape_1 '_' shape_2 '_*_faces_curr_nr_' num2str(i + 1) '_lp_' num2str(use_lp) '_E.mat']));
        end
        filename = fullfile(used_folder, files(end).name);
        not_solvable = 0;

        % Replace _E.mat with _G.mat
        G_filename = strrep(filename, '_E.mat', '_G.mat');
        t_filename = strrep(filename, '_E.mat', '_time.mat');

        try
            mat_G = load(G_filename);
            mat_t = load(t_filename);
            num_faces = sum(mat_G.G(:));
            curr_time = mat_t.tOpt;
        catch
            num_faces = 0;
            curr_time = 420;
        end

        mat = load(filename);
        if calc_mean
            usedEnergy = mat.usedE / (i + 1);
        else
            usedEnergy = mat.usedE;
        end

        if isempty(files)
            disp('File not accessible');
            usedEnergy = NaN;
            curr_time = NaN;
            num_faces = NaN;
        end

        if num_faces == 0
            disp('File not accessible');
            usedEnergy = NaN;
            curr_time = NaN;
            num_faces = NaN;
            not_solvable = 1;
        end
    end
end
