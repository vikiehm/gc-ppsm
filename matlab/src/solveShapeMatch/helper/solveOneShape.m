function [] = solveOneShape(shape_1, shape_2, curr_animal, lambda_mem, num_faces, curr_type, curr_num, output_folder, curr_lambda, curr_dataset, data_folder)
    splited_shape = split(shape_2, '.');
    shape_2 = splited_shape{1};

    splited_shape = split(shape_1, '.');
    shape_1 = splited_shape{1};

    disp(strcat("Start matching ", shape_1, " with ", shape_2));
    high_res_faces = num_faces;

    close all
    [VX_org, FX_all, VY_org, FY_all, FX_org, FY_org, idx_VX, idx_VY, idx_FX, idx_FY, added_verts, VX_all, VY_all, FeatX, FeatY] = get_non_manifold_mesh(shape_1, ...
        shape_2, high_res_faces, data_folder);

    VX_low = VX_org;
    VY_low = VY_org;
    FX_low = FX_org;
    FY_low = FY_org;

    curr_high_res_facesX = size(FX_org, 1);

    close all

    idxFX = 1;
    idxFY = 1;

    % setup and run shape matching
    setup_and_run_shape_matching(VX_low, FX_low, VY_low, FY_low, FeatX, FeatY, ...
        lambda_mem, shape_1, shape_2, curr_high_res_facesX, [], idxFX, idxFY, curr_num, ...
        output_folder, curr_lambda);
end
