function [VX, FX, VY, FY, FX_low, FY_low, idx_VX, idx_VY, idx_FX, idx_FY, added_verts, VX_all, VY_all, FeatX, FeatY] = get_non_manifold_mesh(shape_1, shape_2, num_faces, data_folder)
    extra_info = load(strcat(data_folder, "/", shape_2, "_extra_info.mat"));
    extra_info_X = load(strcat(data_folder, "/", shape_1, "_extra_info.mat"));

    FY_low = double(extra_info.F_low);
    FY = double(extra_info.F);
    VY = extra_info.V_low;
    VY_all = extra_info.V;
    idx_FY = double(extra_info.idx_F);
    idx_VY = extra_info.idx_V;
    added_verts = extra_info.added_verts;
    VX = extra_info_X.V_low;
    VX_all = extra_info_X.V;
    FX_low = double(extra_info_X.F_low);
    FX = double(extra_info_X.F);
    idx_FX = double(extra_info_X.idx_F);
    idx_VX = double(extra_info_X.idx_V);
    FeatX = extra_info_X.features;
    FeatY = extra_info.features;
end
