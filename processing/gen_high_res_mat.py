import glob
import os
import time

import numpy as np
import open3d as o3d
import scipy.io as sio
import torch

def vis_p2p(verts_source, verts_target, faces_source, faces_target, corres, corres_source):
    """
    Visualize point-to-point correspondences between two meshes.
    
    Args:
        verts_source (np.ndarray): Vertices of the source mesh
        verts_target (np.ndarray): Vertices of the target mesh
        faces_source (np.ndarray): Faces of the source mesh
        faces_target (np.ndarray): Faces of the target mesh
        corres (np.ndarray): Correspondence indices between source and target vertices
    """
    import polyscope as ps

    y_shift = np.array([0.7, 0, 0])

        # rotate both shapes by -90 degrees around the z-axis
    rotation_matrix = np.array([[0, 1, 0],
                                [-1, 0, 0],
                                [0, 0, 1]])

    verts_source = verts_source @ rotation_matrix.T
    verts_target = verts_target @ rotation_matrix.T


    # and then rotate shape Y by -90 degrees around the x-axis
    rotation_matrix_y = np.array([[1, 0, 0],
                                   [0, 0, 1],
                                   [0, 1, 0]])
    verts_target = verts_target @ rotation_matrix_y.T
    verts_source = verts_source @ rotation_matrix_y.T

    ps.init()
    
    poly_shape_source = ps.register_surface_mesh("source", verts_source, faces_source)
    poly_shape_target = ps.register_surface_mesh("target", verts_target + y_shift, faces_target)

    color_source = verts_source
    color_source = color_source - np.min(color_source)
    color_source = color_source / np.max(color_source)
    color_source[corres_source==-1] = np.ones(3) * 0.7
    color_target = color_source[corres]
    color_target[corres==-1] = np.ones(3) * 0.7
    poly_shape_source.add_color_quantity("correspondences", color_source, defined_on='vertices', enabled=True)
    poly_shape_target.add_color_quantity("correspondences", color_target, defined_on='vertices', enabled=True)

    #ps.screenshot(f"correspondences.png")
    ps.show()

def all_nn_query(feat_x, feat_y, dim=-2):
    """
    Find correspondences via nearest neighbor query
    Args:
        feat_x: feature vector of shape x. [V1, C].
        feat_y: feature vector of shape y. [V2, C].
        dim: number of dimension
    Returns:
        p2p: point-to-point map (shape y -> shape x). [V2].
    """
    dist = torch.cdist(feat_x, feat_y)  # [V1, V2]
    p2p = dist.topk(feat_x.shape[0], dim=dim, largest=False).indices
    return p2p, dist


def triangle_comb_to_p2p(triangle_comb, size_VX, size_VY):
    """
    Find correspondences via nearest neighbor query
    Args:
        feat_x: feature vector of shape x. [V1, C].
        feat_y: feature vector of shape y. [V2, C].
        dim: number of dimension
    Returns:
        p2p: point-to-point map (shape y -> shape x). [V2].
    """
    p2pX = {}
    p2pY = {}
    p2pXFirst = -1 * np.ones((size_VX, 1))
    p2pYFirst = -1 * np.ones((size_VY, 1))
    for curr_comb in triangle_comb:
        for j in range(3):
            comb_x = curr_comb[j] - 1
            comb_y = curr_comb[j + 3] - 1
            if comb_x not in p2pX:
                p2pX[comb_x] = {comb_y}
            else:
                p2pX[comb_x].add(comb_y)
                p2pXFirst[comb_x] = comb_y
            if comb_y not in p2pY:
                p2pY[comb_y] = {comb_x}
            else:
                p2pY[comb_y].add(comb_x)
                p2pYFirst[comb_y] = comb_x
    return p2pX, p2pY, p2pXFirst, p2pYFirst


def get_best_matching(shape_1, shape_2, used_folder, use_mean=True, used_lambda=0):
    """
    Get best matching between two shapes
    """
    max_iter = 300
    big_nr = np.inf
    # numpy array to save mean curve
    mean_curve = big_nr * np.ones((1, max_iter))
    sum_curve = big_nr * np.ones((1, max_iter))
    mean_runtime = big_nr * np.ones((1, max_iter))
    sum_runtime = big_nr * np.ones((1, max_iter))
    # load mat file
    for i in range(0, max_iter):
        # get all files that have pattern
        files = glob.glob(
            f"{used_folder}*{shape_1}_{shape_2}_*_faces_curr_nr_{i+1}_*_E.mat"
        )
        if len(files) == 0:
            # print("File not accessible")
            mean_curve[0, i] = None
            sum_curve[0, i] = None
            mean_runtime[0, i] = None
            sum_runtime[0, i] = None

        else:
            filename = files[-1]
            # replace _E.mat with _G.mat
            G_filename = filename.replace("_E.mat", "_G.mat")
            t_filename = filename.replace("_E.mat", "_time.mat")
            # get the number of faces
            try:
                mat_G = sio.loadmat(G_filename)
                mat_t = sio.loadmat(t_filename)
                num_faces = mat_G["G"].sum()
                curr_time = mat_t["tOpt"][0, 0]
            except:
                num_faces = 0
                curr_time = 420
            if num_faces == 0:
                # print("File not accessible")
                mean_curve[0, i] = None
                sum_curve[0, i] = None
                mean_runtime[0, i] = None
                sum_runtime[0, i] = None
                continue
            mat = sio.loadmat(filename)
            mean_curve[0, i] = mat["usedE"] / (i + 1)
            sum_curve[0, i] = mat["usedE"] - used_lambda * (i + 1)
            mean_runtime[0, i] = curr_time / (i + 1)
            sum_runtime[0, i] = curr_time
    if use_mean:
        return np.nanargmin(mean_curve) + 1
    else:
        return np.nanargmin(sum_curve) + 1


def get_high_res_corres(p2pFirstX, p2p, VX_high, highToLowX, highToLowY, p2p_high):
    for i in range(VX_high.shape[0]):
        low_X_idx = highToLowX[i]
        if p2pFirstX[low_X_idx] != -1:
            low_Y_idx = p2p[low_X_idx.item()]
            for curr_low_Y_idx in low_Y_idx:
                # get all idx in highToLowY that have value low_Y_idx
                highToLowY_idx = torch.nonzero(highToLowY == curr_low_Y_idx).flatten()
                highToLowY_idx = highToLowY_idx.unique()
                p2p_high[i, highToLowY_idx.long()] = 1
    return p2p_high


if __name__ == "__main__":
    use_mean = True
    save_mode = True
    calc_corres = True
    vis = False
    final_folder = "./high_res_corres/"
    # check if folder exists
    if not os.path.exists(final_folder):
        os.makedirs(final_folder)

    # geodesic distance folder
    dist_folder = "./dists/"

    category = "cat"
    nr_1 = "7"
    nr_2 = "6"
    basename_shape_1 = f"{category}-{nr_1}"
    basename_shape_2 = f"{category}-{nr_2}"
    shape_1 = f"{basename_shape_1}.off"
    shape_2 = f"{basename_shape_2}.off"

    # TODO: ADD path to high resolution features
    featInfoX = sio.loadmat(f"./features/cuts_{category}_shape_{nr_1}.mat")
    featInfoY = sio.loadmat(f"./features/cuts_{category}_shape_{nr_2}.mat")
    featX = torch.from_numpy(featInfoX["featY"])
    featY = torch.from_numpy(featInfoY["featY"])
    used_folder = f"./matlab/src/solveShapeMatch/output/{basename_shape_1}_{basename_shape_2}/"

    file_1 = f"./low_res_meshes/{basename_shape_1}_extra_info.mat"
    file_2 = f"./low_res_meshes/{basename_shape_2}_extra_info.mat"
    # load mat file
    extraInfoX = sio.loadmat(file_1)
    extraInfoY = sio.loadmat(file_2)

    idx_VX = extraInfoX["idx_V"] - 1
    idx_VY = extraInfoY["idx_V"] - 1

    V_extraX = extraInfoX["V"]
    V_extraX = torch.from_numpy(V_extraX)
    V_extraX_low = extraInfoX["V_low"]
    V_extraX_low = torch.from_numpy(V_extraX_low)

    V_extraY = extraInfoY["V"]
    V_extraY = torch.from_numpy(V_extraY)
    V_extraY_low = extraInfoY["V_low"]
    V_extraY_low = torch.from_numpy(V_extraY_low)

    try:
        min_nr = get_best_matching(
            basename_shape_1, basename_shape_2, used_folder, use_mean
        )
        vis_files = glob.glob(
            f"{used_folder}*{basename_shape_1}_{basename_shape_2}_*_faces_curr_nr_{min_nr}_lp_0_Vis_py.mat"
        )
    except Exception as e:
        print("No best matching found")
    try:
        vis_file = vis_files[-1]
    except Exception as e:
        print("No valid file found")
        # load mat file
    Vis_mat = sio.loadmat(vis_file)
    usedComb = Vis_mat["usedComb"]
    VX = Vis_mat["VX"]
    VY = Vis_mat["VY"]

    usedComb = usedComb[
        np.all(usedComb[:, 1:3] <= V_extraX_low.shape[0], 1), :
    ]
    usedComb = usedComb[
        np.all(usedComb[:, 3:6] <= V_extraY_low.shape[0], 1), :
    ]
    p2pX, p2pY, p2pFirstX, p2pFirstY = triangle_comb_to_p2p(
        usedComb, VX.shape[0], VY.shape[0]
    )

    VX_high_path = f"./data/{basename_shape_1}.off"
    VY_high_path = f"./data/{basename_shape_2}.off"

    mesh_X_high = o3d.io.read_triangle_mesh(VX_high_path)
    mesh_Y_high = o3d.io.read_triangle_mesh(VY_high_path)

    VX_high = torch.from_numpy(np.asarray(mesh_X_high.vertices)).float()
    VY_high = torch.from_numpy(np.asarray(mesh_Y_high.vertices)).float()

    FX_high = torch.from_numpy(np.asarray(mesh_X_high.triangles)).float()
    FY_high = torch.from_numpy(np.asarray(mesh_Y_high.triangles)).float()

    nneigh, nneigh_values = all_nn_query(featX, featY)
    nneighAlt, nneighAlt_values = all_nn_query(featY, featX)
    # load dist files for shape 1 and shape 2
    dist_X = torch.cdist(torch.tensor(extraInfoX["V"]), torch.tensor(VX))
    dist_Y = torch.cdist(torch.tensor(extraInfoY["V"]), torch.tensor(VY))

    highToLowX = dist_X.argmin(dim=-1)
    highToLowY = dist_Y.argmin(dim=-1)

    assert not (
        V_extraX.shape < VX_high.shape or V_extraY.shape < VY_high.shape
    )

    p2p_high = torch.zeros((VX_high.shape[0], VY_high.shape[0]))
    p2p_high_alt = torch.zeros((VY_high.shape[0], VX_high.shape[0]))
    # stop time for loop
    start_time = time.time()

    usedXLow = list(p2pX.keys())
    usedYLow = list(p2pY.keys())
    # set all elements of highToLowX to -1 that are not in usedXLow
    usedXHigh = np.isin(highToLowX, usedXLow)
    usedYHigh = np.isin(highToLowY, usedYLow)

    if calc_corres:
        p2p_high = get_high_res_corres(
            p2pFirstX, p2pX, VX_high, highToLowX, highToLowY, p2p_high
        )
        p2p_alt = get_high_res_corres(
            p2pFirstY,
            p2pY,
            VY_high,
            highToLowY,
            highToLowX,
            p2p_high_alt,
        )

        max_val = torch.inf
        nneigh_values[p2p_high == 0] = max_val
        nneighAlt_values[p2p_high_alt == 0] = max_val
        nneigh = torch.argmin(nneigh_values, dim=-1)
        nneigh[nneigh == 0] = -1
        nneighAlt = torch.argmin(nneighAlt_values, dim=-1)
        nneighAlt[nneighAlt == 0] = -1

    # get all elements of color_X that are > 0
    used_indices_X = np.array(usedXHigh) > 0
    used_indices_Y = np.array(usedYHigh) > 0

    # binary array for overlapping region in X and Y
    used_indices_X = used_indices_X.astype(int)
    used_indices_Y = used_indices_Y.astype(int)

    ours_geo_error = nneigh.flatten()
    ours_geo_error_alt = nneighAlt.flatten()
    ours_geo_error_alt[used_indices_Y == 0] = -1
    ours_geo_error[used_indices_X == 0] = -1
    if vis:
        vis_p2p(
            verts_source=VX_high.numpy(),
            verts_target=VY_high.numpy(),
            faces_source=FX_high.numpy(),
            faces_target=FY_high.numpy(),
            corres=nneighAlt.numpy(),
            corres_source=nneigh.numpy(),
        )

    if save_mode:
        # save to matlab file
        sio.savemat(
            f"{final_folder}/{basename_shape_1}-{basename_shape_2}_Vis_partial.mat",
            {
                "VX": VX_high.numpy(),
                "VY": VY_high.numpy(),
                "FX": FX_high.numpy(),
                "FY": FY_high.numpy(),
                "p2p": ours_geo_error.numpy(),
                "used_indices_X": used_indices_X,
                "used_indices_Y": used_indices_Y,
            },
        )
