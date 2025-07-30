from typing import Tuple
import numpy as np
import igl
import open3d as o3d
from pathlib import Path
import torch
from tqdm import tqdm
import scipy.io as sio
from sklearn.neighbors import NearestNeighbors


def triangles_share_edge(triangle1, triangle2):
    common_vertices = set(triangle1) & set(triangle2)
    return len(common_vertices) == 2

def downsample_mesh(
    shape_fp: Path,
    decimate_size: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    mesh = o3d.io.read_triangle_mesh(str(shape_fp))

    non_manifold_verts = mesh.get_non_manifold_vertices()

    extra_info_dict = {}
    added_verts = []
    added_corres = []
    F = np.asarray(mesh.triangles)
    V = np.asarray(mesh.vertices)

    while len(non_manifold_verts) > 0:
        for curr_vert in non_manifold_verts:
            rows, cols = (F == curr_vert).nonzero()
            V = np.concatenate([V, V[curr_vert, :][None]])
            added_verts.append(curr_vert)
            curr_t = F[rows[0], :]
            for i in range(1, len(rows)):
                other_t = F[rows[i], :]
                if triangles_share_edge(other_t, curr_t):
                    F[rows[i], cols[i]] = V.shape[0] - 1
            F[rows[0], cols[0]] = V.shape[0] - 1

        mesh.vertices = o3d.utility.Vector3dVector(V)
        mesh.triangles = o3d.utility.Vector3iVector(F)
        non_manifold_verts = mesh.get_non_manifold_vertices()

    P, _ = igl.orientable_patches(F)
    [FF, _] = igl.orient_outward(V, F, P)
    num_p, P = igl.extract_manifold_patches(F)

    num_p = max(P) + 1
    areas = np.zeros(num_p)

    # find slice with biggest area
    for curr_p in range(num_p):
        areas[curr_p] = igl.doublearea(V, F[P == curr_p, :]).sum()

    # decimate only biggest slice
    used_element = areas.argmax()

    curr_patch = F[P == used_element, :]
    curr_patch, c = igl.bfs_orient(curr_patch)

    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(V)
    mesh.triangles = o3d.utility.Vector3iVector(curr_patch)
    # remove unsused vertices
    mesh.remove_unreferenced_vertices()

    if curr_patch.shape[0] > 1:
        worked = False
        start_decimate_size = decimate_size
        while worked == False:
            worked, V_low, F_low, idx_F, idx_V = igl.decimate(
                V, curr_patch, start_decimate_size
            )
            start_decimate_size += 1
        assert worked
        # start_decimate_size += 1
        meshLow = o3d.geometry.TriangleMesh()
        meshLow.vertices = o3d.utility.Vector3dVector(V_low)
        meshLow.triangles = o3d.utility.Vector3iVector(F_low)
        assert meshLow.is_vertex_manifold()
        assert meshLow.is_edge_manifold()
        assert meshLow.is_orientable()
    else:
        decimate_size = curr_patch.shape[0]
        worked = True

    assert worked or V_low.shape[0] > 0
    basename = Path(shape_fp).stem
    feat_file_name = (
        "cuts_"
        + basename.split("-")[0]
        + "_shape_"
        + basename.split("-")[1]
        + ".mat"
    )
    featHighPath = (
        f"./features/{feat_file_name}"
    )
    featHigh = sio.loadmat(featHighPath)["featY"]

    LowToHighX = (
        NearestNeighbors(n_neighbors=1, algorithm="ball_tree")
        .fit(V)
        .kneighbors(V_low, return_distance=False)
    )
    feat_low = featHigh[LowToHighX.flatten(), :]
    extra_info_dict["features"] = feat_low

    extra_info_dict["added_verts"] = np.array(added_verts) + 1
    extra_info_dict["added_corres"] = added_corres
    extra_info_dict["num_p"] = num_p
    extra_info_dict["P"] = P
    extra_info_dict["F"] = curr_patch + 1
    extra_info_dict["V"] = V
    extra_info_dict["F_low"] = F_low.astype(np.int32) + 1
    extra_info_dict["V_low"] = V_low
    extra_info_dict["idx_V"] = idx_V + 1
    extra_info_dict["idx_F"] = idx_F + 1

    return V_low, F_low, extra_info_dict

if __name__ == "__main__":
    nr_faces = 100

    base_folder = "./data/"
    all_files = list(Path(base_folder).glob("*.off"))
    for _, curr_file in enumerate(tqdm(all_files)):
        test_path = f"{base_folder}/{curr_file.name}"
        stem_name = curr_file.stem
        V, F, extra_info_dict = downsample_mesh(
            test_path, nr_faces
        )
        save_folder = Path(
            "./low_res_meshes"
        )
        save_folder.mkdir(parents=True, exist_ok=True)
        # save to file
        sio.savemat(
            f"{save_folder}/{stem_name}_extra_info.mat",
            extra_info_dict,
        )
        o3d.io.write_triangle_mesh(
            f"{save_folder}/{stem_name}_extra_info.off",
            o3d.geometry.TriangleMesh(
                o3d.utility.Vector3dVector(V), o3d.utility.Vector3iVector(F)
            ),
        )