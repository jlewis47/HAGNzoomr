import cmd
from math import e
from tkinter import Y
import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse

# from shutil import copyfile

HAGN_PATH = "/data52/Horizon-AGN/OUTPUT_DIR/"
# HAGN_PATH = "/data82/dubois/Deuterium/Horizon-noAGN_512_Deuterium/OUTPUT_DIR/"
HAGN_FID_IC_PATH = "/data102/dubois/Lewis/ExtractZoom/TestRun/FiducialIC/"


def mvee(points, tol=0.001):
    """
    Find the minimum volume ellipse.
    Return A, c where the equation for the ellipse given in "center form" is
    (x-c).T * A * (x-c) = 1
    credit https://stackoverflow.com/questions/14016898/port-matlab-bounding-ellipsoid-code-to-python
    """
    points = np.asmatrix(points)
    N, d = points.shape
    Q = np.column_stack((points, np.ones(N))).T
    err = tol + 1.0
    u = np.ones(N) / N
    while err > tol:
        # assert u.sum() == 1 # invariant
        X = Q * np.diag(u) * Q.T
        M = np.diag(Q.T * np.linalg.inv(X) * Q)
        jdx = np.argmax(M)
        step_size = (M[jdx] - d - 1.0) / ((d + 1) * (M[jdx] - 1.0))
        new_u = (1 - step_size) * u
        new_u[jdx] += step_size
        err = np.linalg.norm(new_u - u)
        u = new_u
    c = u * points
    A = np.linalg.inv(points.T * np.diag(u) * points - c.T * c) / d
    return np.asarray(A), np.squeeze(np.asarray(c))


def run_get_music_refmask(output_stt, output_end, x, y, z, r):
    # run get_music_refmask -ini output_stt -inf output_end -xc x -yc y -zc z -rad r
    run_cmd = f"./get_music_refmask -ini {output_stt} -inf {output_end} -xc {x:f} -yc {y:f} -zc {z:f} -rad {r:f}"
    print(run_cmd)
    # os.system(run_cmd)
    with open("get_music_refmask.log", "w") as f:
        subprocess.run(run_cmd.split(), stdout=f, stderr=f)

    return 1


def load_refmask(hid):
    # load refmask
    # pos_final = np.loadtxt(f"./{hid}/music_region_file.txt_final")
    pos = np.loadtxt(f"./{hid}/music_region_file.txt")

    # return pos, pos_final
    return pos


def correct_periodic(pos, med):
    too_large = np.abs(pos - med) > 0.5

    up = pos > 0.5
    dw = pos <= 0.5

    # print(pos)
    # print(pos[too_large])

    pos[too_large * up] -= 1
    pos[too_large * dw] += 1

    return pos


def get_zoom_region(hid, zoom_type="ellipsoid"):

    assert zoom_type in [
        "ellipsoid",
        "sphere",
    ], "type must be either 'ellipsoid' or 'sphere'"

    # get zoom region
    stt_pos = load_refmask(hid)
    # stt_pos, fin_pos = load_refmask(hid)

    med = np.median(stt_pos, axis=0)

    stt_pos = correct_periodic(stt_pos, med)
    # fin_pos = correct_periodic(fin_pos, med)

    # from internet random
    baryctr = np.mean(stt_pos, axis=0)
    coords = stt_pos - baryctr

    # coords = np.asarray(
    #     [
    #         np.random.normal(scale=1, size=500),
    #         np.random.normal(scale=4, size=500),
    #         np.random.normal(scale=2, size=500),
    #     ]
    # ).T

    inertia = np.dot(coords.transpose(), coords)
    e_values, e_vectors = np.linalg.eig(inertia)

    order = np.argsort(e_values)  # [::-1]
    # eval3, eval2, eval1 = e_values[order]
    axis3, axis2, axis1 = e_vectors[:, order].transpose()

    # normalize the eigenvectors
    axis1 /= np.linalg.norm(axis1)
    axis2 /= np.linalg.norm(axis2)
    axis3 /= np.linalg.norm(axis3)

    # now diagonilize tensor to get rotation matrix and eigenvalues
    # 3d recap plot
    # fig = plt.figure()
    # # ax = fig.add_subplot(111, projection="3d")
    # fig, ax = plt.subplots(3, 1, figsize=(8, 8), sharex=True, sharey=True)

    # ax[0].scatter(coords[:, 0], coords[:, 1], alpha=0.01)
    # ax[1].scatter(coords[:, 1], coords[:, 2], alpha=0.01)
    # ax[2].scatter(coords[:, 0], coords[:, 2], alpha=0.01)

    # # print(stt_pos[too_large])

    # ax.quiver(
    #     baryctr[0],
    #     baryctr[1],
    #     baryctr[2],
    #     axis1[0],
    #     axis1[1],
    #     axis1[2],
    #     color="r",
    #     # normalize=True,
    # )
    # ax.quiver(
    #     baryctr[0],
    #     baryctr[1],
    #     baryctr[2],
    #     axis2[0],
    #     axis2[1],
    #     axis2[2],
    #     color="g",
    #     # normalize=True,
    # )
    # ax.quiver(
    #     baryctr[0],
    #     baryctr[1],
    #     baryctr[2],
    #     axis3[0],
    #     axis3[1],
    #     axis3[2],
    #     color="b",
    #     # normalize=True,
    # )

    rot_m = np.asarray([axis1, axis2, axis3])

    rot_coords = np.dot(coords, rot_m)

    # # ax.scatter(rot_coords[:, 0], rot_coords[:, 1], rot_coords[:, 2], alpha=0.01)
    # ax[0].scatter(rot_coords[:, 0], rot_coords[:, 1], alpha=0.01)
    # ax[1].scatter(rot_coords[:, 1], rot_coords[:, 2], alpha=0.01)
    # ax[2].scatter(rot_coords[:, 0], rot_coords[:, 2], alpha=0.01)

    rot_xctr = 0.5 * (rot_coords[:, 0].max() + rot_coords[:, 0].min())
    rot_yctr = 0.5 * (rot_coords[:, 1].max() + rot_coords[:, 1].min())
    rot_zctr = 0.5 * (rot_coords[:, 2].max() + rot_coords[:, 2].min())

    rot_ctr = np.asarray([rot_xctr, rot_yctr, rot_zctr])

    # back to cartesian
    new_ctr = np.dot(rot_ctr, rot_m.T)

    # print(baryctr, new_ctr)

    # ax.scatter([0], [0], [0], s=50, marker="x", c="r")
    # ax.scatter(new_ctr[0], new_ctr[1], new_ctr[2], s=50, marker="+", c="r")
    # ax[0].scatter([0], [0], s=50, marker="x", c="r")
    # ax[0].scatter(new_ctr[0], new_ctr[1], s=50, marker="+", c="r")
    # ax[1].scatter([0], [0], marker="x", c="r")
    # ax[1].scatter(new_ctr[1], new_ctr[2], s=50, marker="+", c="r")
    # ax[2].scatter([0], [0], marker="x", c="r")
    # ax[2].scatter(new_ctr[0], new_ctr[2], s=50, marker="+", c="r")

    # # get max distance from barycentre
    # rmax = np.max(np.abs(coords))
    # print(rmax)
    # print(np.max(np.abs(rot_coords - rot_ctr)))
    # # from new centre
    new_ctr_coords = -coords + new_ctr

    # if zoom_type == "sphere":
    # print(np.linalg.norm(new_ctr_coords, axis=1).shape)
    rmax = np.max(np.linalg.norm(new_ctr_coords, axis=1))
    # print(new_ctr_coords.shape)
    # else:
    # A, c = mvee(new_ctr_coords[::100, :])

    # print(A)
    # print(A, c)

    # print(A.diagonal())
    # amax, bmax, cmax = 1.0 / A.diagonal() ** 0.5

    amax = np.abs(new_ctr_coords[:, 0].max() - new_ctr_coords[:, 0].min()) * 0.5
    bmax = np.abs(new_ctr_coords[:, 1].max() - new_ctr_coords[:, 1].min()) * 0.5
    cmax = np.abs(new_ctr_coords[:, 2].max() - new_ctr_coords[:, 2].min()) * 0.5

    # # first guesses...

    # rr = (
    #     new_ctr_coords[:, 0] ** 2 / amax**2
    #     + new_ctr_coords[:, 1] ** 2 / bmax**2
    #     + new_ctr_coords[:, 2] ** 2 / cmax**2
    # )

    # # print(rr.min(), rr.max())

    if zoom_type == "ellipsoid":

        niter = 0
        last = 0
        step_mult = 1.001
        while np.any(
            new_ctr_coords[:, 0] ** 2 / amax**2
            + new_ctr_coords[:, 1] ** 2 / bmax**2
            + new_ctr_coords[:, 2] ** 2 / cmax**2
            > 1
        ):
            if last == 0:
                amax *= step_mult
                last = 1
            elif last == 1:
                bmax *= step_mult
                last = 2
            elif last == 2:
                cmax *= step_mult
                last = 0

            # print(amax, bmax, cmax)

            niter += 1

        if niter > 0:
            if last == 1:
                amax /= step_mult
            if last == 2:
                bmax /= step_mult
            elif last == 0:
                cmax /= step_mult

        # amax = 0.09
        # bmax = 0.08
        # cmax = 0.1

        rr = (
            new_ctr_coords[:, 0] ** 2 / amax**2
            + new_ctr_coords[:, 1] ** 2 / bmax**2
            + new_ctr_coords[:, 2] ** 2 / cmax**2
        )

        # print(rr.min(), rr.max())

        print(f"After {niter:d} grow iterations, I found:")
        print(amax, bmax, cmax)

    # ax[0].scatter(new_ctr_coords[:, 0], new_ctr_coords[:, 1], alpha=0.01)  # , c=rr)
    # ax[1].scatter(new_ctr_coords[:, 1], new_ctr_coords[:, 2], alpha=0.01)  # , c=rr)
    # ax[2].scatter(new_ctr_coords[:, 0], new_ctr_coords[:, 2], alpha=0.01)  # , c=rr)

    # xyell = Ellipse(
    #     (new_ctr[0], new_ctr[1]), 2 * amax, 2 * bmax, fill=False, edgecolor="r"
    # )
    # xzell = Ellipse(
    #     (new_ctr[0], new_ctr[2]), 2 * amax, 2 * cmax, fill=False, edgecolor="r"
    # )
    # yzell = Ellipse(
    #     (new_ctr[1], new_ctr[2]), 2 * bmax, 2 * cmax, fill=False, edgecolor="r"
    # )

    # ax[0].add_patch(xyell)
    # ax[1].add_patch(yzell)
    # ax[2].add_patch(xzell)

    # xycirl = Circle((new_ctr[0], new_ctr[1]), rmax, fill=False, edgecolor="k")
    # xzcirl = Circle((new_ctr[0], new_ctr[2]), rmax, fill=False, edgecolor="k")
    # yzcirl = Circle((new_ctr[1], new_ctr[2]), rmax, fill=False, edgecolor="k")

    # ax[0].add_patch(xycirl)
    # ax[1].add_patch(yzcirl)
    # ax[2].add_patch(xzcirl)

    # print(new_ctr_coords[0].min(), new_ctr_coords[0].max())

    # fig.savefig("eigen.png")

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection="3d")

    # ax.scatter(new_ctr[0], new_ctr[1], new_ctr[2], s=50, marker="+", c="r")

    # ax.scatter(
    #     new_ctr_coords[:, 0],
    #     new_ctr_coords[:, 1],
    #     new_ctr_coords[:, 2],
    #     alpha=0.01,
    #     # c=rr,
    # )

    # ax.plot(
    #     [new_ctr[0] - amax, new_ctr[0] + amax],
    #     [new_ctr[1], new_ctr[1]],
    #     [new_ctr[2], new_ctr[2]],
    #     c="r",
    # )
    # ax.plot(
    #     [new_ctr[0], new_ctr[0]],
    #     [new_ctr[1] - bmax, new_ctr[1] + bmax],
    #     [new_ctr[2], new_ctr[2]],
    #     c="r",
    # )
    # ax.plot(
    #     [new_ctr[0], new_ctr[0]],
    #     [new_ctr[1], new_ctr[1]],
    #     [new_ctr[2] - cmax, new_ctr[2] + cmax],
    #     c="r",
    # )

    # ax.view_init(-45, -45, -45)

    # fig.savefig("eigen3d.png")
    # # print(rmax)
    # ax.view_init(45, 150, 45)
    # fig.savefig("eigen3d_perp.png")

    # print(np.dot(rot_coords, rot_m.T) / coords)

    # print(f"starting centre: {baryctr}")
    # print(f"new centre: {new_ctr}")
    # print(f"final centre: {bary_fin}")

    # dplcmt = np.linalg.norm(bary_fin - baryctr)
    # print(f"start to end displacement: {dplcmt:f}")

    # rmax = np.max(np.abs(stt_pos - baryctr))
    # rmax = np.max(np.abs(stt_pos - new_ctr))

    # print(rmax, np.max(np.abs(stt_pos - baryctr)))
    # print(new_ctr, baryctr)
    # print(rmax, np.max(np.abs(stt_pos - baryctr)))
    # rmax = r
    # print(median_ctr, baryctr, rmax)
    # print(median_ctr, baryctr, rmax)

    # # reflect back into box if overstep
    # for idim in range(3):
    #     if baryctr[idim] > 1.0:
    #         baryctr[idim] -= 1
    #     elif baryctr[idim] < 0:
    #         baryctr[idim] += 1
    for idim in range(3):
        if new_ctr[idim] > 1.0:
            new_ctr[idim] -= 1
        elif new_ctr[idim] < 0:
            new_ctr[idim] += 1

    # print("")

    # print(baryctr, new_ctr, rot_ctr)

    # return baryctr, rmax

    if zoom_type == "ellipsoid":
        return new_ctr, amax, bmax, cmax
    else:
        return new_ctr, rmax
    # return [], []


def zoom_nml(nml, baryctr, rmax):
    nml["REFINE_PARAMS"]["xzoom"] = baryctr[0]
    nml["REFINE_PARAMS"]["yzoom"] = baryctr[1]
    nml["REFINE_PARAMS"]["zzoom"] = baryctr[2]
    nml["REFINE_PARAMS"]["rzoom"] = rmax


def zoom_ellipsoid_nml(nml, baryctr, amax, bmax, cmax):
    nml["REFINE_PARAMS"]["xzoom"] = baryctr[0]
    nml["REFINE_PARAMS"]["yzoom"] = baryctr[1]
    nml["REFINE_PARAMS"]["zzoom"] = baryctr[2]
    nml["REFINE_PARAMS"]["azoom"] = amax
    nml["REFINE_PARAMS"]["bzoom"] = bmax
    nml["REFINE_PARAMS"]["czoom"] = cmax


def zoom_ic_nml(nml, zoom_path, lvls):

    # fid_fname = os.path.join(fid_path, f"{int(2**(lvls[0]-1)):d}")

    lvl_fnames = []
    for lvl in lvls:

        lvl_res = 2**lvl

        # lvl_fnames.append(os.path.join("./", "ZoomIC", f"{lvl_res:d}Cropped"))
        lvl_fnames.append(os.path.join(zoom_path, f"{lvl_res:d}Cropped"))

    # fnames = [fid_fname] + lvl_fnames

    del nml["init_params"]["initfile"]

    for ifname, fname in enumerate(lvl_fnames):

        nml["init_params"][f"initfile({ifname+1:d})"] = fname


def apply_var_params(nml, params):
    found = np.zeros(len(params.keys()))
    for ikey, (key, val) in enumerate(params.items()):
        for nml_key in nml.keys():
            if key in nml[nml_key]:
                nml[nml_key][key] = val
                found[ikey] = 1

    if np.any(found == 0):
        print(
            "I didn't find all parameters in the nml files... please check their spelling and values"
        )
        print(",".join(np.asarray(list(params.keys()))[found == 0]))
        assert False, "missing parameters in nml file"

    return 1


def extract_grafic_call(fin, fout, baryctr, rmax):
    # print(fin, fout, baryctr, rmax)
    cmd = f"./extract_grafic {fin} {fout} {baryctr[0]:d} {baryctr[1]:d} {baryctr[2]:d} {rmax:d}"

    print(cmd)
    # os.system(cmd)
    with open("extract_grafic.log", "w") as f:
        subprocess.run(cmd.split(), stdout=f, stderr=f)

    return 1


def extract_grafic_cuboid_call(fin, fout, baryctr, ellipsoid):
    (amax, bmax, cmax) = ellipsoid
    # print(fin, fout, baryctr, rmax)
    cmd = f"./extract_grafic_cuboid {fin} {fout} {baryctr[0]:d} {baryctr[1]:d} {baryctr[2]:d} {amax:d}, {bmax:d}, {cmax:d}"

    print(cmd)
    # os.system(cmd)
    with open("extract_grafic.log", "w") as f:
        subprocess.run(cmd.split(), stdout=f, stderr=f)

    return 1


def centre_grafic_call(fin, fout, ctr):
    # print(fin, fout, ctr)
    cmd = f"./centre_grafic {fin} {fout} {ctr[0]:d} {ctr[1]:d} {ctr[2]:d}"

    print(cmd)
    # os.system(cmd)
    with open("centre_grafic.log", "w") as f:
        subprocess.run(cmd.split(), stdout=f, stderr=f)

    return 1


def create_sh(
    tmplt,
    tgt_path,
    pbs_params=None,
    ramses_exec=None,
    nml=None,
    ntasks=None,
    nnodes=None,
    wt=None,
    name=None,
):

    lines = []

    with open(tmplt, "r") as f:
        lines = f.readlines()

    cd_line = np.where(["cd" in l for l in lines])[0][0]

    lines[cd_line] = f"cd {tgt_path}\n"

    if pbs_params is not None:
        last_pbs_line = np.where(["#PBS" in l for l in lines])[0][-1]
        nb_added = 0
        for key, val in pbs_params.items():
            if key in lines:
                continue
            lines.insert(last_pbs_line + nb_added, f"#PBS -{key} {val}\n")
            nb_added += 1

    exec_line = np.where([".nml" in l for l in lines])[0][0]

    cmd_words = lines[exec_line].split(" ")
    nml_pos = np.where([".nml" in w for w in cmd_words])[0][0]

    if nml is not None:
        cmd_words[nml_pos] = nml

    if ramses_exec is not None:
        cmd_words[nml_pos - 1] = "./" + ramses_exec

    if ntasks is not None:
        cmd_words[nml_pos - 2] = str(ntasks)

    # print(cmd_words)

    lines[exec_line] = " ".join(cmd_words)

    if nnodes is not None or wt is not None:
        if wt == None:
            wt = "48:00:00"
        if nnodes == None:
            nnodes = 1

        node_line = np.where(["nodes=" in l for l in lines])[0][0]
        lines[node_line] = f"#PBS -l nodes={nnodes}:ppn=128,walltime={wt}\n"

    if name is not None:
        name_line = np.where(["#PBS -N" in l for l in lines])[0][0]
        lines[name_line] = f"#PBS -N {name}\n"

    with open(os.path.join(tgt_path, "run.sh"), "w") as f:
        f.writelines(lines)
