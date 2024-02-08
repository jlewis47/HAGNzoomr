import cmd
from math import e
import numpy as np
import os
import subprocess

# from shutil import copyfile

HAGN_PATH = "/data52/Horizon-AGN/OUTPUT_DIR/"
HAGN_FID_IC_PATH = "/data102/dubois/Lewis/ExtractZoom/TestRun/FiducialIC/"


def run_get_music_refmask(output_stt, output_end, x, y, z, r):
    # run get_music_refmask -ini output_stt -inf output_end -xc x -yc y -zc z -rad r
    run_cmd = f"./get_music_refmask -ini {output_stt} -inf {output_end} -xc {x:f} -yc {y:f} -zc {z:f} -rad {r:f}"
    print(run_cmd)
    # os.system(run_cmd)
    with open("get_music_refmask.log", "w") as f:
        subprocess.run(run_cmd.split(), stdout=f, stderr=f)

    return 1


def load_refmask():
    # load refmask
    pos = np.loadtxt("music_region_file.txt")

    return pos


def get_zoom_region():
    # get zoom region
    stt_pos = load_refmask()

    median_ctr = np.median(stt_pos, axis=0)

    too_large = np.abs(stt_pos - median_ctr) > 0.5

    up = stt_pos[too_large] > 0.5
    dw = stt_pos[too_large] <= 0.5

    stt_pos[too_large][up] -= 1
    stt_pos[too_large][dw] += 1

    baryctr = np.mean(stt_pos, axis=0)

    rmax = np.max(np.abs(stt_pos - baryctr))

    # reflect back into box if overstep
    for idim in range(3):
        if baryctr[idim] > 1.0:
            baryctr[idim] -= 1
        elif baryctr[idim] < 0:
            baryctr[idim] += 1

    print("")

    return baryctr, rmax


def zoom_nml(nml, baryctr, rmax):
    nml["REFINE_PARAMS"]["xzoom"] = baryctr[0]
    nml["REFINE_PARAMS"]["yzoom"] = baryctr[1]
    nml["REFINE_PARAMS"]["zzoom"] = baryctr[2]
    nml["REFINE_PARAMS"]["rzoom"] = rmax


def zoom_ic_nml(nml, fid_path, zoom_path, lvls):

    fid_fname = os.path.join(fid_path, f"{int(2**(lvls[0]-1)):d}")

    lvl_fnames = []
    for lvl in lvls:

        lvl_res = 2**lvl

        lvl_fnames.append(os.path.join(zoom_path, "ZoomIC", f"{lvl_res:d}Cropped"))

    fnames = [fid_fname] + lvl_fnames

    del nml["init_params"]["initfile"]

    for ifname, fname in enumerate(fnames):

        nml["init_params"][f"initfile({ifname+1:d})"] = fname


def apply_var_params(nml, params):
    for key, val in params.items():
        for nml_key in nml.keys():
            if key in nml[nml_key]:
                nml[nml_key][key] = val

    return 1


def extract_grafic_call(fin, fout, baryctr, rmax):
    # print(fin, fout, baryctr, rmax)
    cmd = f"/home/jlewis/extract_grafic/extract_grafic_notinteractive {fin} {fout} {baryctr[0]:d} {baryctr[1]:d} {baryctr[2]:d} {rmax:d}"

    print(cmd)
    # os.system(cmd)
    with open("extract_grafic.log", "w") as f:
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
        lines[node_line] = f"#PBS -l nodes={nnodes}:ppn=128,walltime={wt}"

    with open(os.path.join(tgt_path, "run.sh"), "w") as f:
        f.writelines(lines)
