import numpy as np
import os
from steps import *

from shutil import copy2, move

from gremlin.read_sim_params import get_nml_params
from f90nml import write

tmplt_nml_path = "/data101/jlewis/sims/dust_fid/lvlmax_22/mh1e11/id292074"
ramses_exec_path = "/home/jlewis/ramses-yomp/bin/ramses_refmask_qhil3d"

ramses_exec = ramses_exec_path.split("/")[-1]

tgt_hid = 33051
tgt_pos = 3.123425e-02, 2.192202e-01, 3.037270e-01
tgt_rvir = 1.093000e-03
tgt_mvir = 5.152000e12
tgt_snap = 197

lvlmax = 20

params = {}
params["lvlmin"] = 7
params["lvlmax"] = lvlmax

# stuff for 200pc/lvlvmax=20 res
params["n_star"] = 1
params["sf_model"] = 0
params["eps_stars"] = 0.1

nml_name = "cosmo.nml"

mass_bin = f"mh1e{int(np.log10(tgt_mvir)):d}"

sim_path = "/data101/jlewis/sims/"
zoom_path = f"{sim_path}/dust_fid/lvlmax_{lvlmax:d}/{mass_bin}/"
zoom_name = f"id{tgt_hid}"

zoom_IC_path = os.path.join(sim_path, "ICs", zoom_name)

nnodes = 1
tpn = 128

r = 1.5 * tgt_rvir
os.makedirs(os.path.join(zoom_path, zoom_name), exist_ok=True)


# get region containing particles
get_refmask = False
# check if file exists for current hid?
if os.path.exists(f"./{tgt_hid:d}"):

    if not os.path.exists(f"./{tgt_hid:d}/music_region_file.txt"):
        get_refmask = True

else:

    os.makedirs(f"./{tgt_hid:d}")
    get_refmask = True

# if not run get_music_refmask
if get_refmask:
    run_get_music_refmask(
        os.path.join(HAGN_PATH, "output_00000"),
        os.path.join(HAGN_PATH, f"output_{tgt_snap:05d}"),
        tgt_pos[0],
        tgt_pos[1],
        tgt_pos[2],
        r,
    )

    # then mv to ./{tgt_hid}
    move("music_region_file.txt", f"./{tgt_hid:d}/.")


# print("ran get_music_refmask")

baryctr, rmax = get_zoom_region(tgt_hid)

print("got zoom region")

# make ics

# zoom ic levels
z_ic_lvls = [n for n in range(8, 12)]
# z_ic_lvls = [n for n in range(8, 10)]

for ilvl, lvl in enumerate(z_ic_lvls):

    print(f"making ICs for level {lvl}...")

    xzoom, yzoom, zzoom = np.int32(baryctr * 2**lvl)

    rzoom = (
        # np.int32(rmax + 2 * (2 + len(z_ic_lvls) - ilvl)) * 2
        np.int32(rmax * 2**lvl + 2 * (3 + len(z_ic_lvls) - ilvl))
    ) * 2  # 4 cells to include all ramses octs w particles

    # print(rzoom, rzoom / 2**lvl, rmax)

    ic_out_path = os.path.join(zoom_IC_path, f"{2**lvl:d}Cropped")
    if not os.path.exists(ic_out_path):
        os.makedirs(ic_out_path)

    if not os.path.exists(os.path.join(ic_out_path, "ic_deltab")):
        extract_grafic_call(
            os.path.join(HAGN_FID_IC_PATH, f"{2**lvl:d}"),
            ic_out_path,
            np.int32([xzoom, yzoom, zzoom]),
            int(rzoom),
        )

print("made ICs")

# nml

nml = get_nml_params(tmplt_nml_path, "cosmo.nml")

zoom_nml(nml, baryctr, rmax)
zoom_ic_nml(nml, HAGN_FID_IC_PATH, os.path.join(zoom_path, zoom_name), z_ic_lvls)
apply_var_params(nml, params)

write(nml, os.path.join(zoom_path, zoom_name, nml_name), force=True)

print("wrote nml")

# submission script
create_sh(
    "run.sh",
    os.path.join(zoom_path, zoom_name),
    ramses_exec=ramses_exec,
    nml=nml_name,
    ntasks=nnodes * tpn,
    nnodes=nnodes,
)

print("wrote submission script")

copy2(ramses_exec_path, os.path.join(zoom_path, zoom_name))

print("copied ramses exec")

copy2("./ramses_swind_Sikey.dat", os.path.join(zoom_path, zoom_name))

print("done!")

print("sim setup at %s" % (os.path.join(zoom_path, zoom_name)))