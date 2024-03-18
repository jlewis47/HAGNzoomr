import numpy as np
import os
from steps import *

from shutil import copy2, move, rmtree

from gremlin.read_sim_params import get_nml_params
from f90nml import write

from hagn.catalogues import make_super_cat, get_cat_hids

# need to add support for halos on edge of box... i.e. recentre to 0.5,0.5,0.5
# can use centre_grafic.f90 from ramses/utils/f90

tmplt_nml_path = "/data101/jlewis/sims/dust_fid/lvlmax_22/mh1e11/id292074"
ramses_exec_path = "/home/jlewis/ramses-yomp-ellipsoid/bin/ramses_refmask_qhil3d"

ramses_exec = ramses_exec_path.split("/")[-1]
tgt_snap = 197
# tgt_snap = 15

supercat = make_super_cat(197, "/data101/jlewis/hagn/super_cats")

# tgt_hid = 74099
tgt_hid = 147479
# tgt_hid = 13310
# tgt_hid = 68373
# tgt_hid = 326421
# tgt_hid = 194228
# tgt_hid = 242704
pties = get_cat_hids(supercat, [tgt_hid])
tgt_pos = np.array([pties["hx"][0], pties["hy"][0], pties["hz"][0]])
tgt_rvir = pties["rvir"][0]
tgt_mvir = pties["mhalo"][0]


overwrite = True


params = {}

##################
# .nml params
##################

lvlmax = 20
lvlmin = 7
params["levelmin"] = lvlmin
params["levelmax"] = lvlmax

# stuff for 200pc/lvlvmax=20 res
params["n_star"] = 2  # cc
params["sf_model"] = 0
params["eps_star"] = 0.1

params["foutput"] = 25

params["n_sink"] = 2  # cc
params["ns_sink"] = 5  # cc

params["ngridmax"] = 750000
params["npartmax"] = 2500000
params["nsinkmax"] = 1000

nml_name = "cosmo.nml"

mass_bin = f"mh1e{int(np.log10(tgt_mvir)):d}"

sim_path = "/data101/jlewis/sims/"
zoom_path = os.path.join(sim_path, f"dust_fid/lvlmax_{lvlmax:d}", mass_bin)
zoom_name = f"id{tgt_hid}"

zoom_IC_path = os.path.join(sim_path, "ICs", zoom_name)

##################
# qsub params
##################

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
    move("music_region_file.txt", f"./{tgt_hid:d}/music_region_file.txt")
    # move("music_region_file.txt_final", f"./{tgt_hid:d}/music_region_file.txt_final")


# print("ran get_music_refmask")

baryctr, amax, bmax, cmax = get_zoom_region(tgt_hid, zoom_type="ellipsoid")
# baryctr, rmax = get_zoom_region(tgt_hid, zoom_type="sphere")
# rmax = np.max([amax, bmax, cmax])
# cmax = 0.09


print(baryctr)
print(f"amax:{amax:f}, bmax:{bmax:f}, cmax:{cmax:f} r:{r:f}, tgt_rvir:{tgt_rvir:f}")

print("got zoom region")


zoom_ctr = np.copy(baryctr)

centre = True

# print([baryctr - 2 * rmax < 0, baryctr + 2 * rmax > 1])

if np.any(
    [
        baryctr[0] + amax * 2 > 1,
        baryctr[0] - amax * 2 < 0,
        baryctr[1] + bmax * 2 > 1,
        baryctr[1] - bmax * 2 < 0,
        baryctr[2] + cmax * 2 > 1,
        baryctr[2] - cmax * 2 < 0,
    ]
):
    centre = True

# make ics


# zoom ic levels
z_ic_lvls = [n for n in range(lvlmin, 12)]
# z_ic_lvls = [n for n in range(lvlmin, 10)]
# z_ic_lvls = [n for n in range(8, 10)]

rmult = 2

for ilvl, lvl in enumerate(z_ic_lvls):

    print(f"making ICs for level {lvl}...")

    xzoom, yzoom, zzoom = np.int32(zoom_ctr * 2**lvl)

    boost = 2 * (2 + len(z_ic_lvls) - ilvl)
    azoom = (
        # np.int32(rmax + 2 * (2 + len(z_ic_lvls) - ilvl)) * 2
        np.int32(amax * 2**lvl + boost)
    ) * rmult  # 4 cells to include all ramses octs w particles
    bzoom = (
        # np.int32(rmax + 2 * (2 + len(z_ic_lvls) - ilvl)) * 2
        np.int32(bmax * 2**lvl + boost)
    ) * rmult  # 4 cells to include all ramses octs w particles
    czoom = (
        # np.int32(rmax + 2 * (2 + len(z_ic_lvls) - ilvl)) * 2
        np.int32(cmax * 2**lvl + boost)
    ) * rmult  # 4 cells to include all ramses octs w particles

    # rzoom = (
    #     np.int32(rmax * 2**lvl + boost)
    # ) * rmult  # 4 cells to include all ramses octs w particles

    # print(rzoom, rzoom / 2**lvl, rmax)
    # print(boost)
    # print(amax, bmax, cmax)
    # print(azoom, bzoom, czoom)

    ic_out_path = os.path.join(zoom_IC_path, f"{2**lvl:d}Cropped")
    if not os.path.exists(ic_out_path):
        os.makedirs(ic_out_path)

    ctr_out_path = os.path.join(zoom_IC_path, f"{2**lvl:d}Ctr")
    if not os.path.exists(ctr_out_path):
        os.makedirs(ctr_out_path)

    if not os.path.exists(os.path.join(ic_out_path, "ic_deltab")) or overwrite:

        src_dir = os.path.join(HAGN_FID_IC_PATH, f"{2**lvl:d}")

        # if ilvl == 0:  # first level covers whole box - just copy as is

        # copy all files in src_dir to ic_out_path, then we might centre
        # then we get the cuboid

        extract_path = src_dir
        if centre:
            centre_grafic_call(src_dir, ctr_out_path, [xzoom, yzoom, zzoom])
            xzoom, yzoom, zzoom = np.int32(np.asarray([0.5, 0.5, 0.5]) * 2**lvl)
            extract_path = ctr_out_path

        else:
            ftocp = os.listdir(src_dir)
            for f in ftocp:
                src_f = os.path.join(src_dir, f)
                copy2(src_f, ic_out_path)

        if ilvl > 0:

            extract_grafic_cuboid_call(
                extract_path,
                ic_out_path,
                np.int32([xzoom, yzoom, zzoom]),
                np.int32([azoom, bzoom, czoom]),
            )

            # extract_grafic_call(
            #     extract_path,
            #     ic_out_path,
            #     np.int32([xzoom, yzoom, zzoom]),
            #     int(rzoom),
            # )

            # remove centered file
            if centre:
                rmtree(ctr_out_path)


print("made ICs")


# nml

nml = get_nml_params(tmplt_nml_path, "cosmo.nml")

if centre:
    zoom_ctr = np.array([0.5, 0.5, 0.5])

# zoom_nml(nml, zoom_ctr, rmax)
zoom_ellipsoid_nml(nml, zoom_ctr, amax, bmax, cmax)
# zoom_ic_nml(nml, HAGN_FID_IC_PATH, zoom_IC_path, z_ic_lvls)
zoom_ic_nml(nml, zoom_IC_path, z_ic_lvls)
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
    name=f"zoom_{zoom_name}",
)

print("wrote submission script")

copy2(ramses_exec_path, os.path.join(zoom_path, zoom_name))

print("copied ramses exec")

copy2("./ramses_swind_Sikey.dat", os.path.join(zoom_path, zoom_name))

print("done!")

print("sim setup at %s" % (os.path.join(zoom_path, zoom_name)))
