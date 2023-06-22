from typing import Optional

import MDAnalysis as mda
import pandas as pd

import simulation_names as sims


def load_simulations(root_path: str) -> dict[str, list[mda.Universe]]:
    return {
        sims.PfHT_MMV12: [
            mda.Universe(
                f"{root_path}/MMV12/replica_{i}/production/PfHT_MMV12.em.{i}.0_200ns.protonly.tpr",
                f"{root_path}/EM_MMV/MMV12/replica_{i}/production/PfHT_MMV12.em.{i}.0_1000ns.skip250.protonly.xtc",
            )
            for i in range(1, 4)
        ],
        sims.W4126_MMV12: [
            mda.Universe(
                f"{root_path}/EM_MMV/W412A_MMV12/replica_{i}/production/W412A_MMV12.em.{i}.0_200ns.protonly.tpr",
                f"{root_path}/EM_MMV/W412A_MMV12/replica_{i}/production/W412A_MMV12.em.{i}.0_1000ns.skip250.protonly.xtc",
            )
            for i in range(1, 4)
        ],
    }


def get_trajectory_lengths(
    sim_unis: dict[str, list[mda.Universe]]
) -> dict[str, list[int]]:
    out = {}
    for simulation, replicas in sim_unis.items():
        out[simulation] = [len(r.trajectory) for r in replicas]
    return out


def load_em_unis() -> tuple[dict[str, list[mda.Universe]], dict[str, list[int]]]:
    """
    will load mda unis of all em analyses so far
    Not so space intensive since I don't save this to memory with mda
    will return a dict with the key being the uni name and the values
    being a list of the replicas starting with 1

    for now, protonly and TPR
    which means that numbering is incorrect, you'll have to add 21

    Note that each simulation is a different length.
    A dict of this information is also returned
    """

    all_uni = {}
    all_len = {}

    ### MMV12, W412A, MMV8, GLUT1 replicas first. All near 1000ns so we can loop
    uni_wt_mmv12 = []
    uni_w412a_mmv12 = []
    uni_wt_mmv8 = []
    uni_glut1_mmv12 = []
    len_wt_mmv12 = []
    len_w412a_mmv12 = []
    len_wt_mmv8 = []
    len_glut1_mmv12 = []
    for i in range(1, 4):
        u = mda.Universe(
            f"../../EM_MMV/MMV12/replica_{i}/production/PfHT_MMV12.em.{i}.0_200ns.protonly.tpr",
            f"../../EM_MMV/MMV12/replica_{i}/production/PfHT_MMV12.em.{i}.0_1000ns.skip250.protonly.xtc",
        )
        uni_wt_mmv12.append(u)
        len_wt_mmv12.append(len(u.trajectory))

        u = mda.Universe(
            f"../../EM_MMV/W412A_MMV12/replica_{i}/production/W412A_MMV12.em.{i}.0_200ns.protonly.tpr",
            f"../../EM_MMV/W412A_MMV12/replica_{i}/production/W412A_MMV12.em.{i}.0_1000ns.skip250.protonly.xtc",
        )
        uni_w412a_mmv12.append(u)
        len_w412a_mmv12.append(len(u.trajectory))

        u = mda.Universe(
            f"../../EM_MMV/MMV8/replica_{i}/production/PfHT_MMV8.em.{i}.0_200ns.protonly.tpr",
            f"../../EM_MMV/MMV8/replica_{i}/production/PfHT_MMV8.em.{i}.0_1000ns.skip250.protonly.xtc",
        )
        uni_wt_mmv8.append(u)
        len_wt_mmv8.append(len(u.trajectory))

        u = mda.Universe(
            f"../../EM_MMV/GLUT1_MMV12/replica_{i}/production/GLUT1_MMV12.em.{i}.0_200ns.protonly.tpr",
            f"../../EM_MMV/GLUT1_MMV12/replica_{i}/production/GLUT1_MMV12.em.{i}.0_1000ns.skip250.protonly.xtc",
        )
        uni_glut1_mmv12.append(u)
        len_glut1_mmv12.append(len(u.trajectory))

    # apo is mixed. No need for a new list name since we finished above, but might as
    # well keep it easy to read if I need to change later
    u1 = mda.Universe(
        "../../EM_MMV/PfHT_apo/replica_1/production/PfHT_apo.em.1.0_200ns.protonly.tpr",
        "../../EM_MMV/PfHT_apo/replica_1/production/PfHT_apo.em.1.0_800ns.skip250.protonly.xtc",
    )
    u2 = mda.Universe(
        "../../EM_MMV/PfHT_apo/replica_2/production/PfHT_apo.em.2.0_200ns.protonly.tpr",
        "../../EM_MMV/PfHT_apo/replica_2/production/PfHT_apo.em.2.0_800ns.skip250.protonly.xtc",
    )
    u3 = mda.Universe(
        "../../EM_MMV/PfHT_apo/replica_3/production/PfHT_apo.em.3.0_200ns.protonly.tpr",
        "../../EM_MMV/PfHT_apo/replica_3/production/PfHT_apo.em.3.0_600ns.skip250.protonly.xtc",
    )
    uni_wt_apo = [u1, u2, u3]
    len_wt_apo = [len(u1.trajectory), len(u2.trajectory), len(u3.trajectory)]

    all_uni[sims.PfHT_MMV12] = uni_wt_mmv12
    all_len[sims.PfHT_MMV12] = len_wt_mmv12
    all_uni[sims.W4126_MMV12] = uni_w412a_mmv12
    all_len[sims.W4126_MMV12] = len_w412a_mmv12
    all_uni[sims.PfHT_MMV8] = uni_wt_mmv8
    all_len[sims.PfHT_MMV8] = len_wt_mmv8
    all_uni[sims.GLUT1_MMV12] = uni_glut1_mmv12
    all_len[sims.GLUT1_MMV12] = len_glut1_mmv12
    all_uni[sims.PfHT_apo] = uni_wt_apo
    all_len[sims.PfHT_apo] = len_wt_apo

    ### W412A replicas next
    # uni_l = []
    return all_uni, all_len


def get_fp_dataframe(
    u: mda.Universe,
    skip: int,
    ligname: str = "resname MMV",
    run_calc: bool = True,
    filename: Optional[str] = None,
) -> pd.DataFrame:
    """
    Here I just run prolif fingerprints on a uni
    One can also save the filename, which will save by default in '../fingerprints_df
    If you don't want to perform the calculation, you can set this to "False" and write
    the filename you wish to load
    """
    import prolif as plf
    import os

    if run_calc:
        ligand = u.select_atoms(ligname)
        protein = u.select_atoms("protein")
        fp = plf.Fingerprint()
        if skip != 1:
            fp.run(u.trajectory[::skip], ligand, protein)
        else:  # MDA does not like if skip == 1 and you don't provide start and stop
            fp.run(u.trajectory, ligand, protein)
        df = fp.to_dataframe()
        df = df.droplevel(
            "ligand", axis=1
        )  ##no need for keeping ligand axis, everything is MMV
        if filename is not None:
            pd.to_pickle(df, f"../fingerprints_df/{filename}.pkl")
        else:
            print("no filename given, returning dataframe, not saving file...")
    else:
        if filename is not None:
            df = pd.read_pickle(f"../fingerprints_df/{filename}.pkl")
        else:
            raise NameError(
                f'please provide a filename to load if not running calculation,\
                            possibly: {os.listdir("../fingerprints_df")}'
            )
            # print(os.listdir('../fingerprints_df'))
    return df
