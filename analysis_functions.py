from typing import Optional
import glob
import MDAnalysis as mda
import pandas as pd
import simulation_names as sims

def load_unis(
        sim_list: list=[
                        sims.PfHT_MMV12,
                        sims.W412A_MMV12,
                        sims.PfHT_MMV8,
                        sims.PfHT_apo,
                        sims.GLUT1_MMV12,
                        sims.PfHT_3361_crystal,
                        sims.PfHT_3361_em                        
                        ] ,  
        protonly_or_wholesys: str='protonly'     
) -> tuple[dict[str, list[mda.Universe]], dict[str, list[int]]]:

    """
    will load mda unis of all em analyses so far
    Not so space intensive since I don't save this to memory with mda
    will return a dict with the key being the uni name and the values
    being a list of the replicas starting with 1


    Accepts a list of lists as [sim_name, n_replicas, path_to_sim],[]...

    for now, protonly and TPR
    which means that numbering is incorrect, you'll have to add 21

    Note that each simulation is a different length. This function will find the longest sim
    of 0_xxx ns long assuming that naming convention is correct 

    A tuple of dicts is returned, one with {sim_name:mda_uni}, and one with {sim_name, sim_length}
    """

    all_uni = {}
    all_lens = {}

    for sim_name, path_name, n_replicas in sim_list:
        uni_list = []
        sim_length_list = []

        for replica in range(1, n_replicas+1):

            ####### get filenames

            ## first tpr - any protonly is OK
            if protonly_or_wholesys == 'protonly':
                tpr_file = glob.glob(f'{path_name}/replica_{replica}/production/*protonly.tpr', recursive=True)
            else:
                tpr_file = glob.glob(f'{path_name}/replica_{replica}/production/*0_200ns.tpr', recursive=True)
            if len(tpr_file) == 0:
                raise FileNotFoundError(f'protonly tpr file missing for {sim_name} rep {replica}')
            else:
                tpr_file = tpr_file[0]

            ## next find longest traj file
            traj_files = glob.glob(f'{path_name}/replica_{replica}/production/*.0_*skip250*{protonly_or_wholesys}.xtc', recursive=True)
            max_len = 0
            longest_trajectory = None
            for traj in traj_files:
                traj_full_name = traj.split('/')[-1] ## full sim name will always be the final name in the path
                traj_len = int(traj_full_name.split('0_')[1].split('ns')[0])  ## naming convention will always have "0_XXXns", so we can split between these 
                
                ## find longest sim length, replace previous until you found maxmimum #. Keep this trajectory as longest_trajectory
                if traj_len > max_len:
                    max_len = traj_len
                    longest_trajectory = traj
            if longest_trajectory == None:
                raise FileNotFoundError(f'Simulation {sim_name} rep {replica} here is not named correctly. \
                                        Are you sure you have "protonly" if protonly= {protonly_or_wholesys}?')
            



            ##### make uni of filenames. Add all reps to this universe list

            u = mda.Universe(tpr_file, longest_trajectory)
            uni_list.append(u)
            sim_length_list.append(len(u.trajectory))


        ### add to dict
        all_uni[sim_name] = uni_list
        all_lens[sim_name] = sim_length_list

    return all_uni, all_lens





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




'''

## older functions that I am not ready to totally remove





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

    all_uni[sims.PfHT_MMV12[0]] = uni_wt_mmv12
    all_len[sims.PfHT_MMV12[0]] = len_wt_mmv12
    all_uni[sims.W412A_MMV12[0]] = uni_w412a_mmv12
    all_len[sims.W412A_MMV12[0]] = len_w412a_mmv12
    all_uni[sims.PfHT_MMV8[0]] = uni_wt_mmv8
    all_len[sims.PfHT_MMV8[0]] = len_wt_mmv8
    all_uni[sims.GLUT1_MMV12[0]] = uni_glut1_mmv12
    all_len[sims.GLUT1_MMV12[0]] = len_glut1_mmv12
    all_uni[sims.PfHT_apo[0]] = uni_wt_apo
    all_len[sims.PfHT_apo[0]] = len_wt_apo

    ### W412A replicas next
    # uni_l = []
    return all_uni, all_len

'''