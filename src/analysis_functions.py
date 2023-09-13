from typing import Optional
import os
import glob

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
import prolif as plf

import src.simulation_metadata as sims


def load_unis_marta(
    sim_list: list[sims.SimulationMetadata],
    filename_xtc: str,
    filename_tpr: Optional[str] = None,
) -> tuple[dict[str, list[mda.Universe]], dict[str, list[int]]]:
    """
    same as below, load mda universes of  each replica
    return tuple of two dicts
    dict 1 = all_unis, key=sim.name, value = list as long as n replicas of unis
    of mda unis
    dict2 = all_lens, key=sim.name, value=list as long as n replicas of unis
    of lengths of unis

    This assumes that the filenames are always the same for each condition
    besides the replica #
    """
    all_uni = {}
    all_lens = {}

    for sim in sim_list:
        uni_list = []
        sim_length_list = []

        for replica in range(1, sim.n_replicas + 1):
            tpr = f"{sim.path}/rep{replica}.prot.tpr"
            xtc = f"{sim.path}/{filename_xtc}.{replica}.skip250.xtc"

            u = mda.Universe(tpr, xtc)
            uni_list.append(u)
            sim_length_list.append(len(u.trajectory))

        all_uni[sim.name] = uni_list
        all_lens[sim.name] = sim_length_list

    return all_uni, all_lens


def load_unis(
    sim_list: list[sims.SimulationMetadata],
    protonly_or_wholesys: str = "protonly",
) -> tuple[dict[str, list[mda.Universe]], dict[str, list[int]]]:
    """
    will load mda unis of all em analyses so far
    Not so space intensive since I don't save this to memory with mda
    will return a dict with the key being the uni name and the values
    being a list of the replicas starting with 1


    Accepts a list of lists as [sim_name, n_replicas, path_to_sim],[]...

    for now, protonly and TPR
    which means that numbering is incorrect, you'll have to add 21

    Note that each simulation is a different length.
    This function will find the longest sim
    of 0_xxx ns long assuming that naming convention is correct

    A tuple of dicts is returned, one with {sim_name:mda_uni},
      and one with {sim_name, sim_length}
    """

    all_uni = {}
    all_lens = {}

    for sim in sim_list:
        uni_list = []
        sim_length_list = []

        for replica in range(1, sim.n_replicas + 1):
            ## first tpr - any protonly is OK
            if protonly_or_wholesys == "protonly":
                tpr_file = glob.glob(
                    f"{sim.path}/replica_{replica}/production/*protonly.tpr",
                    recursive=True,
                )
            else:
                tpr_file = glob.glob(
                    f"{sim.path}/replica_{replica}/production/*0_200ns.tpr",
                    recursive=True,
                )

            if len(tpr_file) == 0:
                raise FileNotFoundError(
                    f"protonly tpr file missing\
                                         for {sim.name} rep {replica}"
                )
            else:
                tpr_file = tpr_file[0]

            ## next find longest traj file
            traj_files = glob.glob(
                f"{sim.path}/replica_{replica}/production/*.0_*skip250*{protonly_or_wholesys}.xtc",
                recursive=True,
            )
            max_len = 0
            longest_trajectory = None
            for traj in traj_files:
                ## full sim name will always be the final name in the path
                ## naming convention will always have "0_XXXns",
                #### so we can split between these
                traj_full_name = traj.split("/")[-1]
                traj_len = int(traj_full_name.split("0_")[1].split("ns")[0])

                ## find longest sim length, replace previous until you found maxmimum #
                # #Keep this trajectory as longest_trajectory
                if traj_len > max_len:
                    max_len = traj_len
                    longest_trajectory = traj
            if longest_trajectory is None:
                raise FileNotFoundError(
                    f'Simulation {sim.name} rep {replica} here is not named correctly. \
                Are you sure you have "protonly" if protonly= {protonly_or_wholesys}?'
                )

            ##### make uni of filenames. Add all reps to this universe list

            u = mda.Universe(tpr_file, longest_trajectory)
            uni_list.append(u)
            sim_length_list.append(len(u.trajectory))

        ### add to dict
        all_uni[sim.name] = uni_list
        all_lens[sim.name] = sim_length_list

    return all_uni, all_lens


#######################################################
##### PROLIF FINGERPRINTING FUNCTIONS ##################
#######################################################


def get_fp_dataframe(
    u: mda.Universe,
    skip: int,
    ligname: str,
    run_calc: bool = True,
    filename: Optional[str] = None,
    filename_path: Optional[str] = "/data/PfHT_MMV/analysis/fingerprints_df",
) -> pd.DataFrame:
    """
    Here I just run prolif fingerprints on a uni
    One can also save the filename, which will save by default in '../fingerprints_df
    If you don't want to perform the calculation, you can set this to "False" and write
    the filename you wish to load
    """

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
            pd.to_pickle(df, f"{filename_path}/{filename}.pkl")
        else:
            print("no filename given, returning dataframe, not saving file...")
    else:
        if filename is not None:
            df = pd.read_pickle(f"{filename_path}/{filename}.pkl")
        else:
            raise NameError(
                f'please provide a filename to load if not running calculation,\
                            possibly:\
                     {os.listdir(f"{filename_path}")}'
            )
            # print(os.listdir('../fingerprints_df'))
    return df


def pct_intxn_per_residue_wide(
    intxn_name: str,
    all_fp_dfs: dict[str, list[pd.DataFrame]],
    sims: list[sims.SimulationMetadata],
    mean_cutoff: float = 0.1,
) -> pd.DataFrame:
    """
    Calculates the percent interaction time for a specified interaction type
    Returns a multilevel dataframe with residue interacting as the index,
      and replica #(python indexed) as the column name, for each condition multi index
      this is in wide format, but we need to be able to keep Nans in place for ensuring
      that zeros are counted

      intxn_name should be HBAcceptor, HBDonor, Hydrophobic, PiStacking generally

      Different use cases:
      Comparing all PfHT simulations:
            Intention is therefore to drop GLUT1 from this dataframe
            and then melt along residues to get a long form dataframe
    """
    all_mean_interactions = {}
    for sim in sims:
        replica_names = [f"replica {i}" for i in range(1, sim.n_replicas + 1)]

        mean_intxn_all_reps = []
        for rep in all_fp_dfs[sim.name]:
            intxn_group_over_time = rep.xs(intxn_name, level="interaction", axis=1)
            mean_intxn = intxn_group_over_time.mean()
            mean_intxn = mean_intxn.loc[mean_intxn > mean_cutoff]
            mean_intxn_all_reps.append(mean_intxn)
        mean_intxn_all_reps = pd.concat(mean_intxn_all_reps, axis=1, keys=replica_names)
        all_mean_interactions[sim.name] = mean_intxn_all_reps

    all_mean_interactions = pd.concat(all_mean_interactions, axis=1)

    return all_mean_interactions


def pct_n_intxn_per_frame_wide(
    intxn_name: str,
    sims: list[sims.SimulationMetadata],
    all_fp_dfs: dict[str, list[pd.DataFrame]],
    percentage: bool = True,
    value_counts: bool = True,
) -> pd.DataFrame:
    """
    In all dataframes, find the sum across columns for a specific interaction type
    this will give you the n interactions per frame (since boolean df)

    If value_counts is set to True:
        will find how many times each column has each n interactions
        index will become n_interactions

    If percentage is set to True:
        find % of total simulation time this number occurs
        by dividing by length of df



    This will return a dataframe in a wide format, which we can then postprocess later
    """
    temp_all = {}
    for sim in sims:
        fp_df_l = all_fp_dfs[sim.name]
        temp_rep = []

        for rep_n, rep_df in enumerate(fp_df_l):
            n_counts = rep_df.xs(intxn_name, level="interaction", axis=1).sum(axis=1)
            if value_counts:
                n_counts = n_counts.value_counts()
                index_name = "n_interactions"
            else:
                index_name = "frame"
            if percentage:
                n_counts = n_counts / len(rep_df)
            n_counts.name = f"replica {rep_n+1}"
            temp_rep.append(n_counts)

        temp_rep = pd.concat(temp_rep, axis=1)
        temp_all[sim.name] = temp_rep

    temp_all = pd.concat(temp_all, axis=1)
    temp_all.index.name = index_name
    return temp_all


def res_intxn_over_time_plf(
    df_l: list[pd.DataFrame],
    df_names: list[str],
    resname: str,
    subtract_num: int = 21,  # PfHT numbering, may need another # for gluts or whatever
) -> tuple[pd.DataFrame, dict]:
    """
    Takes a list of dataframes
    of prolif fingerprint as input

    PfHT_num assumes that you have to subtract 21 from the res number
    provided

    melts dataframe on "Frame", and replaces "False" with np.nan
    so that plotting is easier
    Returns a longform dataframe with colums:
    Frame, interaction, value, and interaction_loc


        ############### FOR DISPLACING ON SEABORN ################
     setting y='original_df' etc will mean that all values overlap
     therefore, we give the y value as a number (based on interaction type, so that
     each y is separated)
     then, we also give each dataframe that went into this function a number,
       which we will
     use to "jitter" the plot a bit, so that you don't have all data points overlapping

    interaction_loc is a handy way to toggle sns scatterplot so that you can show
    distributions side by side
    """

    resname = resname[:3] + str(int(resname[3:]) - subtract_num)

    ### make dataframe longform, replace all False with nan
    intxn_possibilities = []
    new_df_l = []
    for df in df_l:
        df = df[resname]
        intxn_possibilities.extend(df.columns)
        new_df_l.append(df)

    df = pd.concat(new_df_l, axis=1, keys=df_names)
    df = df.reset_index().melt(
        id_vars=["Frame"], var_name=["original_df", "interaction"]
    )
    # df['value']= df['value'].replace(False, np.nan)

    #### SEABORN DISPLACEMENT THINGS ###
    ## for displacing y on seaborn plot, we need to have the locations for the
    #### interaction type
    intxn_name_dict = {}
    intxn_possibilities = set(intxn_possibilities)
    for n, interaction in enumerate(intxn_possibilities):
        intxn_name_dict[interaction] = n + 1
    df["interaction_loc"] = df["interaction"].replace(intxn_name_dict)

    ## also for displacing, need to know how much to displace by,
    # # which we should dictate based
    #### on dataframe name (so each displacement is the same)
    df_name_dict = {}
    for n, df_name in enumerate(df_names):
        df_name_dict[df_name] = n + 1
    df["original_df_name_loc"] = df["original_df"].replace(df_name_dict)

    return df, intxn_name_dict


def process_wide_df(
    df: pd.DataFrame,
    index_col_name: str,
    condition_to_remove: Optional[str] = None,
    index_name: str = "index",
    top_level_melt_name: str = "protein",
    lower_level_melt_name: str = "replica",
    add_21: bool = False,
) -> pd.DataFrame:
    """
    Wide multi index dataframe as input with the following format:
    eventual label as index (ie residue name)
    condition name as top label
    replica name as under label
    This isn't strictly necessary but it is designed this way because
    other functions make dataframes like this
    ie:
    index:              PfHT_MMV12     ... GLUT1_MMV12
            replica 1:   replica 2:         replica 3:
    GLN104     value1       value3            value500
    ASP205     value2       value4

    is a typical format for input
    and output will be like:
           value:     replica:   protein:    residue:
    0      value1     replica1    PfHT_MMV12   GLN104
    1      value3     replica2    PfHT_MMV12   GLN104
    ...
    20     value2     replica3   PfHT_MMV12   ASP205


    will melt array and rename the index to the variable 'index_col_name'
    (this might be residue, n_interactions...). If index is from prolif
    directly, it's probably named "protein"


    """
    ## sometimes there might be values specific only to the
    ### condition you want to remove, so drop these from the dataframe
    if condition_to_remove is not None:
        del df[condition_to_remove]
        df = df.dropna(how="all", axis=0)

    ## for the melting, if you have gathered replicas together, and then
    #### gathered conditions on top of that (as is often the structure of other
    ##### functions of this script), the top_level_melt_name is usually 'protein'
    ###### and lower_level_melt_name is 'residue'
    df = (
        df.fillna(0.0)
        .reset_index()
        .rename({index_name: index_col_name}, axis=1)
        .melt(
            id_vars=[index_col_name],
            var_name=[top_level_melt_name, lower_level_melt_name],
        )
    )

    ### finally, add 21 to res numbers
    ## TODO - conditional if PfHT
    if add_21:
        df["resnr"] = df["residue"].str[3:].astype(int) + 21
        df["residue"] = df["residue"].str[:3] + df["resnr"].astype(str)
        del df["resnr"]

    return df


###################################################
############### MD analysis things #################
####################################################


def calc_rmsd(
    u: mda.Universe,
    ref_name,  ## TODO I don't think this is actually a str
    select="backbone",
    groupselections: list = ["resname MMV"],
    run_calc: bool = True,
    filename: Optional[str] = None,
) -> np.array:
    """
    This will calculate RMSD using MDA
    It's barely different from the mda tool itself, but
    here I give the option to save the np array so you can easily call it later
    The calculation otherwise can take some time


    Will save an np array for each rmsd individually
    """
    ref = mda.Universe(ref_name)

    R = RMSD(
        u,
        ref,
        select=select,
        groupselections=groupselections,
    )
    R.run()

    return R.results.rmsd
