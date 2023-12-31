from dataclasses import dataclass
from typing import Optional

# absolute paths are provided here. I know it's not ideal
### but since some scripts are in diff directories, this
#### is the easiest. B/c of gromacs I will 100% not ever
### be changing the absolute path so I think we are set


@dataclass()
class SimulationMetadata:
    name: str
    path: str
    n_replicas: int
    color: str
    ref_path: Optional[str] = None  # this is a gro file for RMSD calculation
    ligname: Optional[str] = None  # this is ligname in mda Uni file
    protein: Optional[str] = "PfHT1"  # this is for some res# dicts later on


PfHT_bf1_tosh = SimulationMetadata(
    name="PfHT_TOSH",
    path="/data/PfHT_MMV/analysis/mnt_marta",
    n_replicas=5,
    color="blue",
    protein="PfHT1",
    ligname="LIG",
)

PfHT_MMV12 = SimulationMetadata(
    name="PfHT_MMV12",
    path="/data/PfHT_MMV/EM_MMV/MMV12",
    n_replicas=3,
    color="#539C44",
    protein="PfHT1",
    ref_path="/data/PfHT_MMV/EM_MMV/MMV12/replica_1/production/PfHT_MMV12.em.1.start.protonly.gro",
    ligname="MMV",
)
W412A_MMV12 = SimulationMetadata(
    name="W412A_MMV12",
    path="/data/PfHT_MMV/EM_MMV/W412A_MMV12",
    n_replicas=3,
    color="#F2B705",
    protein="PfHT1",
    ref_path="",
)
PfHT_MMV8 = SimulationMetadata(
    name="PfHT_MMV8",
    path="/data/PfHT_MMV/EM_MMV/MMV8",
    n_replicas=3,
    color="#84B3BB",
    protein="PfHT1",
    ref_path="",
)
GLUT1_MMV12 = SimulationMetadata(
    name="GLUT1_MMV12",
    path="/data/PfHT_MMV/EM_MMV/GLUT1_MMV12",
    n_replicas=3,
    color="#BCBCBC",
    protein="GLUT1",
    ligname="MMV",
    ref_path="/data/PfHT_MMV/EM_MMV/GLUT1_MMV12/replica_1/production/GLUT1_MMV12.em.1.start.protonly.gro",
)
PfHT_apo = SimulationMetadata(
    name="PfHT_apo",
    path="/data/PfHT_MMV/EM_MMV/PfHT_apo",
    n_replicas=3,
    color="blue",
    protein="PfHT1",
    ref_path="",
)
PfHT_3361_crystal = SimulationMetadata(
    name="PfHT_3361.crystal",
    path="/data/PfHT_MMV/C3361_sims/crystal_str_6m2l",
    n_replicas=3,
    color="grey",
    protein="PfHT1",
    ref_path="/data/PfHT_MMV/C3361_sims/crystal_str_6m2l/replica_1/production/PfHT_3361.crystal.1.start.protonly.gro",
    ligname="F00",
)
PfHT_3361_em = SimulationMetadata(
    name="PfHT_3361.em",
    path="/data/PfHT_MMV/C3361_sims/EM_3361",
    n_replicas=3,
    color="#ED91B8",
    protein="PfHT1",
    ref_path="/data/PfHT_MMV/C3361_sims/EM_3361/replica_1/production/PfHT_3361.em.1.start.protonly.gro",
    ligname="F00",
)
PfHT_MMV_crystal = SimulationMetadata(
    name="PfHT_MMV_crystal",
    path="/data/PfHT_MMV/crystal_structure_sims/MMV_sims",
    n_replicas=4,
    color="#F28705",
    ligname="MMV",
    protein="PfHT1",
    ref_path="/data/PfHT_MMV/crystal_structure_sims/MMV_sims/replica_1/production/PfHT_MMV.1.start.protonly.gro",
)

PfHT_MMV8_crystal = SimulationMetadata(
    name="PfHT_MMV8_crystal",
    path="/data/PfHT_MMV/crystal_structure_sims/MMV8_sims",
    n_replicas=3,
    color="#84B3BB",
    ligname="MMV",
    protein="PfHT1",
    ref_path="/data/PfHT_MMV/crystal_structure_sims/MMV_sims/replica_1/production/PfHT_MMV.1.start.protonly.gro",
)

W412A_MMV_crystal = SimulationMetadata(
    name="W412A_MMV_crystal",
    path="/data/PfHT_MMV/crystal_structure_sims/W412A_MMV",
    n_replicas=3,
    color="#F2B705",
    ligname="MMV",
    protein="PfHT1",
    ref_path="/data/PfHT_MMV/crystal_structure_sims/MMV_sims/replica_1/production/PfHT_MMV.1.start.protonly.gro",
)

GLUT3_MMV_crystal = SimulationMetadata(
    name="GLUT3_MMV_crystal",
    path="/data/PfHT_MMV/crystal_structure_sims/GLUT3_MMV",
    n_replicas=3,
    color="#BCBCBC",
    ligname="MMV",
    protein="GLUT3",
    ref_path="/data/PfHT_MMV/crystal_structure_sims/GLUT3_MMV/replica_1/production/GLUT3_MMV.1.start.protonly.gro",
)
