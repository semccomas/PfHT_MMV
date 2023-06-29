from dataclasses import dataclass


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


PfHT_MMV12 = SimulationMetadata(
    name="PfHT_MMV12",
    path="/data/PfHT_MMV/EM_MMV/MMV12",
    n_replicas=3,
    color="#539C44",
)
W412A_MMV12 = SimulationMetadata(
    name="W412A_MMV12",
    path="/data/PfHT_MMV/EM_MMV/W412A_MMV12",
    n_replicas=3,
    color="#37456d",
)
PfHT_MMV8 = SimulationMetadata(
    name="PfHT_MMV8",
    path="/data/PfHT_MMV/EM_MMV/MMV8",
    n_replicas=3,
    color="#E27439",
)
GLUT1_MMV12 = SimulationMetadata(
    name="GLUT1_MMV12",
    path="/data/PfHT_MMV/EM_MMV/GLUT1_MMV12",
    n_replicas=3,
    color="grey",
)
PfHT_apo = SimulationMetadata(
    name="PfHT_apo",
    path="/data/PfHT_MMV/EM_MMV/PfHT_apo",
    n_replicas=3,
    color="blue",
)
PfHT_3361_crystal = SimulationMetadata(
    name="PfHT_3361.crystal",
    path="/data/PfHT_MMV/C3361_sims/crystal_str_6m2l",
    n_replicas=3,
    color="green",
)
PfHT_3361_em = SimulationMetadata(
    name="PfHT_3361.em",
    path="/data/PfHT_MMV/C3361_sims/EM_3361",
    n_replicas=3,
    color="red",
)
PfHT_MMV_crystal = SimulationMetadata(
    name="PfHT_MMV_crystal",
    path="/data/PfHT_MMV/crystal_structure_sims/MMV_sims",
    n_replicas=3,
    color="green",
)
