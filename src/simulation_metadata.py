from dataclasses import dataclass


@dataclass()
class SimulationMetadata:
    name: str
    path: str
    n_replicas: int
    color: str


PfHT_MMV12 = SimulationMetadata(
    name="PfHT_MMV12",
    path="../../EM_MMV/MMV12",
    n_replicas=3,
    color="#539C44",
)

W412A_MMV12 = SimulationMetadata(
    name="W412A_MMV12",
    path="../../EM_MMV/W412A_MMV12",
    n_replicas=3,
    color="#37456d",
)
PfHT_MMV8 = SimulationMetadata(
    name="PfHT_MMV8",
    path="../../EM_MMV/MMV8",
    n_replicas=3,
    color="#E27439",
)
GLUT1_MMV12 = SimulationMetadata(
    name="GLUT1_MMV12",
    path="../../EM_MMV/GLUT1_MMV12",
    n_replicas=3,
    color="grey",
)
PfHT_apo = SimulationMetadata(
    name="PfHT_apo",
    path="../../EM_MMV/PfHT_apo",
    n_replicas=3,
    color="blue",
)
PfHT_3361_crystal = SimulationMetadata(
    name="PfHT_3361.crystal",
    path="../../C3361_sims/crystal_str_6m2l",
    n_replicas=3,
    color="green",
)
PfHT_3361_em = SimulationMetadata(
    name="PfHT_3361.em",
    path="../../C3361_sims/EM_3361",
    n_replicas=3,
    color="red",
)
