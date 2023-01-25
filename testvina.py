import subprocess
from mpi4py import MPI
comm = MPI.COMM_WORLD
from vina import Vina
import pickle

with open('ligands_ZINC-in-trials.pkl', 'rb') as f:
    mynewdict = pickle.load(f)

rank = comm.rank
size = comm.size

for num, ligand in enumerate(mynewdict):
    if rank == (num % size):
        v = Vina(sf_name='vina')
        v.set_receptor('1iep_receptor.pdbqt')
        v.compute_vina_maps(center=[15.190, 53.903, 16.917], box_size=[20, 20, 20])
        v.set_ligand_from_string(mynewdict[ligand])
        v.dock()
        v.write_poses(f'output/{num}_{rank}_vina_out_{ligand}', n_poses=1, overwrite=True)
        command = f"grep -i -m 1 'REMARK VINA RESULT:' output/{num}_{rank}_vina_out_{ligand} | awk '{{print $4}}' >> temp_{num}.txt; echo {ligand} >> temp_{num}.txt"
        subprocess.run(command,shell=True)
