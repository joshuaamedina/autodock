from vina import Vina
from mpi4py import MPI
import subprocess
import pickle
import time
import os
import sys
from os.path import exists
from os.path import basename

# Set variables
receptor='1iep_receptor'
docking_type = 'vina'
ligand_library = ''
num_ligands = 10000
model = 'basic'

def sort():
    subprocess.run(["cat results* >> results_merged.txt"], shell=True)
    INPUTFILE = 'results_merged.txt'
    OUTPUTFILE = './output/processed_results.txt'

    result = []

    with open(INPUTFILE) as data:
        line = data.readline()
        while (line):
            filename = basename(line.split()[-1])
            v = data.readline().split()[0]
            result.append(f'{v} {filename}\n')
            line = data.readline()

    with open(OUTPUTFILE, 'w') as data:
        data.writelines(sorted(result, key=lambda x: float(x.split()[1])))

    subprocess.run(["rm results*; mv *map* *.gpf ./output/maps"], shell=True)

def run_docking(ligand, v, filename):
    v.set_ligand_from_string(ligand)
    v.dock()
    v.write_poses(f'./output/{docking_type}/output_{filename}', n_poses=1, overwrite=True)


def prep_maps(receptor):
    if docking_type == 'ad4':
        if exists('{receptor}.gpf'):
            subprocess.run([f"rm {receptor}.gpf"], shell=True)
        subprocess.run([f"python3 /scripts/write-gpf.py --box /configs/config.config /input/receptors/{receptor}.pdbqt"], shell=True)
        subprocess.run([f"/scripts/autogrid4 -p /input/receptors/{receptor}.gpf"], shell=True)

def conversion(receptor):
    if exists('./input/receptors/{receptor}.pdb'):
        subprocess.run([f'./scripts/prepare_receptor -r ./input/receptors/{receptor}.pdb -o ./input/receptors/{receptor}.pdbqt'], shell=True)

def pre_processing():
    prep_maps(receptor)
    conversion(receptor)

def partition(rank,size,N):
    n = N//size + ((N % size) > rank)
    s = rank * (N//size)
    if(N % size) > rank:
        s += rank
    else:
        s += N % size
    return s, s+n

def main():
    # Setup
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    ligands = None
    filenames = None
    # Pre-Processing (only rank 0)
    if rank == 0:
        pre_processing()
        ligandlist = pickle.load(open('./input/ligands_10.pkl', 'rb'))
        filenames = ligandlist[0]
        ligands = ligandlist[1]

    filenames = comm.bcast(filenames, root=0)
    ligands = comm.bcast(ligands, root=0)
    comm.Barrier()

    # Initialize Vina or AD4 configurations
    docking_type = 'vina'
#    ligands = pickle.load(open('./input/ligands_10.pkl', 'rb'))

    if docking_type == 'vina':
        v = Vina(sf_name='vina', cpu=4, verbosity=0)
        v.set_receptor(f'./input/receptors/{receptor}.pdbqt')
        v.compute_vina_maps(center=[15.190, 53.903, 16.917], box_size=[20, 20, 20])
    elif docking_type == 'ad4':
        v = Vina(sf_name='ad4', cpu=0, verbosity=0)
        v.load_maps(map_prefix_filename='1iep_receptor')

    start, end = partition(rank,size,len(ligands))
    # Run docking
    for x in range(start,end):
        filename = filenames[x]
        ligand =  ligands[x]
        run_docking(ligand, v, filename)
        subprocess.run([f"grep -i -m 1 'REMARK VINA RESULT:' ./output/{docking_type}/output_{filename} \
                            | awk '{{print $4}}' >> results_{rank}.txt; echo {filename} \
                            >> results_{rank}.txt"], shell=True)

    # Post-Processing
    comm.Barrier()
    if rank == 0:
        sort()
        end_time = time.time()
    benchmark(start_time, end_time, num_ligands, model, docking_type)

main()
