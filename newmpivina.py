import subprocess
from mpi4py import MPI
comm = MPI.COMM_WORLD

rank = comm.rank
size = comm.size

fp = open("100zinc.txt", 'r')

for line_no, line in enumerate(fp):
    if rank == (line_no % size):
        print(f"RANK({rank}): Line {line_no} = {line}")
        subprocess.run(["vina", "--receptor", f'1iep_receptor.pdbqt',"--ligand",line.strip(),"--config", f'1iep_receptor_vina_box.txt',"--exhaustiveness=32",
                       "--out",#f'output/{(line[:-7])[9:]}_vina_out_{line_no}.pdbqt'])
        command = f"grep -i -m 1 'REMARK VINA RESULT:' output/{line[:-7][9:]}_vina_out_{line_no}.pdbqt | awk '{{print $4}}' >> temp_{line_no}.txt; echo {line.strip()} >> temp_{line_no}.txt"
        subprocess.run(command,shell=True)

fp.close()

