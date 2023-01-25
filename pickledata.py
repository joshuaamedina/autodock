import os
import pickle

ligands = {}

for filename in os.listdir('smallset'):
  f = open(f'smallset/{filename}', 'r')
  lines = f.read()
  ligands[filename] = lines

output_file = open('ligands_ZINC-in-trials.pkl', 'wb')
pickle.dump(ligands, output_file)
