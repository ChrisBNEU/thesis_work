import lmdb
import pickle
from ase import Atoms
from ase.visualize import view
import ase.data 
from ocpmodels.datasets import LmdbDataset
import time

oc20_path = "/work/westgroup/opencatalyst/"
data_dict = pickle.load(open(f"{oc20_path}oc20_data_mapping.pkl","rb"))

mat_proj_id = {
    "Pt":'mp-126', # fcc
    "Cu":'mp-30',  # fcc
    "Ni":'mp-23', # fcc
    "Ru":'mp-33', # hcp
    "Rh":'mp-74', # fcc
    "Ir":'mp-101', # fcc
    "Au":'mp-81', # fcc
    "Pd":'mp-2', # fcc
    "Ag":'mp-124', # fcc
    "Co":'mp-102', # hcp
}
metal_sids = {}
for metal, mpid in mat_proj_id.items():
    count = 0
    for key, value_dict in data_dict.items():
        if value_dict['bulk_mpid'] == mpid:
            count +=1
            print(key, value_dict['ads_symbols'], value_dict['ads_id'], value_dict['miller_index'], )
            # make dictionary of the symbols (e.g. *CH3) miller indices and metal with the "randomXXXXXXX" ID tag
            metal_sids.update({int(key.replace('random','')):(value_dict['ads_symbols'], value_dict['miller_index'], metal)})
    print(f"{count} metal species found for metal {metal}")


# load all of the val and test lmdb datasets
data_lengths = []
data_lmdbs = []
data_10ktr = LmdbDataset({"src": "/work/westgroup/opencatalyst/ocp/data/is2re/100k/train/data.lmdb"})
print("Size of the dataset created:", len(data_10ktr))
data_lengths.append(len(data_10ktr))
data_lmdbs.append(data_10ktr)

dataval_100ktr = LmdbDataset({"src": "/work/westgroup/opencatalyst/ocp/data/is2re/10k/train/data.lmdb"})
print("Size of the dataset created:", len(dataval_100ktr))
data_lengths.append(len(dataval_100ktr))
data_lmdbs.append(dataval_100ktr)

max_size = max(data_lengths)
# helper function for parsing data
def check_data(sid_dict, output_dict, lmdb, num):
    if num >= len(lmdb):
        return
    if int(lmdb[num].sid) in list(sid_dict.keys()):
        print("in")
        lmdb[num]["ads_symbol"] = sid_dict[lmdb[num].sid][0]
        lmdb[num]["miller_index"] = sid_dict[lmdb[num].sid][1]
        lmdb[num]["metal"] = sid_dict[lmdb[num].sid][2]
        output_dict[lmdb[num].sid] = lmdb[num].to_dict()

st = time.time()
# truncate size to 2 
# parse data and pull the info for species identified
metal_dict = {}
for db_num, lmdb_dat in enumerate(data_lmdbs):
    for num in range(0, len(lmdb_dat)):
        # check in domain
        check_data(metal_sids, metal_dict, lmdb_dat, num)

        if num%100==0:
            print(num)
            tim = time.time()

            elapsed = tim - st
            print(f"{num} in db {db_num} took {elapsed} seconds")
        
et = time.time()

tot_elapsed = et - st

print(f"all took {tot_elapsed} seconds") 

with open('metal_species_10_100k.pkl', 'wb') as f:
    pickle.dump(metal_dict, f)