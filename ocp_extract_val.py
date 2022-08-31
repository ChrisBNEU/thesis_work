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
dataval_id = LmdbDataset({"src": "/work/westgroup/opencatalyst/ocp/data/is2re/all/val_id/data.lmdb"})
print("Size of the dataset created:", len(dataval_id))
data_lengths.append(len(dataval_id))
data_lmdbs.append(dataval_id)

dataval_ood_ad = LmdbDataset({"src": "/work/westgroup/opencatalyst/ocp/data/is2re/all/val_ood_ads/data.lmdb"})
print("Size of the dataset created:", len(dataval_ood_ad))
data_lengths.append(len(dataval_ood_ad))
data_lmdbs.append(dataval_ood_ad)

dataval_ood_cat = LmdbDataset({"src": "/work/westgroup/opencatalyst/ocp/data/is2re/all/val_ood_cat/data.lmdb"})
print("Size of the dataset created:", len(dataval_ood_cat))
data_lengths.append(len(dataval_ood_cat))
data_lmdbs.append(dataval_ood_cat)

dataval_ood_both = LmdbDataset({"src": "/work/westgroup/opencatalyst/ocp/data/is2re/all/val_ood_both/data.lmdb"})
print("Size of the dataset created:", len(dataval_ood_both))
data_lengths.append(len(dataval_ood_both))
data_lmdbs.append(dataval_ood_both)

datats_id = LmdbDataset({"src": "/work/westgroup/opencatalyst/ocp/data/is2re/all/test_id/data.lmdb"})
print("Size of the dataset created:", len(datats_id))
data_lengths.append(len(datats_id))
data_lmdbs.append(dataval_id)

datats_ood_ad = LmdbDataset({"src": "/work/westgroup/opencatalyst/ocp/data/is2re/all/test_ood_ads/data.lmdb"})
print("Size of the dataset created:", len(datats_ood_ad))
data_lengths.append(len(datats_ood_ad))
data_lmdbs.append(datats_ood_ad)

datats_ood_cat = LmdbDataset({"src": "/work/westgroup/opencatalyst/ocp/data/is2re/all/test_ood_cat/data.lmdb"})
print("Size of the dataset created:", len(datats_ood_cat))
data_lengths.append(len(datats_ood_cat))
data_lmdbs.append(datats_ood_cat)

datats_ood_both = LmdbDataset({"src": "/work/westgroup/opencatalyst/ocp/data/is2re/all/test_ood_both/data.lmdb"})
print("Size of the dataset created:", len(datats_ood_both))
data_lengths.append(len(datats_ood_both))
data_lmdbs.append(datats_ood_both)

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

with open('metal_species_vt.pkl', 'wb') as f:
    pickle.dump(metal_dict, f)
    
