import lmdb
import pickle
from ase import Atoms
from ase.visualize import view
import ase.data 
from ocpmodels.datasets import LmdbDataset

oc20_path = "/work/westgroup/opencatalyst/"

data_dict = pickle.load(open(f"{oc20_path}oc20_data_mapping.pkl","rb"))
print(len(data_dict))


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
    

dataset = LmdbDataset({"src": "/work/westgroup/opencatalyst/ocp/data/is2re/all/train/data.lmdb"})
print("Size of the dataset created:", len(dataset))


# parse data and pull the info for species identified
metal_dict = {}
for data in dataset:
    if int(data.sid) in list(metal_sids.keys()):
        print(f"{data.sid} match! {metal_sids[data.sid][:]}")
        data["ads_symbol"] = metal_sids[data.sid][0]
        data["miller_index"] = metal_sids[data.sid][1]
        data["metal"] = metal_sids[data.sid][2]
        metal_dict[data.sid] = data.to_dict()


with open('metal_species.pkl', 'wb') as f:
    pickle.dump(metal_dict, f)

