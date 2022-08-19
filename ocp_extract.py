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

for metal, mpid in mat_proj_id:
    
    
count = 0
pt_sids = {}
for key, value_dict in data_dict.items():
    if value_dict['bulk_mpid'] == 'mp-126':
        count +=1
        print(key, value_dict['ads_symbols'], value_dict['ads_id'], value_dict['miller_index'])
        # make dictionary of the symbols (e.g. *CH3) and miller indices with the "randomXXXXXXX" ID tag
        pt_sids.update({int(key.replace('random','')):(value_dict['ads_symbols'], value_dict['miller_index'])})
print(f"{count} Pt species found")


count = 0
cu_sids = {}
for key, value_dict in data_dict.items():
    if value_dict['bulk_mpid'] == 'mp-30':
        count +=1
        print(key, value_dict['ads_symbols'], value_dict['ads_id'], value_dict['miller_index'])
        # make dictionary of the symbols (e.g. *CH3) with the "randomXXXXXXX" ID tag
        cu_sids.update({int(key.replace('random','')):(value_dict['ads_symbols'], value_dict['miller_index'])})
print(f"{count} Cu species found")

ni_sids = {}
for key, value_dict in data_dict.items():
    if value_dict['bulk_mpid'] == 'mp-23':
        count +=1
        print(key, value_dict['ads_symbols'], value_dict['ads_id'], value_dict['miller_index'])
        # make dictionary of the symbols (e.g. *CH3) with the "randomXXXXXXX" ID tag
        ni_sids.update({int(key.replace('random','')):(value_dict['ads_symbols'], value_dict['miller_index'])})
print(f"{count} Ni species found")



dataset = LmdbDataset({"src": "/work/westgroup/opencatalyst/ocp/data/is2re/all/train/data.lmdb"})
print("Size of the dataset created:", len(dataset))


# parse data and pull the info for Cu, Ni, Pt species identified
pt_dict = {}
cu_dict = {}
ni_dict = {}
for data in dataset:
    if int(data.sid) in list(pt_sids.keys()):
        data["miller_index"] = pt_sids[data.sid][1]
        pt_dict[data.sid] = data.to_dict()
    elif int(data.sid) in list(cu_sids.keys()):
        data["miller_index"] = cu_sids[data.sid][1]
        cu_dict[data.sid] = data.to_dict() 
    elif int(data.sid) in list(ni_sids.keys()):
        data["miller_index"] = ni_sids[data.sid][1]
        ni_dict[data.sid] = data.to_dict()


with open('Pt111_species.pkl', 'wb') as f:
    pickle.dump(pt_dict, f)
    
with open('Cu111_species.pkl', 'wb') as f:
    pickle.dump(cu_dict, f)

with open('Ni111_species.pkl', 'wb') as f:
    pickle.dump(ni_dict, f)
