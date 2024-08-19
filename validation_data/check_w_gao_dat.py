#!/usr/bin/env python
# coding: utf-8

# # Check with Gao Data
# verify the rmg calculated thermo by extrapolating the H to 0K against the data from Gao (2016)

# In[1]:


import sys
import os
import re
import numpy
import pandas as pd
import sys
# sys.path.append('/work/westgroup/ChrisB/_04_thesis/rmg_thesis/RMG-Py/')
print(sys.path)
from rmgpy.species import Species

import matplotlib.pyplot as plt
import numpy as np
# Species("OH")


# In[2]:


# load the thermodynamics database
from rmgpy.data.thermo import ThermoDatabase
thermoDatabase = ThermoDatabase()
libraries = ['surfaceThermoPt111']
if sys.platform == "darwin":
    thermoDatabase.load(path="/Users/blais.ch/Documents/_01_code/RMG_env_1/RMG-database/input/thermo", libraries = libraries)
else: 
    thermoDatabase.load(path='/work/westgroup/ChrisB/_04_thesis/rmg_thesis/RMG-database/input/thermo', libraries = libraries)
    
thermoDatabase.load_surface()
pt_library = thermoDatabase.libraries['surfaceThermoPt111']


# In[20]:


class SpeciesDat():
    def __init__(self, adj_list): 
        
        self.rmg_spec = Species().from_adjacency_list(adj_list)
     
        self.new_thermo = {}
        self.Hf_rmg = {}
        
        # gao data structures
        self.Hf_lit = {}
        self.Hf_gao = {}
        self.metal_list = []
        
        # delta Hfs
        self.dH_rmg = {}
        self.dH_lit = {}
        self.dH_gao = {}

    def scale_rmg_thermo(self, metal_from='Pt', metal_to='Rh', 
                         facet_from='111', facet_to='211'):
        
        # need to redefine species every time so we don't stack Hf
        self.rmg_spec.thermo, _, _ = thermoDatabase.get_thermo_data_from_library(self.rmg_spec,pt_library)
        self.rmg_spec.generate_resonance_structures()
        self.spec_thermo = self.rmg_spec.thermo.to_thermo_data()
        self.rmg_spec.thermo = self.spec_thermo
        self.Hf_old = (self.rmg_spec.thermo.H298.value_si)/9.68e4
        new_metal_str = f"{metal_to}({facet_to})"
        self.new_thermo[new_metal_str] = thermoDatabase.correct_binding_energies_extended(self.spec_thermo, 
                                                                           self.rmg_spec, 
                                                                           metal_to_scale_from=metal_from,
                                                                           metal_to_scale_to=metal_to, 
                                                                           facet_to_scale_from=facet_from, 
                                                                           facet_to_scale_to=facet_to)
        self.Hf_rmg[new_metal_str] = self.new_thermo[new_metal_str].H298.value_si/9.68e4
        self.Hf_rmg[new_metal_str]
        
    def scale_all_thermo(self): 
        """
        scale all of the thermo corresponding to each metal in gao data
        """
        if len(self.metal_list) > 0:
            for metal, facet in self.metal_list:
                self.scale_rmg_thermo(metal_from='Pt', metal_to=metal, 
                         facet_from='111', facet_to=facet)
        else: 
            raise Exception("please load gao data first")
                
                
    
    def set_gao_hf(self, hf): 
        self.gao_hf_df = hf
        self.gao_hf_df.columns = ["DFTcal", "predicted"]
        for row in self.gao_hf_df.iterrows(): 
            self.Hf_lit[row[0]] = row[1]['DFTcal']
            self.Hf_gao[row[0]] = row[1]['predicted']
            metal = row[0].split("(")[0]
            facet = row[0].split("(")[1].replace(")", "")
            self.metal_list.append((metal, facet))

            
    def get_dh_rmg(self, surface_to = "Rh(211)"):
        """
        get diff between pt111 database value and scaled value
        """
        dH = self.Hf_rmg[surface_to] - self.Hf_old
        return dH
    
    def get_all_dHs(self, surface_from="Pt(111)"):
        
        Hf_rmg_from = self.Hf_rmg[surface_from]
        Hf_gao_from = self.Hf_gao[surface_from]
        Hf_lit_from = self.Hf_lit[surface_from]
        
        for surface in self.Hf_lit.keys():
            self.dH_rmg[surface] = self.Hf_rmg[surface] - Hf_rmg_from
            self.dH_gao[surface] = self.Hf_gao[surface] - Hf_gao_from
            self.dH_lit[surface] = self.Hf_lit[surface] - Hf_lit_from
            
    
    def get_dh_lit(self, metal_from='Pt', metal_to='Rh', 
                         facet_from='111', facet_to='211'):
        """
        get diff between literature values 
        """
        old_metal_str = f"{metal_from}({facet_from})"
        new_metal_str = f"{metal_to}({facet_to})"
        dH = self.Hf_lit[new_metal_str] - self.Hf_lit[old_metal_str]
        return dH 
    
    
    def plot_dhs(self):

        fig, ax = plt.subplots(1,2, figsize=(15, 7.5))

        ########################
        # plot rmg
        ########################
        x_rmg = list(self.dH_rmg.values())
        y = list(self.dH_lit.values())
        labels_rmg = list(self.dH_rmg.keys())

        ax[0].set_title(f"parity for predicted vs dft $\Delta H$ from Pt for {spec}*")
        ax[0].set_xlabel(r"$\Delta H_{f}^{298K}$ for rmg species from Pt(111) (eV)")
        ax[0].set_ylabel(r"$\Delta E_{B}^{0K}$ for lit species from Pt(111) (eV)")
        ax[0].plot(x_rmg,y, ".")

        # plot parity line
        buff = 0.1
        par_x = np.array([min(x_rmg) - buff, max(x_rmg) + buff])
        par_y = np.array([min(x_rmg) - buff, max(x_rmg) + buff])

        ax[0].plot(par_x, par_y, "-.", color="k")

        # fill between for 0.3 eV error on RMG
        err_y_plus = par_y + 0.3
        err_y_minus = par_y - 0.3

        ax[0].fill_between(par_x, err_y_plus, err_y_minus, color="k", alpha = 0.1)

        # label the metals
        for i, label in enumerate(labels_rmg):
            ax[0].annotate(label, (x_rmg[i], y[i]))


        ##############################
        # plot gao
        ##############################
        x_gao = list(self.dH_gao.values())
        labels_gao = list(self.dH_gao.keys())

        ax[1].set_title(f"parity for predicted vs dft $\Delta H$ from Pt for {spec}*")
        ax[1].set_xlabel(r"$\Delta H_{f}^{298K}$ for rmg species from Pt(111) (eV)")
        ax[1].set_ylabel(r"$\Delta E_{B}^{0K}$ for lit species from Pt(111) (eV)")
        ax[1].plot(x_gao,y, ".")

        # plot parity line
        buff = 0.1
        par_x = np.array([min(x_gao) - buff, max(x_gao) + buff])
        par_y = np.array([min(x_gao) - buff, max(x_gao) + buff])


        ax[1].plot(par_x, par_y, "-.", color="k")

        # fill between for 0.3 eV error on RMG
        # err_x = par_x
        err_y_plus = par_y + 0.3
        # err_x_minus = par_x - 0.3
        err_y_minus = par_y - 0.3

        ax[1].fill_between(par_x, err_y_plus, err_y_minus, color="k", alpha = 0.1)

        # label the metals
        for i, label in enumerate(labels_gao):
            ax[1].annotate(label, (x_gao[i], y[i]))
    
    


# In[21]:


pt_bes = thermoDatabase.surface['metal'].find_binding_energies('Pt111')
rh_bes = thermoDatabase.surface['metal'].find_binding_energies('Rh211')


# In[22]:


(rh_bes['C'].value_si - pt_bes['C'].value_si)/9.68e4


# In[23]:


gao_species = {}

# # C
# spec = 'C'
# adj = """
# 1 C u0 p0 c0 {2,Q}
# 2 X u0 p0 c0 {1,Q}
# """
# gao_species[spec] = SpeciesDat(adj)
# gao_species[spec].set_gao_hf(
#     pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=20, usecols="C:E"))

# gao_species[spec].scale_all_thermo()
# gao_species[spec].get_all_dHs()

# # display(spec,gao_species[spec].rmg_spec)
# gao_species[spec].plot_dhs()


# # In[24]:


# pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=20, usecols="H:J")


# In[25]:


# gao_species = {}

# # C
# spec = 'C'
# adj = """
# 1 C u0 p0 c0 {2,Q}
# 2 X u0 p0 c0 {1,Q}
# """
# gao_species[spec] = SpeciesDat(adj)
# gao_species[spec].set_gao_hf(
#     pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=20, usecols="C:E"))

# surface_from = "Pt(111)"
# gao_species[spec].scale_all_thermo()
# gao_species[spec].get_all_dHs()
# # display(spec,gao_species[spec].rmg_spec)
# gao_species[spec].plot_dhs()
gao_calcs_path = "/work/westgroup/ChrisB/_04_thesis/Thesis_repo/validation_data/gao_calcs.xlsx"

# CH
spec = 'CH'
adj = """
1 C u0 p0 c0 {2,S} {3,T}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0 {1,T}
"""
gao_species[spec] = SpeciesDat(adj)
gao_species[spec].set_gao_hf(
    pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=20, usecols="H:J"))

gao_species[spec].scale_all_thermo()
gao_species[spec].get_all_dHs()

# display(spec,gao_species[spec].rmg_spec)
gao_species[spec].plot_dhs()

# CH2
spec = 'CH2'
adj = """
1 C u0 p0 c0 {2,S} {3,S} {4,D}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 X u0 p0 c0 {1,D}
"""
gao_species[spec] = SpeciesDat(adj)
gao_species[spec] = SpeciesDat(adj)
gao_species[spec].set_gao_hf(
    pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=21, usecols="M:O"))

gao_species[spec].scale_all_thermo()
gao_species[spec].get_all_dHs()

# display(spec,gao_species[spec].rmg_spec)
gao_species[spec].plot_dhs()

# CH3
spec = 'CH3'
adj = """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0 {1,S}
"""
gao_species[spec] = SpeciesDat(adj)
gao_species[spec].set_gao_hf(
    pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=21, usecols="R:T"))

gao_species[spec].scale_all_thermo()
gao_species[spec].get_all_dHs()

# display(spec,gao_species[spec].rmg_spec)
gao_species[spec].plot_dhs()

# CO
spec = 'CO'
adj = """
1 C u0 p0 c0 {2,D} {3,D} 
2 O u0 p2 c0 {1,D}
3 X u0 p0 c0 {1,D}
"""
gao_species[spec] = SpeciesDat(adj)
gao_species[spec].set_gao_hf(
    pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=20, usecols="W:Y"))

gao_species[spec].scale_all_thermo()
gao_species[spec].get_all_dHs()

# display(spec,gao_species[spec].rmg_spec)
gao_species[spec].plot_dhs()

# COOH
spec = 'COOH'
adj = """
1 C u0 p0 c0 {2,D} {3,S} {5,S}
2 O u0 p2 c0 {1,D}
3 O u0 p2 c0 {1,S} {4,S}
4 H u0 p0 c0 {3,S}
5 X u0 p0 c0 {1,S}
"""
gao_species[spec] = SpeciesDat(adj)
gao_species[spec].set_gao_hf(
    pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=21, usecols="AB:AD"))

gao_species[spec].scale_all_thermo()
gao_species[spec].get_all_dHs()

# display(spec,gao_species[spec].rmg_spec)
gao_species[spec].plot_dhs()

# CHO
spec = 'CHO'
adj = """
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 O u0 p2 c0 {1,D}
3 H u0 p0 c0 {1,S}
4 X u0 p0 c0 {1,S}
"""
gao_species[spec] = SpeciesDat(adj)
gao_species[spec].set_gao_hf(
    pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=21, usecols="AG:AI"))

gao_species[spec].scale_all_thermo()
gao_species[spec].get_all_dHs()

# display(spec,gao_species[spec].rmg_spec)
gao_species[spec].plot_dhs()

# COH
spec = 'COH'
adj = """
1 C u0 p0 c0 {2,S} {4,T}
2 O u0 p2 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
4 X u0 p0 c0 {1,T}
"""
gao_species[spec] = SpeciesDat(adj)
gao_species[spec].set_gao_hf(
    pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=20, usecols="AL:AN"))

gao_species[spec].scale_all_thermo()
gao_species[spec].get_all_dHs()

# display(spec,gao_species[spec].rmg_spec)
gao_species[spec].plot_dhs()

# CHOH
spec = 'CHOH'
adj = """
1 C u0 p0 c0 {2,S} {4,S} {5,D}
2 O u0 p2 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0 {1,D}
"""
gao_species[spec] = SpeciesDat(adj)
gao_species[spec].set_gao_hf(
    pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=18, usecols="AQ:AS"))

gao_species[spec].scale_all_thermo()
gao_species[spec].get_all_dHs()

# display(spec,gao_species[spec].rmg_spec)
gao_species[spec].plot_dhs()

# N
spec = 'N'
adj = """
1 N u0 p1 c0 {2,T} 
2 X u0 p0 c0 {1,T}
"""
gao_species[spec] = SpeciesDat(adj)
gao_species[spec].set_gao_hf(
    pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=29, usecols="BB:BD"))

gao_species[spec].scale_all_thermo()
gao_species[spec].get_all_dHs()

# display(spec,gao_species[spec].rmg_spec)
gao_species[spec].plot_dhs()

# NH
spec = 'NH'
adj = """
1 N u0 p1 c0 {2,S} {3,D} 
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0 {1,D}
"""
gao_species[spec] = SpeciesDat(adj)
gao_species[spec].set_gao_hf(
    pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=29, usecols="BG:BI"))

gao_species[spec].scale_all_thermo()
gao_species[spec].get_all_dHs()

# display(spec,gao_species[spec].rmg_spec)
gao_species[spec].plot_dhs()

# NNH2 
spec = 'NNH2'
adj = """
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 N u0 p1 c0 {1,S} {5,D} 
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0 {2,D}
"""
gao_species[spec] = SpeciesDat(adj)
gao_species[spec].set_gao_hf(
    pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=22, usecols="BL:BN"))

gao_species[spec].scale_all_thermo()
gao_species[spec].get_all_dHs()

# display(spec,gao_species[spec].rmg_spec)
gao_species[spec].plot_dhs()

# NH2
spec = 'NH2'
adj = """
1 N u0 p1 c0 {2,S} {3,S} {4,S} 
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 X u0 p0 c0 {1,S}
"""
gao_species[spec] = SpeciesDat(adj)
gao_species[spec].set_gao_hf(
    pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=20, usecols=""))

gao_species[spec].scale_all_thermo()
gao_species[spec].get_all_dHs()

# display(spec,gao_species[spec].rmg_spec)
gao_species[spec].plot_dhs()

# O
spec = 'O'
adj = """
1 O u0 p2 c0 {2,D}
2 X u0 p0 c0 {1,D}
"""
gao_species[spec] = SpeciesDat(adj)
gao_species[spec].set_gao_hf(
    pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=20, usecols=""))

gao_species[spec].scale_all_thermo()
gao_species[spec].get_all_dHs()

# display(spec,gao_species[spec].rmg_spec)
gao_species[spec].plot_dhs()

# OH
spec = 'OH'
adj = """
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0 {1,S}
"""
gao_species[spec] = SpeciesDat(adj)
gao_species[spec].set_gao_hf(
    pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=20, usecols=""))

gao_species[spec].scale_all_thermo()
gao_species[spec].get_all_dHs()

# display(spec,gao_species[spec].rmg_spec)
gao_species[spec].plot_dhs()

# OOH
spec = 'OOH'
adj = """
1 O u0 p2 c0 {2,S} {4,S}
2 O u0 p2 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
4 X u0 p0 c0 {1,S}
"""
gao_species[spec] = SpeciesDat(adj)
gao_species[spec].set_gao_hf(
    pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=20, usecols=""))

gao_species[spec].scale_all_thermo()
gao_species[spec].get_all_dHs()

#display(spec,gao_species[spec].rmg_spec)
gao_species[spec].plot_dhs()

# OCH3
spec = 'OCH3'
adj = """
1 O u0 p2 c0 {2,S} {6,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 X u0 p0 c0 {1,S}
"""
gao_species[spec] = SpeciesDat(adj)
gao_species[spec].set_gao_hf(
    pd.read_excel(gao_calcs_path, sheet_name = "Fig. 7", index_col=0, skiprows=[0], nrows=20, usecols=""))

gao_species[spec].scale_all_thermo()
gao_species[spec].get_all_dHs()

#display(spec,gao_species[spec].rmg_spec)
gao_species[spec].plot_dhs()


# ## Read in gao data






# %%
