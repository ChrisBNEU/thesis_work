#!/usr/bin/env python
# coding: utf-8
def get_pref_site(species, metal_atoms, cn_nums, bound_atom, bonds):
    """
    gets the preferred site for a single adatom

    :param species: the adsorbate species
    :param metal_atoms: the number of metal atoms available on the surface
    :param cn_nums: the coordination numbers for the different sites on the surface
    :return: the preferred site
    """
    # gets the preferred site for a single adatom
    # get max metal atoms available for the facet (e.g. 111 = 3 for hollow site) 
    max_ma = max(metal_atoms.values())

    if bound_atom.symbol == 'C':
        pref_site = bonds
    elif bound_atom.symbol == 'O': 
        pref_site = bonds + 1
    elif bound_atom.symbol == 'N': 
        pref_site = bonds + 1
    elif bound_atom.symbol == 'H':
        pref_site = 1

    # if the preference is higher than the largest num of metal atoms, 
    # then just go to the site with most available metal atoms
    pref_site = max_ma if pref_site > max_ma else pref_site

    # get all possible sites on the surface
    poss_sites = []
    for poss_site, ma_num in metal_atoms.items():
        if pref_site == ma_num: 
            poss_sites.append(poss_site)

    # if multiple suitable sites are found, use one with lowest coordination number 
    # site 1:
    if len(poss_sites) == 1:
        site = poss_sites[0]

    elif len(poss_sites1) > 1:
        ma_site = [k for k, v in metal_atoms1.items() if v == pref_site1]

        # get the coordination numbers corresponding to the max ma sites
        cns_ma = {k:cn_nums1[k] for k in ma_site1}
        site = min(cns_ma, key=cns_ma.get)

    # no sites found, see if it is larger than max sites.
    elif len(poss_sites) == 0:
        raise ValueError(f"no site could be found for bound atom {bound_atom.symbol} in species {species.label}")
    
    return site

def get_sites(species, metal_atoms1, cn_nums1, metal_atoms2, cn_nums2):
    """
    gets the adatoms, bond orders, and sites for a given species. 

    :param species: the adsorbate species
    :param metal_atoms: the number of metal atoms available on the surface
    :param cn_nums: the coordination numbers for the different sites on the surface
    :return: a list of tuples with the bound atom, site, bond order, and max bond order
    """
    sites = {}
    max_bond_orders = {'C': 4., 'O': 2., 'N': 3., 'H': 1.}

    for atom in species.molecule[0].atoms: 
        if atom.is_surface_site():
            # vdw check
            if len(atom.bonds) == 0:
                return None
            else:
                bound_atom = list(atom.bonds.keys())[0]
                # what do I do for bonds? see where it is used and if I can use total bonds for getting preferred site. 
                bonds = list(atom.bonds.values())[0].get_order_num()
                # max_bond_order = max_bond_orders[bound_atom.symbol]

                # I am hardcoding this for now. if there is an oh group, 
                # then there effectively 2 functional groups in the advanced lsr equation
                count_coh = 0
                count_cr = 0
                if bound_atom.symbol == 'C':
                    for atom, bond in bound_atom.bonds.items():
                        # check for OH, can have 5 bonds to C
                        if atom.is_oxygen() and bond.is_single():
                            count_oh += 1
                        else: 
                            groups_xm.append(4.)
                            count_r += bond.get_order_num()
                            
                if count_coh > 0:
                    xm1 = 5
                    xm2 = 4
                    x1 = count_coh
                    x2 = count_cr

                else: 
                    xm1 = max_bond_orders[bound_atom.symbol]
                    xm2 = 0
                    x1 = xm1 - bonds
                    x2 = 0

                # calculate alphas 
                alpha1 = (xm1 - x1)/(xm1 + 1)
                alpha2 = (x2)/(xm2 + 1)
                alpha = alpha1 - alpha2

                # if the site is specified, use it for metal one (i.e. the database we are scaling from)
                if len(atom.site) > 0:
                    site1 = atom.site
                else: 
                    site1 = self.get_pref_site(species, metal_atoms1, cn_nums1, bound_atom, bonds)
                
                site2 = self.get_pref_site(species, metal_atoms2, cn_nums2, bound_atom, bonds)

            sites[bound_atom] = {"site1": site1, "site2": site2, "alpha" : alpha}

    return sites

def correct_binding_energies_extended(thermo, species, metal_to_scale_from=None,
                                    metal_to_scale_to=None, facet_to_scale_from=None, facet_to_scale_to=None):
    """
    Uses the binding energy correction proposed by gao, which allows for scaling
    from one metal and facet (e.g. Pt111) to a completely different metal and 
    facet (e.g. Cu111). 

    :param thermo: the thermo data to correct
    :param species: the species to correct
    :param metal_to_scale_from: the metal to scale from (e.g. Pt)
    :param metal_to_scale_to: the metal to scale to (e.g. Cu)
    :param facet_to_scale_from: the facet to scale from (e.g. 111)
    :param facet_to_scale_to: the facet to scale to (e.g. 111)
    :return: the corrected thermo data

    """

    if metal_to_scale_from == metal_to_scale_to and facet_to_scale_from == facet_to_scale_to:
        return thermo
    elif metal_to_scale_from is None or metal_to_scale_to is None or facet_to_scale_from is None or facet_to_scale_to is None:
        raise ValueError("If you are scaling, you must specify both the metal and the facet to scale from and to.")

    # get the required attributes for every facet and metal involved
    cn1_dict = self.surface['site_properties'].get_all_coordination_numbers_on_facet(facet_to_scale_from)
    cn2_dict = self.surface['site_properties'].get_all_coordination_numbers_on_facet(facet_to_scale_to)

    ma1_dict = self.surface['site_properties'].get_all_metal_atoms_on_facet(facet_to_scale_from)
    ma2_dict = self.surface['site_properties'].get_all_metal_atoms_on_facet(facet_to_scale_to)

    psi1 = self.surface['metal_properties'].get_psi(metal_to_scale_from)
    psi2 = self.surface['metal_properties'].get_psi(metal_to_scale_to)

    # determine the preferred sites for the species on each surface
    # only call once so we can link the sites
    surf_sites = self.get_sites(species, ma1_dict, cn1_dict, ma2_dict, cn2_dict)

    # check for vdw
    if not surf_sites or len(surf_sites) == 0: 
        print(f"species {species.label} has no sites")
        print(thermo.H298.value_si/9.68e4)
        return thermo

    # print for logging
    metal1_str = f"{metal_to_scale_from}({facet_to_scale_from})"
    metal2_str = f"{metal_to_scale_to}({facet_to_scale_to})"

    theta = 0
    E_ads_new = 0

    # scale the enthalpy of formation at 298. if there are multiple bound atoms
    # then we will add like we do in lsrs
    BE_diff = 0
    comments = []
    for atom, site_info in surf_sites.items():
        site1 = site_info['site1']
        site2 = site_info['site2']
        alpha = site_info['alpha']
        cn1 = cn1_dict[site1]
        cn2 = cn2_dict[site2]

        # get relative difference in binding energies. assuming relative enthalpy diff from
        # Hf0k to Hf298 remains the same for both species. 
        BE_diff += 0.1*alpha*(psi2-psi1) + 0.2*(1-alpha)*(cn2-cn1)

        comments.append(f"{atom.symbol} scaled from {metal1_str} at {site1} to {metal2_str} at site {site2} with alpha={alpha:.2f}")

    # update the enthalpy of formation to it's new value
    H298_orig = (thermo.H298.value_si)/9.68e4
    H298_new = H298_orig + BE_diff
    thermo.comment += " Binding energy corrected by Gao relations ({})".format(', '.join(comments))

    # adjust the H298
    thermo.H298.value_si = H298_new*9.68e4

    return thermo

# test the function. load model each time 
curr_dir = os.path.dirname(os.path.abspath(__file__))
folder = os.path.join(curr_dir, 'test_mech', 'chemkin')

chemkin_path = os.path.join(folder, 'chem_annotated-gas.inp')
chemkin_surf_path = os.path.join(folder, 'chem_annotated-surface.inp')
dictionary_path = os.path.join(folder, 'species_dictionary.txt')

# load_chemkin_file
species, reactions = load_chemkin_file(chemkin_surf_path, dictionary_path)
# identify a surface species to use for testing
surf_species = species[5]
surf_species.generate_resonance_structures()
# display(surf_species)
spec_thermo = surf_species.thermo.to_thermo_data()
surf_species.thermo = spec_thermo
thermo_list = self.get_all_thermo_data(surf_species)
E_ad_old = (spec_thermo.H298.value_si)/9.68e4

new_thermo = correct_binding_energies_extended(spec_thermo, surf_species, metal_to_scale_from='Pt',
                                       metal_to_scale_to='Rh', facet_to_scale_from='111', facet_to_scale_to='211')

E_ad_new = (new_thermo.H298.value_si)/9.68e4
print("old E_ad = ", E_ad_old)
print("new E_ad = ", E_ad_new)
print(new_thermo.comment)



