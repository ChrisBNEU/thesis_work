#!/usr/bin/env python
# encoding: utf-8


# for now, focusing on adding in the normal, non=defect sites for 111 surfaces:
# 111 top site
# 111 bridge site
# 111 fcc hollow site
# 111 hcp hollow site

name = "Metal properties"
shortDesc = ""
longDesc = """

subsequent use for advanced scaling relations from Gao et al. in "Determining the 
adsorption energies of small molecules with the intrinsic properties of adsorbates 
and substrates", 2020, doi: 10.1038/s41467-020-14969-8

valence_electrons: number of valence electrons for the metal atom
electronegativity: pauling electronegativity of the metal atom (I assume in eV, 
not sure if it matters since it is a relative scale)
phi: 
"""

entry(
    index = 1,
    label = "Sc",
    metal = "Sc",
    valence_electrons = 3,
    electronegativity = 1.36,
    phi = 6.62,
    beta = 1,
)

entry(
    index = 2,
    label = "Ti",
    metal = "Ti",
    valence_electrons = 4,
    electronegativity = 1.54,
    phi = 10.39,
    beta = 1,
)

entry(
    index = 3,
    label = "V",
    metal = "V",
    valence_electrons = 5,
    electronegativity = 1.63,
    phi = 15.34,
    beta = 1,
)

entry(
    index = 4,
    label = "Cr",
    metal = "Cr",
    valence_electrons = 6,
    electronegativity = 1.66,
    phi = 21.69,
    beta = 1,
)

entry(
    index = 5,
    label = "Mn",
    metal = "Mn",
    valence_electrons = 7,
    electronegativity = 1.55,
    phi = 31.61,
    beta = 1,
)
# copilot, please suggest multiple lines of code to fill in the rest of the database entries
# for the metals in the list below

entry(
    index = 6,
    label = "Fe",
    metal = "Fe",
    valence_electrons = 8,
    electronegativity = 1.83,
    phi = 34.97,
    beta = 1,
)

entry(
    index = 7,
    label = "Co",
    metal = "Co",
    valence_electrons = 9,
    electronegativity = 1.88,
    phi = 43.09,
    beta = 1,
)

entry(
    index = 8,
    label = "Ni",
    metal = "Ni",
    valence_electrons = 10,
    electronegativity = 1.91,
    phi = 52.36,
    beta = 1,
)

entry(
    index = 9,
    label = "Cu",
    metal = "Cu",
    valence_electrons = 11,
    electronegativity = 1.90,
    phi = 63.68,
    beta = 1,
)

entry(
    index = 10,
    label = "Zn",
    metal = "Zn",
    valence_electrons = 12,
    electronegativity = 1.65,
    phi = 87.27,
    beta = 1,
)

entry(
    index = 11,
    label = "Y",
    metal = "Y",
    valence_electrons = 3,
    electronegativity = 1.22,
    phi = 7.38,
    beta = 1,
)

entry(
    index = 12,
    label = "Zr",
    metal = "Zr",
    valence_electrons = 4,
    electronegativity = 1.33,
    phi = 12.03,
    beta = 1,
)

entry(
    index = 13,
    label = "Nb",
    metal = "Nb",
    valence_electrons = 5,
    electronegativity = 1.60,
    phi = 15.63,
    beta = 1,
)

entry(
    index = 14,
    label = "Mo",
    metal = "Mo",
    valence_electrons = 6,
    electronegativity = 1.16,
    phi = 31.03,
    beta = 1,
)

entry(
    index = 15,
    label = "Tc",
    metal = "Tc",
    valence_electrons = 7,
    electronegativity = 1.90,
    phi = 25.79,
    beta = 1,
)

entry(
    index = 16,
    label = "Ru",
    metal = "Ru",
    valence_electrons = 8,
    electronegativity = 2.20,
    phi = 29.09,
    beta = 1,
)

entry(
    index = 17,
    label = "Rh",
    metal = "Rh",
    valence_electrons = 9,
    electronegativity = 2.28,
    phi = 35.53,
    beta = 1,
)

entry(
    index = 18,
    label = "Pd",
    metal = "Pd",
    valence_electrons = 10,
    electronegativity = 2.20,
    phi = 45.45,
    beta = 1,
)

entry(
    index = 19,
    label = "Ag",
    metal = "Ag",
    valence_electrons = 11,
    electronegativity = 1.93,
    phi = 87.10,
    beta = 1,
)

entry(
    index = 20,
    label = "Cd",
    metal = "Cd",  
    valence_electrons = 12,
    electronegativity = 1.69,
    phi = 85.21,
    beta = 1,
)

entry(
    index = 21,
    label = "La",
    metal = "La",
    valence_electrons = 3,
    electronegativity = 1.10,
    phi = 8.18,
    beta = 1,
)

entry(
    index = 22,
    label = "Hf",
    metal = "Hf",
    valence_electrons = 4,
    electronegativity = 1.30,
    phi = 12.31,
    beta = 1,
)

entry(
    index = 23,
    label = "Ta",
    metal = "Ta",
    valence_electrons = 5,
    electronegativity = 1.50,
    phi = 16.67,
    beta = 1,
)

entry(
    index = 24,
    label = "W",
    metal = "W",
    valence_electrons = 6,
    electronegativity = 2.36,
    phi = 15.25,
    beta = 1,
)

entry(
    index = 25,
    label = "Re",
    metal = "Re",
    valence_electrons = 7,
    electronegativity = 1.90,
    phi = 25.79,
    beta = 1,
)

entry(
    index = 26,
    label = "Os",
    metal = "Os",
    valence_electrons = 8,
    electronegativity = 2.20,
    phi = 29.09,
    beta = 1,
)

entry(
    index = 27,
    label = "Ir",
    metal = "Ir",
    valence_electrons = 9,
    electronegativity = 2.20,
    phi = 36.82,
    beta = 1,
)

entry(
    index = 28,
    label = "Pt",
    metal = "Pt",
    valence_electrons = 10,
    electronegativity = 2.28,
    phi = 43.86,
    beta = 1,
)

entry(
    index = 29,
    label = "Au",
    metal = "Au",
    valence_electrons = 11,
    electronegativity = 2.54,
    phi = 75.92,
    beta = 1,
)

entry(
    index = 30,
    label = "Hg",
    metal = "Hg",
    valence_electrons = 12,
    electronegativity = 2.00,
    phi = 72,
    beta = 1,
)















metal_dict = {
    "Sc": [3, 1.36, 6.62],
    "Ti": [4, 1.54, 10.39], 
    "V": [5, 1.63, 15.34], 
    "Cr": [6, 1.66, 21.69], 
    "Mn": [7, 1.55, 31.61], 
    "Fe": [8, 1.83, 34.97], 
    "Co": [9, 1.88, 43.09], 
    "Ni": [10, 1.91, 52.36], 
    "Cu": [11, 1.90, 63.68], 
    "Zn": [12, 1.65, 87.27], 
    "Y": [3, 1.22, 7.38], 
    "Zr": [4, 1.33, 12.03], 
    "Nb": [5, 1.60, 15.63], 
    "Mo": [6, 1.16, 31.03], 
    "Tc": [7, 1.90, 25.79], 
    "Ru": [8, 2.20, 29.09], 
    "Rh": [9, 2.28, 35.53], 
    "Pd": [10, 2.20, 45.45], 
    "Ag": [11, 1.93, 87.10],
    "Cd": [12, 1.69, 85.21], 
    "La": [3, 1.10, 8.18], 
    "Hf": [4, 1.30, 12.31], 
    "Ta": [5, 1.50, 16.67], 
    "W": [6, 2.36, 15.25], 
    "Re": [7, 1.90, 25.79], 
    "Os": [8, 2.20, 29.09],
    "Ir": [9, 2.20, 36.82], 
    "Pt": [10, 2.28, 43.86],
    "Au": [11, 2.54, 75.92], 
    "Hg": [12, 2.00, 72],
}






