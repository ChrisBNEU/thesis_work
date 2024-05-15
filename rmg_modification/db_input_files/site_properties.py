#!/usr/bin/env python
# encoding: utf-8


# for now, focusing on adding in the normal, non=defect sites for 111 surfaces:
# 111 top site
# 111 bridge site
# 111 fcc hollow site
# 111 hcp hollow site

name = "Metal Binding Energies"
shortDesc = ""
longDesc = """

generalized coordination numbers first postulated by Calle-Vallejo et al. in 
"Fast Prediction of Adsorption Properties for Platinum Nanocatalysts with Generalized 
Coordination Numbers", 2014 doi: 10.1002/anie.201402958

subsequent use for advanced scaling relations from Gao et al. in "Determining the 
adsorption energies of small molecules with the intrinsic properties of adsorbates 
and substrates", 2020, doi: 10.1038/s41467-020-14969-8
"""
entry(
    index = 1,
    label = "111_top",
    facet = "111",
    metal = "Pt",
    site = "Terrace"
    coordination_number =7.5,
    shortDesc = """fcc""",
    longDesc = 
"""
111 surface on-top site. 
""",
)

entry(
    index = 1,
    label = "111_bridge",
    facet = "111",
    metal = "Pt",
    site = "Terrace"
    coordination_number =7.33,
    shortDesc = """fcc""",
    longDesc = 
"""
111 surface bridge site.
""",
)

entry(
    index = 0,
    label = "111_fcc_hollow",
    facet = "111",
    metal = "Pt",
    site = ""
    coordination_number = 7.5,
    shortDesc = """fcc""",
    longDesc = 
"""
111 surface FCC hollow site. 
""",
)

entry(
    index = 1,
    label = "111_hcp_hollow",
    facet = "111",
    metal = "Pt",
    site = "Terrace"
    coordination_number = 6.955,
    shortDesc = """fcc""",
    longDesc = 
"""
111 surface HCP hollow site.
""",
)






