import numpy as np
import os,sys,pickle
from catmap import ReactionModel
from catmap import analyze
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

model = ReactionModel(setup_file='c3.mkm')
model.output_variables.append('production_rate')
model.output_variables.append('free_energy')
model.output_variables.append('reverse_rate')
model.output_variables.append('forward_rate')
model.output_variables.append('rxn_order')
model.run()

'''
mm = analyze.MatrixMap(model)
mm.plot_variable = 'rate_control'
mm.log_scale = False
mm.min = -2
mm.max = 2
mm.plot(save='rate_control.pdf')

exit()
'''

vm = analyze.VectorMap(model)
vm.plot_variable = 'production_rate'
vm.log_scale = True
vm.min = 1e-10
vm.max = 1e1
vm.threshold = 1e-25
vm.include_labels = ['CO2_g']
vm.descriptor_labels = ['T', 'P']
vm.subplots_adjust_kwargs = {'left':0.2,'right':0.8,'bottom':0.15}
vm.plot(save='CO2.pdf')

model.descriptor_dict = {'Pd':[425,1]}
ma = analyze.MechanismAnalysis(model)
label_size = 8
ma.kwarg_dict = {
        'non-assisted':{'label_positions':'top','label_args':{'color':'b','rotation':45,'ha':'left','size':label_size}},
	'O-assisted':{'label_positions':'bot','label_args':{'color':'r','rotation':45,'ha':'right','size':label_size}}
}
ma.subplots_adjust_kwargs = {'top': 0.8, 'bottom':0.22, 'left':0.15}
ma.energy_type = 'free_energy' #can also be free_energy
ma.include_labels = True
ma.pressure_correction = False
#ma.coverage_correction = False
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylabel('$\Delta G$ / eV',fontsize=12)
ax.tick_params(labelsize=12)
colors = ['b','r']
lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='-') for c in colors]
labels = ['direct','O*-assisted']
plt.legend(lines,labels)
ma.plot_variant_colors=['b','r']
ma.plot(ax=ax,save='FED.pdf')
print ma.data_dict
