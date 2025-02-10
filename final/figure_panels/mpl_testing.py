#%%
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl


fig, axs = plt.subplot_mosaic([['x'],['y'],['z']],
                              gridspec_kw={'height_ratios':[1, 1,1 ]},
                              constrained_layout=True)

pc = axs['x'].pcolormesh(np.random.randn(20, 20))
fig.colorbar(pc, ax=axs['x'], location='bottom', pad=0.25)
plt.show()


#%%

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def do_stuff(cell): #just so the plots show up
    ax = plt.subplot(cell)
    ax.plot()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
plt.subplots_adjust(hspace=0.0)
#make outer gridspec
outer = gridspec.GridSpec(2, 1, height_ratios = [1, 2]) 
#make nested gridspecs
gs1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec = outer[0], hspace = 1)
gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec = outer[1], hspace = 0)
for cell in gs1:
    do_stuff(cell)
for cell in gs2:
    do_stuff(cell)
plt.show()



#%%

nrows = 3

fig = plt.figure()
fig.set_figheight(nrows * 3)
sfigs = fig.subfigures(2, 1, height_ratios = [1,2])
axstop = sfigs[0].subplots(1, 1)
axsbot = sfigs[1].subplots(2, 1)

sfigs[1].subplots_adjust(hspace = 0)

#%%


n_rows = 3 #Scale figure height based on number of plots

fig.set_figheight(n_rows * 3)

ax1 = fig.add_subplot(2,1,1)

fig2 = fig.add_subfigure(1,2)

ax2 = fig2.add_subplot(2,1,1)
ax3 = fig2.add_subplot(2,1,2, sharex = ax2)

fig.subplots_adjust(hspace = 0.25)



#%%

nrows = 3
fig = plt.figure()

fig.set_figheight(nrows * 3)
fig.set_figwidth(5)

# fig.set_figheight(8)

sfigs = fig.subfigures(2,1, height_ratios = [1,2], hspace = 0)

sfigs[0].subplots(1,2,sharey = True)

(ax1, ax2) = sfigs[1].subplots(nrows = 2)



labels = ['B)', 'C)']

axlist = [ax1, ax2]

for ax, label in zip(axlist, labels):
    ax.text(-0.15, 1.1, label, transform=ax.transAxes, va='bottom', ha='left')

sfigs[0].get_axes()[0].text(-0.15, 2.75, 'A)', transform=ax1.transAxes)
