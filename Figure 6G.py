import os
import pickle
import matplotlib.pyplot as plt
import commot as ct

prefix = os.getcwd()

with open(f'{prefix}/Data/Figure 6G.pkl', 'rb') as f:
    data = pickle.load(f)

for s in ['A4','B4','C4','D4']:
    input = data[s]
    input.layers['counts'] = input.X
    ct.tl.communication_direction(input, database_name='cellchat', pathway_name='IL1', k=5)
    plt.figure(num=1, figsize=(8, 8),dpi=300)
    ct.pl.plot_cell_communication(input, database_name='cellchat', pathway_name='IL1', plot_method='stream', background_legend=True,
                                  scale=0.00003, ndsize=8, grid_density=0.4, summary='sender', background='image', clustering='Histology', cmap='Alphabet',
                                  normalize_v = True, normalize_v_quantile=0.995)
    plt.savefig(f'{prefix}/Results/Figure 6G.'+s+'.CCI.png', format='png', dpi=300)
