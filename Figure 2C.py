import os
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
from scipy.sparse import issparse
import pickle

def plot_category(adata, group, background_alpha, col_map, subset=None, legend=False, output=None, dpi=300, figsize=(10,10)):
    col_map = {each:col_map[each] for each in col_map if each in adata.obs[group].tolist()}
    col_order = list(col_map.keys())
    color_values = [col_map[t] for t in col_order]
    temp_color_cmap_pa = mp.colors.ListedColormap(color_values)
    color_vector = [col_order.index(each) for each in adata.obs[group]]
    
    if subset is None:
        pixel_row = adata.obs['pixel_row']
        pixel_col = adata.obs['pixel_col']
    else:
        index = adata.obs[group].isin(subset)
        pixel_row = adata.obs['pixel_row'].loc[index]
        pixel_col = adata.obs['pixel_col'].loc[index]
        score = [score[i] for i in range(len(index)) if index[i]]
    
    fig,ax = plt.subplots(figsize=figsize)
    ax.imshow(list(adata.uns['spatial'].values())[0]['images']['hires'], alpha = background_alpha)
    img_ax = ax.hexbin(pixel_col, pixel_row, C = color_vector,
                    gridsize = adata.uns['hexagon']['gridsize'],
                    edgecolors = 'face', linewidth = -0.15,
                    cmap=temp_color_cmap_pa,
                    extent = adata.uns['hexagon']['pixel_extent'],
                    alpha = 1)
    ax.axis('off')
    ax.set_frame_on(False)
    if legend:
        legend_handles = [Patch(color=color, label=category) for category, color in col_map.items()]
        legend = ax.legend(handles=legend_handles, loc='center', title = group, bbox_to_anchor=(1.14, 0.5), fontsize='large')
    plt.savefig(output,dpi=dpi)
    plt.clf()

def plot_continuous(adata, variable, background_alpha, cmap = 'magma', Vmin=None, Vmax=None, group=None, subset = None, legend=True, output=None, figsize = (10, 10)):
    if variable in adata.var_names:
        gene_index = adata.var_names.get_loc(variable)
        score = adata.X[:, gene_index]
        if issparse(score):
            score = score.toarray().ravel().tolist()
        else:
            score = np.array(score).ravel().tolist()
    else:
        score = adata.obs[variable].tolist()
    
    if subset is None:
        pixel_row = adata.obs['pixel_row']
        pixel_col = adata.obs['pixel_col']
    else:
        index = adata.obs[group].isin(subset)
        pixel_row = adata.obs['pixel_row'].loc[index]
        pixel_col = adata.obs['pixel_col'].loc[index]
        score = [score[i] for i in range(len(index)) if index[i]]
    
    fig,ax = plt.subplots(figsize=figsize)
    ax.imshow(list(adata.uns['spatial'].values())[0]['images']['hires'], alpha = background_alpha)
    img_ax = ax.hexbin(pixel_col, pixel_row, C = score,
                    gridsize = adata.uns['hexagon']['gridsize'],
                    edgecolors = 'face', linewidth = 0,
                    cmap=cmap,
                    extent = adata.uns['hexagon']['pixel_extent'],
                    vmin = Vmin, vmax = Vmax,
                    alpha = 1)
    ax.axis('off')
    ax.set_frame_on(False)
    if legend:
        cbar_ax = fig.add_axes([0.9, 0.4, 0.03, 0.3]) 
        fig.colorbar(img_ax, cax=cbar_ax, location='right')
    plt.savefig(output)
    plt.clf()

prefix = os.getcwd()

with open(f'{prefix}/Data/Figure 2C.pkl', 'rb') as f:
    data = pickle.load(f)

col_map_histology = data['col_map_histology']
col_map_clone = data['col_map_clone']

for s in ['P3_AIS', 'P3_LUAD', 'P20_AAH', 'P20_LUAD', 'P23_AIS1', 'P23_AIS2', 'P23_LUAD']:
    input = data[s]
    plot_category(adata = input, group = 'Histology', background_alpha=1, col_map = col_map_histology, 
                  output = f'{prefix}/Results/Figure 2C.{s}.Histology.png', figsize = (10, 10))
    plot_category(adata = input[~input.obs.Clone.isna()], group = 'Clone',background_alpha=0.3, col_map = col_map_clone, 
                  output = f'{prefix}/Results/Figure 2C.{s}.Clone.png', figsize = (10, 10))
    plot_continuous(adata = input[~input.obs.Pseudotime.isna()], variable = 'Pseudotime',background_alpha=0.3, Vmin=input.obs.Pseudotime.min(), Vmax=input.obs.Pseudotime.max(), 
                    output = f'{prefix}/Results/Figure 2C.{s}.Pseudotime.png')
    plot_continuous(adata = input[~input.obs.CytoTRACE_Score.isna()], variable = 'CytoTRACE_Score',background_alpha=0.3, cmap = 'turbo', Vmin=input.obs.CytoTRACE_Score.min(), Vmax=input.obs.CytoTRACE_Score.max(), 
                    output = f'{prefix}/Results/Figure 2C.{s}.CytoTRACE.png')

