import os
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pickle
import numpy as np
from PIL import Image

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

prefix = os.getcwd()

with open(f'{prefix}/Data/Figure 6E 6F 7D.pkl', 'rb') as f:
    data = pickle.load(f)

col_map_histology = data['col_map_histology']
score_colors = data['score_colors']

# Figure 6E
for s in ['A2', 'B2', 'A4', 'B4', 'C4', 'D4']:
    input = data[s]
    plot_category(adata = input, group = 'Histology', background_alpha=1, col_map = col_map_histology, 
                  output = f'{prefix}/Results/Figure 6E.{s}.Histology.png', figsize = (10, 10))

# Figure 6F
for s in ['A2', 'B2', 'A4', 'B4', 'C4', 'D4']:
    stacked_scores=data[str(s)+'_stacked_scores']
    Inflammatory_score=stacked_scores[...,0]
    KAC_score=stacked_scores[...,1]
    max_scores = np.nanmax(stacked_scores, axis=-1)
    label_image=data[str(s)+'_label_image']
    q70=[np.quantile(Inflammatory_score[~np.isnan(Inflammatory_score)],0.7),np.quantile(KAC_score[~np.isnan(KAC_score)],0.7)]
    q90=[np.quantile(Inflammatory_score[~np.isnan(Inflammatory_score)],0.9),np.quantile(KAC_score[~np.isnan(KAC_score)],0.9)]
    for i in range(stacked_scores.shape[0]):
        for j in range(stacked_scores.shape[1]):
            max_score = max_scores[i, j]
            if np.isnan(max_score):
                label_image[i, j] = [0, 0, 0] 
            elif (Inflammatory_score[i,j] > q90[0]) & (KAC_score[i,j] > q90[1]):
                label_image[i, j] = [252,127,0] 
            elif (Inflammatory_score[i,j] < q70[0]) & (KAC_score[i,j] < q70[1]):
                label_image[i, j] = [0, 0, 0]
            else:
                max_index = np.where(stacked_scores[i, j] == max_score)[0][0]
                label_image[i, j] = score_colors[[*score_colors.keys()][max_index]]
    label_image = Image.fromarray(label_image)
    label_image.save(f'{prefix}/Results/Figure 6F.'+str(s)+'.overlay.png')


# Figure 7D
for s in ['IgG_3m', 'antiIL1b_3m_1', 'antiIL1b_3m_2', 'IgG_7m', 'antiIL1b_7m_1', 'antiIL1b_7m_2']:
    input = data[s]
    plot_category(adata = input, group = 'Histology', background_alpha=1, col_map = col_map_histology, 
                  output = f'{prefix}/Results/Figure 7D.{s}.Histology.png', figsize = (10, 10))
