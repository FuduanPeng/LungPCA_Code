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

with open(f'{prefix}/Data/Figure 3K 3N.pkl', 'rb') as f:
    data = pickle.load(f)

col_map_histology = data['col_map_histology']
col_map_clone = data['col_map_clone']
MP_colors = data['MP_colors']
col_map_RPII_clone = data['col_map_RPII_clone']

# Figure 3K
for s in ['P12_AIS', 'P12_LUAD', 'P24_AAH', 'P24_LUAD']:
    input = data[s]
    plot_category(adata = input, group = 'Histology', background_alpha=1, col_map = col_map_histology, 
                  output = f'{prefix}/Results/Figure 3K.{s}.Histology.png', figsize = (10, 10))
    plot_category(adata = input[~input.obs.Clone.isna()], group = 'Clone',background_alpha=0.3, col_map = col_map_clone, 
                  output = f'{prefix}/Results/Figure 3K.{s}.Clone.png', figsize = (10, 10))
    
    stacked_MPs = data[str(s)+'_stacked_MPs']
    max_scores = np.nanmax(stacked_MPs, axis=-1)
    label_image = data[str(s)+'_label_image']
    for i in range(stacked_MPs.shape[0]):
      for j in range(stacked_MPs.shape[1]):
          max_score = max_scores[i, j]
          if np.isnan(max_score):
              label_image[i, j] = [0, 0, 0]  
          else:
              max_index = np.where(stacked_MPs[i, j] == max_score)[0][0]
              label_image[i, j] = MP_colors[[*MP_colors.keys()][max_index]]
    label_image = Image.fromarray(label_image)
    label_image.save(f'{prefix}/Results/Figure 3K.'+str(s)+'.MPs.png')

# Figure 3N
P12_AIS_RPII_clone = data['P12_AIS_RPII_clone']
plot_category(adata = P12_AIS_RPII_clone,group = 'Clone_Type',background_alpha=1, col_map = col_map_RPII_clone, 
              output = f'{prefix}/Results/Figure 3N.RPII_clone.png', figsize = (10, 10))
