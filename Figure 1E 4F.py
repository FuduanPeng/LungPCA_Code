import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import pickle


def mkdir(path):
    dirname = os.path.dirname(path)
    if dirname != '':
        os.makedirs(dirname, exist_ok=True)

def load_pickle(filename, verbose=True):
    with open(filename, 'rb') as file:
        data = pickle.load(file)
    if verbose:
        print(f'Pickle loaded from {filename}')
    return data

def load_image(filename, verbose=True):
    img = Image.open(filename)
    img = np.array(img)
    if img.ndim == 3 and img.shape[-1] == 4:
        img = img[..., :3]  # remove alpha channel
    if verbose:
        print(f'Image loaded from {filename}')
    return img

def save_image(img, filename):
    mkdir(filename)
    Image.fromarray(img).save(filename)
    print(filename)

def plot_super(x, outfile, underground=None, truncate=None):
    x = x.copy()
    mask = np.isfinite(x)
    if truncate is not None:
        x -= np.nanmean(x)
        x /= np.nanstd(x) + 1e-12
        x = np.clip(x, truncate[0], truncate[1])
    x -= np.nanmin(x)
    x /= np.nanmax(x) + 1e-12
    cmap = plt.get_cmap('turbo')
    # cmap = cmap_turbo_truncated
    if underground is not None:
        under = underground.mean(-1, keepdims=True)
        under -= under.min()
        under /= under.max() + 1e-12
    img = cmap(x)[..., :3]
    if underground is not None:
        img = img * 0.5 + under * 0.5
    img[~mask] = 0
    img = (img * 255).astype(np.uint8)
    save_image(img, outfile)


# Figure 1E
prefix = os.getcwd()

sample_list = ['P24_LUAD','P8_LUAD','P11_AAH','P6_LUAD','P3_AIS','P15_MIA']
gene_list = ['AGER','CEACAM5','SFTPC','MS4A1','APOE','COL14A1'] 

for i in range(len(sample_list)):
  mask = load_image(f'{prefix}/Data/{sample_list[i]}.png') > 0
  cnts = load_pickle(f'{prefix}/Data/{sample_list[i]}.{gene_list[i]}.pickle')
  cnts[~mask] = np.nan
  plot_super(cnts, f'{prefix}/Results/Figure 1E.{sample_list[i]}.{gene_list[i]}.png')


# Figure 4F
prefix = os.getcwd()

for s in ['P24_AAH','P24_LUAD']:
   for gn in ['IL1R1','IL1B']:
      mask = load_image(f'{prefix}/Data/{s}.png') > 0
      cnts = load_pickle(f'{prefix}/Data/{s}.{gn}.pickle')
      cnts[~mask] = np.nan
      plot_super(cnts, f'{prefix}/Results/Figure 4F.{s}.{gn}.png')

