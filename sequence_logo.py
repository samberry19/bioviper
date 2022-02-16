import numpy as np
import matplotlib.pyplot as plt
import logomaker as lm

def sequence_logo(aln, s, color_scheme='chemistry', vpad=.025, width=.9, figsize=(15,2),font='Arial',
                    spacing=1, anchor=25, rotation=90, fmt='%d', fontsize=12):

    counts_mat = lm.alignment_to_matrix([''.join(np.array(aln[k])[s]).upper() for k in range(len(aln))], to_type='counts')
    prob_mat = lm.transform_matrix(counts_mat, from_type='counts', to_type='probability', pseudocount=0.01)
    info_mat = lm.transform_matrix(prob_mat, from_type='probability', to_type='information')

    lm_logo=lm.Logo(info_mat, figsize=figsize, color_scheme=color_scheme, vpad=vpad, width=width, font_name=font)

    # style using Axes methods
    lm_logo.ax.set_ylabel('information (bits)')

    lm_logo.style_xticks(spacing=spacing, anchor=anchor, rotation=rotation, fmt=fmt, fontsize=fontsize)

    return lm_logo
