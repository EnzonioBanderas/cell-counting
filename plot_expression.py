import seaborn as sns
import os
from glob import glob
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from statsmodels.stats.oneway import anova_oneway
nIterBootstrap = 10000
from statsmodels.stats.multitest import multipletests
from scipy.stats import linregress



# Define
data_folder = 'Immunostainings_Pax5_mergedSelectionMaria' #################################################################################################
expression_selection_acronym_list = ['SNr', 'SNc', 'VTA', 'PBG', 'MRN', 'PB', 'PAG', 'CUN', 'PPN', 'IPN']
expression_selection_name_list = ['substantia nigra pars reticula', 'substantia nigra pars compacta',
                                     'ventral tegmental area', 'PBG', 'midbrain reticular formation',
                                     'parabrachial complex', 'periaqueductal gray', 'CUN',
                                     'pendunculopontine nucleus', 'IPN'] # PBG, CUN and IPN are not in human atlases. Ventral tegmental area both in subcortical and AAN atlases. 'parabrachial pigmented nucleus' likely a duplicate of 'parabrachial complex'.
expression_selection_sharptrack_acronym_list = ['SNr', 'SNc',
                                     'VTA', 'IPN', 'PAG', 'MRN',
                                     'PPN', 'PBG', 'CUN',
                                     'PB']



## Concatenate n_slice tables per subject into one n_slice table containing all data
n_slice_table_path_list = glob(os.path.join(data_folder, '*', 'processed', 'n_slice_table.csv'))
n_slice_table_list = list()
for pername_table_path in n_slice_table_path_list:
    n_slice_table = pd.read_csv(pername_table_path)

    genotype = pername_table_path.split(os.sep)[1].split('_')[1]
    subject = pername_table_path.split(os.sep)[1].split('_')[2]

    n_slice_table['genotype'] = genotype
    n_slice_table['subject'] = subject

    n_slice_table_list.append(n_slice_table)
n_slice_table_all = pd.concat(n_slice_table_list, ignore_index=True)
n_slice_table_all.to_csv('n_slice_table_all.csv')
n_slice_table_all = n_slice_table_all[np.isin(n_slice_table_all['acronym'], expression_selection_sharptrack_acronym_list)]
n_slice_table_all.to_csv('n_slice_table_expSel.csv')

## Concatenate expression tables per subject into one expression table containing all data
expression_table_path_list = glob(os.path.join(data_folder, '*', 'processed', 'pername_table.csv'))
expression_table_list = list()
for pername_table_path in expression_table_path_list:
    expression_table = pd.read_csv(pername_table_path)

    genotype = pername_table_path.split(os.sep)[1].split('_')[1]
    subject = pername_table_path.split(os.sep)[1].split('_')[2]

    expression_table['genotype'] = genotype
    expression_table['subject'] = subject

    expression_table_list.append(expression_table)
expression_table_all = pd.concat(expression_table_list, ignore_index=True)

## Adjust expression_table_all
# Calculate approximate roi_count per mm^2, assuming uniform resolution of 10 micrometer/pixel
expression_table_all['roi_count_permm'] = expression_table_all['roi_count_perpixel'] * 10000



customHueOrder = ['2A#7', '2B#8', '2B#10',
             '2A#2', '2A#5', '2B#6',
             '2A#9', '2B#1', '2B#3']
# Figure 11 mouse new significant
customPalette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
    (1.0, 0.4980392156862745, 0.054901960784313725),
    (1.0, 0.4980392156862745, 0.054901960784313725),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
# fig11 = plt.figure()
sns.set_theme(style="whitegrid")
g = sns.catplot(
    data=n_slice_table_all, kind="bar",
    x="acronym", y="n_slice", hue="subject",
    alpha=.6, height=6,
    order=expression_selection_sharptrack_acronym_list,
    hue_order=customHueOrder,
    palette=customPalette
)
g.despine(left=True)
g.set_axis_labels("", "n_slice")
g.legend.set_title("")
plt.savefig('n_slice_coloredByGenotype.png')

# Figure 12 mouse new significant
# fig12 = plt.figure()
sns.set_theme(style="whitegrid")
g = sns.catplot(
    data=n_slice_table_all, kind="bar",
    x="subject", y="n_slice", hue="acronym",
    palette="dark", alpha=.6, height=6,
    hue_order = expression_selection_sharptrack_acronym_list
)
g.despine(left=True)
g.set_axis_labels("", "n_slice")
g.legend.set_title("")
plt.savefig('n_slice_groupedByMouse.png')

# Figure 13 mouse new significant
# fig13 = plt.figure()
sns.set_theme(style="whitegrid")
g = sns.catplot(
    data=n_slice_table_all, kind="bar",
    x="acronym", y="n_slice", hue="subject",
    palette="dark", alpha=.6, height=6,
    order=expression_selection_sharptrack_acronym_list,
    hue_order=customHueOrder
)
g.despine(left=True)
g.set_axis_labels("", "n_slice")
g.legend.set_title("")
plt.savefig('n_slice.png')

# Figure 14 mouse new significant
n_slice_table_all_merged = pd.merge(left=n_slice_table_all, right=expression_table_all,
         left_on=['subject', 'name'], right_on=['subject', 'name'])
# fig14 = plt.figure()
sns.set_theme(style="whitegrid")
g = sns.catplot(
    data=n_slice_table_all_merged, kind="bar",
    x="acronym_x", y="roi_count_permm", hue="subject",
    palette="dark", alpha=.6, height=6,
    order=expression_selection_sharptrack_acronym_list,
    hue_order=customHueOrder
)
g.despine(left=True)
g.set_axis_labels("", "roi_count_permm")
g.legend.set_title("")
plt.savefig('roi_count_permm.png')