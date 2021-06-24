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
data_folder = 'Immunostainings_Pax5_mergedSelectionMaria'
expression_selection_acronym_list = ['SNr', 'SNc', 'VTA', 'PBG', 'MRN', 'PB', 'PAG', 'CUN', 'PPN', 'IPN']
expression_selection_name_list = ['substantia nigra pars reticula', 'substantia nigra pars compacta',
                                     'ventral tegmental area', 'PBG', 'midbrain reticular formation',
                                     'parabrachial complex', 'periaqueductal gray', 'CUN',
                                     'pendunculopontine nucleus', 'IPN'] # PBG, CUN and IPN are not in human atlases. Ventral tegmental area both in subcortical and AAN atlases. 'parabrachial pigmented nucleus' likely a duplicate of 'parabrachial complex'.
expression_selection_sharptrack_name_list = ['Substantia nigra reticular part', 'Substantia nigra compact part',
                                     'Ventral tegmental area', 'Parabigeminal nucleus', 'Midbrain reticular nucleus',
                                     'Parabrachial nucleus', 'Periaqueductal gray', 'Cuneiform nucleus',
                                     'Pedunculopontine nucleus', 'Interpeduncular nucleus']
expression_selection_sharptrack_acronym_list = ['SNr', 'SNc',
                                     'VTA', 'IPN', 'PAG', 'MRN',
                                     'PPN', 'PBG', 'CUN',
                                     'PB']

# Load MRI volume tables
mouse_volume_table = pd.read_csv(
    os.path.join('..', 'mep-scripts', 'Data', 'Mouse', 'Analysis', 'all_volumes_mouse.csv'))
mouse_volume_pername_table = pd.read_csv(
    os.path.join('..', 'mep-scripts', 'Data', 'Mouse', 'Analysis', 'pername_volumes_mouse.csv'))
mouse_structure_table = pd.read_csv(
    os.path.join('..', 'mep-scripts', 'Data', 'Mouse', 'Reference', 'structure_graph_mc.csv'))
human_volume_pername_table = pd.read_csv(
    os.path.join('..', 'mep-scripts', 'Data', 'Human', 'Analysis', 'pername_volumes_human.csv'))
human_volume_table = pd.read_csv(
    os.path.join('..', 'mep-scripts', 'Data', 'Human', 'Analysis', 'all_volumes_human.csv'))

# Adjust tables for combining with expression tables
mouse_volume_table['genotype'] = 'Pax5++'
mouse_volume_table.loc[mouse_volume_table['Genotype'] == 'KO', 'genotype'] = 'Pax5-pR31Q-'
mouse_volume_table['VolumeRootNormalized'] = mouse_volume_table['Volume'] / mouse_volume_table['rootVolume']
mouse_volume_pername_table = pd.merge(left=mouse_volume_pername_table.reset_index(),
                                      right=mouse_volume_table.reset_index()[['name', 'acronym']].drop_duplicates(),
                                      left_on='name',
                                      right_on='name')

# Define other selections
significant_mouse_volume_selection = np.array(mouse_volume_pername_table.loc[mouse_volume_pername_table['pValFDR_BrainNormalized']<0.05, 'acronym'])

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




## Concatenate expression tables per subject into one expression-section table containing all data
expression_table_path_list = glob(os.path.join(data_folder, '*', 'processed', 'pernameSection_table.csv'))
expression_table_list = list()
for pername_table_path in expression_table_path_list:
    expression_table = pd.read_csv(pername_table_path)

    genotype = pername_table_path.split(os.sep)[1].split('_')[1]
    subject = pername_table_path.split(os.sep)[1].split('_')[2]

    expression_table['genotype'] = genotype
    expression_table['subject'] = subject

    expression_table_list.append(expression_table)
expressionSection_table_all = pd.concat(expression_table_list, ignore_index=True)

## Add IPN info by combining lower level IPN structure information
expression_table_IPN = expressionSection_table_all[[namestr.__contains__('Interpeduncular nucleus') for namestr in expressionSection_table_all['name']]]
expression_table_IPN = expression_table_IPN[expression_table_IPN['acronym'] != 'IPN']
IPN_roi_count = expression_table_IPN.groupby(['subject', 'section'])['roi_count'].sum()
IPN_pixel_count = expression_table_IPN.groupby(['subject', 'section'])['pixel_count'].sum()
IPN_roi_count_perpixel = IPN_roi_count / IPN_pixel_count
IPN_signal_sum =  expression_table_IPN.groupby(['subject', 'section'])['signal_sum'].sum()
IPN_signal_mean = IPN_signal_sum / IPN_pixel_count

# Loop through subjects and assign calculated IPN values to expression_table_all
acronym_logical = expressionSection_table_all['acronym'] == 'IPN'
for uSubject in np.unique(expressionSection_table_all['subject']):
    for uSection in np.unique(IPN_roi_count[uSubject].reset_index()['section']):
        subject_logical = expressionSection_table_all['subject'] == uSubject
        section_logical = expressionSection_table_all['section'] == uSection
        section_subject_logical = np.logical_and(subject_logical, section_logical)
        expressionSection_table_all.loc[np.logical_and(acronym_logical, section_subject_logical), 'roi_count'] = IPN_roi_count[uSubject][uSection]
        expressionSection_table_all.loc[np.logical_and(acronym_logical, section_subject_logical), 'pixel_count'] = IPN_pixel_count[uSubject][uSection]
        expressionSection_table_all.loc[np.logical_and(acronym_logical, section_subject_logical), 'roi_count_perpixel'] = IPN_roi_count_perpixel[uSubject][uSection]
        expressionSection_table_all.loc[np.logical_and(acronym_logical, section_subject_logical), 'signal_sum'] = IPN_signal_sum[uSubject][uSection]
        expressionSection_table_all.loc[np.logical_and(acronym_logical, section_subject_logical), 'signal_mean'] = IPN_signal_mean[uSubject][uSection]

## Adjust expression_table_all
# Calculate approximate roi_count per mm^2, assuming uniform resolution of 10 micrometer/pixel
expressionSection_table_all['roi_count_permm'] = expressionSection_table_all['roi_count_perpixel'] * 10000

#

# Normalize roi count with total volume per mouse and structure
# expressionSection_table_all = expressionSection_table_all[expressionSection_table_all['roi_count'] != 0]
expressionSection_table_all = expressionSection_table_all.reset_index()
expressionSection_table_all['voxel_count'] = expressionSection_table_all['pixel_count'] * 40
expressionSection_table_all['voxel_volume'] = expressionSection_table_all['voxel_count'] / 4000 # per square micrometer
expressionSection_table_all = expressionSection_table_all.join(expressionSection_table_all.groupby(['subject', 'name'])['voxel_volume'].sum(), on=['subject', 'name'], rsuffix='_sum')
expressionSection_table_all['roi_count_permum3'] = expressionSection_table_all['roi_count'] / expressionSection_table_all['voxel_volume_sum']
expressionSection_table_all.sort_values(by='roi_count_permum3', ascending=False)

# Sum expressionSection_table_all over sections and compute measures
expressionMouse_table_all = expressionSection_table_all.groupby(['subject', 'genotype', 'name', 'acronym'])[['roi_count', 'voxel_volume']].sum().reset_index()
expressionMouse_table_all['roi_count_permum3'] = expressionMouse_table_all['roi_count'] / expressionMouse_table_all['voxel_volume']

# Write expressionSection_table_all to csv
expressionSection_table_all.to_csv(os.path.join('pernameSection_Pax5_counting.csv'))



## Concatenate roi tables per subject into one roi table containing all data
for data_string in [data_folder]:
    roi_table_path_list = glob(os.path.join(data_string, '*', 'processed', 'roi_table_all.csv'))
    roi_table_list = list()
    for roi_table_path in roi_table_path_list:
        roi_table = pd.read_csv(roi_table_path)

        genotype = roi_table_path.split(os.sep)[1].split('_')[1]
        subject = roi_table_path.split(os.sep)[1].split('_')[2]

        roi_table['genotype'] = genotype
        roi_table['subject'] = subject

        roi_table_list.append(roi_table)
    roi_table_all = pd.concat(roi_table_list, ignore_index=True)
    if data_string == 'Immunostainings_Pax5':
        roi_table_all.groupby('subject')['filename'].nunique().reset_index().to_csv('slice_count_persubject_table.csv')
    else:
        roi_table_all.groupby('subject')['filename'].nunique().reset_index().to_csv('slice_count_persubject_extra_table.csv')

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


## Add IPN info by combining lower level IPN structure information
expression_table_IPN = expression_table_all[[namestr.__contains__('Interpeduncular nucleus') for namestr in expression_table_all['name']]]
expression_table_IPN = expression_table_IPN[expression_table_IPN['acronym'] != 'IPN']
IPN_roi_count = expression_table_IPN.groupby('subject')['roi_count'].sum()
IPN_pixel_count = expression_table_IPN.groupby('subject')['pixel_count'].sum()
IPN_roi_count_perpixel = IPN_roi_count / IPN_pixel_count
IPN_signal_sum =  expression_table_IPN.groupby('subject')['signal_sum'].sum()
IPN_signal_mean = IPN_signal_sum / IPN_pixel_count

# Loop through subjects and assign calculated IPN values to expression_table_all
acronym_logical = expression_table_all['acronym'] == 'IPN'
for uSubject in np.unique(expression_table_all['subject']):
    subject_logical = expression_table_all['subject'] == uSubject
    expression_table_all.loc[np.logical_and(acronym_logical, subject_logical), 'roi_count'] = IPN_roi_count[uSubject]
    expression_table_all.loc[np.logical_and(acronym_logical, subject_logical), 'pixel_count'] = IPN_pixel_count[uSubject]
    expression_table_all.loc[np.logical_and(acronym_logical, subject_logical), 'roi_count_perpixel'] = IPN_roi_count_perpixel[uSubject]
    expression_table_all.loc[np.logical_and(acronym_logical, subject_logical), 'signal_sum'] = IPN_signal_sum[uSubject]
    expression_table_all.loc[np.logical_and(acronym_logical, subject_logical), 'signal_mean'] = IPN_signal_mean[uSubject]



## Adjust expression_table_all
# Calculate approximate roi_count per mm^2, assuming uniform resolution of 10 micrometer/pixel
expression_table_all['roi_count_permm'] = expression_table_all['roi_count_perpixel'] * 10000



# # Get top most expressive structures
# pername_series_meanroicountpermm = expression_table_all.groupby(['acronym', 'name'])['roi_count_permm'].mean()
# acronym_top_list = list(pername_series_meanroicountpermm.sort_values()[-1:-20:-1].index)




# Normalize roi count with total volume per mouse and structure
expression_table_all['voxel_count'] = 0
for iRow in range(expression_table_all.shape[0]):
    expression_table_all.loc[iRow, 'voxel_count'] = expression_table_all['pixel_count'].iloc[iRow] * 40
expression_table_all['voxel_volume'] = expression_table_all['voxel_count'] / 4000 # per square micrometer
expression_table_all['roi_count_permum3'] = expression_table_all['roi_count'] / expression_table_all['voxel_volume']
expression_table_all.sort_values(by='roi_count_permum3', ascending=False)

# Write expression_table_all to csv
expression_table_all.to_csv(os.path.join('pername_Pax5_counting.csv'))



# Selection
expression_table_selection = expression_table_all[np.isin(expression_table_all['acronym'], expression_selection_acronym_list)]
expression_table_selection.to_csv(os.path.join('pername_Pax5_counting_selection.csv'))
mouse_volume_table_selection = mouse_volume_table[np.isin(mouse_volume_table['acronym'], expression_selection_acronym_list)]
mouse_volume_table_selection.to_csv(os.path.join('volume_table_expSel.csv'))
expressionSection_table_selection = expressionSection_table_all[np.isin(expressionSection_table_all['acronym'], expression_selection_acronym_list)]
expressionSection_table_selection.to_csv(os.path.join('pernameSection_Pax5_counting_selection.csv'))
expressionMouseSelection_table_all = expressionMouse_table_all[np.isin(expressionMouse_table_all['acronym'], expression_selection_acronym_list)]
expressionMouseSelection_table_all.to_csv(os.path.join('pernameMouseSelection_Pax5_counting_selection.csv'))



## Plotting
# Figure 1: Expression of expSel structures
fig1 = plt.figure()
customPalette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
ax1 = sns.boxplot(x="acronym",
                  y="roi_count_permum3",
                  hue="genotype",
                  data=expression_table_selection,
                  hue_order=['Pax5++', 'Pax5+-', 'Pax5-pR31Q-'],
                  order=expression_selection_acronym_list,
                  palette=customPalette)
# sns.swarmplot(x="acronym",
#               y="roi_count_permum3",
#               hue="genotype",
#               data=expression_table_selection,
#               hue_order=['Pax5++', 'Pax5+-', 'Pax5-pR31Q-'],
#               order=expression_selection_acronym_list,
#               color='.25',
#               dodge=0.4)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=60)
ax1.set(xlabel='structure acronym', ylabel='cell count per mum^3')
handles, labels = ax1.get_legend_handles_labels()
l = plt.legend(handles[0:3], labels[0:3])
ax1.yaxis.grid(False) # Hide the horizontal gridlines
ax1.xaxis.grid(True) # Show the vertical gridlines
plt.show()
plt.savefig('histo_selection_boxplot_mouse.png')
expression_table_selection.to_csv('histo_selection_boxplot_mouse.csv')



# Figure 2: Mouse Volumes of expSel structures
fig2 = plt.figure()
customPalette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
ax2 = sns.barplot(x="acronym",
                  y="VolumeNormalized",
                  hue="genotype",
                  data=mouse_volume_table_selection,
                  hue_order=['Pax5++', 'Pax5-pR31Q-'],
                  order=expression_selection_acronym_list,
                  ci=68,
                  palette=customPalette)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=60)
ax2.set(xlabel='structure acronym', ylabel='Normalized volume')
plt.show()



# Figure 3 mouse old significant
fig3 = plt.figure()
mri_name_list = ['Substantia nigra, reticular part', 'Substantia nigra, compact part', 'Lobules IV-V',
                 'Midbrain reticular nucleus', 'Folium-tuber vermis (VII)']
mri_sig_acronym_list = ['SNr', 'SNc', 'MRN', 'CUL4, 5', 'FOTU']
expression_table_selection = expression_table_all[np.isin(expression_table_all['name'], mri_name_list)]
mouse_volume_table_selection = mouse_volume_table[np.isin(mouse_volume_table['name'], mri_name_list)]

# Calculate effect size for names where there is at least one sample for KO and WT
# S_p = np.sqrt(((nWT - 1) * np.power(std_WT, 2) + (nKO - 1) * np.power(std_KO, 2)) / (N - 2))
# S_p_AN = np.sqrt(((nWT - 1) * np.power(std_WT_AN, 2) + (nKO - 1) * np.power(std_KO_AN, 2)) / (N - 2))
# S_p_RN = np.sqrt(((nWT - 1) * np.power(std_WT_RN, 2) + (nKO - 1) * np.power(std_KO_RN, 2)) / (N - 2))
# cohenD = (mean_WT - mean_KO) / S_p

customPalette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
# sns.set_theme(style="whitegrid")
ax3 = sns.boxplot(x="acronym",
                  y="VolumeNormalized",
                  hue="genotype",
                  data=mouse_volume_table_selection,
                  hue_order=['Pax5++', 'Pax5-pR31Q-'],
                  order=mri_sig_acronym_list,
                  palette=customPalette)
sns.swarmplot(x="acronym",
              y="VolumeNormalized",
              hue="genotype",
              data=mouse_volume_table_selection,
              hue_order=['Pax5++', 'Pax5-pR31Q-'],
              order=mri_sig_acronym_list,
              dodge=0.4,
              color='.25')
ax3.set_xticklabels(ax3.get_xticklabels(), rotation=60)
ax3.set(xlabel='structure acronym', ylabel='Normalized volume')
ax3.yaxis.grid(False) # Hide the horizontal gridlines
ax3.xaxis.grid(True) # Show the vertical gridlines
plt.show()
handles, labels = ax3.get_legend_handles_labels()
l = plt.legend(handles[0:2], labels[0:2])
plt.savefig('mri_sig_boxplot_mouse.png')
mouse_volume_table_selection.to_csv('mri_sig_boxplot_mouse.csv')



fig4 = plt.figure()
name_selection_list = ['substantia nigra pars reticula',
                       'substantia nigra pars compacta',
                       'midbrain reticular formation',
                       'Left_I_IV', 'Right_I_IV',
                       'Left_V', 'Right_V',
                       'Vermis_VIIb']
customPalette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
customPalette2 = [
    (0.25, 0.25, 0.25),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
human_volume_table_selection = human_volume_table[np.isin(human_volume_table['name'], name_selection_list)]
human_volume_table_selection['genotype'] = 'control'
human_volume_table_selection.loc[human_volume_table_selection['subject'] == 'patient', 'genotype'] = 'patient'
ax4 = sns.boxplot(x="name",
                  y="VolumeNormalized",
                  hue="genotype",
                  data=human_volume_table_selection,
                  hue_order=['control', 'patient'],
                  order=name_selection_list,
                  palette=customPalette)
ax4 = sns.swarmplot(x="name",
                  y="VolumeNormalized",
                  hue="genotype",
                  data=human_volume_table_selection,
                  hue_order=['control', 'patient'],
                  order=name_selection_list,
                  palette=customPalette2,
                    dodge=0.4)
ax4.set_xticklabels(ax4.get_xticklabels(), rotation=60)
ax4.set(xlabel='structure name', ylabel='Normalized volume')
handles, labels = ax4.get_legend_handles_labels()
l = plt.legend(handles[0:2], ['controls', 'patient'])
ax4.set_xticklabels(['SNr', 'SNc', 'MRN', 'Left_I_IV', 'Right_I_IV', 'Left_V', 'Right_V', 'Vermis_VII'])
plt.show()
plt.savefig('mri_sig_boxplot_human.png')
human_volume_table_selection.to_csv('mri_sig_boxplot_human.csv')




# Human SUIT volume overview plot
fig5 = plt.figure()
customPalette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
customPalette2 = [
    (0.25, 0.25, 0.25),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
human_volume_table_selection = human_volume_table[np.isin(human_volume_table['atlas'], ['Lobules-SUIT'])]
human_volume_table_selection['genotype'] = 'control'
human_volume_table_selection.loc[human_volume_table_selection['subject'] == 'patient', 'genotype'] = 'patient'
plotAx_human_SUIT = sns.boxplot(x="name",
                                y="VolumeNormalized",
                                hue="genotype",
                                data=human_volume_table_selection,
                                hue_order=['control', 'patient'],
                                palette=customPalette)
plotAx_human_SUIT = sns.swarmplot(x="name",
                                  y="VolumeNormalized",
                                  hue="genotype",
                                  data=human_volume_table_selection,
                                  hue_order=['control', 'patient'],
                                  palette=customPalette2,
                                  dodge=0.4)
# ax5.set_xticklabels(ax5.get_xticklabels(), rotation=60)
plotAx_human_SUIT.set(xlabel='structure name', ylabel='Normalized volume')
handles, labels = plotAx_human_SUIT.get_legend_handles_labels()
l = plt.legend(handles[0:2], ['controls', 'patient'])
# ax4.set_xticklabels(['SNr', 'SNc', 'MRN', 'Left_I_IV', 'Right_I_IV', 'Left_V', 'Right_V', 'Vermis_VII'])
mng = plt.get_current_fig_manager()
# mng.window.state('zoomed') #works fine on Windows!
# mng.frame.Maximize(True)
# mng.window.showMaximized
plotAx_human_SUIT.set_xticklabels(['L_I-IV', 'R_I-IV', 'L_V', 'R_V',
                     'L_VI', 'V_VI', 'R_VI',
                     'L_CrusI', 'V_CrusI', 'R_CrusI',
                     'L_CrusII', 'V_CrusII', 'R_CrusII',
                     'L_VIIb', 'V_VIIb', 'R_VIIb',
                     'L_VIIIa', 'V_VIIIa', 'R_VIIIa',
                     'L_VIIIb', 'V_VIIIb', 'R_VIIIb',
                     'L_IX', 'V_IX', 'R_IX',
                     'L_X', 'V_X', 'R_X',
                     'L_Dentate', 'R_Dentate', 'L_Interposed', 'R_Interposed',
                     'L_Fastigial', 'R_Fastigial',
                     'CE', 'L_CE', 'R_CE', 'PO_CE', 'AN_CE',
                     'L_PO_CE', 'R_PO_CE', 'L_AN_CE', 'R_AN_CE',
                     'V_CE'])
plotAx_human_SUIT.set_xticklabels(plotAx_human_SUIT.get_xticklabels(), rotation=60)
plt.show()
mng.full_screen_toggle()
plt.savefig('mri_CE_boxplot_human.png')
human_volume_table_selection.to_csv('mri_CE_boxplot_human.csv')



# Human subcortical volume overview plot
fig6 = plt.figure()
customPalette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
customPalette2 = [
    (0.25, 0.25, 0.25),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
human_volume_table_selection = human_volume_table[np.isin(human_volume_table['atlas'], ['subcortical'])]
human_volume_table_selection['genotype'] = 'control'
human_volume_table_selection.loc[human_volume_table_selection['subject'] == 'patient', 'genotype'] = 'patient'
plotAx_human_subcortical = sns.boxplot(x="name",
                                y="VolumeNormalized",
                                hue="genotype",
                                data=human_volume_table_selection,
                                hue_order=['control', 'patient'],
                                palette=customPalette)
plotAx_human_subcortical = sns.swarmplot(x="name",
                                  y="VolumeNormalized",
                                  hue="genotype",
                                  data=human_volume_table_selection,
                                  palette=customPalette2,
                                  dodge=0.4)
plotAx_human_subcortical.set(xlabel='structure name', ylabel='Normalized volume')
handles, labels = plotAx_human_subcortical.get_legend_handles_labels()
l = plt.legend(handles[0:2], ['controls', 'patient'])
# ax5.set_xticklabels(['SNr', 'SNc', 'MRN', 'Left_I_IV', 'Right_I_IV', 'Left_V', 'Right_V', 'Vermis_VII'])
mng = plt.get_current_fig_manager()
# mng.window.state('zoomed') #works fine on Windows!
# mng.frame.Maximize(True)
# mng.window.showMaximized
plotAx_human_subcortical.set_xticklabels(['PU', 'CU', 'NA', 'EA',
                     'GPext', 'GPint', 'SNc',
                     'RN', 'SNr', 'PN',
                     'VTA', 'VP', 'HN',
                     'Hypothal', 'MN', 'SubthalN'])
plotAx_human_subcortical.set_xticklabels(plotAx_human_subcortical.get_xticklabels(), rotation=60)
plt.show()
mng.full_screen_toggle()
plt.savefig('mri_subcortical_boxplot_human.png')
human_volume_table_selection.to_csv('mri_subcortical_boxplot_human.csv')



# Human CerebrA-mask volume overview plot
fig7 = plt.figure()
customPalette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
customPalette2 = [
    (0.25, 0.25, 0.25),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
human_volume_table_selection = human_volume_table[np.isin(human_volume_table['atlas'], ['CerebrA', 'mask'])] ###
human_volume_table_selection['genotype'] = 'control'
human_volume_table_selection.loc[human_volume_table_selection['subject'] == 'patient', 'genotype'] = 'patient'
plotAx_human_CerebrA = sns.boxplot(x="name",
                                y="VolumeNormalized",
                                hue="genotype",
                                data=human_volume_table_selection,
                                hue_order=['control', 'patient'],
                                palette=customPalette) ###
plotAx_human_CerebrA = sns.swarmplot(x="name",
                                  y="VolumeNormalized",
                                  hue="genotype",
                                  data=human_volume_table_selection,
                                  palette=customPalette2,
                                  dodge=0.4) ###
plotAx_human_CerebrA.set(xlabel='structure name', ylabel='Normalized volume') ###
handles, labels = plotAx_human_CerebrA.get_legend_handles_labels() ###
l = plt.legend(handles[0:2], ['controls', 'patient'])
# ax5.set_xticklabels(['SNr', 'SNc', 'MRN', 'Left_I_IV', 'Right_I_IV', 'Left_V', 'Right_V', 'Vermis_VII'])
mng = plt.get_current_fig_manager()
# mng.window.state('zoomed') #works fine on Windows!
# mng.frame.Maximize(True)
# mng.window.showMaximized
# plotAx_human_CerebrA.set_xticklabels(['PU', 'CU', 'NA', 'EA',
#                      'GPext', 'GPint', 'SNc',
#                      'RN', 'SNr', 'PN',
#                      'VTA', 'VP', 'HN',
#                      'Hypothal', 'MN', 'SubthalN']) ###
plotAx_human_CerebrA.set_xticklabels(plotAx_human_CerebrA.get_xticklabels(), rotation=60) ###
plt.show()
mng.full_screen_toggle()
plt.savefig('mri_CerebrA_boxplot_human.png') ###
human_volume_table_selection.to_csv('mri_CerebrA_boxplot_human.csv') ###



# Human AAN volume overview plot
fig8 = plt.figure()
customPalette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
customPalette2 = [
    (0.25, 0.25, 0.25),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
human_volume_table_selection = human_volume_table[np.isin(human_volume_table['atlas'], ['AAN'])] ###
human_volume_table_selection['genotype'] = 'control'
human_volume_table_selection.loc[human_volume_table_selection['subject'] == 'patient', 'genotype'] = 'patient'
plotAx_human_AAN = sns.boxplot(x="name",
                                y="VolumeNormalized",
                                hue="genotype",
                                data=human_volume_table_selection,
                                hue_order=['control', 'patient'],
                                palette=customPalette) ###
plotAx_human_AAN = sns.swarmplot(x="name",
                                  y="VolumeNormalized",
                                  hue="genotype",
                                  data=human_volume_table_selection,
                                  palette=customPalette2,
                                  dodge=0.4) ###
plotAx_human_AAN.set(xlabel='structure name', ylabel='Normalized volume') ###
handles, labels = plotAx_human_AAN.get_legend_handles_labels() ###
l = plt.legend(handles[0:2], ['controls', 'patient'])
# ax5.set_xticklabels(['SNr', 'SNc', 'MRN', 'Left_I_IV', 'Right_I_IV', 'Left_V', 'Right_V', 'Vermis_VII'])
mng = plt.get_current_fig_manager()
# mng.window.state('zoomed') #works fine on Windows!
# mng.frame.Maximize(True)
# mng.window.showMaximized
# plotAx_human_AAN.set_xticklabels(['DR', 'LC', 'MRN', 'MR',
#                      'PAG', 'PB', 'PO', 'PPN', 'VTA']) ######
plotAx_human_AAN.set_xticklabels(plotAx_human_AAN.get_xticklabels(), rotation=60) ###
plt.show()
mng.full_screen_toggle()
plt.savefig('mri_AAN_boxplot_human.png') ###
human_volume_table_selection.to_csv('mri_AAN_boxplot_human.csv') ###



# Human expSel volume overview plot
fig9 = plt.figure()
expression_human_selection_name_list = expression_selection_name_list.copy()
expression_human_selection_name_list.append('parabrachial pigmented nucleus')
customPalette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
customPalette2 = [
    (0.25, 0.25, 0.25),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
human_volume_table_selection = human_volume_table[np.isin(human_volume_table['name'], expression_human_selection_name_list)] ###
human_volume_table_selection['genotype'] = 'control'
human_volume_table_selection.loc[human_volume_table_selection['subject'] == 'patient', 'genotype'] = 'patient'
VTA_AAN_logical = np.logical_and(human_volume_table_selection['name'] == 'ventral tegmental area', human_volume_table_selection['atlas'] == 'AAN')
human_volume_table_selection.loc[VTA_AAN_logical, 'name'] = 'ventral tegmental area AAN'
plotAx_human_expSel = sns.boxplot(x="name",
                                y="VolumeNormalized",
                                hue="genotype",
                                data=human_volume_table_selection,
                                hue_order=['control', 'patient'],
                                palette=customPalette) ###
plotAx_human_expSel = sns.swarmplot(x="name",
                                  y="VolumeNormalized",
                                  hue="genotype",
                                  data=human_volume_table_selection,
                                  palette=customPalette2,
                                  dodge=0.4) ###
plotAx_human_expSel.set(xlabel='structure name', ylabel='Normalized volume') ###
handles, labels = plotAx_human_expSel.get_legend_handles_labels() ###
l = plt.legend(handles[0:2], ['controls', 'patient'])
# ax5.set_xticklabels(['SNr', 'SNc', 'MRN', 'Left_I_IV', 'Right_I_IV', 'Left_V', 'Right_V', 'Vermis_VII'])
mng = plt.get_current_fig_manager()
# mng.window.state('zoomed') #works fine on Windows!
# mng.frame.Maximize(True)
# mng.window.showMaximized
plotAx_human_expSel.set_xticklabels(['SNc', 'SNr', 'PB', 'VTA',
                     'MRN_AAN', 'PAG_AAN', 'PB_AAN', 'PPN_AAN', 'VTA_AAN']) ######
plotAx_human_expSel.set_xticklabels(plotAx_human_expSel.get_xticklabels(), rotation=60) ###
plt.show()
mng.full_screen_toggle()
plt.savefig('mri_expSel_boxplot_human.png') ###
human_volume_table_selection.to_csv('mri_expSel_boxplot_human.csv') ###



# Figure 10 mouse new significant
mouse_volume_table_selection = mouse_volume_table[np.isin(mouse_volume_table['acronym'], significant_mouse_volume_selection)]
mouse_volume_table_selection.to_csv(os.path.join('significant_mouse_volume_table.csv'))
fig10 = plt.figure()

# Calculate effect size for names where there is at least one sample for KO and WT
# S_p = np.sqrt(((nWT - 1) * np.power(std_WT, 2) + (nKO - 1) * np.power(std_KO, 2)) / (N - 2))
# S_p_AN = np.sqrt(((nWT - 1) * np.power(std_WT_AN, 2) + (nKO - 1) * np.power(std_KO_AN, 2)) / (N - 2))
# S_p_RN = np.sqrt(((nWT - 1) * np.power(std_WT_RN, 2) + (nKO - 1) * np.power(std_KO_RN, 2)) / (N - 2))
# cohenD = (mean_WT - mean_KO) / S_p

customPalette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
# sns.set_theme(style="whitegrid")
ax10 = sns.boxplot(x="acronym",
                  y="VolumeNormalized",
                  hue="genotype",
                  data=mouse_volume_table_selection,
                  hue_order=['Pax5++', 'Pax5-pR31Q-'],
                  order=significant_mouse_volume_selection,
                  palette=customPalette)
sns.swarmplot(x="acronym",
              y="VolumeNormalized",
              hue="genotype",
              data=mouse_volume_table_selection,
              hue_order=['Pax5++', 'Pax5-pR31Q-'],
              order=significant_mouse_volume_selection,
              dodge=0.4,
              color='.25')
ax10.set_xticklabels(ax10.get_xticklabels(), rotation=60)
ax10.set(xlabel='structure acronym', ylabel='Normalized volume')
ax10.yaxis.grid(False) # Hide the horizontal gridlines
ax10.xaxis.grid(True) # Show the vertical gridlines
plt.show()
handles, labels = ax10.get_legend_handles_labels()
l = plt.legend(handles[0:2], labels[0:2])
plt.savefig('mri_sigNew_boxplot_mouse.png')


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
fig11 = plt.figure()
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
fig12 = plt.figure()
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
fig13 = plt.figure()
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
fig14 = plt.figure()
sns.set_theme(style="whitegrid")
g = sns.catplot(
    data=n_slice_table_all_merged, kind="bar",
    x="acronym_x", y="roi_count_permum3", hue="subject",
    palette="dark", alpha=.6, height=6,
    order=expression_selection_sharptrack_acronym_list,
    hue_order=customHueOrder
)
g.despine(left=True)
g.set_axis_labels("", "roi_count_permum3")
g.legend.set_title("")
plt.savefig('roi_count_permum3.png')

# Figure 15: per section boxplots
fig15 = plt.figure()
customPalette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
ax15 = sns.boxplot(x="acronym",
                  y="roi_count_permum3",
                  hue="genotype",
                  data=expressionSection_table_selection,
                  hue_order=['Pax5++', 'Pax5+-', 'Pax5-pR31Q-'],
                  order=expression_selection_acronym_list,
                  palette=customPalette)
# sns.swarmplot(x="acronym",
#               y="roi_count_permum3",
#               hue="genotype",
#               data=expressionSection_table_selection,
#               hue_order=['Pax5++', 'Pax5+-', 'Pax5-pR31Q-'],
#               order=expression_selection_acronym_list,
#               color='.25',
#               dodge=0.4)
ax15.set_xticklabels(ax15.get_xticklabels(), rotation=60)
ax15.set(xlabel='structure acronym', ylabel='cell count per mum^3')
handles, labels = ax15.get_legend_handles_labels()
l = plt.legend(handles[0:3], labels[0:3])
ax15.yaxis.grid(False) # Hide the horizontal gridlines
ax15.xaxis.grid(True) # Show the vertical gridlines
plt.show()
plt.savefig('histo_selectionSection_boxplot_mouse.png')
expressionSection_table_selection.to_csv('histo_selectionSection_boxplot_mouse.csv')

# Figure 16: per mouse boxplots?
fig16 = plt.figure()
customPalette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
    (1.0, 0.4980392156862745, 0.054901960784313725)
]
ax16 = sns.boxplot(x="acronym",
                  y="roi_count_permum3",
                  hue="genotype",
                  data=expressionMouseSelection_table_all,
                  hue_order=['Pax5++', 'Pax5+-', 'Pax5-pR31Q-'],
                  order=expression_selection_acronym_list,
                  palette=customPalette)
# sns.swarmplot(x="acronym",
#               y="roi_count_permum3",
#               hue="genotype",
#               data=expressionMouseSelection_table_all,
#               hue_order=['Pax5++', 'Pax5+-', 'Pax5-pR31Q-'],
#               order=expression_selection_acronym_list,
#               color='.25',
#               dodge=0.4)
ax16.set_xticklabels(ax16.get_xticklabels(), rotation=60)
ax16.set(xlabel='structure acronym', ylabel='cell count per mum^3')
handles, labels = ax16.get_legend_handles_labels()
l = plt.legend(handles[0:3], labels[0:3])
ax16.yaxis.grid(False) # Hide the horizontal gridlines
ax16.xaxis.grid(True) # Show the vertical gridlines
plt.show()
plt.savefig('histo_MouseSelection_boxplot_mouse.png')
expressionMouseSelection_table_all.to_csv('histo_MouseSelection_boxplot_mouse.csv')

# ################
# #                   data=human_volume_table_selection,
#                   hue_order=['control', 'patient'],
#                   palette=customPalette2,
#                     dodge=0.4)
# ax5.set_xticklabels(ax5.get_xticklabels(), rotation=60)
# ax5.set(xlabel='structure name', ylabel='Normalized volume')
# handles, labels = ax5.get_legend_handles_labels()
# l = plt.legend(handles[0:2], ['controls', 'patient'])
# # ax5.set_xticklabels(['SNr', 'SNc', 'MRN', 'Left_I_IV', 'Right_I_IV', 'Left_V', 'Right_V', 'Vermis_VII'])
# ax5.set_xticklabels(['L_I-IV', 'R_I-IV', 'L_V', 'R_V',
#                      'L_VI', 'V_VI', 'R_VI',
#                      'L_CrusI', 'V_CrusI', 'R_CrusI',
#                      'L_CrusII', 'V_CrusII', 'R_CrusII',
#                      'L_VIIb', 'V_VIIb', 'R_VIIb',
#                      'L_VIIIa', 'V_VIIIa', 'R_VIIIa',
#                      'L_VIIIb', 'V_VIIIb', 'R_VIIIb',
#                      'L_IX', 'V_IX', 'R_IX',
#                      'L_X', 'V_X', 'R_X',
#                      'L_Dentate', 'R_Dentate', 'L_Interposed', 'R_Interposed',
#                      'L_Fastigial', 'R_Fastigial',
#                      'CE', 'L_CE', 'R_CE', 'PO_CE', 'AN_CE',
#                      'L_PO_CE', 'R_PO_CE', 'L_AN_CE', 'R_AN_CE',
#                      'V_CE'])
# mng = plt.get_current_fig_manager()
# # mng.window.state('zoomed') #works fine on Windows!
# # mng.frame.Maximize(True)
# # mng.window.showMaximized
# plt.show()
# mng.full_screen_toggle()
# plt.savefig('mri_subcortical_boxplot_human.png')
# human_volume_table_selection.to_csv('mri_subcortical_boxplot_human.csv')
# ################




# # pername expression processing to through pername
# pername_pergen_table = pername_table_all.groupby(['genotype', 'name'])['roi_count_permm'].mean()
# pername_pergen_table = pername_pergen_table.reset_index().rename(columns={'roi_count_permm': 'roi_count_permm_mean'})
# # pername merging with mri data
# mouse_volume_pergen_table = mouse_volume_pername_table.groupby(['genotype', 'name'])['VolumeNormalized'].mean()
# mouse_volume_pergen_table = mouse_volume_pergen_table.reset_index().rename(columns={'VolumeNormalized': 'VolumeNormalized_mean'})
# mouse_volume_pergen_table = mouse_volume_pergen_table[np.logical_not(np.isnan(mouse_volume_pergen_table['VolumeNormalized_mean']))]
# pergen_merged_table = pd.merge(left=pername_pergen_table,
#                                right=mouse_volume_pergen_table,
#                                left_on=['genotype', 'name'],
#                                right_on=['genotype', 'name'])
#
# # scatter plot average for different genotypes
# mouse_volume_table[mouse_volume_table['name'] == 'Accessory abducens nucleus']

name_uniq = np.unique(np.array(expression_table_all['name'].astype('category')))
nName = len(name_uniq)
pername_expression_table_list = list()
for nameStruct in name_uniq:
    pername_table_nameStruct = expression_table_all[expression_table_all['name'] == nameStruct]
    acronymStruct = pername_table_nameStruct.iloc[0]['acronym']
    pername_table_nameStruct_WT = pername_table_nameStruct.loc[pername_table_nameStruct['genotype'] == 'Pax5++']
    pername_table_nameStruct_het = pername_table_nameStruct.loc[pername_table_nameStruct['genotype'] == 'Pax5+-']
    pername_table_nameStruct_KO = pername_table_nameStruct.loc[pername_table_nameStruct['genotype'] == 'Pax5-pR31Q-']
    stats_res = anova_oneway(data=[pername_table_nameStruct_WT['roi_count_permum3'],
                                   pername_table_nameStruct_het['roi_count_permum3'],
                                   pername_table_nameStruct_KO['roi_count_permum3']],
                             use_var='unequal')
    t_stat = stats_res[0]
    p_val = stats_res[1]

    mean_WT = np.mean(pername_table_nameStruct_WT['roi_count_permum3'])
    mean_het = np.mean(pername_table_nameStruct_het['roi_count_permum3'])
    mean_KO = np.mean(pername_table_nameStruct_KO['roi_count_permum3'])
    std_WT = np.std(pername_table_nameStruct_WT['roi_count_permum3'])
    std_het = np.std(pername_table_nameStruct_het['roi_count_permum3'])
    std_KO = np.std(pername_table_nameStruct_KO['roi_count_permum3'])
    nWT = len(pername_table_nameStruct_WT['roi_count_permum3'])
    nhet = len(pername_table_nameStruct_het['roi_count_permum3'])
    nKO = len(pername_table_nameStruct_KO['roi_count_permum3'])
    mean_count_WT = np.mean(pername_table_nameStruct_WT['roi_count'])
    mean_count_het = np.mean(pername_table_nameStruct_het['roi_count'])
    mean_count_KO = np.mean(pername_table_nameStruct_KO['roi_count'])
    N = nWT + nKO
    # print(f'nWT{nWT}, nKO={nKO}, N={N}')
    S_p = np.sqrt(((nWT - 1) * np.power(std_WT, 2) + (nKO - 1) * np.power(std_KO, 2)) / (N - 2))
    cohenD = (mean_WT - mean_KO) / S_p

    S_a = np.sqrt((np.power(std_WT, 2) + np.power(std_KO, 2)) / 2)
    hedge_correction = (N - 3) / (N - 2.25)
    cohenD_ac = ((mean_WT - mean_KO) / S_a) * hedge_correction
    # print(f'cohenD_RN = {cohenD_RN}, cohenD_RN_a = {cohenD_RN_a}, '
    #       f'cohenD_RN_ac = {cohenD_RN_ac}, cohenD_ac = {cohenD_ac}')

    cohenD_BS = np.empty(nIterBootstrap)
    cohenD_RN_BS = np.empty(nIterBootstrap)
    for iBS in range(nIterBootstrap):
        WT_BS = np.random.normal(mean_WT, std_WT, nWT)
        KO_BS = np.random.normal(mean_KO, std_KO, nKO)
        mean_WT_BS = np.mean(WT_BS)
        mean_KO_BS = np.mean(KO_BS)
        std_WT_BS = np.std(WT_BS)
        std_KO_BS = np.std(KO_BS)
        S_p = np.sqrt(((nWT - 1) * np.power(std_WT_BS, 2) + (nKO - 1) * np.power(std_KO_BS, 2)) / (nKO + nWT - 2))
        S_a = np.sqrt((np.power(std_WT, 2) + np.power(std_KO, 2)) / 2)
        cohenD_BS[iBS] = ((mean_WT_BS - mean_KO_BS) / S_a) * hedge_correction

    cohenD_ac_CI = [np.quantile(cohenD_BS, .025), np.quantile(cohenD_BS, .975)]
    # cohenD_check = [np.quantile(cohenD_BS, .5), np.median(cohenD_BS)]
    # print(cohenD_check)

    pername_expression_table_list.append(pd.DataFrame({'name': [nameStruct],
                                                       'acronym': [acronymStruct],
                                                       'cohenD': [cohenD_ac],
                                                       'cohenD_CI': [cohenD_ac_CI],
                                                       't_stat': [t_stat],
                                                       'pVal': [p_val],
                                                       'WT_mean': [mean_WT],
                                                       'WT_std': [std_WT],
                                                       'WT_mean_count': [mean_count_WT],
                                                       'het_mean': [mean_het],
                                                       'het_std': [std_het],
                                                       'het_mean_count': [mean_count_het],
                                                       'KO_mean': [mean_KO],
                                                       'KO_std': [std_KO],
                                                       'KO_mean_count': [mean_count_KO]}))

    # nameStruct_filename = "".join([c for c in nameStruct if c.isalpha() or c.isdigit() or c == ' ']).rstrip()
    # mouse_table_nameStruct.to_csv(os.path.join(analysis_path, 'perstructure', nameStruct_filename+'_volumes_mouse.csv'))

pername_expression_table = pd.concat(pername_expression_table_list, ignore_index=True)
pername_expression_table = pername_expression_table[np.logical_not(np.isnan(pername_expression_table['pVal']))]
pername_expression_table['pValBon'] = multipletests(pername_expression_table['pVal'], method='bonferroni')[1]
pername_expression_table['pValFDR'] = multipletests(pername_expression_table['pVal'], method='fdr_bh')[1]
pername_expression_table = pername_expression_table.sort_values(by='pVal')
pername_expression_table.loc[pername_expression_table['name'] == 'root', 'name'] = 'nonroot'
pername_expression_table.loc[pername_expression_table['acronym'] == 'root', 'acronym'] = 'nonroot'
pername_expression_table.loc[pername_expression_table['acronym'] == 'CUL4 5', 'acronym'] = 'CUL4, 5'
pername_expression_table.loc[pername_expression_table['acronym'] == 'Mmd', 'acronym'] = 'MMd'
pername_expression_table.to_csv('pername_expression_table.csv')

pername_merged_table = pd.merge(left=pername_expression_table,
         right=mouse_volume_pername_table,
         left_on='acronym',
         right_on='acronym',
         how='inner') # outer for testing purposes
pername_merged_table = pername_merged_table.drop(columns=['name_y', 'pValBon_x', 'pValBon_y', 'pValBon_BrainNormalized',
                                                          'WT_mean_AllenNormalized', 'WT_std_AllenNormalized',
                                                          'KO_mean_AllenNormalized', 'KO_std_AllenNormalized'])
pername_merged_table = pername_merged_table.rename(columns={'cohenD_x': 'cohenD_histo',
                                                            'cohenD_y': 'cohenD_mri',
                                                            'cohenD_CI_x': 'cohenD_CI_histo',
                                                            'cohenD_CI_y': 'cohenD_CI_mri',
                                                            'WT_mean_x': 'WT_mean_histo',
                                                            'WT_mean_y': 'WT_mean_mri',
                                                            'WT_mean_count': 'WT_mean_count_histo',
                                                            'WT_std_x': 'WT_std_histo',
                                                            'WT_std_y': 'WT_std_mri',
                                                            'het_mean': 'het_mean_histo',
                                                            'het_mean_count': 'het_mean_count_histo',
                                                            'het_std': 'het_std_histo',
                                                            'KO_mean_x': 'KO_mean_histo',
                                                            'KO_mean_y': 'KO_mean_mri',
                                                            'KO_mean_count': 'KO_mean_count_histo',
                                                            'KO_std_x': 'KO_std_histo',
                                                            'KO_std_y': 'KO_std_mri',
                                                            't_stat_x': 't_stat_histo',
                                                            't_stat_y': 't_stat_mri',
                                                            'pVal_x': 'pVal_histo',
                                                            'pVal_y': 'pVal_mri',
                                                            'pValFDR_x': 'pValFDR_histo',
                                                            'pValFDR_y': 'pValFDR_mri',
                                                            'pVal_BrainNormalized': 'pVal_norm_mri',
                                                            'pValFDR_BrainNormalized': 'pValFDR_norm_mri',
                                                            'WT_mean_BrainNormalized': 'WT_mean_norm_mri',
                                                            'WT_std_BrainNormalized': 'WT_std_norm_mri',
                                                            'KO_mean_BrainNormalized': 'KO_mean_norm_mri',
                                                            'KO_std_BrainNormalized': 'KO_std_norm_mri',
                                                            'cohenD_BrainNormalized': 'cohenD_norm_mri',
                                                            'cohenD_BrainNormalized_CI': 'cohenD_CI_norm_mri',
                                                            'name_x': 'name'})
pername_merged_table.to_csv('pername_merged_table.csv')



#
fig10 = plt.figure()
sns.regplot(data=pername_merged_table,
                x='cohenD_histo', y='cohenD_mri')
plt.show()
sns.regplot(data=pername_merged_table,
                x='cohenD_histo', y='cohenD_norm_mri')
plt.show()

notnan = np.logical_and(np.logical_not(np.isnan(pername_merged_table['cohenD_histo'])),
                        np.logical_not(np.isnan(pername_merged_table['cohenD_mri'])))
notnan2 = np.logical_and(np.logical_not(np.isnan(pername_merged_table['cohenD_histo'])),
                        np.logical_not(np.isnan(pername_merged_table['cohenD_norm_mri'])))
slope, intercept, r_value, p_value, std_err = linregress(pername_merged_table.loc[notnan, 'cohenD_histo'],
                                                         pername_merged_table.loc[notnan, 'cohenD_mri'])
slope2, intercept2, r2_value, p2_value, std2_err = linregress(pername_merged_table.loc[notnan2, 'cohenD_histo'],
                                                         pername_merged_table.loc[notnan2, 'cohenD_norm_mri'])


# pername_expression_table = pername_expression_table.reindex(columns = ['name',
#                                                              'cohenD', 'cohenD_BrainNormalized',
#                                                              'cohenD_CI', 'cohenD_BrainNormalized_CI',
#                                                              't_stat', 'pVal', 'pVal_BrainNormalized',
#                                                              'pValBon', 'pValBon_BrainNormalized',
#                                                              'pValFDR', 'pValFDR_BrainNormalized',
#                                                              'WT_mean', 'WT_std',
#                                                              'KO_mean', 'KO_std',
#                                                              'WT_mean_AllenNormalized', 'WT_std_AllenNormalized',
#                                                              'KO_mean_AllenNormalized', 'KO_std_AllenNormalized',
#                                                              'WT_mean_BrainNormalized', 'WT_std_BrainNormalized',
#                                                              'KO_mean_BrainNormalized', 'KO_std_BrainNormalized'])
# mouse_table_pername.to_csv(pername_table_path_list[iIncludeInTest - 1])
