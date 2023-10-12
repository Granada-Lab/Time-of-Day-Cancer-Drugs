########################################################################################################################
# SCRIPT CLUSTER AND CORRELATION ANALYSIS OF CIRCADIAN RHYTMICITIES IN TRIPLE NEGATIVE BREAST CANCER (TNBC) ############
# Collaboration with The Granada Lab from Charite, Berlin. ###### Code written by Jeff DIDIER, Sebastien DE LANDTSHEER #
########################################################################################################################
#########################
# ## Importing libraries
#########################
import os

import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt

from adjustText import adjust_text
from matplotlib.patches import Patch
from sklearn.decomposition import PCA
from functions.plot_functions import plot_LDA
from functions.data_functions import remove_columns_with_zeros, binarize_with_kmeans

###########################
# ## 1. DATA PREPROCESSING
###########################

# import circadian data
data_bmal1 = pd.read_excel("../Data_fig_6/TNBC_Circadian_Values_v5.xlsx", sheet_name="BMAL1")
data_per2 = pd.read_excel("../Data_fig_6/TNBC_Circadian_Values_v5.xlsx", sheet_name="PER2")
data_bmal1_per2 = pd.read_excel("../Data_fig_6/TNBC_Circadian_Values_v5.xlsx", sheet_name="BMAL1_PER2")
data_scores = pd.read_excel("../Data_fig_6/TNBC_Circadian_Values_v5.xlsx", sheet_name="Scores")
data_growth = pd.read_excel("../Data_fig_6/TNBC_Circadian_Values_v5.xlsx", sheet_name="Growth")
data_tod = pd.read_excel("../Data_fig_6/TNBC_Circadian_Values_v5.xlsx", sheet_name="ToD")

# ## at the moment we exclude MCF10A : Is this still the case? Yes, because it was not in CCLE

# import CCLE data
data_exp = pd.read_csv("../Data_fig_6/CCLE_expression_from_CCLE_2022_11_22.csv")
data_info = pd.read_csv("../Data_fig_6/sample_info.csv").set_index('CCLE_Name')

# removing the hyphens in CCLE names and retrieving the DepMap IDs, setting as index, and removing MCF10A
for sets in [data_bmal1, data_per2, data_bmal1_per2, data_scores, data_growth]:
    sets['cellline'] = [x.replace('-', '') for x in sets.loc[:, 'cellline']]
    sets['DepMap_ID'] = [data_info[data_info['stripped_cell_line_name'] == x].loc[:,
                         'DepMap_ID'].to_numpy()[0] for x in sets.loc[:, 'cellline']]
    sets.set_index('DepMap_ID', inplace=True)
    sets.sort_index(inplace=True)
    sets.drop(sets[sets['cellline'] == 'MCF10A'].index, inplace=True)
    sets.drop('cellline', axis=1, inplace=True)

# ## Expression dataset
# filter data: only cell lines we need
data_exp_filtered = data_exp.loc[data_exp['Unnamed: 0'].isin(data_bmal1.index)].set_index('Unnamed: 0')
# remove genes for which there is more  than one "0.00"
data_exp_filtered = remove_columns_with_zeros(data_exp_filtered, 2)  # 19221 -> 14901 (4320 genes removed)
# remove the number from the gene name
colnames = data_exp_filtered.columns
colnames_correct = [x.split()[0] for x in colnames]
data_exp_filtered.columns = colnames_correct

# expression data - normalizing the expression (min-max) [ xscaled = (x - min(x)) / (max(x) - min(x)) ]
data_exp_norm = (data_exp_filtered-data_exp_filtered.min())/(data_exp_filtered.max()-data_exp_filtered.min())

# SUBSETTING EXPRESSION DATA TO A SHORT AND LONGER CORE-CLOCK GENES LIST (and alternative names)
short_ccg_list = pd.read_excel("../Data_fig_6/circadian_clock_gene_lists_v2.xlsx", sheet_name="genelist_short",
                               skiprows=1, header=None, names=['Gene', 'AlternativeName'])  # short gene list (16 genes)
long_ccg_list = pd.read_excel("../Data_fig_6/circadian_clock_gene_lists_v2.xlsx", sheet_name="genelist_long",
                              skiprows=1, header=None, names=['Gene', 'AlternativeName'])  # long gene list (60 genes)

# check if all genes are in the CCLE data
for num, ccg in enumerate([short_ccg_list, long_ccg_list]):
    missing = []
    for gene in ccg['Gene']:
        if gene not in data_exp_norm.columns:
            print(f'{gene} not in CCLE, checking for alternative name.')
            alt_gene = list(ccg[ccg['Gene'] == gene]['AlternativeName'])
            if pd.isna(alt_gene):
                print(f'{gene} has no alternative name.\n')
                missing.append(gene)
            elif alt_gene[0] not in data_exp_norm.columns or alt_gene[0].replace('‐', '') not in data_exp_norm:
                print(f'Alternative name of {gene} <{alt_gene[0]} or {alt_gene[0].replace("‐", "")}> not in CCLE.\n')
                missing.append(gene)
            else:
                print(f'Alternative name of {gene} <{alt_gene[0]}> found in CCLE.\n')
                ccg.loc[ccg['Gene'] == gene, 'Gene'] = alt_gene[0]
    print(f'{missing} not in CCLE data <{"short list" if num == 0 else "long list"}>.\n')

# 13 first-name genes of the long lists were not found, but their alternative names were found and replaced accordingly
# Still, 7 genes of the long list['BTRCP', 'CUL4', 'GATAD2a', 'NCoR', 'SUV39H', 'TNFA', 'TRAP150'] are not among the
# cleaned and normalized CCLE data, they either have no alternatives or their alternative names was also not in CCLE
# we will therefore remove them from the lists

to_remove = ['BTRCP', 'CUL4', 'GATAD2a', 'NCoR', 'SUV39H', 'TNFA', 'TRAP150']
short_ccg_clean = [ele for ele in short_ccg_list['Gene'] if ele not in to_remove]
long_ccg_clean = [ele for ele in long_ccg_list['Gene'] if ele not in to_remove]

# short and long ccg dataframes
short_ccg_norm = data_exp_norm[short_ccg_clean]
long_ccg_norm = data_exp_norm[long_ccg_clean]

X_short_unscaled = data_exp_filtered[short_ccg_clean].sort_index()
X_long_unscaled = data_exp_filtered[long_ccg_clean].sort_index()

X_short = short_ccg_norm.sort_index()  # final short core clock genes expression (16 genes)
X_long = long_ccg_norm.sort_index()  # final long core clock genes expression (53 genes)
X_full = data_exp_norm.sort_index()  # Full gene set, most likely not needed

# ## our targets considered for analysis
targets = list(data_bmal1_per2.columns) + list(data_scores.columns) + list(data_growth.columns)

# adding binarized oscillator strength in two ways based on median and kmeans binarization of the different scores
scores = ['meanscore', 'mra_score', 'cvphdiff_score', 'ridgelength_score']
y_oscillators_median = data_scores[scores] > data_scores[scores].median()
y_oscillators_median.columns = ['meanscore_med', 'mra_score_med', 'cvphdiff_score_med', 'ridgelength_score_med']
y_oscillators_kmeans = data_scores[scores] < binarize_with_kmeans(data_scores[scores], random_state=42)
y_oscillators_kmeans.columns = ['meanscore_km', 'mra_score_km', 'cvphdiff_score_km', 'ridgelength_score_km']

# add those columns to our targets
targets_final = pd.Series(targets +
                          list(y_oscillators_median.columns) +
                          list(y_oscillators_kmeans.columns)).unique()

# add the kmeans, median, and growth rate values to the 3 final sets BMAL1, PER2, BMAL1_PER2
# we only need subtype once, but it is present in the gene-set and the growth rate, let's remove it from growth rate
data_growth.drop('subtype', axis=1, inplace=True)
data_scores.drop('subtype', axis=1, inplace=True)
sets = [data_bmal1, data_per2, data_bmal1_per2]
for i, df in enumerate(sets):
    sets[i] = pd.concat([df, data_growth, data_scores, y_oscillators_kmeans, y_oscillators_median],
                        axis=1).reset_index(drop=False)
    sets[i].set_index('DepMap_ID', inplace=True)
# assign the concatenated sets to their variables, rowname = DepMap_ID, first column = subtype, final 8 columns are bool
final_bmal1 = sets[0]
final_per2 = sets[1]
final_bmal1_per2 = sets[2]

# ## regarding time of day drug sensitivity data
# clear the ToD data set , keep maxrange_spline for each of the drug, subtype, and data availability
for col in data_tod.columns[2:]:
    if 'Unnamed' not in col:
        # store the drug name if encountered
        drug_name_tmp = col
    if 'maxrange_spline' not in data_tod[col].values:
        # drop the col if maxrange spline is not included
        data_tod.drop(col, axis=1, inplace=True)
    else:
        # bypass pycharm fals warning
        drug_name_tmp = globals().get('drug_name_tmp')
        # rename the specific column with the drug name to map maxrange spline with correct drug
        data_tod.rename(columns={col: drug_name_tmp}, inplace=True)
# drop first row (nans) and re index
data_tod = data_tod.iloc[1:, :].reset_index(drop=True)

# get the tod data in the exact same order as the others
data_tod['cellline'] = [x.replace('-', '') for x in data_tod.loc[:, 'cellline']]
data_tod['DepMap_ID'] = [data_info[data_info['stripped_cell_line_name'] == x].loc[:,
                         'DepMap_ID'].to_numpy()[0] for x in data_tod.loc[:, 'cellline']]
data_tod.set_index('DepMap_ID', inplace=True)
data_tod.sort_index(inplace=True)
data_tod.drop(data_tod[data_tod['cellline'] == 'MCF10A'].index, inplace=True)
data_tod.drop('cellline', axis=1, inplace=True)

# we can basically drop the drug Olaparib as there is one single value available
data_tod.drop('Olaparib', axis=1, inplace=True)


#####################################################
# ## Version 5 Time of day drug sensitivity analysis
#####################################################
formats = ['png', 'svg']

if os.path.isdir('./results') is False:
    os.mkdir('./results')

# ## Figure 6 b) and d) ################################################################################################
# check correlation of drug sensitivity and short CCLE expression (restrain to available drug sensititvity data)

corr_matrix = pd.DataFrame(index=data_tod.columns[1:],
                           columns=X_short_unscaled.columns,
                           dtype=float)
p_val_matrix = pd.DataFrame(index=data_tod.columns[1:],
                            columns=X_short_unscaled.columns,
                            dtype=float)
for para in X_short_unscaled.columns:
    for tod_drug in data_tod.columns[1:]:
        spearman_rank_stats = \
            stats.spearmanr(X_short_unscaled[para][data_tod[tod_drug].notnull()],
                            data_tod[tod_drug][data_tod[tod_drug].notnull()])
        corr_matrix.loc[tod_drug, para] = spearman_rank_stats.correlation
        p_val_matrix.loc[tod_drug, para] = spearman_rank_stats.pvalue
# save corr and pvals as excel
with pd.ExcelWriter('./results/6_b.xlsx') as writer:
    corr_matrix.to_excel(writer, sheet_name='Correlation_coefficient')
    p_val_matrix.to_excel(writer, sheet_name='P_values')
# figures
fig = sns.clustermap(corr_matrix, cmap=plt.get_cmap('coolwarm'),
                     figsize=(10, 10), metric='correlation',
                     dendrogram_ratio=(0.1, 0.4), cbar_pos=(0.05, 0.82, 0.02, 0.15), annot=True, fmt='.2f',
                     vmin=-1, center=0, vmax=1, linewidth=.003,
                     annot_kws={'ha': "center",
                                'va': "center",
                                'fontsize': 8,
                                'color': 'k'})
# annotate only interesting correlations
for text in fig.ax_heatmap.texts:
    text.set_visible(True if float(text.get_text()) > 0.5 or float(text.get_text()) < - 0.5 else False)
# modify re-arranged labels to include sample size
dendro_labels = fig.ax_heatmap.yaxis.get_majorticklabels()
dendro_labels = [lab.get_text() for lab in dendro_labels]
samples_per_drug = []
for tod_drug in dendro_labels:
    samples_per_drug.append(len(data_tod[tod_drug][data_tod[tod_drug].notnull()]))
y_labels = [x + f' ({s})' for x, s in zip(dendro_labels, samples_per_drug)]
# set the labels plus sample size
fig.ax_heatmap.set_yticklabels(y_labels, fontsize=10)
# add frame lines around heatmap
fig.ax_heatmap.axhline(y=0, color='k', linewidth=2)
fig.ax_heatmap.axhline(y=corr_matrix.shape[0], color='k', linewidth=2)
fig.ax_heatmap.axvline(x=0, color='k', linewidth=2)
fig.ax_heatmap.axvline(x=corr_matrix.shape[1], color='k', linewidth=2)
fig.ax_col_dendrogram.set_title("Spearman Rank correlation of CCLE core clock genes short list and drug "
                                "sensitivity\nSample size in parenthesis")
# adapt colorbar ticks
cbar = fig.ax_heatmap.collections[0].colorbar
cbar.ax.tick_params(labelsize=8)
cbar.ax.set_yticks([-1, 0, 1], [-1, 0, 1])
# complementary subfigures: which ToD best correlates with all genes, and which gene best correlates with all ToD
fig.fig.subplots_adjust(bottom=0.5)
# restore color bar posiion
x0, y0, _w, _h = fig.cbar_pos
fig.ax_cbar.set_position([x0*1.5, y0, _w*1.5 , _h])
# add axes (per drugs 20% of the plot width if long ccle, else 30%)
ax1 = fig.fig.add_axes([0.103, 0.1, 0.788*0.3, 0.3])
drug_mean_ranked = corr_matrix.abs().mean(axis=1).sort_values(ascending=False).index.to_list()
abs_cumul_per_drug = corr_matrix.abs().T[drug_mean_ranked].melt(var_name=['drug'],
                                                                value_name='Spearman rank correlations')
sns.boxplot(data=abs_cumul_per_drug, x='drug', y='Spearman rank correlations', palette='Greens_r', showmeans=True,
            meanprops={"marker": "o",
                       "markerfacecolor": "white",
                       "markeredgecolor": "black",
                       "markersize": "5",
                       "alpha": 0.8,
                       "zorder": 5})
sns.swarmplot(data=abs_cumul_per_drug, x='drug', y='Spearman rank correlations', size=3.5, palette='Reds_r',
              linewidth=.3, edgecolor='gray')
ax1.set_xticklabels(drug_mean_ranked, rotation=90)
ax1.set_xlabel(None)
# add axes (per genes 80% of the plot width if long ccle, else 70%)
ax2 = fig.fig.add_axes([0.3394+0.05, 0.1, (0.788*0.7)-0.05, 0.3])
gene_mean_ranked = corr_matrix.abs().mean(axis=0).sort_values(ascending=False).index.to_list()
abs_cumul_per_gene = corr_matrix.abs()[gene_mean_ranked].melt(var_name=['gene'],
                                                              value_name='Spearman rank correlations')
sns.boxplot(data=abs_cumul_per_gene, x='gene', y='Spearman rank correlations', palette='Purples_r', showmeans=True,
            meanprops={"marker": "o",
                       "markerfacecolor": "white",
                       "markeredgecolor": "black",
                       "markersize": "5",
                       "alpha": 0.8,
                       "zorder": 5})
sns.swarmplot(data=abs_cumul_per_gene, x='gene', y='Spearman rank correlations', size=3.5, palette='Oranges_r',
              linewidth=.3, edgecolor='gray')
ax2.set_xticklabels(gene_mean_ranked, rotation=90)
ax2.set_xlabel(None)
ax2.set_ylabel(None)
for file in formats:
    plt.savefig(f'./results/6_b_d.{file}', bbox_inches='tight',
                dpi=300)
plt.close()

# Short CCLE: Paclitaxel, Cisplatin and Adavosertib show most correlations above 0.5
# (5 genes, P: CRY2, DBP, PER3, RORB, RORC; A: DBP, PER1, PER2, PER3, RORB; C: NR1D1, RORB, CRY2, PER2, DBP)
# Gene RORB have the most correlated drugs
# (RORB: Torin2, Alpelisib, Paclitaxel, Adavosertib, Cisplatin)
# Paclitaxel and PER3 have the overall highest correlation (-0.88)

# ## Figure S6 a) ######################################################################################################
# simple clustermap of the scaled ccle core clock gene expressions for each cell line, correlation metrics

# only the 9 cell lines that we effectively used, we take 5FU here cause all drugs have the same samples tested
df = X_short[data_tod['5FU'].notnull()]
fig = sns.clustermap(df, cmap=plt.get_cmap('coolwarm'), figsize=(8, 10),
                     metric='euclidean', method='complete', dendrogram_ratio=(0.2, 0.2),
                     cbar_pos=(0.05, 0.88, 0.02, 0.10), annot=False, vmin=0, center=0.5, vmax=1, linewidth=.003)
fig.ax_heatmap.set_ylabel(None)
dendro_labels = fig.ax_heatmap.yaxis.get_majorticklabels()
dendro_labels = [lab.get_text() for lab in dendro_labels]
cellline_labels = [data_info['cell_line_name'][data_info['DepMap_ID'] == x][0] for x in dendro_labels]
fig.ax_heatmap.set_yticklabels(cellline_labels, fontsize=10, rotation=0)
# add frame lines around heatmap
fig.ax_heatmap.axhline(y=0, color='k', linewidth=2)
fig.ax_heatmap.axhline(y=df.shape[0], color='k', linewidth=2)
fig.ax_heatmap.axvline(x=0, color='k', linewidth=2)
fig.ax_heatmap.axvline(x=df.shape[1], color='k', linewidth=2)
fig.ax_col_dendrogram.set_title("Clustermap of the short CCG expression in the 9 breast cancer cell lines\n"
                                "euclidean distance")
# adapt colorbar ticks
cbar = fig.ax_heatmap.collections[0].colorbar
cbar.ax.tick_params(labelsize=8)
cbar.ax.set_yticks([0, 0.5, 1], [0, 0.5, 1])
for file in formats:
    plt.savefig(f'./results/S6_a.{file}', bbox_inches='tight',
                dpi=300)
plt.close()

# ## Figure S6 i) ######################################################################################################
# Deeper PCA plots for Paclitaxel
# Predict drug sensitivity using CCLE short and long core clock genes (with LDA, we need to binarize the ToD)
drugs = data_tod.columns[1:]
# get a mask to restore the nans afer binarization
nan_mask = data_tod[drugs].isna()
data_tod_med = data_tod[drugs] > data_tod[drugs].median()
# restore nans
data_tod_med[nan_mask] = np.nan
data_tod_med.columns = data_tod_med.columns + '_med'
data_tod_kmeans = binarize_with_kmeans(data_tod[drugs], random_state=42)
data_tod_kmeans.replace({0: False, 1: True}, inplace=True)
data_tod_kmeans.columns = data_tod_kmeans.columns + '_km'
# add subtype to both
bin_tod_df = pd.concat([data_tod['subtype'], data_tod_med, data_tod_kmeans], axis=1)

# need the med and kms here because the ToD truths are different fror each binarized drug form
best_drugs = ['Paclitaxel_med']
remaining_drugs = ['5FU_med', 'Torin2_med', 'Alpelisib_med', 'Doxorubicin_med', 'Adavosertib_med', 'Alisertib_med',
                   'Cisplatin_med', '5FU_km', 'Torin2_km', 'Alpelisib_km', 'Doxorubicin_km', 'Paclitaxel_km',
                   'Adavosertib_km', 'Alisertib_km', 'Cisplatin_km']

width = 0.25  # for the barplots
for drug in best_drugs + remaining_drugs:
    df = X_short[bin_tod_df[drug].notnull()]
    truth, _ = bin_tod_df[bin_tod_df[drug].notnull()][drug].factorize()  # this truth changes with each drug
    # plot PCA loadings (comp 1 & 2), barplot, and elbow curve
    pca = PCA(n_components=4, random_state=42)  # we now we only need 4 here
    pca.fit(df)
    PCs = pca.fit_transform(df)
    PCdf = pd.DataFrame(data=PCs, columns=["PC" + str(i) for i in range(1, PCs.shape[1] + 1)])
    PCdf.index = df.index
    # Match PC names to loadings
    pc_loadings = dict(zip(PCdf.columns, pca.components_))
    # Matrix of corr coefficients between pcs and features
    loadings_df = pd.DataFrame.from_dict(pc_loadings)
    loadings_df['feature_names'] = df.columns
    loadings_df = loadings_df.set_index('feature_names')
    # plot the loading plot with the scatterplot
    xs = pca.components_[0] * PCdf['PC1'].max()**2
    ys = pca.components_[1] * PCdf['PC2'].max()**2
    # loadings plot
    g = sns.lmplot(x='PC1', y='PC2', data=pd.concat([PCdf,
                                                     pd.DataFrame(truth, index=df.index,
                                                                  columns=['ToD sensitivity'])],
                                                    axis=1),
                   fit_reg=False, hue='ToD sensitivity', palette={0: 'b', 1: 'r'})
    sns.move_legend(g, "upper right", bbox_to_anchor=(1.02, .98), frameon=True)
    texts = []
    for i, varnames in enumerate(df.columns):
        plt.arrow(0, 0, xs[i], ys[i], color='grey', linestyle='--', linewidth=0.5, head_width=0.01)
        plt.scatter(xs[i], ys[i], s=75, marker='*', color='k', alpha=0.5)
        texts.append(plt.text(xs[i], ys[i], varnames, fontsize=9))
    adjust_text(texts, only_move={'points': '', 'text': 'y'}, force_text=0.05, autoalign='y')
    plt.xlabel(f'PC1 ({round(pca.explained_variance_[0]*100, 2)}%)')
    plt.ylabel(f'PC2 ({round(pca.explained_variance_[1]*100, 2)}%)')
    plt.title(f'PCA 2D Biplot on short core clock genes, colorized by {drug}\n'
              f'Loadings were upscaled to fit PCs')
    for file in formats:
        plt.savefig(f'./results/S6_i_{drug}_biplot.{file}',
                    bbox_inches='tight', dpi=300)
    plt.close()
    # ranked loadings

# plotting the ranked loadings of the first three principals per drug binarized (as most drugs have same samples, some
# of the rankings will be identical
drugs = data_tod.columns[1:]
with pd.ExcelWriter(f'./results/S6_i_loadings.xlsx') as writer:
    for drug in drugs:
        df = X_short[data_tod[drug].notnull()]
        # plot PCA loadings (comp 1 , 2)
        pca = PCA(n_components=2, random_state=42)  # we now we only need 3 here
        pca.fit(df)
        PCs = pca.fit_transform(df)
        PCdf = pd.DataFrame(data=PCs, columns=["PC" + str(i) for i in range(1, PCs.shape[1] + 1)])
        PCdf.index = df.index
        # Match PC names to loadings
        pc_loadings = dict(zip(PCdf.columns, pca.components_))
        # Matrix of corr coefficients between pcs and features
        loadings_df = pd.DataFrame.from_dict(pc_loadings)
        loadings_df['feature_names'] = df.columns
        loadings_df = loadings_df.set_index('feature_names')
        # save loadings, in fact shorts and longs are the same for all drugs, except Cisplatin which has 1 sample less
        pd.concat([loadings_df["PC1"][loadings_df["PC1"].abs().sort_values(ascending=False).index],
                   # to keep the same order, we will use the loadings of PC2 but sorted like PC1
                   loadings_df["PC2"][loadings_df["PC1"].abs().sort_values(ascending=False).index]],
                   axis=1).to_excel(writer, sheet_name=drug)
        # figure
        f, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
        ax1.bar(np.arange(len(loadings_df.index)), loadings_df["PC1"].abs().sort_values(ascending=False), width=0.5,
                label="PC1", edgecolor="black",
                hatch=['////' if i < 0 else None for i in loadings_df["PC1"][
                    loadings_df["PC1"].abs().sort_values(ascending=False).index]])
        ax2.bar(np.arange(len(loadings_df.index)), loadings_df["PC2"].abs().sort_values(ascending=False), width=0.5,
                label="PC2", edgecolor="black",
                hatch=['////' if i < 0 else None for i in loadings_df["PC2"][
                    loadings_df["PC2"].abs().sort_values(ascending=False).index]])
        for ax in [ax1, ax2]:
            ax.set_xticks(np.arange(len(loadings_df.index)))
            ax.tick_params(axis="x", labelsize=7, labelrotation=90)
            lgd = ax.legend(fontsize=12, title_fontsize=18)
            handles, labs = lgd.axes.get_legend_handles_labels()
            handles.append(Patch(facecolor='white', edgecolor='black', hatch='////'))
            labs.append('Negative')
            lgd._legend_box = None
            lgd._init_legend_box(handles, labs)
            lgd._set_loc(lgd._loc)
            lgd.set_title(lgd.get_title().get_text())
            for text in lgd.get_texts():
                text.set_weight('bold')
            for num, ha in enumerate(lgd.legendHandles):
                if num < len(lgd.legendHandles) - 1:
                    ha.set_hatch(None)
        ax1.set_title(f"PC1 ({round(pca.explained_variance_[0]*100, 2)}%)")
        ax1.set_xticklabels(loadings_df["PC1"].abs().sort_values(ascending=False).index)
        ax2.set_title(f"PC2 ({round(pca.explained_variance_[1]*100, 2)}%)")
        ax2.set_xticklabels(loadings_df["PC2"].abs().sort_values(ascending=False).index)

        plt.suptitle(f"{drug}: Ranked PC loadings in short CCG expression")
        plt.tight_layout()
        for file in formats:
            plt.savefig(f'./results/S6_i_loadings_{drug}.{file}',
                        bbox_inches='tight', dpi=300)
        plt.close()

# ## Figure 6e-f, S6b-h) ##############################################################################################
# Predict drug sensitivity using CCLE short and long core clock genes (with LDA, we need to binarize the ToD)

# LDA with both subtypes and cell line name annotated, plus saving all loadings
short_loading_df = pd.DataFrame()
long_loading_df = pd.DataFrame()
for drug_bin in bin_tod_df.columns[1:]:
    if 'med' in drug_bin:
        df = X_short[bin_tod_df[drug_bin].notnull()]
        # plot LDA with both annotations
        lda_norm_loads = plot_LDA(df,
                                  pd.concat([df,
                                             bin_tod_df[drug_bin][bin_tod_df[drug_bin].notnull()].astype(int),
                                             bin_tod_df['subtype'][bin_tod_df[drug_bin].notnull()]], axis=1), drug_bin,
                                  lda_output_target=drug_bin, annot_subtype=True, annot_cellline=True,
                                  data_info=data_info, autoalign='y', only_move={'points': 'y', 'text': 'xy'},
                                  force_text=1, return_loading=True, title=f'LDA expression analysis in short '
                                                                           f'CCLE CCG set targeting {drug_bin}')
        for file in formats:
            plt.savefig(f'./results/6_e_{drug_bin}.{file}' if 'Paclitaxel' in drug_bin else
                        f'./results/S6_b-h_{drug_bin}.{file}',
                        bbox_inches='tight', dpi=300)
        plt.close()
        short_loading_df = pd.concat([short_loading_df, lda_norm_loads], axis=0)

# save all together
with pd.ExcelWriter(f'./results/LDA_loadings_all.xlsx') as writer:
    short_loading_df.to_excel(writer, sheet_name='short')


# Read the file in case you don't run the above chunk
# short_loading_df = pd.read_excel(r'./results/LDA_loadings_all.xlsx', sheet_name='short',
# index_col='Unnamed: 0')

# generate the cumulative LDA loadings of short and long genes for each drug binarized by MEDIAN only
# only medians
df = short_loading_df.loc[[x for x in short_loading_df.index if '_med' in x]]
# multiply by 100 to showcase pourcentage of loading
df = df * 100
# get the rank of best gene
gene_mean_ranked = df.mean(axis=0).sort_values(ascending=False).index.to_list()
abs_cumul_per_gene = df[gene_mean_ranked].melt(var_name=['gene'], value_name='Normalized LDA loading [%]')
plt.figure(figsize=(10, 10))
fig = sns.boxplot(data=abs_cumul_per_gene, x='gene', y='Normalized LDA loading [%]', palette='Purples_r',
                  showmeans=True,
                  meanprops={"marker": "o",
                             "markerfacecolor": "white",
                             "markeredgecolor": "black",
                             "markersize": "5",
                             "alpha": 0.8,
                             "zorder": 5})
sns.swarmplot(data=abs_cumul_per_gene, x='gene', y='Normalized LDA loading [%]',
              size=4.5,
              palette='Reds_r',
              linewidth=.3, edgecolor='gray', ax=fig)
fig.set_xticklabels(gene_mean_ranked, rotation=90)
fig.set_xlabel(None)
fig.set_title(f"Distribution of short clock gene set LDA loadings across drugs")
for file in formats:
    plt.savefig(f'./results/6_f.{file}',
                bbox_inches='tight', dpi=300)
plt.close()

# ## Figure 6 c) #######################################################################################################
# gene correlations with ToD
version_5_best_genes = dict({'Paclitaxel': ['DBP', 'PER3']})

for drug, genes in version_5_best_genes.items():
    for gene in genes:
        data = pd.concat([X_short_unscaled[gene][data_tod[drug].notnull()],
                         data_tod[drug][data_tod[drug].notnull()]], axis=1)
        g = sns.jointplot(data=data.astype(float),
                          x=gene, y=drug, kind='reg', scatter=False, color='grey')
        data_scatter = pd.concat([data, final_bmal1_per2['subtype'][data_tod[drug].notnull()]], axis=1)
        scat = sns.scatterplot(data=data_scatter,
                               x=gene, y=drug, hue='subtype', ax=g.ax_joint, s=75)
        g.ax_joint.legend(bbox_to_anchor=(0.5, -0.3), loc='lower center', borderaxespad=0, ncol=3)
        s, p = stats.spearmanr(X_short_unscaled[gene][data_tod[drug].notnull()],
                               data_tod[drug][data_tod[drug].notnull()])
        if p < 0.05:
            print(f'Significance found: {drug} vs {gene}, p={round(p, 3)} !')
        g.ax_joint.annotate(f'Spearman $\\rho$ = {s:.3f},\n$\\rho$\u00b2 = {s**2:.3f}\np-value = {p:.3f}',
                            xy=(0.05, 0.95), xycoords='axes fraction',
                            ha='left', va='center',
                            bbox={'boxstyle': 'round', 'fc': 'powderblue', 'ec': 'navy', 'alpha': 0.3})
        g.set_axis_labels(xlabel=gene + ' log2(TPM)', ylabel=drug, size=15)
        g.fig.suptitle(f'Joint linear regression plot {gene} vs {drug}')
        g.fig.tight_layout()
        plt.show()
        for file in formats:
            plt.savefig(f'./results/6_c_{gene}_vs_{drug}.{file}',
                        bbox_inches='tight', dpi=300)
        plt.close()

# if using top genes of Paclitaxel (version 5)
# Significance found: Paclitaxel vs DBP, p=0.03 !
# Significance found: Paclitaxel vs PER3, p=0.002 !

########################################################################################################################
# END OF SCRIPT ########################################################################################################
########################################################################################################################
