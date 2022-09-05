# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% tags=["parameters"]
index_dir = "Data/Index/HEIndex0"
groups_file = "Data/Index/SrcFiles/IndexGroups.csv"
groups_cols = ["Age(weeks)", "Genotype", "Diet", "Sex"]
stranded = False
gff_file = "Data/Annotations/genomic.gff"
ue_regions_file = "Data/HyperEditing/ClusterScreening.SE_0.05_0.6_30_0.6_0.1_0.8_0.2/Original/UERegions/UEmergeSd20.stranded.A2G.bed"
MISMATCH = "A2G"
code_dir = "/private6/Projects/Nothobranchius_furzeri/Code"

# %%
# # %cd ..
# %cd {code_dir}

# %% [markdown]
# # Imports

# %%
from pathlib import Path
import sys

from sklearn.decomposition import PCA, SparsePCA, KernelPCA
from sklearn.preprocessing import StandardScaler
from more_itertools import powerset
from pybedtools import BedTool
import plotly.express as px
import pandas as pd
import numpy as np

sys.path.append(str(Path(".")))
from General.consts import STRANDED_MISMATCH_TYPES, UNSTRANDED_MISMATCH_TYPES
from EDA.pandas_utils import reorder_df_by_wanted_cols
from EditingUtils.pybedtools_utils import flip_strand

# %% [markdown]
# # Data

# %% [markdown]
# ## Getting Index data

# %%
index_dir = Path(index_dir)

# %%
mismatches = STRANDED_MISMATCH_TYPES if stranded else UNSTRANDED_MISMATCH_TYPES

# %%
groups_df = pd.read_csv(groups_file)
groups_df.head()


# %% papermill={"duration": 0.144608, "end_time": "2021-12-21T16:38:25.656889", "exception": false, "start_time": "2021-12-21T16:38:25.512281", "status": "completed"} tags=[]
def get_per_region_per_sample_dfs(index_dir, groups_df, mismatches):
    """
    Return a list of $SAMPLE/StrandDerivingCountsPerRegion.csv dataframes.
    """
    
    region_cols = ["GenomicRegion", "Start", "End"]
    id_cols = ["Group", "Sample"]
    # editing_indices_cols = [f"{mismatch}EditingIndex" for mismatch in mismatches]  # TODO check mismatches thing
    # wanted_cols = region_cols + id_cols + editing_indices_cols
    
    per_region_per_sample_dfs = []
        
    samples = groups_df["Sample"]
    sub_dirs = [Path(index_dir, f"{sample}") for sample in samples]
    for sub_dir in sub_dirs:
        per_region_per_sample_path = list(sub_dir.glob("**/StrandDerivingCountsPerRegion.csv"))[0]
        per_region_per_sample_df = pd.read_csv(per_region_per_sample_path)
        
        editing_indices_cols = [col for col in per_region_per_sample_df.columns if col.endswith("EditingIndex")]
        wanted_cols = region_cols + id_cols + editing_indices_cols
        
        per_region_per_sample_df = per_region_per_sample_df.filter(wanted_cols)
        for col in ["GenomicRegion", "Group", "Sample"]:
            # per_region_per_sample_df[col] = per_region_per_sample_df[col].astype("category")  # TODO verify it's OK
            per_region_per_sample_df[col] = per_region_per_sample_df[col].astype("str")
        per_region_per_sample_dfs.append(per_region_per_sample_df)
        
    return per_region_per_sample_dfs


# %%
per_region_per_sample_dfs = get_per_region_per_sample_dfs(index_dir, groups_df, mismatches)

per_region_per_sample_dfs[0].head(3)


# %%
def merge_per_region_per_sample_dfs(per_region_per_sample_dfs):
    id_vars = ["GenomicRegion", "Start", "End", "Group"]
    # value_vars = [f"{mismatch}EditingIndex" for mismatch in mismatches] # TODO check mismatches thing
    value_vars = [col for col in per_region_per_sample_dfs[0].columns if col.endswith("EditingIndex")]
    var_name = "Mismatch"

    new_per_region_per_sample_dfs = []

    for df in per_region_per_sample_dfs:
        sample = df["Sample"].values[0]
        value_name = sample
        df = pd.melt(df, id_vars=id_vars, value_vars=value_vars, var_name=var_name, value_name=value_name)
        df.loc[:, "Mismatch"] = df.loc[:, "Mismatch"].str[:3]
        new_per_region_per_sample_dfs.append(df)

    merged_df = new_per_region_per_sample_dfs[0].drop("Group", axis=1)
    
    for df in new_per_region_per_sample_dfs[1:]:
        df = df.drop("Group", axis=1)
        merged_df = merged_df.merge(df, 
                                    on=["GenomicRegion", "Start", "End", "Mismatch"],
                                    how="outer", 
                                   )
    
    merged_df = merged_df.fillna(0)
    
    no_zeroes_merged_df = merged_df.copy()
    cols = no_zeroes_merged_df.columns[4:]
    no_zeroes_merged_df[cols] = no_zeroes_merged_df[no_zeroes_merged_df[cols] > 0][cols]
    no_zeroes_merged_df = no_zeroes_merged_df.dropna()
    
    return merged_df, no_zeroes_merged_df


# %%
merged_df, no_zeroes_merged_df = merge_per_region_per_sample_dfs(per_region_per_sample_dfs)

# outpath = ...
# merged_df.to_csv(outpath,index=False)

# %%
merged_df

# %%
no_zeroes_merged_df

# %%
a2g_merged_df = merged_df.loc[merged_df["Mismatch"]==MISMATCH]
a2g_merged_df

# %%
a2g_no_zeroes_merged_df = no_zeroes_merged_df.loc[no_zeroes_merged_df["Mismatch"]==MISMATCH]
a2g_no_zeroes_merged_df


# %% [markdown]
# ## Preparing PCA data

# %%
def long_merged_df_to_wide(merged_df):
    
    long_merged_df = merged_df.copy()
    region_col = long_merged_df["GenomicRegion"].astype(str) + "_" + long_merged_df["Start"].astype(str) + "_" + \
                 long_merged_df["End"].astype(str) + "_" + long_merged_df["Mismatch"].astype(str)
    long_merged_df.insert(0, "Region", region_col)
    cols_to_drop = ["GenomicRegion", "Start", "End", "Mismatch"]
    long_merged_df = long_merged_df.drop(cols_to_drop, axis=1)
    long_merged_df = long_merged_df.set_index("Region")

    wide_merged_df = long_merged_df.transpose()
    wide_merged_df = wide_merged_df.reset_index().rename(columns={"index": "Sample"})
    wide_merged_df = wide_merged_df.rename_axis(None, axis = 1)
    wide_merged_df = wide_merged_df.merge(groups_df.loc[:, ["Sample", "Group"]], on="Sample")
    
    expanded_group_col_df = wide_merged_df["Group"].str.split("-", expand=True)
    col_nums = [x for x in expanded_group_col_df.columns]
    expanded_group_col_df = expanded_group_col_df.rename(columns={col_num: group for col_num, group
                                                                  in zip(col_nums, groups_cols)})
    wide_merged_df = pd.concat([wide_merged_df, expanded_group_col_df], axis=1)
    
    first_cols = ["Sample", "Group"] + groups_cols
    wide_merged_df = reorder_df_by_wanted_cols(wide_merged_df, first_cols)
    
    return wide_merged_df


# %%
wide_merged_df = long_merged_df_to_wide(merged_df)

wide_merged_df.head()

# %%
wide_no_zeroes_merged_df = long_merged_df_to_wide(no_zeroes_merged_df)

wide_no_zeroes_merged_df.head()

# %%
id_cols = ["Sample", "Group"] + groups_cols
a2g_cols = [col for col in wide_merged_df.columns if col.endswith(f"_{MISMATCH}")]
cols = id_cols + a2g_cols
a2g_wide_merged_df = wide_merged_df.filter(cols)
a2g_wide_merged_df.head(3)

# %%
id_cols = ["Sample", "Group"] + groups_cols
a2g_cols = [col for col in wide_no_zeroes_merged_df.columns if col.endswith(f"_{MISMATCH}")]
cols = id_cols + a2g_cols
a2g_wide_no_zeroes_merged_df = wide_no_zeroes_merged_df.filter(cols)
a2g_wide_no_zeroes_merged_df.head(3)

# %%
len(a2g_wide_no_zeroes_merged_df)


# %% [markdown]
# # Functions

# %% [markdown]
# ## PCA functions

# %% [markdown]
# One PCA creation & plot, samples are colored according to `group_col`

# %%
def plot_explained_variance(pca):
    # https://plotly.com/python/pca-visualization/#plotting-explained-variance
    cumulative_explained_variance = np.cumsum(pca.explained_variance_ratio_)
    fig = px.area(
        x=range(1, cumulative_explained_variance.shape[0] + 1),
        y=cumulative_explained_variance,
        labels={"x": "# Components", "y": "Explained Variance"}
    )
    fig.update_layout(template="simple_white")
    fig.show()


# %%
def most_important_feature_of_each_pc(pca, features_names):
    # https://stackoverflow.com/a/50845697/10249633 - part 2 of the answer
    
    # number of components
    n_pcs = pca.components_.shape[0]

    # get the index of the most important feature on EACH component
    # LIST COMPREHENSION HERE
    most_important = [np.abs(pca.components_[i]).argmax() for i in range(n_pcs)]

    # get the names
    most_important_names = [features_names[most_important[i]] for i in range(n_pcs)]

    # LIST COMPREHENSION HERE AGAIN
    dic = {'PC{}'.format(i): most_important_names[i] for i in range(n_pcs)}

    # build the dataframe
    df = pd.DataFrame(dic.items()).rename(columns={0: "PC", 1: "Feature"})
    
    return df


# %%
def pca_2d_scatter(df, index_of_first_feature=6, scale=True, group_col="max_group", print_components=True,
                   main_title="", save_to_out_path="", pca_function=PCA, n_components=2,
                   sample_col="Sample", pca_df_out_path=None, return_compoponents=False, outliers=None,
                   get_most_important_features=False
                  ):
    if outliers:
        df = df.loc[~df[sample_col].isin(outliers)]
    
    X = df.iloc[:, index_of_first_feature:]
    features_names = X.columns
    
    if scale:
        scaler = StandardScaler()
        scaler.fit(X)
        X = scaler.transform(X)
    
    pca = pca_function(n_components=n_components, random_state=1892)  # our lord & savior J.R.R.T was born in 1892
    
    components = pca.fit_transform(X)

    components = pd.DataFrame({sample_col: df[sample_col],
                               group_col: df[group_col],
                               "PC1": components[:, 0],
                               "PC2": components[:, 1]
                              })

    if pca_df_out_path:
        components.to_csv(pca_df_out_path, index=False)
    
    x = "PC1"
    y = "PC2"
    
    if print_components:
        print("Feature importance (rows = PCs, cols = features):\n")
        print(pca.components_, "\n")
    
    if pca_function == PCA:
        cumulative_explained_variance_first_2_pcs = pca.explained_variance_ratio_[:2].sum() * 100
        total_var_str = f"Total explained variance: {cumulative_explained_variance_first_2_pcs:.2f}%"
        main_title = f"{main_title}<br>{total_var_str}" if main_title else total_var_str
        
        if group_col:
            fig = px.scatter(components, x=x, y=y, color=components[group_col], 
                             title=main_title,
#                              hover_data=components[sample_col],
                             hover_name=sample_col,
#                              labels={"color": group_col}
                            )
        else:
            fig = px.scatter(components, x=x, y=t, 
                             title=main_title,
                            )
        fig.update_xaxes(title_text=f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f} %)")
        fig.update_yaxes(title_text=f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f} %)")
    
    else:
        if group_col:
            fig = px.scatter(components, x=x, y=y, color=components[group_col], 
#                              hover_data=components[sample_col],
                             hover_name=sample_col,
#                              labels={"color": group_col}
                            )
        else:
            fig = px.scatter(components, x=x, y=y)
        fig.update_xaxes(title_text=f"PC1")
        fig.update_yaxes(title_text=f"PC2")
    
    if main_title:
        fig.update_layout(title_text=main_title)
    fig.update_layout(title_x=0.5, template="simple_white")
    
    if save_to_out_path:
        fig.write_image(save_to_out_path)
    
#     if return_compoponents:
#         return components, fig
    
    fig.show()
    
    if pca_function == PCA:
        if n_components is None:
            plot_explained_variance(pca)
        if get_most_important_features:
            return most_important_feature_of_each_pc(pca, features_names)


# %% [markdown]
# Run PCA with custom group comrised of certain base variables (e.g. `group_col="Sex-Age"`)

# %%
def custom_group_values(df, group_cols):
    # s = pd.Series([df[group].to_list() for group in group_cols])
    # t = pd.Series([*vars] for vars in zip(*s)).str.join("-")
    # return t
    return pd.Series([*vars] for vars in zip(*pd.Series([df[group].to_list() for group in group_cols]))).str.join("-")


def custom_group_pca(df, index_of_first_feature, group_cols, print_components=False, 
                     main_title=None, pca_df_out_path=None, pca_function=PCA, n_components=2,
                     get_most_important_features=False
                    ):
    new_df = df.copy()
    custom_group_name = "-".join(group_cols)
    if custom_group_name in new_df.columns:
        return pca_2d_scatter(
            new_df,
            index_of_first_feature=index_of_first_feature,
            group_col=custom_group_name,
            print_components=print_components,
            main_title=main_title if main_title else custom_group_name,
            pca_df_out_path=pca_df_out_path,
            pca_function=pca_function,
            n_components=n_components,
            get_most_important_features=get_most_important_features
    )
    else:
        custom_group_vals = custom_group_values(new_df, group_cols)
        new_df.insert(index_of_first_feature, custom_group_name, custom_group_vals)
        return pca_2d_scatter(
            new_df,
            index_of_first_feature=index_of_first_feature+1,
            group_col=custom_group_name,
            print_components=print_components,
            main_title=main_title if main_title else custom_group_name,
            pca_df_out_path=pca_df_out_path,
            pca_function=pca_function,
            n_components=n_components,
            get_most_important_features=get_most_important_features
        )


# %% [markdown]
# Run PCA and color it according to all possible subgroups

# %%
def powerset_pca(df, all_group_cols, index_of_first_feature=6, main_title=None, pca_df_out_path=None,
                pca_function=PCA, n_components=2, get_most_important_features=False):
    most_important_features_got = False
    for group_cols in powerset(all_group_cols):
        if len(group_cols) == 0:  # the empty group
            continue
        elif get_most_important_features and not most_important_features_got:
            most_important_features = custom_group_pca(df, index_of_first_feature, group_cols, 
                                                         main_title=main_title, pca_df_out_path=pca_df_out_path,
                                                         pca_function=pca_function, n_components=n_components,
                                                         get_most_important_features=get_most_important_features
                                                        )
            most_important_features_got = True
        else:
            custom_group_pca(df, index_of_first_feature, group_cols, 
                             main_title=main_title, pca_df_out_path=pca_df_out_path,
                             pca_function=pca_function, n_components=n_components,
                             get_most_important_features=get_most_important_features
                            )
            pca_df_out_path = None # even if it wasn't None, it's enough to save the df for only one coloring 
                                   # as the PCA df is exactly the same
        if get_most_important_features:
            return most_important_features


# %% [markdown] tags=[]
# ## Iris set tests

# %%
df = px.data.iris()
features = ["sepal_width", "sepal_length", "petal_width", "petal_length"]
pca = PCA()
components = pca.fit_transform(df[features])
components.shape

# %%
# how much variance is explained by each PC
pca.explained_variance_ratio_

# %%
# rows = PCs, sorted by decreasing order of explained variance ratio
# cols = features scores in each PC
pca.components_

# %%
abs(pca.components_)

# %%
exp_var_cumul = np.cumsum(pca.explained_variance_ratio_)
exp_var_cumul

# %%
px.area(
    x=range(1, exp_var_cumul.shape[0] + 1),
    y=exp_var_cumul,
    labels={"x": "# Components", "y": "Explained Variance"}
)

# %% [markdown] tags=[]
# # PCA Plots

# %% [markdown]
# ## All regions

# %% tags=[]
df = a2g_wide_merged_df
groups = groups_cols
main_title = "All regions"
# pca_df_out_path = ...
powerset_pca(df, groups, main_title=main_title)

# %% [markdown] tags=[]
# ## All regions, sparse PCA

# %% tags=[]
pca_2d_scatter(a2g_wide_merged_df, 
               index_of_first_feature=6,
               group_col="Sex",
               print_components=False,
               pca_function=SparsePCA,
               n_components=2
              )

# %% tags=[]
pca_2d_scatter(a2g_wide_merged_df, 
               index_of_first_feature=6,
               group_col="Sex",
               print_components=False,
               pca_function=SparsePCA,
               n_components=2,
               outliers=["sample18", "sample26"]
              )

# %% [markdown] tags=[]
# ## All regions, kernel PCA

# %% tags=[]
df = a2g_wide_merged_df
# groups = groups_cols
groups = ["Sex", "Diet"]
main_title = "All regions, kernel PCA"
pca_function = KernelPCA
powerset_pca(df, groups, main_title=main_title, pca_function=KernelPCA)

# %% [markdown]
# ## All regions, Sex-Diet subgroups

# %%
sex_diet_subgroups = [("Male", "Full"), 
                      ("Male", "Fasted"), 
                      ("Female", "Fasted")]
index_of_first_feature = 6
group_cols = ["Age(weeks)", "Genotype"]
for sex, diet in sex_diet_subgroups:
    df = a2g_wide_merged_df.loc[(a2g_wide_merged_df["Sex"]==sex) & 
                                (a2g_wide_merged_df["Diet"]==diet)]
    df = df.reset_index(drop=True)
    main_title = f"{sex}-{diet}, all regions"
    custom_group_pca(df, 
                     index_of_first_feature, 
                     group_cols, 
                     main_title=main_title
                    )

# %% [markdown] tags=[]
# ## Non-zero regions

# %% tags=[]
df = a2g_wide_no_zeroes_merged_df
groups = groups_cols
main_title = "Non-zero regions"
# pca_df_out_path = ...
powerset_pca(df, groups, main_title=main_title)

# %% [markdown]
# ## Non-zero regions, Sex-Diet subgroups

# %%
sex_diet_subgroups = [("Male", "Full"), 
                      ("Male", "Fasted"), 
                      ("Female", "Fasted")]
index_of_first_feature = 6
group_cols = ["Age(weeks)", "Genotype"]
for sex, diet in sex_diet_subgroups:
    df = a2g_wide_no_zeroes_merged_df.loc[(a2g_wide_no_zeroes_merged_df["Sex"]==sex) & 
                                          (a2g_wide_no_zeroes_merged_df["Diet"]==diet)]
    df = df.reset_index(drop=True)
    main_title = f"{sex}-{diet}, non-zero regions"
    custom_group_pca(df, 
                     index_of_first_feature, 
                     group_cols, 
                     main_title=main_title
                    )

# %% [markdown]
# # Intersections

# %% [markdown] tags=[]
# ## Intersecting gff with indexed regions

# %%
gff_bed = BedTool(gff_file)
gff_bed.head()

# %%
# original HE UE regions
ue_regions_bed6 = BedTool(ue_regions_file)
if MISMATCH in ["C2G", "C2T", "G2T", "T2G", "T2C", "T2A"]: # done for technical reasons (legacy code)
    ue_regions_bed6 = flip_strand(ue_regions_bed6)
ue_regions_bed6.head()

# %%
# regions used for the Index (not all of the original UE regions were included due to conflicting strands)
a2g_indexed_bed3 = BedTool.from_dataframe(a2g_merged_df.loc[:, ["GenomicRegion", "Start", "End"]]).sort()
a2g_indexed_bed3.head()

# %%
# # regions used for the Index are strandless, 
# # so we intersect them with their original ancestors to restore the strand information
# a2g_indexed_bed6 = ue_regions_bed6.intersect(a2g_indexed_bed3, u=True).sort()
# a2g_indexed_bed6.head()

# regions used for the Index are strandless, 
# so we intersect them with their original ancestors to restore the strand information.
# since the indexed regions may have been extanded, we use -wa and -wb for that purpose.
# finally, we filter duplicates. the duplicates arrise from the use of -wb which yields intersection with different
# features in ue_regions_bed6; however, their strads can't conflict and so we may rely on the first intsance.
a2g_indexed_bed6 =  a2g_indexed_bed3.intersect(ue_regions_bed6, wa=True, wb=True).sort()
from pybedtools import Interval
features = []
for feature in a2g_indexed_bed6:
    new_feature = Interval(feature[0], int(feature[1]), int(feature[2]), strand=feature[8])
    features.append(new_feature)
a2g_indexed_bed6 = BedTool(features).merge(s=True, c=[4, 5, 6], o="first")
a2g_indexed_bed6.head()

# %%
# gff with indexed regions
# indexed_gff_bed = gff_bed.intersect(a2g_indexed_bed6, u=True)
indexed_gff_bed = gff_bed.intersect(a2g_indexed_bed6, s=True)
print(len(indexed_gff_bed), "\n")
indexed_gff_bed.head()

# %%
set(f[2] for f in indexed_gff_bed)

# %%
indexed_genes_gff_bed  = indexed_gff_bed.filter(lambda f: f[2]=="gene").saveas()
print(len(indexed_genes_gff_bed), "\n")
indexed_genes_gff_bed.head()

# %%
indexed_trna_lnc_gff_bed = indexed_gff_bed.filter(lambda f: f[2] in ["tRNA", "lnc_RNA"]).saveas()
print(len(indexed_trna_lnc_gff_bed), "\n")
indexed_trna_lnc_gff_bed.head()

# %% [markdown]
# ## Intersecting gff with PCA's most-important features

# %% [markdown]
# ### Non-zero regions PCA

# %% tags=[]
df = a2g_wide_no_zeroes_merged_df
groups = ["Sex"]
main_title = "Non-zero regions"
most_import_features_df = powerset_pca(df, 
                                       groups, 
                                       main_title=main_title, 
                                       n_components=None, 
                                       get_most_important_features=True
                                      )
most_import_features_df.head()

# %%
unique_features_df = (most_import_features_df.
                      loc[:, "Feature"].
                      str[:-4].
                      str.rsplit("_", n=2, expand=True).
                      drop_duplicates().
                      reset_index(drop=True)
                     )
unique_features_df.insert(3, "name", ".")
unique_features_df.insert(4, "score", unique_features_df.index + 1)  # corresponds to PC, e.g. score=1 <--> PC=1
unique_features_df.insert(5, "strand", ".")
unique_features_df

# %%
most_import_features_bed6 = BedTool.from_dataframe(unique_features_df)
most_import_features_bed6.head()

# %%
wa_names = [f"{f}_feature" for f in most_import_features_bed6.to_dataframe().columns.to_list()]
wb_names = [f"{f}_gff" for f in indexed_gff_bed.to_dataframe().columns.to_list()]
wb_names

# %%
df = most_import_features_bed6.intersect(indexed_gff_bed, wa=True, wb=True).to_dataframe(names=wa_names+wb_names)

df

# %%
no_zeroes_merged_df.loc[no_zeroes_merged_df["Mismatch"]==MISMATCH].head()

# %%
# outpath = ...

# (df.merge(no_zeroes_merged_df.loc[no_zeroes_merged_df["Mismatch"]=="A2G"],
#          left_on=["chrom_feature", "start_feature", "end_feature"],
#          right_on=["GenomicRegion", "Start", "End"]).
#       to_csv(outpath, index=False)
# )

# %%
df2 = (df.
         loc[df["feature_gff"]=="gene", ["score_feature", "attributes_gff"]].
         reset_index(drop=True)
        )
df2["gene"] = (df2["attributes_gff"].str.split(";", n=1, expand=True).iloc[:, 0])
df2 = df2.drop(columns=["attributes_gff"], axis=1).rename(columns={"score_feature": "PC"})

df2
