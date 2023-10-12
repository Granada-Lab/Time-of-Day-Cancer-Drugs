import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import scipy.stats as ss
from sklearn.model_selection import train_test_split
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_score
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis


def reformat_drugs(dataframe):
    """reshapes a CCLE pandas dataframe from 'one line per datapoint' to a more convenient
    'one line per sample' format, meaning the response of a given cell line to different drugs
    will be placed on the same line in different columns."""

    drug_names = dataframe["Compound"].unique()
    dataframe["Compound"].value_counts()
    # concatenate the drug info with one line per cell line
    merged = pd.DataFrame()
    for thisDrug in drug_names:
        dataframe_spec = dataframe.loc[dataframe["Compound"] == thisDrug]
        dataframe_spec_clean = dataframe_spec.drop(
            columns=[
                "Primary Cell Line Name",
                "Compound",
                "Target",
                "Activity SD",
                "Num Data",
                "FitType",
            ]
        )
        dataframe_spec_clean.columns = [
            "CCLE Cell Line Name",
            thisDrug + "_dr_doses",
            thisDrug + "_dr_responses",
            thisDrug + "_EC50",
            thisDrug + "_IC50",
            thisDrug + "_Amax",
            thisDrug + "_ActArea",
        ]

        if merged.empty:
            merged = dataframe_spec_clean.copy()
        else:
            merged = pd.merge(
                merged,
                dataframe_spec_clean,
                how="left",
                on="CCLE Cell Line Name",
                sort=False,
                suffixes=("_x", "_y"),
                copy=True,
            )
    merged_dataframe = merged.set_index("CCLE Cell Line Name")

    actarea_cols = [x for x in merged_dataframe.columns if 'ActArea' in x]

    return merged_dataframe[actarea_cols]


def remove_columns_with_zeros(df, max_zeros):
    """
    Removes columns that contain more than a certain number of zeros
    Made by ChatGPT
    :param df: some DataFrame
    :param max_zeros: an arbitrary threshold
    :return: a much cleaner DataFrame
    """
    # Calculating the number of zeros in each column
    num_zeros = (df == 0).sum(axis=0)
    # Selecting columns with less than or equal to max_zeros zeros
    df = df.loc[:, num_zeros <= max_zeros]

    return df


def binarize_with_kmeans(df, n_clusters=2, keep_order=True, random_state=None):
    """
    Binarize a DataFrame using k-means clustering independently for each column.
    df : pandas DataFrame
    n_clusters : int, optional
        The number of clusters to use for k-means clustering. Default is 2.

    Returns:
    --------
    pandas DataFrame
        The binarized DataFrame with the same shape as the input DataFrame.
    """
    df_binary = pd.DataFrame(index=df.index)

    for column in df.columns:
        kmeans = KMeans(n_clusters=n_clusters, random_state=random_state)
        labels = kmeans.fit_predict(df[column][df[column].notnull()].values.reshape(-1, 1))
        centroids = kmeans.cluster_centers_
        is_inverted = centroids[0] > centroids[1]
        if is_inverted and keep_order:
            labels = np.logical_not(labels).astype(int)
        if df[column].isnull().values.any():
            # suppress SettingWithCopyWarning warning
            pd.options.mode.chained_assignment = None
            nan_mask = df[column].isna()
            df_binary[column] = nan_mask
            # put nan back where they were
            df_binary[column][nan_mask] = np.nan
            # on inverted nan location, add the labels
            df_binary[column][~nan_mask] = labels
        else:
            df_binary[column] = labels
    return df_binary


# Corrected Cramer's V correlation between categorical features
def cramers_corrected_stat(x, y):
    """
    Function to calculate corrected Cramers V statistic for categorical-categorical association. Uses correction
    from Bergsma and Wicher, Journal of the Korean Statistical Society 42 (2013): 323-328.

    Parameters
    ----------
    x : np.array
        array of first vector or column to analyze Cramers V correlation with the second
    y : np.array
        array of second vector or column to analyze Cramers V correlation with the first

    Returns
    -------
    result : float
        float value of the corrected Cramers V correlation coefficient between x and y
    """
    result = -1
    if len(np.unique(x)) == 1:
        print("First variable is constant")
    elif len(np.unique(y)) == 1:
        print("Target variable is constant")
    else:
        conf_matrix = pd.crosstab(x, y)
        if conf_matrix.shape[0] == 2:
            correct = False
        else:
            correct = True
        chi_2 = ss.chi2_contingency(conf_matrix, correction=correct)[0]
        n = sum(conf_matrix.sum())
        phi2 = chi_2 / n
        r, k = conf_matrix.shape
        phi2corr = max(0, phi2 - ((k - 1) * (r - 1)) / (n - 1))
        r_corr = r - ((r - 1) ** 2) / (n - 1)
        k_corr = k - ((k - 1) ** 2) / (n - 1)
        result = np.sqrt(phi2corr / min((k_corr - 1), (r_corr - 1)))
    return round(result, 6)


def LDA_loocv(data=None, y=None, target=None):
    """
    Function to perform leave one out cross validation using LDA prediction.

    Parameters
    ----------
    data : pandas.core.frame.DataFrame
        Feature matrix
    y : pandas.core.frame.DataFrame
        Target matrix
    target : string
        LDA-specific output target to control on what feature the supervised LDA is fitted, as it would not work for
        continuous features

    Returns
    -------
    acc : float
        float value rounded to 3 decimals, accuracy of the leave one out LDA cross-validation
    acc_sd : float
        standard deviation of accuracy during loocv
    bsl : float
        brier score loss
    bsl_sd : float
        precision recall standard deviation
    """
    # define cross validation method
    cv = LeaveOneOut()
    # initialize LDA model
    comp = len(y[target].astype(int).unique()) - 1
    LDA = LinearDiscriminantAnalysis(n_components=comp)
    # Use loocv to evaluate model
    # return the accuracy score
    scores = cross_val_score(estimator=LDA, X=data, y=y[target].astype(int),
                             scoring='accuracy', cv=cv, n_jobs=-1)
    acc = np.mean(scores)
    acc_sd = np.std(scores)
    # return brier score loss
    scores = cross_val_score(estimator=LDA, X=data, y=y[target].astype(int),
                             scoring='neg_brier_score', cv=cv, n_jobs=-1)
    bsl = np.mean(scores)
    bsl_sd = np.std(scores)
    return acc, acc_sd, bsl, bsl_sd
