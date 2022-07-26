import seaborn as sns
from matplotlib import pyplot as plt
import scipy
import scipy.cluster.hierarchy as sch
import adjustText
import numpy as np
import pandas as pd

def cluster_corr(corr_array, inplace=False):
    """
    Rearranges the correlation matrix, corr_array, so that groups of highly
    correlated variables are next to eachother

    Parameters
    ----------
    corr_array : pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix

    Returns
    -------
    pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix with the columns and rows rearranged
    """
    pairwise_distances = sch.distance.pdist(corr_array)
    linkage = sch.linkage(pairwise_distances, method='complete')
    cluster_distance_threshold = pairwise_distances.max() / 2
    idx_to_cluster_array = sch.fcluster(linkage, cluster_distance_threshold,
                                        criterion='distance')
    idx = np.argsort(idx_to_cluster_array)


def quick_volcano_plot(results):
    plt.figure(figsize=(10, 5))
    # make scatterplot
    p1 = sns.scatterplot(data=results, x="coefficients", y="neglog10p",
                         hue="neglog10p", size="neglog10p",
                         sizes=(20, 200), legend=False, palette="magma")
    # label data points
    texts = [plt.text(results.coefficients[i], results.neglog10p[i],
                      results.biomarker[i]) for i in \
             range(len(results.coefficients)) \
             if results.neglog10p[i] > -np.log10(0.05)]
    adjustText.adjust_text(texts)

    p1.axhline(y=-1 * np.log10(0.05), ls="--", color="red")
    p1.axvline(x=0, color="black")
    plt.title('Volcano plot')
    # Set x-axis label
    plt.xlabel('Coefficients')
    # Set y-axis label
    plt.ylabel('-10log(FDR)')
    plt.show()

def volcano_plot(results, color='coolwarm', y="neglog10p", x="coefficients",
                 p_lim=0.05, plot_title='Volcano Plot', x_label='Coefficients',
                 y_label='-10logFDR'):
    # coolwarm, magma
    plt.figure(figsize=(10, 5))
    # make scatterplot
    p1 = sns.scatterplot(data=results, x=x, y=y,
                         hue=y, size=y,
                         sizes=(20, 200), legend=False, palette=color)
    # label data points
    texts = [plt.text(results.coefficients[i], results[y][i] ,
                      results.biomarker[i]) for i in \
             range(len(results.coefficients)) \
             if results[y][i] > -1 * np.log10(p_lim)]
    adjustText.adjust_text(texts)

    p1.axhline(y=-1 * np.log10(p_lim), ls="--", color="red")
    p1.axvline(x=0, color="black")
    plt.title(plot_title)
    # Set x-axis label
    plt.xlabel(x_label)
    # Set y-axis label
    plt.ylabel(y_label)
    plt.show()