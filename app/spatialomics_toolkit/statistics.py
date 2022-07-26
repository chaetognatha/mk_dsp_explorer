import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection as fdr
from lifelines import CoxPHFitter
import multiprocessing as mp

def split_df(df, a, b, by='colname'):
    """
    returns two dataframes, either by slicing on column names a to b or
    column indeces a to b, a or b can be empty to indicate that we want to
    slice from beginning of dataframe or until the end of dataframe
    :param df: pandas dataframe
    :param by: colname or index
    :param a: first column to extract
    :param b: last column to extract
    :return: tuple of pandas data frames
    """
    if by=='index':
        df1 = df.iloc[a:b]
    elif by == 'colname':
        df1 = df.loc[a:b]
    else: print('by param accepts colname or index!')
    df2 = df.drop(df1.columns.tolist(), axis=1)
    return df1, df2

def univariate_cox_regression(df, biomarker, time_col, event_col, patient_id_col):
    """
    perform cox regression on a single biomarker
    :param df: pandas dataframe
    :param biomarker: pandas series / column
    :param time_col: see survival
    :param event_col: see survival
    :param patient_id_col: see survival
    :return: coefficient and p-value if significant
    """
    results = survival(df, time_col, event_col, patient_id_col)


def lmm(df, biomarkers, FE, RE, fcol='', filter=''):
    if fcol and filter:
        df = df[list(df[fcol] == filter)]
    # model fitting loop
    p_values = []
    coefficients = []
    names = []
    for b in biomarkers:
        model = smf.mixedlm(f"Q('{b}') ~ {FE}", data=df,
                            groups=df[RE])
        fit = model.fit()
        coefficients.append(fit.params.values[1])
        p_values.append(fit.pvalues.values[1])
        names.append(b)

    results = pd.DataFrame({"biomarker": names,
                            "coefficients": coefficients,
                            "p-values": p_values,
                            "FDR": fdr(p_values, alpha=0.05)[1]})
    results['neglog10p'] = - np.log10(results['FDR'])
    results['neglog10p-unadj'] = - np.log10(results['p-values'])
    return results

def survival(df, time, status, patient_id, a=0.05, l1_ratio=0.9):
    """
    Function for returning a lifelines clustered Cox regression model
    :param df: pandas dataframe (beware of collinearity)
    :param time: time
    :param status: status 1/0 or True/False
    :param patient_id: this column will be weighted to control for repeated measures
    :param a: alpha - constant that multiplies the penalty 
    :param l1_ratio: ratio of l1 to l2 penalty, if it is 1 it is a Lasso
    :return: fitted model object that can be used to print summary and plot
    """
    cph = CoxPHFitter(alpha=a, l1_ratio=l1_ratio)
    cph.fit(df,
            time,
            event_col=status,
            cluster_col=patient_id
            )
    return cph

