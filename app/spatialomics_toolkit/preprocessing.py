import pandas as pd
import numpy as np

"""
Author: Mattis Knulst
Description: This module contains functions that are useful in 
pre-processing tables containing data from GeoMX

"""


def clean_import(filehandle, sep='\t', select_name=None, select_index=None):
    """
    Imports a T/CSV and cleans column names from problematic symbols.
    :param filehandle: path to file containing a table
    :param sep: table separation symbol, default=TAB '\t'
    :param select_name: list of column names to include
    :param select_index: tuple of indeces for slicing columns based on index
    :return: pandas dataframe with formatted column names
    """
    raw = pd.read_csv(filehandle, sep=sep)
    # clean column names
    raw.columns = raw.columns.str.replace(r'[#,@,&,\s,\,,\.,-]', '', regex=True)
    if select_index and select_name:
        select_name.extend(raw.columns.tolist()[select_index[0]:select_index[1]])
        raw = raw.loc[:, select_name]
    elif select_name:
        raw = raw.loc[:, select_name]
    elif select_index:
        raw = raw.iloc[:, select_index[0]:select_index[1]]
    return raw


def scrub_data(df, controls=None, remove_na=False,
               reset_index=True, na_subset=None, bin_col=None):
    """
    Scrubs a dataframe and optionally removes controls, NA,
    and/or resets the index of the dataframe
    :param df: pandas dataframe
    :param controls: list of column names to drop, default: do nothing
    :param remove_na: remove NA from column if na_subset is string or list or
     from entire df, default is to not remove NA
    :param reset_index: default True (makes df forget about what was removed)
    :param na_subset: column name or list of columns to look for NA in
    :param bin_col: list of columns that should be reformatted to sparse
    :return: neatly scrubbed pandas dataframe
    """
    if controls:
        for c in controls:
            try:
                df.drop(c, axis=1, inplace=True)
            except:
                print(f'Could not drop {c} check spelling!')
    if remove_na:
        if na_subset:
            df.dropna(subset=na_subset, inplace=True)
        else:
            df.dropna(inplace=True)
    if bin_col:
        for b in bin_col:
            df=pd.get_dummies(data=df, columns=b, drop_first=True)
    if reset_index:
        df.reset_index(drop=True, inplace=True)
    return df



