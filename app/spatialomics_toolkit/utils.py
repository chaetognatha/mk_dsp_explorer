from typing import Dict, Tuple

from dash import Dash, html, dash_table
import pandas as pd
import plotly.express as px


def dataframe_to_dash_datable(
    df: pd.DataFrame,
    separator: str = "_",
    properties: Dict = {},
    column_properties: Dict = {},
) -> Tuple:
    column_dicts, column_renaming = [], []
    for column in df.columns:
        if isinstance(column, (list, tuple)):
            column_str = [str(x) for x in column if pd.notnull(x)]
            if "" in column_str:
                n_empty = column_str.count("")
                column_str.remove("")
                column_str = [""]*n_empty + column_str
            column_id = f"{separator}".join([x for x in column_str if x != ""])
        else:
            column_id = str(column)
        column_renaming.append(column_id)
        column_dict = {"id": column_id, "name": column_str, **properties}
        if column_id in column_properties:
            column_dict.update(column_properties[column_id])
        column_dicts.append(column_dict)
    df.columns = column_renaming
    df = df.to_dict("records")
    return df, column_dicts