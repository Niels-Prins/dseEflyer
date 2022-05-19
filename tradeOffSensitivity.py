#%%
# Import libraries
import pandas as pd
import numpy as np
from math import *
import itertools as it
from tqdm import tqdm

range_list = np.arange(-0.07, 0.08, 0.01)

#%%
# Loop through all sheets
for i in range(2, 9):
    combi_dict = {}
    value_dict = {}
    winnerDict = {}
    # Import Excel file
    P = pd.read_excel(
        "C:/Users/niels/OneDrive/E-Flyer/3_Midterm/Trade-off table.xlsx", sheet_name=i
    )

    P = P.drop(P.index[1], axis=0).rename(columns={"Unnamed: 0": "DesignOptions"})
    idx = P.columns.get_loc("Score")
    lastRow = P.index[P["DesignOptions"] == "Definitions"].to_list()[0] - 2
    P = P.drop(P.columns[idx + 1 :], axis=1).drop(P.index[lastRow:], axis=0)

    initWeight = P.loc[[0]].drop(["DesignOptions", "Score"], axis=1).to_dict("list")
    for key in initWeight:
        initWeight[key] = initWeight[key][0]
        value = float(initWeight[key])
        value_dict[key] = [value + percentage for percentage in range_list]

    P = P.drop(P.index[0], axis=0).reset_index(drop=True)

    designOptionsDict = P["DesignOptions"].to_dict()
    for key in designOptionsDict.values():
        scores_weights[key] = 0
        scores_dropout[key] = 0

    keys, values = zip(*value_dict.items())
    combi_dict = [
        {key: value for (key, value) in (zip(keys, v))} for v in it.product(*values)
    ]

    combi_dict = [x for x in combi_dict if np.isclose(sum(x.values()), 1)]

    #%%
    for weightcombi in tqdm(combi_dict):
        scores_weights = score_calc(weightcombi, scores_weights, designOptionsDict)

    total_combi_weight = sum(scores_weights.values())

    #Criterion Dropout
    for key in initWeight:
        temp = initWeight.copy()
        temp[key] = 0
        scores_dropout = score_calc(temp, scores_dropout, designOptionsDict)
        
    total_dropout = sum(scores_dropout.values())

    for key in scores_weights:
        print(
            f"{key} & {round((scores_weights[key]/total_combi_weight) * 100, 2)}\% & {round((scores_dropout[key]/total_dropout) * 100, 2)}\% \\\ \hline"
        )

 # %%
