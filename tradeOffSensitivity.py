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
    lastRow = P.index[P["DesignOptions"] == "Definitions"].to_list()[0] - 1
    P = P.drop(P.columns[idx + 1 :], axis=1).drop(P.index[lastRow:], axis=0)

    initWeight = P.loc[[0]].drop(["DesignOptions", "Score"], axis=1).to_dict("list")

    #%%
    for key in initWeight:
        value = float(initWeight[key][0])
        value_dict[key] = [value + percentage for percentage in range_list]

    P = P.drop(P.index[0], axis=0).reset_index(drop=True)

    test = P["DesignOptions"].to_dict()
    for key in test.values():
        winnerDict[key] = 0

    keys, values = zip(*value_dict.items())
    combi_dict = [
        {key: value for (key, value) in (zip(keys, v))} for v in it.product(*values)
    ]

    combi_dict = [x for x in combi_dict if np.isclose(sum(x.values()),1)]

    #%%
    for key in tqdm(combi_dict):
        values = pd.DataFrame(key, index=[0]).astype(float)
        values = list(key.values())
        Q = P.copy().drop(["DesignOptions", "Score"], axis=1).astype(float)
        Q = (
            Q.apply(lambda x: x * values, axis=1)
            .reset_index(drop=True)
            .div(sum(values))
        )
        idx = Q.sum(axis=1).idxmax()
        winnerDict[test[idx]] += 1
    
    total = sum(winnerDict.values())
    for key in winnerDict:
        print(f"{key} & {round((winnerDict[key]/total) * 100, 2)} \\ \hline %")

# %%
