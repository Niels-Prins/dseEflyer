import pandas as pd
import numpy as np

A = [1.00, 10.10]
B = [2.00, 20.20]

test = pd.DataFrame([A,B], columns = ["A", "B"])
print(test)
print(test.loc[test["B"] == 20.2, "A"].iloc[0])