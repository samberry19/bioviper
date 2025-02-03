import pandas as pd
import numpy as np

def selector(df, dic):

    '''A little tool for selecting from pandas dataframes by passing a dictionary, e.g.
            selector(df, {"color":"red", "shape":["square", "circle"]})
       
        For advanced usage, you can pass a function and it will return where True, e.g.
            selector(df, ["name": lambda name: "Sam" in name])
           
        You can also use this to select things greater than or less than a value, e.g.
            selector(df, ["enrichment": lambda enr: enr > 1])'''
   
    X = df.copy()

    for key,val in dic.items():
       
        # If you pass a tuple, list or numpy array
        if isinstance(val, (tuple, list, np.ndarray)):
            where = np.any(np.array([X[key]==v for v in val]),axis=0)
            X = X.loc[where]
           
        # If you pass a function
        elif isinstance(val, type(lambda x: x+1)):
            X = X.loc[X[key].apply(val)]
            
        elif isinstance(val, str) and '>' in val:
            
            if val.count('>')==2:
                X = X.loc[(X[key]<float(val.split('>')[0])) & (X[key]>float(val.split('>')[2]))]
            else:
                X = X.loc[X[key]>float(val.split('>')[1])]
                      
        elif isinstance(val, str) and '<' in val:
            
            if val.count('<')==2:
                X = X.loc[(X[key]>float(val.split('<')[0])) & (X[key]<float(val.split('<')[2]))]
            else:
                X = X.loc[X[key]<float(val.split('<')[1])]
           
        # Otherwise we assume it's a single value
        else:
            X = X.loc[X[key]==val]

    return X