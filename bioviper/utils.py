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

def get_variant(reference, mutations, start_index=1):

    '''Take in a reference sequence and a list of mutations and return the variant sequence'''

    if mutations[0]=='WT':
        return reference.copy()
    
    else:
        variant = reference.copy()
        for mut in mutations:
            pos = int(mut[1:-1])-start_index
            variant[pos] = mut[-1]
        return variant

def call_variant(aa_seq, wt, start_index=1):
    
    '''Take in a protein variant and compare it to a reference sequence'''
    
    # if the input is a string, convert to a numpy array
    if isinstance(aa_seq, str):
        aa_seq = np.array(list(aa_seq))
    if isinstance(wt, str):
        wt = np.array(list(wt))
    
    mut_pos = np.where(aa_seq!=wt)[0]
        
    if len(mut_pos)==0:
        return 'WT'
    else:
        return ','.join([wt[i]+str(i+start_index)+aa_seq[i] for i in mut_pos])