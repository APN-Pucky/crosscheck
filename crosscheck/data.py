import pandas as pd
import matplotlib.pyplot as plt

def rel_diff(comp_df,*val_dfs, column=None,plot=True,figsize=(20,5),logy=True,labels=None):
    """
    Compute the relative difference between two dataframes
    """
    # TODO check that momenta are the same and valid
    rets = []
    for i,val_df in enumerate(val_dfs):
        ret = abs(((comp_df[column]-val_df[column])/comp_df[column]))
        if plot:
            ret.plot(figsize=figsize,logy=logy,label=(None if labels is None else ("rel. " +labels[i+1]+ " vs. " + labels[0])))
            if labels is not None:
                plt.legend()
        rets.append(ret)
    return rets

def ratio(comp_df,*val_dfs, column=None,plot=True,figsize=(20,5),logy=True,labels=None):
    """
    Compute the relative difference between two dataframes
    """
    # TODO check that momenta are the same and valid
    rets = []
    for i,val_df in enumerate(val_dfs):
        ret = ((val_df[column])/comp_df[column])
        if plot:
            ret.plot(figsize=figsize,logy=logy,label=(None if labels is None else (labels[i+1]+ "/" + labels[0])))
            if labels is not None:
                plt.legend()
        rets.append(ret)
    return rets
