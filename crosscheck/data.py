import pandas as pd
import matplotlib.pyplot as plt

def rel_diff(val_df, comp_df,column=None,plot=True,figsize=(20,5),logy=True):
    """
    Compute the relative difference between two dataframes
    """
    # TODO check that momenta are the same and valid
    ret = abs(((comp_df[column]-val_df[column])/comp_df[column]))
    if plot:
        ret.plot(figsize=figsize,logy=logy) # LO
        plt.show()
    return ret

