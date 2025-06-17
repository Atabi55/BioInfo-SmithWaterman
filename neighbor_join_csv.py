import pandas as pd
import numpy as np
from ete3 import Tree



def neighbor_join():
    #read csv file for input
    df = pd.read_csv( 'taxa.dat', index_col=0, header=0, sep=r'\s+' )

    #convert "-" into NaN
    df = df.map(lambda x: pd.NA if x == "-" else x)
    #df = df.map(lambda x: np.NaN if x == "-" else x)

    #merge dataframe with its mirror to fill the other half
    df_mirrored = df.T.combine_first(df)

    #fill diagonal with 0
    for i in range(len(df_mirrored)):
        df_mirrored.iat[i, i] = 0

    df_mirrored = df_mirrored.apply(pd.to_numeric, errors='coerce')
    #set datatype of whole dataframe to Float64
    df_mirrored = df_mirrored.astype("Float64")


    tree_string = "" #will store the final result tree in string form

    if len(df_mirrored) == 2:
        table = df_mirrored
        tree_string = f"({table.columns[0]}:{table.iloc[0][1]/2},{table.columns[1]}:{table.iloc[0][1]/2});"
        return tree_string

    i = len(df_mirrored) - 2
    while i:
        ri = df_mirrored.sum(axis=1).tolist()
        ri = [x / (len(df_mirrored) - 2) for x in ri]
        print(ri)

        #get indices of upper triangle elements in dataframe
        rows, cols = np.triu_indices(len(df_mirrored), k=1)

        for r, c in zip(rows, cols):
            Q = df_mirrored
            Q.iat[r, c] = Q.iat[r, c] - ri[r] - ri[c]

        print(Q)

        n = len(Q)
        rows, cols = np.tril_indices(n, k=-1)

        # Set those positions to NaN
        for r, c in zip(rows, cols):
            Q.iat[r, c] = np.nan


        # merge dataframe with its mirror to fill the other half
        Q = Q.T.combine_first(Q)
        print(Q)
        break #debug, to stop loop


        #now we have



        i = i - 1


    """
    t2 = Tree('(A:1,(B:1,(C:1,D:1):0.5):0.5);')
    t2 = Tree('((((C:11,D:17):7.25),((B:6.75,E:14.25):4.75)),A:4.75);')
    t2.show()
    """

    return tree_string


def draw_tree():
    #inputs: list of species, list of distance between species
    pass


if __name__ == "__main__":
    tree_tuple_string = neighbor_join()
    print(tree_tuple_string)
    t3 = Tree(tree_tuple_string)
    t3.show()