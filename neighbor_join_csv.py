import pandas as pd
import numpy as np
from ete3 import Tree, TreeStyle



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

    # Initialize labels for Newick string
    labels = {name: name for name in df_mirrored.columns}

    D = df_mirrored.copy()

    while len(D) > 2:
        n = len(D)
        # Calculate r values
        r = D.sum(axis=1) / (n - 2)

        # Compute Q matrix
        Q = pd.DataFrame(np.zeros(D.shape), index=D.index, columns=D.columns)
        for i in D.index:
            for j in D.columns:
                if i != j:
                    Q.loc[i, j] = D.loc[i, j] - r[i] - r[j]
                else:
                    Q.loc[i, j] = np.nan

        # Find pair with minimum Q value
        min_idx = Q.stack().idxmin()
        i, j = min_idx

        # Calculate branch lengths
        dij = D.loc[i, j]
        li = 0.5 * dij + 0.5 * (r[i] - r[j])
        lj = dij - li

        # Create new label in Newick format
        new_label = f"({labels[i]}:{li:.5f},{labels[j]}:{lj:.5f})"

        # Update distance matrix
        new_row = {}
        for k in D.index:
            if k not in [i, j]:
                dik = D.loc[i, k]
                djk = D.loc[j, k]
                new_row[k] = 0.5 * (dik + djk - dij)

        # Remove i and j, add new node
        D = D.drop(index=[i, j], columns=[i, j])
        D.loc[new_label] = pd.Series(new_row)
        D[new_label] = pd.Series(new_row)
        D.loc[new_label, new_label] = 0

        # Update labels
        labels[new_label] = new_label
        labels.pop(i)
        labels.pop(j)

    # Final join
    i, j = D.index
    dij = D.loc[i, j]
    newick = f"({labels[i]}:{dij/2:.5f},{labels[j]}:{dij/2:.5f})"
    return newick



if __name__ == "__main__":
    tree_tuple_string = neighbor_join() + ";"
    print(tree_tuple_string)#prints the tree in Newick format
    t3 = Tree(tree_tuple_string)
    ts = TreeStyle()
    ts.show_branch_length = True
    t3.render("output_tree.png", w=600, units="px", tree_style=ts)