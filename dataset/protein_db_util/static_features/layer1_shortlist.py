import pandas as pd
import numpy as np
from io import StringIO

csv_data = """ProteinID,PocketID,Ligand,DockingScore,PocketProbability,CenterX,CenterY,CenterZ,SASpoints,GridEdgeLength,ABCadj,Corrected_PLIP_contacts,InteractionTypes,Domain,Protein_Pocket_Ligand
CIMG_09093,4,CIMG_09093_prepared_Beauvericin,-10.0,0.444,-5.6583,-42.1991,-20.4558,64,28.0,0.0,12,HydrogenBond;Hydrophobic,NONE,CIMG_09093_4_CIMG_09093_prepared_Beauvericin
CIMG_09093,3,CIMG_09093_prepared_Beauvericin,-8.9,0.831,6.2066,21.3999,-3.143,142,28.0,0.0,5,HydrogenBond;Hydrophobic,NONE,CIMG_09093_3_CIMG_09093_prepared_Beauvericin
CIMG_09093,4,CIMG_09093_prepared_Tacrolimus,-9.4,0.444,-5.6583,-42.1991,-20.4558,64,28.0,0.0,11,HydrogenBond;Hydrophobic,NONE,CIMG_09093_4_CIMG_09093_prepared_Tacrolimus
CIMG_09093,3,CIMG_09093_prepared_Tacrolimus,-8.7,0.831,6.2066,21.3999,-3.143,142,28.0,0.0,9,HydrogenBond;Hydrophobic,NONE,CIMG_09093_3_CIMG_09093_prepared_Tacrolimus
CIMG_09093,4,CIMG_09093_prepared_Milbemycin,-9.3,0.444,-5.6583,-42.1991,-20.4558,64,24.0,0.0,9,HydrogenBond;Hydrophobic,NONE,CIMG_09093_4_CIMG_09093_prepared_Milbemycin
CIMG_09093,3,CIMG_09093_prepared_Milbemycin,-9.6,0.831,6.2066,21.3999,-3.143,142,28.0,0.0,5,HydrogenBond;Hydrophobic;SaltBridge,NONE,CIMG_09093_3_CIMG_09093_prepared_Milbemycin
CIMG_00780,5,CIMG_00780_prepared_Milbemycin,-10.4,0.366,4.605,1.5681,-12.1955,123,28.0,0.0,5,HydrogenBond;Hydrophobic,NONE,CIMG_00780_5_CIMG_00780_prepared_Milbemycin
CIMG_00780,1,CIMG_00780_prepared_Milbemycin,-9.6,0.913,5.2501,-35.776,-19.6413,153,30.0,0.05,7,Hydrophobic,TMD,CIMG_00780_1_CIMG_00780_prepared_Milbemycin
CIMG_00780,4,CIMG_00780_prepared_Milbemycin,-10.0,0.391,7.3607,-20.5703,-4.3144,84,26.0,0.0,8,HydrogenBond;Hydrophobic,NONE,CIMG_00780_4_CIMG_00780_prepared_Milbemycin
CIMG_00780,3,CIMG_00780_prepared_Milbemycin,-10.6,0.424,-0.5827,28.9579,2.6713,100,26.0,0.0,4,HydrogenBond;Hydrophobic,NONE,CIMG_00780_3_CIMG_00780_prepared_Milbemycin
CIMG_00780,8,CIMG_00780_prepared_Milbemycin,-9.5,0.211,-6.978,16.6939,22.4535,63,24.0,0.0,6,HydrogenBond;Hydrophobic,NONE,CIMG_00780_8_CIMG_00780_prepared_Milbemycin
CIMG_00780,2,CIMG_00780_prepared_Milbemycin,-10.3,0.578,-7.4362,-12.9353,-14.6505,93,26.0,0.0,6,HydrogenBond;Hydrophobic,NONE,CIMG_00780_2_CIMG_00780_prepared_Milbemycin
CIMG_00780,7,CIMG_00780_prepared_Milbemycin,-9.1,0.257,-10.9808,-27.4943,-3.8876,56,24.0,0.05,4,Hydrophobic,TMD,CIMG_00780_7_CIMG_00780_prepared_Milbemycin
CIMG_00780,6,CIMG_00780_prepared_Milbemycin,-8.9,0.334,18.6047,-16.3728,-15.4714,50,22.0,0.0,6,Hydrophobic,NONE,CIMG_00780_6_CIMG_00780_prepared_Milbemycin
CIMG_00780,1,CIMG_00780_prepared_Beauvericin,-8.6,0.913,5.2501,-35.776,-19.6413,153,30.0,0.05,9,Hydrophobic;PiStacking,TMD,CIMG_00780_1_CIMG_00780_prepared_Beauvericin
CIMG_00780,3,CIMG_00780_prepared_Beauvericin,-11.5,0.424,-0.5827,28.9579,2.6713,100,28.0,0.0,10,Hydrophobic;PiStacking;SaltBridge,NONE,CIMG_00780_3_CIMG_00780_prepared_Beauvericin
CIMG_00780,8,CIMG_00780_prepared_Beauvericin,-8.2,0.211,-6.978,16.6939,22.4535,63,28.0,0.0,7,HydrogenBond;Hydrophobic,NONE,CIMG_00780_8_CIMG_00780_prepared_Beauvericin
CIMG_00780,7,CIMG_00780_prepared_Beauvericin,-8.5,0.257,-10.9808,-27.4943,-3.8876,56,28.0,0.05,9,Hydrophobic;PiStacking,TMD,CIMG_00780_7_CIMG_00780_prepared_Beauvericin
CIMG_00780,6,CIMG_00780_prepared_Beauvericin,-8.4,0.334,18.6047,-16.3728,-15.4714,50,28.0,0.0,10,Hydrophobic,NONE,CIMG_00780_6_CIMG_00780_prepared_Beauvericin
CIMG_00780,3,CIMG_00780_prepared_Tacrolimus,-9.7,0.424,-0.5827,28.9579,2.6713,100,28.0,0.0,10,HydrogenBond;Hydrophobic,NONE,CIMG_00780_3_CIMG_00780_prepared_Tacrolimus
CIMG_00780,2,CIMG_00780_prepared_Tacrolimus,-10.1,0.578,-7.4362,-12.9353,-14.6505,93,28.0,0.0,9,HydrogenBond;Hydrophobic;SaltBridge,NONE,CIMG_00780_2_CIMG_00780_prepared_Tacrolimus
CIMG_00780,7,CIMG_00780_prepared_Tacrolimus,-8.4,0.257,-10.9808,-27.4943,-3.8876,56,28.0,0.05,9,Hydrophobic,TMD,CIMG_00780_7_CIMG_00780_prepared_Tacrolimus
CIMG_00780,6,CIMG_00780_prepared_Tacrolimus,-8.2,0.334,18.6047,-16.3728,-15.4714,50,28.0,0.0,6,HydrogenBond;Hydrophobic,NONE,CIMG_00780_6_CIMG_00780_prepared_Tacrolimus
"""

df = pd.read_csv(StringIO(csv_data))

# Pareto front computation
def pareto_front(group):
    dominated = []
    for i, row_i in group.iterrows():
        dom = False
        for j, row_j in group.iterrows():
            if (row_j['DockingScore'] <= row_i['DockingScore']) and (row_j['PocketProbability'] >= row_i['PocketProbability']):
                if (row_j['DockingScore'] < row_i['DockingScore']) or (row_j['PocketProbability'] > row_i['PocketProbability']):
                    dom = True
                    break
        dominated.append(dom)
    group['Pareto'] = ~np.array(dominated)
    return group

df = df.groupby('ProteinID', group_keys=False).apply(pareto_front)

# select top 3 per protein based on non-dominated; if fewer than 3 pareto, fill with best docking
shortlists = []
for prot, sub in df.groupby('ProteinID'):
    pareto = sub[sub['Pareto']]
    if len(pareto) >= 3:
        shortlist = pareto.sort_values(['DockingScore', 'PocketProbability']).head(3)
    else:
        remaining = sub[~sub['Pareto']].sort_values(['DockingScore', 'PocketProbability'])
        shortlist = pd.concat([pareto, remaining]).head(3)
    shortlists.append(shortlist)

short_df = pd.concat(shortlists).reset_index(drop=True)
display_dataframe_to_user("MD_Shortlist_Pareto", short_df[['ProteinID', 'PocketID', 'Ligand', 'DockingScore', 'PocketProbability', 'Corrected_PLIP_contacts', 'InteractionTypes', 'Domain']])

