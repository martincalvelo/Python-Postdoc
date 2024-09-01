import pandas as pd
import MDAnalysis as mda
import os

# Execute with this bash script. Iterates for all folders of all simulations and create a dataframe (stored as csv) merging the results of analysis_ABlina_database.py
##########################################################################################################################

##!/bin/bash
#route="/home/usc/qf/mcs/lustre/XUNTA_POSTDOC/PRODUCTION"
#route2="/home/usc/qf/mcs/lustre/XUNTA_POSTDOC/scripts"
#for memb in POPC POPC-POPE_8-5 POPC-POPE-POPS_8-5-1 POPC-POPE-POPS-CHOL_8-5-1-4 POPC-POPE-POPS-CHOL-DPSM_8-5-1-4-2 POPC-POPE-POPS-CHOL-DPSM-GM1_8-5-1-4-2-1 POPC-POPE-POPS-CHOL-DPSM-GM1_8-5-1-4-2-2 POPC-POPE-POPS-CHOL-DPSM-GM1_8-5-1-4-2-4 POPC-GM1_20-1 POPC-GM1_8-1 POPC-CHOL-GM1_15-4-1 POPC-CHOL-GM1_13-6-1 POPC-CHOL-GM1_11-8-1 POPC-POPE-POPS-CHOL-DPSM-GM3_8-5-1-4-2-1 POPC-POPE-POPS-CHOL-DPSM-GM3_8-5-1-4-2-2 POPC-POPE-POPS-CHOL-DPSM-GM3_8-5-1-4-2-4 POPC-GM3_20-1 POPC-GM3_8-1; do
#        for pep in "1ba4" "1iyt" "2otk" "5oqv" "1z0q"; do 
#                for i in 16 32; do
#                        cd $route/$memb/$pep/$i/production/
#                                python3 $route2/group_csvs_ablina.py
#                        cd $route/$memb/$pep/$i/production_rep2/
#                                python3 $route2/group_csvs_ablina.py
#                        cd $route/$memb/$pep/$i/production_rep3/
#                                python3 $route2/group_csvs_ablina.py
#                done
#        done
#done

##########################################################################################################################
#Dictionary for results
Results_av = { "Contacts_mPOPC" : {},
          "Contacts_mPOPE" : {},
          "Contacts_mPOPS" : {},
          "Contacts_mCHOL" : {},
          "Contacts_mDPSM" : {},
          "Contacts_mGM3" : {},
          "Contacts_mGM1" : {},
          "Contacts_mGB3" : {},
          "Contacts_mtotal" : {},
          "Contacts_aPOPC" : {},
          "Contacts_aPOPE" : {},
          "Contacts_aPOPS" : {},
          "Contacts_aCHOL" : {},
          "Contacts_aDPSM" : {},
          "Contacts_aGM3" : {},
          "Contacts_aGM1" : {},
          "Contacts_aGB3" : {},
          "Contacts_atotal" : {},
          "Aggregation" : {},
          "Num_clusters" : {},
          "Av_distance" : {},
          "Min_distance" : {},
          "AB_in_membrane" : {},
          "Area_per_lipid" : {},
          "Deposition" : {} }        

##########################################################################################################################

#Read previous dataframes
df_top = pd.read_csv("info_top.csv", header=0) 
df_contacts= pd.read_csv("Results_vs_time_contacts-aggreg.csv", header=0)
df_apl= pd.read_csv("area_per_lipid.csv", header=0)

#Create a dataframe with the average results of contacts and aggregation, since we have 3 replicas per system
#First, contacts with lipids
list_lipids=['POPC','POPE','POPS','CHOL','DPSM','GM1','GM3','GB3']
for i in list_lipids:
    name_key="Contacts_m"+str(i)
    Results_av[ name_key ] = round(df_contacts[ name_key ].mean(), 3)
Results_av[ "Contacts_mtotal" ] = round(df_contacts[ "Contacts_mtotal" ].mean(), 3)
for i in list_lipids:
    name_key="Contacts_a"+str(i)
    Results_av[ name_key ] = round(df_contacts[ name_key ].mean(), 3)
Results_av[ "Contacts_atotal" ] = round(df_contacts[ "Contacts_atotal" ].mean(), 3)
#Now, the other keys
other_keys= ["Aggregation", "Num_clusters", "Av_distance", "Min_distance", "AB_in_membrane", "Deposition"]
for name_key in other_keys:
    Results_av [ name_key ] = round(df_contacts[ name_key ].mean(), 3)
Results_av [ "Area_per_lipid" ] = round(df_apl["Area_per_lipid"].mean() , 3)
df_av = pd.DataFrame(Results_av, index=[0])

#Merge the df with the averaged results with the df with the topology information
df = pd.concat([df_top, df_av], axis=1)

#Write output. If the csv with the results already exits, new results will be appended to th existing file
route_output='/home/usc/qf/mcs/lustre/XUNTA_POSTDOC/RESULTS/'
output_file = route_output+"results_aggregation-and-contacts.csv"
if os.path.exists(output_file) is True:
    df.to_csv(output_file, mode='a', header=False, index=False)
else:
    df.to_csv(output_file, mode='w', header=True, index=False)
