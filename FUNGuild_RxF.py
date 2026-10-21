import os
import pandas as pd
os.getcwd()
os.listdir()

genus = pd.read_table('C:\\Users\\Kat\\VS_Code\\Bioinformatics\\export-R\\genus-freq\\genus_table.tsv', sep='\t', skiprows=1) #ignore the first row
genus.columns.values[0] = 'taxonomy'
genus.head(10)
genus.to_csv("genus_table.tsv", sep='\t')

#git clone https://github.com/UMNFuN/FUNGuild
#cd FUNGuild
#python Guilds_v1.1.py -otu /Users/Kat/VS_Code/genus_table.tsv -db fungi -m -u