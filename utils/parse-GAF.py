import pandas as pd

gaf_file_path = "/Users/damianofassina/Downloads/goa_human.gaf"

# Read the GAF file into a pandas DataFrame
gaf_df = pd.read_csv(gaf_file_path, sep='\t', comment='!', header=None)

gaf_df.shape

# Define the column names based on the GAF 2.1 format
gaf_df.columns = [
    "DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", 
    "DB:Reference", "Evidence_Code", "With (or) From", "Aspect", 
    "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type", 
    "Taxon", "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID"
]

# Display the DataFrame
gaf_df.head()

gaf_df.groupby('Qualifier').count()

# I can append to the ALL_gene.csv the columns that I want from the GAF dataframe
# and then include them as properties in the graph node.
# DB_Object_Name, DB_Object_Type, Qualifier, GO_ID.
# Also aspects is important!! P, F or C! (Process, function or component)
# Name and Type sono unici, mentre Qualifier and GO_ID can have multiple values.
# There are hundreds of rows with Symbol = INSR, so Qualifier can have a few dozens
# values and GO_ID even hundreds. 
# It would be: <SYMBOL> <Qualifier> <GO_ID>. For example:
# INSR is active in GO:0003723, which should be searched in the Gene Ontology browser

genes = pd.read_csv('../data/ALL_gene.csv')

dict_to_merge = {'Name': [],
                 'Type': [],
                 'Qualifier-GO_ID-Aspect': []
                 }

list_counter = 0
for symbol in genes['Gene_name']:
    # Get only the rows which correspond to the same gene symbol:
    curr_gaf_df = gaf_df[gaf_df['DB_Object_Symbol'] == symbol ]
    n_rows = curr_gaf_df.shape[0]
    curr_gaf_df.groupby('DB_Object_Name').count()
    name = curr_gaf_df['DB_Object_Name'].sample().item()
    dict_to_merge['Name'].append(name)
    obj_type = curr_gaf_df['DB_Object_Type'].sample().item()
    dict_to_merge['Type'].append(obj_type)
    
    # Get all the 'qualifier' attributes for each row:
    qualis = curr_gaf_df['Qualifier'].to_list()
    # Get all the GO_ID attributes for each row:
    go_ids = curr_gaf_df['GO_ID'].to_list()
    # Get all the 'Aspect' attributes for each row:
    aspects = curr_gaf_df['Aspect'].to_list()
    
    dict_to_merge['Qualifier-GO_ID-Aspect'].append([])
    for i in range(len(qualis)):
        curr_D = {}
        curr_D['GO_ID'] = go_ids[i]
        curr_D['Qualifier'] = qualis[i]
        curr_D['Aspect'] = aspects[i]
        dict_to_merge['Qualifier-GO_ID-Aspect'][list_counter].append(curr_D)
    # print(len(dict_to_merge['Qualifier-GO_ID-Aspect'][list_counter]))
    list_counter += 1
    
df_to_merge = pd.DataFrame(dict_to_merge)

df_merged = pd.concat([genes.reset_index(drop=True),df_to_merge.reset_index(drop=True)], axis=1)

df_merged.head(3)

df_merged.drop('id',axis=1).to_csv(f'../data/ALL_gene_rich.csv', header=['Symbol','Name','Type','Qualifier-GO_ID-Aspect'])

# This is how the df should look like before saving it again to ALL_genes.csv

# gene_Symbol, Name,             Type,    Qualifier,GO_ID,Aspect
#   INSR     , Insulin receptor, protein, [ {'GO_ID':0001540,
#                                            'Qualifier': 'is_active_in',
#                                            'Aspect': 'P'},
#                                            {'GO_ID': 000.... } 
#                                         ]