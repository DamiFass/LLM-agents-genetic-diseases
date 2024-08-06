from lxml import etree
import pandas as pd

# Compounds symbol-name dict:
compound_dict = {
    'C00002': 'ATP',
    'C00022': 'Pyruvate',
    'C00238': 'K+',
    'C00031': 'Glucose',
    'C00162': 'Fatty acid',
    'C00076': 'Ca2+'
}

# Load the KGML file
file_path = '../data/KGML/hsa04930.xml'
tree = etree.parse(file_path)
root = tree.getroot()

# Extract gene-pathway relationships and associated diseases
entries = root.findall(".//entry")
relations = root.findall(".//relation")

# Get a dict with entry types as keys and occurences as values:
D = {}
for entry in entries:
    D[entry.get('type')] = D.get(entry.get('type'),0) + 1
    
entries_types = list(D.keys())
D_names = {}
for name in entries_types:
    new_key = name + '_names'
    D_names[new_key] = []

dict_entry_id = {}
entries_types_names = list(D_names.keys())
for entry in entries:
    index = entries_types.index(entry.get('type'))
    # Get the graphics name (i.e. the symbol, not the entry code)
    if index == 1: # if it's a compound
        entry_name = entry.find('graphics').get('name')
        D_names[ entries_types_names[index] ].append(compound_dict[entry_name])
        dict_entry_id[entry.get('id')] = compound_dict[entry_name]
    else:
        entry_name = entry.find('graphics').get('name')
        firts_comma_index = entry_name.find(',')
        # Saving only the firs symbol, as it can have multiple ones (e.g. PKLR, PK1, PKL, )
        D_names[ entries_types_names[index] ].append(entry_name[:firts_comma_index])
        dict_entry_id[entry.get('id')] = entry_name[:firts_comma_index]
        
# Make sure there's no repetitions in D_names:
for key in D_names:
    D_names[key] = list(set(D_names[key]))

D
D_names.keys()
D_names
dict_entry_id

rel_types = {}
rel_counter = 0
for relation in relations:
    rel_counter += 1
    
    # rel_types[relation.get('type')] = rel_types.get(relation.get('type'),0) + 1
    if relation.get('type') == 'PCrel':
        name = 'PCrel'
    if len(relation.findall('subtype')) > 1:
        name = ''
        for subtype in relation.findall('subtype'):
            name = name + subtype.get('name') + ' ' 
    else:
        for subtype in relation.findall('subtype'):
            name = subtype.get('name')
            
    # rel_types[name] = rel_types.get( name , 0 ) + 1
    
    if name not in list(rel_types.keys()):
        rel_types[name] = {'entry1': [], 'entry2': []}

    rel_types[name]['entry1'].append( dict_entry_id[relation.get('entry1')] )
    rel_types[name]['entry2'].append( dict_entry_id[relation.get('entry2')] )
        
rel_types['activation']
print(f'Relations counter: {rel_counter}')


for key in D_names:
    all_nodes = []
    if key != 'map_names':
        for elem in D_names[key]:
            all_nodes.append(elem)
        all_nodes_df = pd.DataFrame(all_nodes)
        all_nodes_df.to_csv(f'../data/ALL_{key[:-6]}.csv', header=['id',f'{key[:-6]}_name'])


    for key in rel_types:
        current_df = pd.DataFrame(rel_types[key])
        current_df.to_csv(f'../data/{key}.csv', header=['id','entry1','entry2'])
