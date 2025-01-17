import re
from collections import Counter
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import MolFromSmiles, MolFromSmarts, AllChem
import numpy as np
import pandas as pd
from tqdm import tqdm

def get_functional_group_position(smiles, smarts):

    molecule = Chem.MolFromSmiles(smiles)
    pattern = Chem.MolFromSmarts(smarts)
    matches = molecule.GetSubstructMatches(pattern)
    # Return the position of the first atom in the first match
    return matches[0][0] if matches else None


def determine_alpha_beta_BAs(smiles):

    smarts_pattern = '[#6]12-[#6]-[#6]-[#6]3-[#6](-[#6]-1-[#6]-[#6]-[#6]-2)-[#6]-[#6]-[#6]1-[#6]-3-[#6]-[#6]-[#6]-[#6]-1'
    smarts_pattern_withMet = "[#6]-[#6]12-[#6]-[#6]-[#6]3-[#6](-[#6]-[#6]-[#6]4-[#6]-3(-[#6]-[#6]-[#6]-[#6]-4)-[#6])-[#6]-1-[#6]-[#6]-[#6]-2"


    mol = MolFromSmiles(Chem.MolToSmiles(MolFromSmiles(smiles)))

    #remove C=C doublebonds
    len_db_all = len(mol.GetSubstructMatches(Chem.MolFromSmarts('C=C')))

    if not len(mol.GetSubstructMatches(Chem.MolFromSmarts(smarts_pattern))) == 1:
        if len_db_all > 0:
            rxn = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1][C:2]')

            # Iterate until there are no more double bonds
            while mol.HasSubstructMatch(Chem.MolFromSmarts('C=C')):
                products = rxn.RunReactants((mol,))
                mol = products[0][0]
                mol.UpdatePropertyCache(strict=False)

    mol = AllChem.AddHs(mol)


    pattern = MolFromSmarts(smarts_pattern)
    pattern_withMet = MolFromSmarts(smarts_pattern_withMet)


    matches = mol.GetSubstructMatches(pattern)
    matches_pattern_withMet = mol.GetSubstructMatches(pattern_withMet)

    full_match = matches[0]
    full_match_pattern_withMet_2 = matches_pattern_withMet[0]

    # Get the atom indices at position 0 and 12 in the matched substructure
    atom_0 = full_match[0]
    atom_12 = full_match[12]

    # Find the neighbors of atom_0 and atom_12
    methyl_1_id = [a.GetIdx() for a in mol.GetAtomWithIdx(atom_0).GetNeighbors() if a.GetIdx() not in full_match][0]
    methyl_2_id = [a.GetIdx() for a in mol.GetAtomWithIdx(atom_12).GetNeighbors() if a.GetIdx() not in full_match][0]

    # Generate 3D coordinates
    Chem.SanitizeMol(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0)
    AllChem.MMFFOptimizeMolecule(mol)

    # Get the coordinates of methyl_1_id and methyl_2_id
    conf = mol.GetConformer()
    methyl_1_coords = np.array(conf.GetAtomPosition(methyl_1_id))
    methyl_2_coords = np.array(conf.GetAtomPosition(methyl_2_id))

    beta_oxygen_neighbors = []
    alpha_oxygen_neighbors = []

    # Iterate over each atom in the matches
    for idx, match in enumerate(full_match_pattern_withMet_2):
        atom = mol.GetAtomWithIdx(match)

        # Get the neighbors of the atom
        neighbors = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetIdx() not in full_match_pattern_withMet_2]

        # Skip if there are not exactly 2 neighbors
        if len(neighbors) != 2:
            continue

        # Get the coordinates of the neighbors
        neighbor_1_coords = np.array(conf.GetAtomPosition(neighbors[0].GetIdx()))
        neighbor_2_coords = np.array(conf.GetAtomPosition(neighbors[1].GetIdx()))

        # Calculate distances to methyl1 and methyl2 for each neighbor
        distances_1 = [np.linalg.norm(methyl_1_coords - neighbor_1_coords),
                       np.linalg.norm(methyl_2_coords - neighbor_1_coords)]
        distances_2 = [np.linalg.norm(methyl_1_coords - neighbor_2_coords),
                       np.linalg.norm(methyl_2_coords - neighbor_2_coords)]

        # Determine the winning and losing neighbor
        if min(distances_1) < min(distances_2):
            if neighbors[0].GetAtomicNum() == 8:
                beta_oxygen_neighbors.append(idx)
            if neighbors[1].GetAtomicNum() == 8:
                alpha_oxygen_neighbors.append(idx)
        else:
            if neighbors[0].GetAtomicNum() == 8:
                alpha_oxygen_neighbors.append(idx)
            if neighbors[1].GetAtomicNum() == 8:
                beta_oxygen_neighbors.append(idx)

    beta_oxygen_neighbors = [id for id in beta_oxygen_neighbors if any(atom.GetAtomicNum() == 1 for atom in mol.GetAtomWithIdx(matches_pattern_withMet[0][id]).GetNeighbors())]
    alpha_oxygen_neighbors = [id for id in alpha_oxygen_neighbors if any(atom.GetAtomicNum() == 1 for atom in mol.GetAtomWithIdx(matches_pattern_withMet[0][id]).GetNeighbors())]

    output = {}

    for idx in beta_oxygen_neighbors:
        output[f'[{idx}*]O_beta'] = output.get(f'[{idx}*]O_beta', 0) + 1

    for idx in alpha_oxygen_neighbors:
        output[f'[{idx}*]O_alpha'] = output.get(f'[{idx}*]O_alpha', 0) + 1

    return output


def split_molecule(smiles, split_position):
    molecule = Chem.MolFromSmiles(smiles)

    # Identify the bonds to break
    bonds = [bond for bond in molecule.GetBonds()
             if bond.GetBeginAtomIdx() == split_position
             or bond.GetEndAtomIdx() == split_position]
    fragments = Chem.FragmentOnBonds(molecule, [bond.GetIdx() for bond in bonds])

    # Find the fragment that contains [123*]
    for fragment in Chem.GetMolFrags(fragments, asMols=True):
        if fragment.HasSubstructMatch(Chem.MolFromSmarts('[123*]')):
            return Chem.MolToSmiles(fragment)
    return None


def get_smiles_of_smarts_on_smiles(smiles):
    rm_dbs = False
    on_smarts = "[#6]-[#6]12-[#6]-[#6]-[#6]3-[#6](-[#6]-[#6]-[#6]4-[#6]-3(-[#6]-[#6]-[#6]-[#6]-4)-[#6])-[#6]-1-[#6]-[#6]-[#6]-2"
    smiles_unster = smiles.replace('@', '')
    mol = Chem.MolFromSmiles(smiles_unster)
    len_db_all = len(mol.GetSubstructMatches(Chem.MolFromSmarts('C=C')))

    if not len(mol.GetSubstructMatches(Chem.MolFromSmarts(on_smarts))) == 1:
        rm_dbs = True
        if len_db_all > 0:
            rxn = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1][C:2]')

            # Iterate until there are no more double bonds
            while mol.HasSubstructMatch(Chem.MolFromSmarts('C=C')):
                products = rxn.RunReactants((mol,))
                mol = products[0][0]

    chem_subst_Mol = Chem.ReplaceCore(mol, Chem.MolFromSmarts(on_smarts), labelByIndex=True)
    substs = Chem.MolToSmiles(chem_subst_Mol)
    on_smarts_smiles = substs.split('.')
    if rm_dbs == True:
        len_db_NOT_core = len(chem_subst_Mol.GetSubstructMatches(Chem.MolFromSmarts('C=C')))
        on_smarts_smiles.extend(["[*]C=C"] * (len_db_all - len_db_NOT_core))

    return on_smarts_smiles


def interpret_subst_as_BA(smiles_subst):
    subst_not_tail_s = [string for string in smiles_subst if not string.startswith('[18*]')]
    subst_not_tail_s_2 = [re.sub(r'\[\d+(\*)\]', '[\\1]', string) for string in subst_not_tail_s]

    subst_not_tail_s = subst_not_tail_s + subst_not_tail_s_2
    counts_on_core = Counter(subst_not_tail_s)

    tail_s = [string for string in smiles_subst if string.startswith('[18*]')][0]
    tail_s = tail_s.replace('[18*]', '[123*]')

    carboxylic_acid_smarts = 'C(=O)[O;h1]'
    amide_smarts = 'C(=O)[N;H1,H2]'
    ester_smarts = 'C(=O)[O]'

    positions = {}
    positions['Carboxylic Acid'] = get_functional_group_position(tail_s, carboxylic_acid_smarts)
    positions['Amide'] = get_functional_group_position(tail_s, amide_smarts)
    positions['Ester'] = get_functional_group_position(tail_s, ester_smarts)

    molecule = Chem.MolFromSmiles(tail_s)
    pattern = Chem.MolFromSmarts('[123*]')
    match = molecule.GetSubstructMatches(pattern)
    position_123 = match[0][0] if match else None

    OH_subst = Chem.MolFromSmarts('[OH]')
    # Filter out None positions and get the functional group closest to [123*]
    if not all(value is None for value in positions.values()):
        positions = {k: abs(v - position_123) for k, v in positions.items() if v is not None}
        closest_group = min(positions, key=positions.get)
        closest_position = positions[closest_group]
        tail_s = split_molecule(tail_s, closest_position)
        tail_unc = Chem.MolFromSmiles(tail_s)
    else:
        tail_unc = Chem.MolFromSmiles(tail_s)
        closest_group = 'None'
    return {'tail_has': closest_group,
            'core_counts': counts_on_core,
            'tail_OH_counts': len(tail_unc.GetSubstructMatches(OH_subst))}


def wrap_bile_acid_annotation(smiles):

    try:
        output = interpret_subst_as_BA(get_smiles_of_smarts_on_smiles(smiles))
        try:
            output['core_counts'].update(determine_alpha_beta_BAs(smiles))
        except:
            print(smiles)
        output['SMILES'] = smiles
    except:
        output = {'SMILES': smiles,
                  'tail_has': None,
                  'core_counts': None,
                  'tail_OH_counts': None}
    return output


def process_BA_smiles_csv(csv_input_path, csv_output_path):

    df = pd.read_csv(csv_input_path)
    unique_smiles = df['SMILES'].unique()

    results = [wrap_bile_acid_annotation(smiles) for smiles in tqdm(unique_smiles)]

    results_df = pd.DataFrame(results)
    df = pd.merge(df, results_df, on='SMILES')
    core_counts_df = df['core_counts'].apply(lambda x: pd.Series(x, dtype='int'))
    df = pd.concat([df, core_counts_df], axis=1)
    df = df.drop(columns='core_counts')

    df.to_csv(csv_output_path, index=False)

    return df


def mol_with_atom_index(mol):
    IPythonConsole.ipython_useSVG = True

    atoms = mol.GetNumAtoms()
    for idx in range(atoms):
        mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx()))
    return mol


if __name__ == '__main__':


    df = process_BA_smiles_csv(csv_input_path='C:/PostDoc/Classify Bile Acids/BILELIB19_Names_ok.csv',
                               csv_output_path='C:/PostDoc/Classify Bile Acids/BILELIB19_Names_ok_labaled.csv')

    print(df)