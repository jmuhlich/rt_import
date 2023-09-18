import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MolStandardize

# I think the log level is actually global across rdkit...
Chem.inchi.logger.setLevel(3)

standardizer = MolStandardize.Standardizer()
remover = MolStandardize.fragment.FragmentRemover()

def mol_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.MolFromInchi(Chem.MolToInchi(mol))
    mol = remover.remove(mol)
    mol = standardizer.standardize(mol)
    return mol

def make_fingerprint(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2)

rt = pd.read_csv("rt-smallmolecules-20230908.csv")
rt["Mol"] = rt["smiles"].map(mol_from_smiles)
database = list(zip(rt["lincs_id"], rt["Mol"]))

q = pd.read_csv("OKL_Compounds_20230831_sde_out.csv")
q["Mol"] = q["smiles"].map(mol_from_smiles)
queries = list(q[["name", "Mol"]].to_records(index=False))

fps_db = [make_fingerprint(rec[1]) for rec in database]
fps_q = [make_fingerprint(rec[1]) for rec in queries]

matches = []
for fp1 in fps_q:
    best_similarity = 0
    best_names = []
    for i, fp2 in enumerate(fps_db):
        similarity = DataStructs.FingerprintSimilarity(
            fp1, fp2, DataStructs.TanimotoSimilarity
        )
        if similarity > best_similarity:
            best_similarity = similarity
            best_names = [database[i][0]]
        elif similarity == best_similarity:
            best_names.append(database[i][0])
    matches.append((best_names, best_similarity))

results = pd.DataFrame({
    "Query": [x[0] for x in queries],
    "Matches": [",".join(x[0]) for x in matches],
    "Similarity": [x[1] for x in matches],
})
m = results[results["Similarity"] > 0.9]
if len(m):
    print(m.to_string(index=False))
else:
    print("No matches found")
