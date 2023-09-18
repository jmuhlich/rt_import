import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

rt = pd.read_csv("rt-smallmolecules-20221121.csv")
database = list(zip(rt["lincs_id"], rt["smiles"].map(AllChem.MolFromSmiles)))

q = pd.read_excel("OTAVA cpd list.xlsx")
q["Mol"] = q["SMILES"].apply(AllChem.MolFromSmiles)
queries = list(q[["COMPOUND ID", "Mol"]].to_records(index=False))

def make_fingerprint(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)

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
m = results[results["Similarity"] > 0.95]
print(m.to_string(index=False))
