import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, MolStandardize
from tqdm import tqdm

RDLogger.DisableLog("rdApp.*")

standardizer = MolStandardize.Standardizer()
remover = MolStandardize.fragment.FragmentRemover()
enumerator = MolStandardize.rdMolStandardize.TautomerEnumerator()
enumerator.SetRemoveSp3Stereo(False)
enumerator.SetRemoveBondStereo(False)


def mol_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = remover.remove(mol)
    mol = standardizer.standardize(mol)
    mol = enumerator.Canonicalize(mol)
    return mol

def make_fingerprint(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2)

def similarity(fp1, fp2):
    return DataStructs.FingerprintSimilarity(
        fp1, fp2, DataStructs.TanimotoSimilarity
    )

rt = (
    pd.read_csv("rt-smallmolecules-20230922.csv")
    .set_index("lincs_id", verify_integrity=True)
)
tqdm.pandas(desc="Parse RT smiles")
rt["Mol"] = rt["smiles"].progress_map(mol_from_smiles)
rt["Fp"] = rt["Mol"].map(make_fingerprint)

lib = (
    pd.read_csv("OKL-v3-vendor-data.csv")
    .set_index("name", verify_integrity=True)
)
tqdm.pandas(desc="Parse new library smiles")
lib["Mol"] = lib["smiles"].progress_map(mol_from_smiles)
lib["Fp"] = lib["Mol"].map(make_fingerprint)

sim = (
    pd.merge(
        rt["Fp"].reset_index(),
        lib["Fp"].reset_index(),
        how="cross",
        suffixes=("_rt", "_lib"),
    )
    .set_index(["lincs_id", "name"])
)
sim["Similarity"] = sim.apply(
    lambda x: similarity(x["Fp_rt"], x["Fp_lib"]), axis=1
)

matches = (
    sim.loc[sim["Similarity"] > 0.9, "Similarity"]
    .rename_axis(index={"name": "library match"})
)
print()
if len(matches):
    print("Matches found:\n")
    print(matches.to_string())
else:
    print("No matches found")
