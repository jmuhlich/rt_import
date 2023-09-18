import csv
import pandas as pd
import rdkit.Chem.AllChem

df = pd.read_excel("OTAVA cpd list.xlsx")
df = df.rename(columns={
    "SALTDATA": "salt",
    "CHEMICAL NAME": "name",
    "SMILES": "smiles",
    "COMPOUND ID": "provider_catalog_id",
    "QUANTITY mg": "order_amount",
    "PROVIDER": "provider",
    "LOT NUMBER": "provider_batch_id",
    "Order date": "order_date",
    "Concentration": "concentration",
    "Solvent": "solvent",
})
df["Mol"] = df["smiles"].map(rdkit.Chem.AllChem.MolFromSmiles)
df["inchi"] = df["Mol"].map(rdkit.Chem.AllChem.MolToInchi)
df["inchi_key"] = df["Mol"].map(rdkit.Chem.AllChem.MolToInchiKey)
df["salt"] = df["salt"].replace("HCl", "2")
df["salt"] = df["salt"].replace("", None)
df["provider_catalog_id"] = df["provider_catalog_id"].str.replace(r"^P", "", regex=True)
df["order_amount"] = df["order_amount"].astype(str) + " mg"
df["concentration"] = df["concentration"].str.replace(r" mM$", "", regex=True).astype(float)
df["solvent"] = df["solvent"].replace("DMSO", "35")

c = df[["lincs_id", "name", "smiles", "inchi", "inchi_key"]].copy()
c["type"] = "small_molecule"
c["date_entered"] = "2022-11-21"
c["curated_by"] = "Corey Ma"
c.to_csv("to_load_canonicals.csv", index=False)

b = df[[
    "lincs_id",
    "salt",
    "provider",
    "provider_catalog_id",
    "provider_batch_id",
    "order_date",
    "order_amount",
    "concentration",
    "solvent",
]].copy()
b["type"] = "small_molecule_batch"
b["center_batch_id"] = 1
b.to_csv("to_load_batches.csv", index=False)
