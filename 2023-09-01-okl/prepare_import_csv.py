import pandas as pd

df = pd.read_csv("OKL-v3-curated.csv")
df["canonical_exists"] = df["canonical_exists"].fillna(False)
df["salt"] = df["salt"].replace("mesylate", "3")
df["salt"] = df["salt"].replace("hydrochloride hydrate", "7")
df["salt"] = df["salt"].replace("dihydrochloride", "41")

c = df.loc[
    ~df["canonical_exists"],
    ["lincs_id", "name", "cas_number", "chembl_id", "smiles", "inchi", "inchi_key"]
].copy()
c["loaded"] = ""
c["type"] = "small_molecule"
c["date_entered"] = "2023-09-29"
c["curated_by"] = "Jeremy Muhlich"
c.to_csv("to_load_canonicals.csv", index=False)

b = df[[
    "lincs_id",
    "center_batch_id",
    "provider",
    "provider_catalog_id",
    "provider_batch_id",
    "salt",
    "molecular_weight",
    "order_amount",
    "order_date",
    "order_number",
    "requestor_name",
    "flag_for_aliquot",
    "reserve_for_qc",
    "admin_comments",
]].copy()
b["loaded"] = ""
b["type"] = "small_molecule_batch"
b["concentration"] = 10
b["solvent"] = 35 # DMSO
b.to_csv("to_load_batches.csv", index=False)
