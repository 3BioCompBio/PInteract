
# 🧬 PInteract: 

**PInteract** is a tool for the identification of π-involving interactions in protein structures and complexes, including protein–protein, protein–DNA, protein-RNA, and protein–ligand systems. Based on geometric criteria, it can detect:
- Individual π-interaction types, including cation–π, amino-π, His-π, sulfur–π and π-π;
- Clusters or chains of π interactions;
- Stair motifs, which are recurrent motifs at protein-DNA/RNA interfaces and combine π-π stacking, cation/amino/His-π and H-bond interactions.


📄 _For full details, please refer to our publication:_  
**D. Li, F. Pucci, M. Rooman** [PInteract: detecting aromatic-involving interactions in proteins and protein-nucleic acid complexes](https://www.google.com/) _(submitted)_

---

## 📦 Installation

To install **PInteract**, first clone the repository and set up the alias:

```bash
cd /your_path/
git clone https://github.com/3BioCompBio/PInteract.git 
alias PInteract="/your_path/PInteract/exec/PInteract"
```

Then go to the PInteract directory and compile:

```bash
cd PInteract
make
```

## ▶️ Usage

To run PInteract to identify the π interactions, do:

* Copy your `.pdb` files into a directory (e.g., `my_pdbs/`).

* Run:
```bash
PInteract
```

* When prompted, provide the directory name (e.g., `my_pdbs`) containing your structures.

## ⚙️ Hyper-parameters
Default parameters are stored in the files `parameters_default` and `parameters`. You may modify them by modifying the file named `parameters`. Do not modify the `parameters_default` file.

Key parameters:

| Parameter     | Description (see our article for details on the calculation of the distance and on the cylindrical model)                           |
| ------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `CATDmax`     | Maximum distance between the closest atoms in the functional groups for cation–π, amino–π and His–π interactions (in Ångström)      |
| `CATangmax1`  | The radius of the cylinder's basis is CATangmax1 times the radius of the aromatic ring for cation–π, amino–π and His–π interactions |
| `PiPiDmax`    | Maximum distance the closest atoms in the functional groups for π–π interactions                                                    |
| `PiPiangmax1` | The radius of the cylinder's basis is PiPiangmax1 times the radius of the aromatic ring for π–π interactions                        |
| `SPiDmax`     | Maximum distance the closest atoms in the functional groups for sulfur–π interactions                                               |
| `SPiangmax1`  | The radius of the cylinder's basis is SPiangmax1 times the radius of the aromatic ring for sulfur–π interactions                    |
| `redundancy`  | 1 = report all interactions for all chains; 0 = suppress duplicates in identical chains                                             |
| ------------- | ----------------------------------------------------------------------------------------------------------------------------------- |

## 📤 Output
PInteract generates three output files in the directory containing the pdb files (ensure that the directory has write permissions enabled):

PInteract.csv: Table of individual π interactions; each row corresponds to a single interaction.

PInteract1.csv: Table of π-interaction chains and of stair motifs. Even though a π-chain may involve more than three residues, and a stair motif more than one stair of three interactions, each row records a triplet to maintain uniformity.

PInteract.txt: Human-readable summary, including all information in the CSV files, with comments prefixed by #.


---

## 🔗 Citation

Please cite our publication when using **PInteract** in your research:

```bibtex
@article{,
  title     = {PInteract: detecting aromatic-involving interactions in proteins and protein-nucleic acid complexes},
  author    = {Li, D. and Pucci, F. and Rooman, M.},
  journal   = {submitted},
  year      = {2025}
}
