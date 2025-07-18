
# 🧬 PInteract: 

**PInteract** is a tool for the identification of π-involving interactions in protein structures and complexes, including protein–protein, protein–DNA/RNA, and protein–ligand systems. Based on geometric criteria, it can detect:
- Individual π-interaction types, including cation–π, amino-π, His-π, sulfur–π and π-π;
- Clusters of π interactions (i.e., π-chains);
- Specific recurrent spatial motifs.


📄 _For full details, please refer to our publication:_  
**D. Li, F. Pucci, M. Rooman** [PInteract: detecting aromatic-involving interactions in proteins and protein-nucleic acid complexes](https://www.google.com/) _(submitted)_

---

## 📦 Installation

To install **PInteract**, first clone the repository and set up the alias:

```bash
cd /your/path/
git clone https://github.com/3BioCompBio/PInteract.git 
alias PInteract="/your/path/PInteract/exec/PInteract"
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
Default parameters are stored in the file `parameters_default`. You may override them by creating a new file named `parameters` in the same directory. Do not modify `parameters_default` directly.

Key parameters:

| Parameter     | Description                                                                                |
| ------------- | ------------------------------------------------------------------------------------------ |
| `CATDmax`     | Max distance for cation–π, amino–π, and His–π interactions (in Ångström) |
| `CATangmax1`  | Angular factor defining the interaction cylinder (see paper)                               |
| `PiPiDmax`    | Max distance for π–π interactions                                              |
| `PiPiangmax1` | Angular factor for π–π                                                         |
| `SPiDmax`     | Max distance for sulfur–π interactions                                               |
| `SPiangmax1`  | Angular factor for sulfur–π                                                          |
| `redundancy`  | 1 = report all interactions for all chains; 0 = suppress duplicates in identical chains    |

## 📤 Output
PInteract generates three output files:

PInteract.csv: Table of individual π interactions; each row corresponds to a single interaction.

PInteract1.csv: Table of π-interaction chains (triads). Even though a π-chain may involve more than three residues, each row records a triplet to maintain uniformity.

PInteract.txt: Human-readable summary, including all information in the CSV files, with comments prefixed by \#.


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
