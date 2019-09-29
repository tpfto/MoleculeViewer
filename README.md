# MoleculeViewer
### A package for visualizing molecules.

[![DOI](https://zenodo.org/badge/104575485.svg)](https://zenodo.org/badge/latestdoi/104575485)

![cucurbituril rendered using MoleculeViewer](https://user-images.githubusercontent.com/13274842/65828934-ed138600-e2d2-11e9-9b07-32fb147249d1.png)

MoleculeViewer is a *Mathematica* package that renders molecules in a manner resembling physical molecular models. In particular, multiple bonds are depicted as out-of-plane bonds.

The package supports additional features like stylized molecule depictions (ball-and-stick, spacefilling, etc.), multiple atom legends, custom coloring of atoms, atom highlighting and tooltips.

The package also provides a number of auxiliary functions that use the services of [Open Babel](http://openbabel.org/), [Imago OCR](https://lifescience.opensource.epam.com/imago/index.html), [JME](http://www.molinspiration.com/jme/), [Accelrys JDraw](http://download.accelrys.com/freeware/accelrys_draw/ReleaseNotes_AccelrysDraw_4.1.pdf), [JChemPaint](https://jchempaint.github.io/), [ChEMBL Beaker](https://github.com/chembl/chembl_beaker), [ChemSpider](http://www.chemspider.com/), [NCI CACTUS](https://cactus.nci.nih.gov/chemical/structure), [PubChem](https://pubchem.ncbi.nlm.nih.gov/), [RCSB PDB](https://www.rcsb.org/), and [ZINC](https://zinc.docking.org/).

Download the paclet from the [releases](https://github.com/tpfto/MoleculeViewer/releases) page, and install it by evaluating the following in *Mathematica*:

    Needs["PacletManager`"]
    PacletInstall["/path/to/paclet/MoleculeViewer-1.0.paclet"]

where `"/path/to/paclet/MoleculeViewer-1.0.paclet"` should be replaced with the actual path to the downloaded paclet file. 

See the file `molviewer.nb` for more detailed information on the package. A molecule gallery `gallery.nb` is also provided.
