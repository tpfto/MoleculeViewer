# MoleculeViewer
### A package for visualizing molecules.

MoleculeViewer is a *Mathematica* package that renders molecules in a manner resembling physical molecular models. In particular, multiple bonds are depicted as out-of-plane bonds.

The package supports features like multiple atom legends, custom coloring of atoms, atom highlighting and tooltips.

The package also provides a number of auxiliary functions that use the services of [Open Babel](http://openbabel.org/), [JME](http://www.molinspiration.com/jme/), [PubChem](https://pubchem.ncbi.nlm.nih.gov/), and [ChemSpider](http://www.chemspider.com/).

Download the paclet from the [releases](https://github.com/tpfto/MoleculeViewer/releases) page, and install it by evaluating the following in *Mathematica*:

    Needs["PacletManager`"]
    PacletInstall["/path/to/paclet/MoleculeViewer-1.0.paclet"]

where `"/path/to/paclet/MoleculeViewer-1.0.paclet"` should be replaced with the actual path to the downloaded paclet file. 

See the file `molviewer.nb` for more detailed information on the package. A molecule gallery `gallery.nb` is also provided.