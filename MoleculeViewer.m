
(* :Title: Molecule Viewer *)

(* :Author: J. M. *)

(* :Summary:

     This package implements an improved molecule viewer derived from Mathematica's default viewer,
     as well as providing a number of useful utility functions for obtaining and visualizing molecular models.

 *)

(* :Copyright:

     © 2017-2019 by J. M. (pleasureoffiguring(AT)gmail(DOT)com)
     This work is free. It comes without any warranty, to the extent permitted by applicable law.
     You can redistribute and/or modify it under the terms of the MIT License.

 *)

(* :Package Version: 3.0 *)

(* :Mathematica Version: 8.0 *)

(* :History:

     3.0 - now under MIT License, added PlotStyle for molecule depictions ("BallAndStick", etc.), added new multiple bond style,
             added support for changing default units of length, added data access functions for PDB and ZINC,
             added Imago OCR support, added support for structure drawing with Accelrys JDraw and JChemPaint,
             added ChEMBL Beaker support functions, added molecule manipulation functions,
             limited support of Entity objects and Molecule objects
     2.1 - removed support for deprecated ChemSpider SOAP API
     2.0 - improved handling of double bonds, added PubChem support, switched to new ChemSpider API, better handling of files and URLs
     1.1 - fixed a few reported bugs, added explicit input dialog for $ChemSpiderAPIToken, support suboptions for RunOpenBabel
     1.0 - initial release

*)

(* :Keywords:
     3D, applet, Beaker, chemical structures, chemistry, ChemSpider, ChEMBL, graphics, Imago OCR, Java, JChemPaint,
     JDraw, JME, models, molecules, OCR, Open Babel, PubChem *)

(* :Limitations:
     currently limited handling of imported chemical file formats
     needs Imago OCR and Open Babel to be installed for some functionality
*)

(* :Requirements:

     Imago OCR (https://lifescience.opensource.epam.com/imago/index.html) for chemical structure OCR,
     Open Babel (http://openbabel.org/) for 3D coordinate generation,
     and whichever of JME (http://www.molinspiration.com/jme/), JChemPaint (http://jchempaint.github.io) or
     Accelrys JDraw (http://accelrys.com/products/informatics/cheminformatics/draw/jdraw.php), for structure drawing, to be installed.

     ChemSpider API key (https://developer.rsc.org/) needed for ChemSpider functionality.

*)

BeginPackage["MoleculeViewer`"]

Unprotect[MoleculeViewer, MoleculeRecognize, MakeMOL, MakeInChIKey,
               GetCACTUS, GetChemicalData, GetChemSpider, GetPDB, GetProteinData, GetPubChem, GetZINC,
               MoleculeDraw, GetMolecule, SetMolecule, JChemPaint, JDraw, JME,
               ChEMBLMOLTo3D, ChEMBLSMILESTo3D, ChEMBLInChITo3D, ChEMBLInChIToMOL, ChEMBLSMILESToMOL,
               ChEMBLAddHydrogens, ChEMBLSanitize, ChEMBLMakeLineNotation, ChEMBLMakeInChIKey, ChEMBLDescriptors,
               RunImago, RunOpenBabel];

(* usage messages *)

MoleculeViewer::usage = "\!\(\*RowBox[{\"MoleculeViewer\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) displays a molecular model corresponding to \!\(\*StyleBox[\"expr\", \"TI\"]\). \!\(\*StyleBox[\"expr\", \"TI\"]\) can be a file name, a URL, a chemical name, a SMILES or an InChI string, or a list containing structural information."

MoleculeDraw::usage = "\!\(\*RowBox[{\"MoleculeDraw\", \"[\", \"]\"}]\) opens a blank instance of the default structure drawing applet.\n\!\(\*RowBox[{\"MoleculeDraw\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) opens the default structure drawing applet and sets the molecule displayed to \!\(\*StyleBox[\"expr\", \"TI\"]\).\n\!\(\*RowBox[{\"MoleculeDraw\", \"[\", RowBox[{StyleBox[\"\\\"Applet\\\"\",ShowStringCharacters->True], \"\[Rule]\", StyleBox[\"\\\"\\!\\(\\*StyleBox[\\\"applet\\\",\\\"TI\\\"]\\)\\\"\", ShowStringCharacters->True]}], \"]\"}]\) launches a specific structure drawing applet.\n\!\(\*RowBox[{\"MoleculeDraw\", \"[\", RowBox[{StyleBox[\"\\\"\\!\\(\\*StyleBox[\\\"applet\\\",\\\"TI\\\"]\\)\\\"\", ShowStringCharacters->True]}], \"]\"}]\) represents an operator form of MoleculeDraw that can be applied to an expression."

MoleculeRecognize::usage = "\!\(\*RowBox[{\"MoleculeRecognize\", \"[\", StyleBox[\"image\", \"TI\"], \"]\"}]\) recognizes molecular information in \!\(\*StyleBox[\"image\", \"TI\"]\) to generate structural information suitable for MoleculeViewer."

GetCACTUS::usage = "\!\(\*RowBox[{\"GetCACTUS\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) uses the NCI/CADD CACTUS Chemical Identifier Resolver to generate structural information suitable for MoleculeViewer."

GetChemicalData::usage = "\!\(\*RowBox[{\"GetChemicalData\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) uses the built-in ChemicalData function to generate structural information suitable for MoleculeViewer.\n\!\(\*RowBox[{\"GetChemicalData\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"],\",\", StyleBox[\"\\\"SMILES\\\"\",ShowStringCharacters->True]}], \"]\"}]\) or \!\(\*RowBox[{\"GetChemicalData\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"],\",\", StyleBox[\"\\\"InChI\\\"\",ShowStringCharacters->True]}], \"]\"}]\) will generate the corresponding line notation."

GetChemSpider::usage = "\!\(\*RowBox[{\"GetChemSpider\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) uses the ChemSpider API to generate structural information suitable for MoleculeViewer.\n\!\(\*RowBox[{\"ChemSpider\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"],\",\", StyleBox[\"\\\"SMILES\\\"\",ShowStringCharacters->True]}], \"]\"}]\) or \!\(\*RowBox[{\"GetChemSpider\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"],\",\", StyleBox[\"\\\"InChI\\\"\",ShowStringCharacters->True]}], \"]\"}]\) will generate the corresponding line notation."

GetMolecule::usage = "\!\(\*RowBox[{\"GetMolecule\", \"[\", \"]\"}]\) extracts structural information from the molecule currently displayed in the default structure drawing applet.\n\!\(\*RowBox[{\"GetMolecule\", \"[\", RowBox[{StyleBox[\"\\\"\\!\\(\\*StyleBox[\\\"applet\\\",\\\"TI\\\"]\\)\\\"\", ShowStringCharacters->True]}], \"]\"}]\) or \!\(\*RowBox[{\"GetMolecule\", \"[\", RowBox[{StyleBox[\"\\\"Applet\\\"\",ShowStringCharacters->True], \"->\", StyleBox[\"\\\"\\!\\(\\*StyleBox[\\\"applet\\\",\\\"TI\\\"]\\)\\\"\", ShowStringCharacters->True]}], \"]\"}]\) extracts structural information from a specific applet."

GetPDB::usage = "\!\(\*RowBox[{\"GetPDB\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) uses the RCSB PDB API to generate structural information suitable for MoleculeViewer."

GetProteinData::usage = "\!\(\*RowBox[{\"GetProteinData\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) uses the built-in ProteinData function to generate structural information suitable for MoleculeViewer."

GetPubChem::usage = "\!\(\*RowBox[{\"GetPubChem\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) uses the PubChem API to generate structural information suitable for MoleculeViewer.\n\!\(\*RowBox[{\"GetPubChem\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"],\",\", StyleBox[\"\\\"SMILES\\\"\",ShowStringCharacters->True]}], \"]\"}]\) or \!\(\*RowBox[{\"GetPubChem\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"],\",\", StyleBox[\"\\\"InChI\\\"\",ShowStringCharacters->True]}], \"]\"}]\) will generate the corresponding line notation."

GetZINC::usage = "\!\(\*RowBox[{\"GetZINC\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) uses the ZINC API to generate structural information suitable for MoleculeViewer.\n\!\(\*RowBox[{\"GetZINC\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"],\",\", StyleBox[\"\\\"SMILES\\\"\",ShowStringCharacters->True]}], \"]\"}]\) or \!\(\*RowBox[{\"GetZINC\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"],\",\", StyleBox[\"\\\"InChI\\\"\",ShowStringCharacters->True]}], \"]\"}]\) will generate the corresponding line notation."

MakeInChIKey::usage = "\!\(\*RowBox[{\"MakeInChIKey\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) constructs the InChI key corresponding to \!\(\*StyleBox[\"expr\", \"TI\"]\)."

MakeMOL::usage = "\!\(\*RowBox[{\"MakeMOL\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) constructs an MDL MOL file corresponding to \!\(\*StyleBox[\"expr\", \"TI\"]\)."

MoleculeCentroid::usage = "\!\(\*RowBox[{\"MoleculeCentroid\", \"[\", StyleBox[\"mol\", \"TI\"], \"]\"}]\) computes the center of mass of the molecule \!\(\*StyleBox[\"mol\", \"TI\"]\)."

MoleculeJoin::usage = "\!\(\*RowBox[{\"MoleculeJoin\", \"[\", RowBox[{SubscriptBox[StyleBox[\"mol\", \"TI\"], StyleBox[\"1\", \"TR\"]], \",\", SubscriptBox[StyleBox[\"mol\", \"TI\"], StyleBox[\"2\", \"TR\"]], \",\", StyleBox[\"\[Ellipsis]\", \"TR\"]}], \"]\"}]\) concatenates two or more molecules."

MoleculeNormalize::usage = "\!\(\*RowBox[{\"MoleculeNormalize\", \"[\", StyleBox[\"mol\", \"TI\"], \"]\"}]\) reorients the molecule \!\(\*StyleBox[\"mol\", \"TI\"]\) so that its center of mass is at the origin, and it lies along its plane of best fit."

MoleculeSplit::usage = "\!\(\*RowBox[{\"MoleculeSplit\", \"[\", StyleBox[\"mol\", \"TI\"], \"]\"}]\) splits the molecule \!\(\*StyleBox[\"mol\", \"TI\"]\) into its connected components."

MoleculeTransform::usage = "\!\(\*RowBox[{\"MoleculeTransform\", \"[\", RowBox[{StyleBox[\"mol\", \"TI\"], \",\", StyleBox[\"tfun\", \"TI\"]}], \"]\"}]\) transforms the molecule \!\(\*StyleBox[\"mol\", \"TI\"]\) with the transformation function \!\(\*StyleBox[\"tfun\", \"TI\"]\)."

RunImago::usage = "\!\(\*RowBox[{\"RunImago\", \"[\", StyleBox[\"image\", \"TI\"], \"]\"}]\) is an interface to Imago OCR, for generating an MDL MOL file for further use."

RunOpenBabel::usage = "\!\(\*RowBox[{\"RunOpenBabel\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) is an interface to Open Babel, for generating structural information suitable for MoleculeViewer."

SetDrawingApplet::usage = "\!\(\*RowBox[{\"SetDrawingApplet\", \"[\", RowBox[{StyleBox[\"\\\"\\!\\(\\*StyleBox[\\\"applet\\\",\\\"TI\\\"]\\)\\\"\", ShowStringCharacters->True]}], \"]\"}]\) sets the default Java chemical drawing applet used by MoleculeDraw, GetMolecule, and SetMolecule."

SetMolecule::usage = "\!\(\*RowBox[{\"SetMolecule\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) sets the chemical stucture in the applet opened by MoleculeDraw to \!\(\*StyleBox[\"expr\", \"TI\"]\).\n\!\(\*RowBox[{\"SetMolecule\", \"[\", RowBox[{StyleBox[\"\\\"\\!\\(\\*StyleBox[\\\"applet\\\",\\\"TI\\\"]\\)\\\"\", ShowStringCharacters->True]}], \"]\"}]\) represents an operator form of SetMolecule that can be applied to an expression."

BondRotate::usage = "\!\(\*RowBox[{\"BondRotate\", \"[\", RowBox[{StyleBox[\"mol\", \"TI\"], \",\", StyleBox[\"bond\", \"TI\"], \",\", StyleBox[\"\[Theta]\", \"TR\"]}], \"]\"}]\) rotates the molecule \!\(\*StyleBox[\"mol\", \"TI\"]\) counterclockwise around \!\(\*StyleBox[\"bond\", \"TI\"]\) by \!\(\*StyleBox[\"\[Theta]\", \"TR\"]\) radians."

ShiftIsolatedAtoms::usage = "\!\(\*RowBox[{\"ShiftIsolatedAtoms\", \"[\", StyleBox[\"mol\", \"TI\"], \"]\"}]\) moves the isolated atoms in molecule \!\(\*StyleBox[\"mol\", \"TI\"]\) closer to the connected components."

ChEMBLAddHydrogens::usage = "\!\(\*RowBox[{\"ChEMBLAddHydrogens\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) uses the ChEMBL Beaker API to add hydrogen atoms to \!\(\*StyleBox[\"expr\", \"TI\"]\)."

ChEMBLDescriptors::usage = "\!\(\*RowBox[{\"ChEMBLDescriptors\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) uses the ChEMBL Beaker API to generate chemical descriptors for \!\(\*StyleBox[\"expr\", \"TI\"]\)."

ChEMBLInChITo3D::usage = "\!\(\*RowBox[{\"ChEMBLInChITo3D\", \"[\", StyleBox[\"inchi\", \"TI\"], \"]\"}]\) uses the ChEMBL Beaker API to generate structural information from the InChI string \!\(\*StyleBox[\"inchi\", \"TI\"]\)."

ChEMBLInChIToMOL::usage = "\!\(\*RowBox[{\"ChEMBLInChIToMOL\", \"[\", StyleBox[\"inchi\", \"TI\"], \"]\"}]\) uses the ChEMBL Beaker API to generate an MDL MOL file from the InChI string \!\(\*StyleBox[\"inchi\", \"TI\"]\)."

ChEMBLMakeInChIKey::usage = "\!\(\*RowBox[{\"ChEMBLMakeInChIKey\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) uses the ChEMBL Beaker API to generate an InChI key from \!\(\*StyleBox[\"expr\", \"TI\"]\)."

ChEMBLMakeLineNotation::usage = "\!\(\*RowBox[{\"ChEMBLMakeLineNotation\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"\\\"\\!\\(\\*StyleBox[\\\"type\\\",\\\"TI\\\"]\\)\\\"\", ShowStringCharacters->True]}], \"]\"}]\) uses the ChEMBL Beaker API to generate line notation from \!\(\*StyleBox[\"expr\", \"TI\"]\), where \!\(\*RowBox[{StyleBox[\"\\\"\\!\\(\\*StyleBox[\\\"type\\\",\\\"TI\\\"]\\)\\\"\", ShowStringCharacters->True]}]\) can be either of \"SMILES\" or \"InChI\"."

ChEMBLMOLTo3D::usage = "\!\(\*RowBox[{\"ChEMBLMOLTo3D\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) uses the ChEMBL Beaker API to generate structural information from \!\(\*StyleBox[\"expr\", \"TI\"]\) in MDL MOL format."

ChEMBLSanitize::usage = "\!\(\*RowBox[{\"ChEMBLSanitize\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) uses the ChEMBL Beaker API to sanitize \!\(\*StyleBox[\"expr\", \"TI\"]\)."

ChEMBLSMILESTo3D::usage = "\!\(\*RowBox[{\"ChEMBLSMILESTo3D\", \"[\", StyleBox[\"smiles\", \"TI\"], \"]\"}]\) uses the ChEMBL Beaker API to generate structural information from the SMILES string \!\(\*StyleBox[\"smiles\", \"TI\"]\)."

ChEMBLSMILESToMOL::usage = "\!\(\*RowBox[{\"ChEMBLSMILESToMOL\", \"[\", StyleBox[\"smiles\", \"TI\"], \"]\"}]\) uses the ChEMBL Beaker API to generate an MDL MOL file from the SMILES string \!\(\*StyleBox[\"smiles\", \"TI\"]\)."

$ChemSpiderAPIToken::usage = "$ChemSpiderAPIToken is the API token used by the ChemSpider API."

$DefaultLengthUnit::usage = "$DefaultLengthUnit is the default unit of length used by MoleculeViewer and related functions."

$DrawingApplets::usage = "$DrawingApplets gives the list of chemical structure drawing Java applets currently supported by MoleculeViewer."

$OpenBabelVersion::usage = "$OpenBabelVersion displays the version of Open Babel available."

OPSINNameToStructure;
$JME; JME; $JDraw; JDraw; $JChemPaint; JChemPaint;
GetJMEMOL; SetJMEMOL; GetJMESMILES; LaunchJME; (* for future removal; retained for backward compatibility *)

(* declare symbols for older versions *)
Highlighted; MemoryConstraint; PlotLegends;
RunProcess; ProcessEnvironment; ProcessDirectory;

(* error/warning messages *)

MoleculeViewer::badchem = "Invalid structural information given to MoleculeViewer.";
MoleculeViewer::badfile = "Valid structural information could not be generated from the file.";
MoleculeViewer::badstrng = "Valid structural information could not be generated from the input string.";
MoleculeViewer::ftconv = "The file `1` is not in a format supported by Import. Conversion with OpenBabel will be attempted.";
MoleculeViewer::is2d = "Two-dimensional coordinates detected; embedding to 3D.";
MoleculeViewer::nbnd = "Incomplete or missing bond information for PlotStyle \[Rule] `1`; using default PlotStyle instead.";
MoleculeViewer::nstl = "PlotStyle \[Rule] `1` is not supported; using default PlotStyle instead.";

$DefaultLengthUnit::badunit = "Invalid unit of length `1` assigned to $DefaultLengthUnit.";

GetChemSpider::token = "ChemSpider API key in $ChemSpiderAPIToken not detected or invalid. Please obtain one from https://developer.rsc.org/.";
GetChemSpider::obsmet = "Method \[Rule] \"SOAP\" is now obsolete. Please obtain a new API key from https://developer.rsc.org/.";
GetChemSpider::noconv = "Cannot convert input of type `1` into `2`.";

MoleculeDraw::nopre = "No preprocessing method supported.";
MoleculeDraw::noapp = GetMolecule::noapp = SetMolecule::noapp = "The `1` applet could not be loaded.";
GetMolecule::nopost = "No postprocessing method supported.";
GetMolecule::nsupo = "`1` output is not supported for the `2` applet.";
SetMolecule::nsupi = "`1` input is not natively supported for the `2` applet. Conversion will be attempted.";

RunOpenBabel::nobab = "Open Babel not installed; please download from http://openbabel.org/ or check Environment[\"PATH\"].";
RunOpenBabel::bdmtd = "Value of option Method \[Rule] `1` is not Automatic, \"Conformer\", \"Generate\", or \"Minimize\". The default method `2` will be used.";
RunOpenBabel::bdff = "Value of option \"ForceField\" \[Rule] `1` is not Automatic, \"GAFF\", \"Ghemical\", \"MMFF94\", \"MMFF94s\", or \"UFF\".";

MoleculeViewer::mem = RunOpenBabel::mem = "Memory allocated by Open Babel exceeded `1` bytes, and was then aborted. Increasing the value of the MemoryConstraint option may yield a result."; 
MoleculeViewer::time = RunOpenBabel::time = "Time spent by Open Babel exceeded `1` seconds, and was then aborted. Increasing the value of the TimeConstraint option may yield a result.";

RunImago::noimg = "Imago OCR not installed; please download from https://lifescience.opensource.epam.com/download/imago.html or check Environment[\"PATH\"].";
RunImago::mem = "Memory allocated by Imago OCR exceeded `1` bytes, and was then aborted. Increasing the value of the MemoryConstraint option may yield a result."; 
RunImago::time = "Time spent by Imago OCR exceeded `1` seconds, and was then aborted. Increasing the value of the TimeConstraint option may yield a result.";

MoleculeRecognize::noocr = "No suitable OCR method found.";
MoleculeRecognize::badimg = "Valid structural information could not be generated from the image.";

ChEMBLMOLTo3D::bdmtd = ChEMBLSMILESTo3D::bdmtd = "Value of option Method \[Rule] `1` is not Automatic, \"ETKDG\", \"KDG\", \"MMFF\", or \"UFF\". The default method `2` will be used.";

MakeMOL::ncnv = "Unable to convert to MOL file.";

BondRotate::mult = "Cannot rotate molecule about the `1` bond `2`.";
BondRotate::nobnd = "Bond `1` cannot be found in the input molecule.";
BondRotate::ring = "Cannot rotate molecule about the bond `1` in a ring.";
MoleculeTransform::notort = "Matrix `1` is not an orthogonal matrix.";
MoleculeTransform::notrgd = "Matrix `1` does not correspond to a rigid transformation; stereochemistry might be affected.";

Begin["`Private`"]

Needs["Utilities`URLTools`"]

(* some variables for internal use *)
$debug = False;
$adjustNormals = True;
$getAllPDBComponents = True;

(* length unit for coordinates and labels *)
If[$VersionNumber >= 9.,
    Quantity[$DefaultLengthUnit = If[$VersionNumber >= 12., "Angstroms", "Picometers"]];
    $ValidLengthUnitQ := Quiet[TrueQ[UnitDimensions[$DefaultLengthUnit] === {{"LengthUnit", 1}}], {Quantity::unkunit}];
    conversionFactor[uIn_, uOut_] := QuantityMagnitude[UnitConvert[Quantity[uIn], uOut]];
    $conversionFactor := If[$ValidLengthUnitQ, conversionFactor["Picometers", $DefaultLengthUnit],
                                       Message[$DefaultLengthUnit::badunit, $DefaultLengthUnit]; Return[$Failed]];
    unitLabel[u_] := With[{s = QuantityUnits`Private`QuantityLabel[u, "Format" -> "UnitLabel"]}, If[StringQ[s], StringTrim[s, "\""], DisplayForm[s]]],
    Needs["Units`"];
    $DefaultLengthUnit = Units`Pico Units`Meter;
    $ValidLengthUnitQ := Quiet[Variables[Units`SI[$DefaultLengthUnit]] === {Units`Meter}];
    conversionFactor[uIn_, uOut_] := Units`Convert[uIn, uOut]/uOut;
    $conversionFactor := If[$ValidLengthUnitQ, conversionFactor[Units`Pico Units`Meter, $DefaultLengthUnit],
                                       Message[$DefaultLengthUnit::badunit, $DefaultLengthUnit]; Return[$Failed]];
    unitLabel[u_] := Which[u === Units`Pico Units`Meter, "pm", u === Units`Angstrom, "\[Angstrom]", True, u]
];

(* utility functions *)

If[$VersionNumber < 10.2, Nothing := Unevaluated[]];

$atoms = Array[ElementData[#, "Symbol"] &, 118] ~Join~ {"D", "T"};
validAtomQ[s_] := StringQ[s] && MemberQ[$atoms, s]

(* CPK-type colors extracted from ElementData[elem, "IconColor"] *)
$colorRules = Thread[$atoms -> ((Apply[RGBColor, IntegerDigits[FromDigits[#, 16], 256, 3]/255.] & /@
                                ImportString[Uncompress["1:eJw9Ur25HDEIfIEbcQdIiB+1YadOQEKp23G51sLeOwU3HywwzPDT//76/e/H11fsc39/\
                                                                        oPVNH7SW+olE4zD4g85xcnrQhOc9aJhAYFUgYsvs2Mb6INu390ikElsSwZ4tY4HWLLscW\
                                                                        PqiowCJQo70/iCeEd0epL6lpqmye3YR82g5Lfh5D/Kz+qppfF+yV5uzYvfPCgUwIxaDCcWA\
                                                                        YAMlWqpQWdkKDq8upjltmrWZWb/ZYhVm2pLzZJx9Jvt92aSmvGjM1I/QqDirfBQXCovckq\
                                                                        8PpdU4C1Zuicvcc0bXyZaoDYWZWQB2zc3v9/eVkntqdvbrEa/sTHIFTH7BWF5uEAQs3+DtR\
                                                                        83eLBGqpH7Xwnh92y5UaLhT1u7tVKpttnsAyeVMH6nBmrpGViyUGPmdhwBm50t9YO7rSLM\
                                                                        Utz0OQjk4pO/ybUDPPTRwVUwmWq+LFb8vu7i65ZYWtKz8MArNmBI3Se2l02EuDTbER7Urdr\
                                                                        Fy3S8DGrURCU2uG49xSiFzWnWxV5QzMoaqY5aDrXHNAJAJyZ66yyl0XT1QzhhE1nJM2nkHwv\
                                                                        feq3NoW1TsRb00WHy8dut8NU+k94iyYm7iYmB9RN3Bo6RWbI6mtdtCLzUskCT7OfTD5Ux/li9nL\
                                                                        sXzok5WqEXdgffm5b63piNr7TQaNW01HOXqbG28mjbAVfzgoL0ovtFG+w/xJTNp"], "List", "Numeric" -> False]) ~Join~
                                {RGBColor[0.6925, 0.715, 0.828], RGBColor[0.6925, 0.8092, 0.715]})];

(* D and T colors: {Blend[{ColorData["Atoms", "H"], RGBColor[0, 0, 0.753]}, 1/4], Blend[{ColorData["Atoms", "H"], RGBColor[0, 0.628, 0]}, 1/4]} *)

(* Pyykkö's covalent radii; https://doi.org/10.1002/chem.200800987 ; entries for deuterium and tritium were slightly reduced from hydrogen *)
$elementCovalentRadii = {32, 46, 133, 102, 85, 75, 71, 63, 64, 67, 155, 139, 126, 116, 111, 
                                       103, 99, 96, 196, 171, 148, 136, 134, 122, 119, 116, 111, 110, 112, 
                                       118, 124, 121, 121, 116, 114, 117, 210, 185, 163, 154, 147, 138, 
                                       128, 125, 125, 120, 128, 136, 142, 140, 140, 136, 133, 131, 232, 
                                       196, 180, 163, 176, 174, 173, 172, 168, 169, 168, 167, 166, 165, 
                                       164, 170, 162, 152, 146, 137, 131, 129, 122, 123, 124, 133, 144, 
                                       144, 151, 145, 147, 142, 223, 201, 186, 175, 169, 170, 171, 172, 
                                       166, 166, 168, 168, 165, 167, 173, 176, 161, 157, 149, 143, 141, 
                                       134, 129, 128, 121, 122, 136, 143, 162, 175, 165, 157, (* D and T *) 30, 28};
$sizeRules = Thread[$atoms -> N[$elementCovalentRadii/2]];

(* Alvarez's van der Waals radii; https://doi.org/10.1039/C3DT50599E ; missing entries given a value of 250 (sodium radius) *)
$vdWRadii = N[{120, 143, 212, 198, 191, 177, 166, 150, 146, 158, 250, 251, 225, 219, 190,
                         189, 182, 183, 273, 262, 258, 246, 242, 245, 245, 244, 240, 240, 238, 239,
                         232, 229, 188, 182, 186, 225, 321, 284, 275, 252, 256, 245, 244, 246, 244,
                         215, 253, 249, 243, 242, 247, 199, 204, 206, 348, 303, 298, 288, 292, 295,
                         250, 290, 287, 283, 279, 287, 281, 283, 279, 280, 274, 263, 253, 257, 249,
                         248, 241, 229, 232, 245, 247, 260, 254, 250, 250, 250, 250, 250, 280, 293,
                         288, 271, 282, 281, 283, 305, 340, 305, 270, 250, 250, 250, 250, 250, 250,
                         250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, (* D and T *) 120, 120}];
$vdWRules = Thread[$atoms -> $vdWRadii];

(* for making translucent colors *)
If[$VersionNumber >= 10.,

    translucentColorQ[col_] := ColorQ[col] && With[{l = Length[col]}, Switch[Head[col],
                                                                            GrayLevel, l == 2,
                                                                            RGBColor | Hue, l == 4,
                                                                            CMYKColor, l == 5,
                                                                            XYZColor | LABColor | LUVColor | LCHColor, l == 4,
                                                                            _, False] && 0 <= Last[col] < 1],
                                       
    translucentColorQ[col_] := With[{l = Length[col]},
                                                    Switch[Head[col],
                                                               GrayLevel, l == 2,
                                                               RGBColor | Hue, l == 4,
                                                               CMYKColor, l == 5,
                                                               _, False] && 0 <= Last[col] < 1]
]

makeTranslucent[dir_] := If[translucentColorQ[dir] || Head[dir] === Opacity, dir, Opacity[0.6, dir]]

(* Hill system order for atoms *)
hillOrder["C", _] = 1; hillOrder[_, "C"] = -1;
hillOrder["H", _] = 1; hillOrder[_, "H"] = -1;
hillOrder["D", _] = 1; hillOrder[_, "D"] = -1;
hillOrder["T", _] = 1; hillOrder[_, "T"] = -1;
hillOrder[s1_, s2_] := Order[s1, s2];

(* impose time and memory constraints on an evaluation *)
SetAttributes[constrainedEvaluate, HoldFirst];
constrainedEvaluate[expr_, fun_, mem_, time_] := Module[{m = mem, t = time},
                 Switch[m,
                            _Quantity, m = Check[QuantityMagnitude[UnitConvert[m, "Bytes"]], Return[$Failed, Module], Quantity::compat],
                            _?NumericQ | Infinity, m = Round[Abs[m]],
                            _, Return[$Failed, Module]];
                 Switch[t,
                            _Quantity, t = Check[QuantityMagnitude[UnitConvert[t, "Seconds"]], Return[$Failed, Module], Quantity::compat],
                            _?NumericQ | Infinity, t = Round[Abs[t]],
                            _, Return[$Failed, Module]];
                 TimeConstrained[MemoryConstrained[expr, m, Message[fun::mem, m]; $Failed], t, Message[fun::time, t]; $Failed]]

(* percent-encoding for URLs *)
percentEncode[s_String] := If[$VersionNumber >= 10., StringReplace[URLEncode[s], "+" -> "%20"], 
                                             StringReplace[s, RegularExpression["[^\\w/\\?=:._~]"] :> 
                                                                      StringJoin["%", ToUpperCase[IntegerString[First[ToCharacterCode["$0"]], 16]]]]]

(* find normal of best-fit plane *)
getNormal[mat_?MatrixQ, sv : (True | False) : False] := With[{svd = SingularValueDecomposition[Standardize[mat, Mean, 1 &], {-1}, Tolerance -> 0]}, 
     If[TrueQ[sv], FlattenAt[Flatten[Rest[svd], {1, 3}], 1], Flatten[Last[svd]]]]

(* validate input format for viewer *)
validateChem[{vertexTypes_, edgeRules_, edgeTypes_, atomPositions_}] :=
           With[{dim = Dimensions[atomPositions], e = Length[edgeRules], v = Length[vertexTypes]},
                   (e == Length[edgeTypes] && v == First[dim]) &&
                   VectorQ[vertexTypes, StringQ] && Complement[vertexTypes, $atoms] === {} && 
                   MatchQ[edgeRules, {(_Integer?Positive -> _Integer?Positive) ..} | {}] && Max[List @@@ edgeRules] <= v &&
                   MatchQ[edgeTypes, {("Single" | "Double" | "Triple" | "Aromatic") ..} | {}] &&
                   ((MatrixQ[atomPositions, NumberQ] ||
                      (MatrixQ[atomPositions, QuantityUnits`NumberQuantityQ] &&
                       MatrixQ[QuantityUnit[atomPositions], UnitDimensions[#] === {{"LengthUnit", 1}} &])) &&
                    MatchQ[Last[dim], 2 | 3] && Positive[Norm[atomPositions, Infinity]])]

validateChem[__] := False

(* for XYZ, PDB, CIF formats *)
inferBonds[vertexTypes_List, atomPositions_List?MatrixQ] := inferBonds[vertexTypes, atomPositions, 40. $conversionFactor, 25. $conversionFactor, 8]

inferBonds[vertexTypes_List, atomPositions_List?MatrixQ, minDistance_, bondTolerance_, maxNeighbors_Integer] /; Length[vertexTypes] == Length[atomPositions] :=
       Module[{pos = atomPositions, atomP, curAtom, maxD, near, nearD, neighbors, nF, rads},
                   If[MatrixQ[pos, QuantityUnits`NumberQuantityQ], pos = QuantityMagnitude[UnitConvert[pos, $DefaultLengthUnit]]];
                   pos = SetPrecision[pos, MachinePrecision];
                   rads = N[vertexTypes /. Thread[$atoms -> ($elementCovalentRadii $conversionFactor)]];
                   maxD = Max[rads] + bondTolerance;
  
                   nF = Nearest[pos -> If[$VersionNumber >= 11.1, "Index", Automatic]];
                   neighbors = Function[idx,
                                                   curAtom = rads[[idx]]; atomP = pos[[idx]];
                                                   near = DeleteCases[nF[atomP, {maxNeighbors, curAtom + maxD}], idx];
                                                   nearD = SquaredEuclideanDistance[#, atomP] & /@ pos[[near]];
                                                   Pick[near, MapThread[minDistance < #1 < #2 &,
                                                                                    {nearD, (curAtom + rads[[near]] + bondTolerance)^2}]]];
  
                   Union[Sort /@ Flatten[Array[Thread[# -> neighbors[#]] &, Length[vertexTypes]]]]]

(* regular expressions for sanity-checking InChI and InChIKey *)
inchirx = RegularExpression["^(InChI=1)(S*)/([a-zA-Z0-9.])+/(.+)"];
inkeyrx = RegularExpression["^(InChIKey=)?([A-Z]{14})-([A-Z]{10})-([A-Z])$"];

(* SMILES sanity check *)
smilesElems = With[{syms = Drop[$atoms, -2], pos = Transpose[{{1, 5, 6, 7, 8, 9, 15, 16, 19, 23, 39, 53, 74, 92}}]}, Delete[syms, pos] ~Join~ Extract[syms, pos]];
validSymbols = smilesElems ~Join~ {"c", "n", "o", "p", "s", "as", "se"} ~Join~ {"-", "=", "#", ":", "(", ")", "[", "]", "%", ".", "/", "\\", "@", "+", DigitCharacter};

validSMILESQ[s_] := StringQ[s] && ! StringFreeQ[s, LetterCharacter] && s === StringJoin[StringCases[s, validSymbols]]

(* format numbers for Tooltip *)
numDisplay[x_?NumberQ] := If[0.001 < Abs[x] < 999., x, ScientificForm[x, ExponentStep -> If[Abs[x] > 1, 3, 1]]]

(* incircle of a polygon, for rendering aromatic rings *)
incirc[pts_?MatrixQ, scale_: 1] := Block[{cen, cp, nrm, psh, rad, tf},
   cen = Mean[pts]; psh = Standardize[pts, Mean, 1 &];
   rad = scale Mean[Norm /@ psh] Cos[Pi/Length[pts]];
   nrm = Flatten[Last[SingularValueDecomposition[psh, {-1}, Tolerance -> 0]]];
   tf = Quiet[Composition[TranslationTransform[cen], RotationTransform[{{0, 0, 1}, nrm}]], RotationTransform::spln];
   cp = tf[rad {{1, 0, 0}, {1, 1, 0}, {-1, 1, 0}, {-1, 0, 0}, {-1, -1, 0}, {1, -1, 0}, {1, 0, 0}}];
   BSplineCurve[cp, SplineDegree -> 2, SplineKnots -> {0, 0, 0, 1/4, 1/2, 1/2, 3/4, 1, 1, 1}, SplineWeights -> {1, 1/2, 1/2, 1, 1/2, 1/2, 1}]]

(* construct molecular formula *)
buildFormula[alist : {__String}] := Row[Subscript @@@ Sort[Tally[alist], hillOrder[#1[[1]], #2[[1]]] == 1 &] /. Subscript[s_, 1] :> s]

(* for computing molecular weight *)
SetAttributes[atomicMass, Listable];
atomicMass[s_String] := With[{u = If[$VersionNumber >= 9., Quantity["AtomicMassUnit"], 1]}, Switch[s, "D", 2. u, "T", 3. u, _, ElementData[s, "AtomicMass"]]]

buildMolarMass[alist : {__String}] := With[{s = Sum[at[[2]] atomicMass[at[[1]]], {at, Tally[alist]}]}, 
                                                                If[$VersionNumber >= 9., UnitConvert[s Quantity["AvogadroConstant"], "Grams"/"Moles"], ToString[s] <> " g/mol"]]

(* molecule manipulation utilities *)

(* centroid of a molecule *)
MoleculeCentroid[{atoms : {__String}, edgeRules : ({} | None), edgeTypes : ({} | None), coords_?MatrixQ}, opts___] :=
            With[{wts = If[$VersionNumber >= 9., QuantityMagnitude, Identity][atomicMass[atoms]]}, wts.coords/Total[wts]]

MoleculeCentroid[{atoms : {__String}, edgeRules : {__Rule}, edgeTypes_List, coords_?MatrixQ}, opts : OptionsPattern[{"IncludeIsolatedAtoms" -> True}]] :=
            Module[{al = atoms, pts = coords, pos, wts},
                         If[! TrueQ[OptionValue["IncludeIsolatedAtoms"]],
                             pos = Transpose[{Complement[Range[Length[atoms]], Union[Flatten[List @@@ edgeRules]]]}];
                             al = Delete[al, pos]; pts = Delete[pts, pos]];
                         wts = If[$VersionNumber >= 9., QuantityMagnitude, Identity][atomicMass[al]];
                         wts.pts/Total[wts]]

(* center and reorient a molecule *)
MoleculeNormalize[mol : {{__String}, {__Rule}, _List, _?MatrixQ}] :=
            With[{cen = MoleculeCentroid[mol], nrm = getNormal[QuantityMagnitude[Last[mol]]]}, 
                    Quiet[MapAt[Composition[RotationTransform[{nrm, {0, 0, 1}}], TranslationTransform[-cen],
                                                           If[MatrixQ[#, QuantityUnits`NumberQuantityQ],
                                                               QuantityMagnitude[UnitConvert[#, $DefaultLengthUnit]], #] &], mol, 4], RotationTransform::spln]]

(* split molecule into connected components *)
MoleculeSplit[mol : {{__String}, {} | None, {} | None, _?MatrixQ}, _] := mol

MoleculeSplit[{atoms : {__String}, edgeRules : {__Rule}, edgeTypes_List, coords_?MatrixQ}] :=
            Module[{mcs, molExtract, skel},
                        skel = Graph[Range[Length[atoms]], edgeRules, DirectedEdges -> False, EdgeWeight -> edgeTypes, GraphLayout -> None];
                        molExtract[g_Graph] := With[{id = VertexList[g], ed = EdgeList[g]},
                                                                    {atoms[[id]], Composition[Sort, Rule] @@@ (ed /. MapIndexed[#1 -> First[#2] &, id]), 
                                                                      PropertyValue[{skel, #}, EdgeWeight] & /@ ed, coords[[id]]}];
                        molExtract /@ ConnectedGraphComponents[skel]]

(* join two or more molecules *)
MoleculeJoin[mols : (({{__String}, {__Rule} | {} | None, {__String} | {} | None, _?MatrixQ} | {} | $Failed | _Missing) ..)] :=
            Module[{molList = DeleteCases[{mols}, {} | $Failed | _Missing, {1}], atoms},
                        atoms = molList[[All, 1]];
                        {Flatten[atoms],
                          Flatten[MapThread[Replace[#1, k_Integer :> k + #2, {2}] &,
                                                      {molList[[All, 2]], Prepend[Accumulate[Length /@ Most[atoms]], 0]}]], 
                          Flatten[molList[[All, 3]]],
                          Flatten[If[MatrixQ[#, QuantityUnits`NumberQuantityQ], QuantityMagnitude[UnitConvert[#, $DefaultLengthUnit]], #] & /@
                                     molList[[All, 4]], 1]}]

(* bring distant isolated atoms closer *)
ShiftIsolatedAtoms[mol : {{__String}, {} | None, {} | None, _?MatrixQ}, _] := mol

ShiftIsolatedAtoms[{atoms : {__String}, edgeRules : {__Rule}, edgeTypes_List, coords_?MatrixQ}, h_: 1.6] :=
      Module[{al = atoms, pts = coords, cent, dist, isd, isf, isol, pos, wts},
                 pos = Transpose[{Complement[Range[Length[atoms]], Union[Flatten[List @@@ edgeRules]]]}];
                 al = Delete[al, pos]; pts = Delete[pts, pos];
                 wts =  If[$VersionNumber >= 9., QuantityMagnitude, Identity][atomicMass[al]];
                 cent = wts.pts/Total[wts]; dist = Max[Norm[cent - #] & /@ pts];
                 isol = Extract[coords, pos]; isd = Norm[cent - #] & /@ isol; isf = Clip[isd, {0, h dist}]/isd;
                 {atoms, edgeRules, edgeTypes, ReplacePart[coords, Thread[pos -> (1 - isf) ConstantArray[cent, Length[isol]] + isf isol]]}]

(* rotate about a single bond *)
BondRotate[{atoms : {__String}, edgeRules : {__Rule}, edgeTypes_List, coords_?MatrixQ}, bond : (_Integer -> _Integer), th_?NumericQ] :=
       Module[{b0, b1, bpos, btyp, cc, ed, p0, sel},
                   bpos = Position[edgeRules, bond | Reverse[bond]];

                   If[bpos =!= {}, bpos = First[bpos], Message[BondRotate::nobnd, bond]; Return[$Failed, Module]];
                   If[(btyp = Extract[edgeTypes, bpos]) =!= "Single",
                       Message[BondRotate::mult, ToLowerCase[btyp], bond]; Return[$Failed, Module]];
                   If[Length[edgeRules] == 1, Return[{atoms, edgeRules, edgeTypes, coords}, Module]];

                   {b0, b1} = List @@ bond; ed = Delete[edgeRules, bpos];
                   cc = ConnectedComponents[Graph[ed, DirectedEdges -> False, GraphLayout -> None]];
                   If[Length[cc] > 1, sel = Flatten[Select[cc, MemberQ[#, b1] &]],
                       Message[BondRotate::ring, bond]; Return[$Failed, Module]];
                   p0 = coords[[b0]];
                   {atoms, edgeRules, edgeTypes, 
                     ReplacePart[coords, Thread[sel -> RotationTransform[th, coords[[b1]] - p0, p0][coords[[sel]]]]]}]

(* apply a geometric transform to a molecule *)

MoleculeTransform[mol : {{__String}, {__Rule} | {} | None, {__String} | {} | None, _?MatrixQ}, tfun_TransformationFunction] := Module[{pos, tmat},
            pos = Last[mol];
            If[MatrixQ[pos, QuantityUnits`NumberQuantityQ], 
                pos = QuantityMagnitude[UnitConvert[pos, $DefaultLengthUnit]]];
            tmat = Drop[TransformationMatrix[tfun], -1, -1];

            If[Last[Dimensions[pos]] != Length[tmat], Message[MoleculeTransform::icdims, tmat]; Return[$Failed, Module]];
            If[! OrthogonalMatrixQ[tmat], Message[MoleculeTransform::notort, tmat]; Return[$Failed, Module]];
            If[Det[tmat] != 1, Message[MoleculeTransform::notrgd, tmat]];

            ReplacePart[mol, {-1} -> tfun[pos]]]

MoleculeTransform[mol_, vec_?VectorQ] := With[{pos = Last[mol]},
            If[Length[vec] == Last[Dimensions[pos]],
                MoleculeTransform[mol, TranslationTransform[vec]], Message[MoleculeTransform::icdims, vec]; $Failed]]

MoleculeTransform[mol_, mat_?MatrixQ] := With[{pos = Last[mol]},
            If[(Equal @@ Dimensions[mat]) && Length[mat] == Last[Dimensions[pos]],
               MoleculeTransform[mol, AffineTransform[mat]], Message[MoleculeTransform::icdims, mat]; $Failed]]

MoleculeTransform[mol_, {mat_?MatrixQ, vec_?VectorQ}] := 
            With[{pos = Last[mol], tf = Check[AffineTransform[{mat, vec}], $Failed, AffineTransform::idimr]},
                    If[tf =!= $Failed && Length[vec] == Last[Dimensions[pos]], MoleculeTransform[mol, tf],
                        Message[MoleculeTransform::icdims, Drop[TransformationMatrix[tf], -1, -1]]; $Failed]]

(* basic importing functions *)
molFromFile[mol_String, type_String : "MOL"] := Block[{res},
     If[$VersionNumber < 12., 
         res = Check[Import[mol, {type, {"VertexTypes", "EdgeRules", "EdgeTypes", "VertexCoordinates"}}],
                            Message[MoleculeViewer::badfile]; Return[$Failed, Block]];
         res = MapAt[(# $conversionFactor) &, res, {4}];
         If[validateChem[res], res, Message[MoleculeViewer::badfile]; $Failed],
         res = Check[Import[mol, type], Message[MoleculeViewer::badfile]; Return[$Failed, Block]];
         processMoleculeObject[res]]]

molFromString[mol_String] := Block[{res},
     If[$VersionNumber < 12.,
         res = Check[ImportString[mol, {"MOL", {"VertexTypes", "EdgeRules", "EdgeTypes", "VertexCoordinates"}}], Return[$Failed, Block]];
         res = MapAt[(# $conversionFactor) &, res, {4}];
         If[validateChem[res], res, $Failed],
         processMoleculeObject[ImportString[mol, "MOL"]]]]

processChemicalFile[fName_String, type : (_String | Automatic) : Automatic] := Module[{atoms, bonds, ext, fileName, het, ida, idx, pos, prop, prts, res},

           If[FileExistsQ[fName], fileName = fName,
               fileName = FindFile[fName];
               If[fileName === $Failed, Message[MoleculeViewer::fnfnd, fName]; Return[$Failed, Module]]];

           ext = type;
           If[ext === Automatic, ext = FileFormat[fileName];
               If[ext === "", Message[MoleculeViewer::fftype, fName]; Return[$Failed, Module]]];
           (* special handling for CIF *)
           If[StringMatchQ[ext, "CIF", IgnoreCase -> True], ext = "MMCIF"]; If[TrueQ[$debug], Print[ext]];

           prop = Quiet[Import[fileName, {ext, "Elements"}]];
           If[prop === $Failed || FreeQ[prop, "VertexTypes" | "EdgeRules" | "EdgeTypes" | "VertexCoordinates"],
               Message[MoleculeViewer::ftconv, fName]; Return[openBabelConvert[fileName, FileExtension[fileName], True, Automatic, 60], Module]];

           If[FreeQ[prop, "EdgeRules" | "EdgeTypes"],

               (* xyz, CIF or PDB format *)
               If[TrueQ[$debug], Print["XYZ/PDB/CIF"]];
               prts = {"VertexTypes", "VertexCoordinates"};
               If[StringMatchQ[ext, "PDB", IgnoreCase -> True],
                   If[TrueQ[$getAllPDBComponents], AppendTo[prts, "ResidueIndex"], prts = Join[prts, {"ResidueIndex", "AdditionalIndex"}]]];
               res = Check[Import[fileName, {ext, prts}], Message[MoleculeViewer::badfile]; Return[$Failed, Module]];
               {atoms, pos, idx, ida} = PadRight[res, 4, Missing["NotAvailable"]];
               pos = If[MatrixQ[pos, NumberQ], pos $conversionFactor, QuantityMagnitude[UnitConvert[pos, $DefaultLengthUnit]]];
               If[ListQ[idx], idx = Transpose[{Sort[Flatten[idx]]}]; (* separate out main and additional atoms for PDB *)
                   If[TrueQ[$getAllPDBComponents], het = {Delete[atoms, idx], Delete[pos, idx]},
                       If[ListQ[ida], ida = Transpose[{Sort[Flatten[ida]]}];
                           het = {Extract[atoms, ida], Extract[pos, ida]},
                           het = {{}, {}}]];
                   If[het =!= {{}, {}},
                       res = {Extract[atoms, idx], Extract[pos, idx]};
                       If[VectorQ[atoms, StringQ] && MatrixQ[pos], {res, het}, Message[MoleculeViewer::badfile]; $Failed],
                       If[VectorQ[atoms, StringQ] && MatrixQ[pos], {{atoms, pos}, Missing["NotAvailable"]}, Message[MoleculeViewer::badfile]; $Failed]],
                   bonds = inferBonds[atoms, pos];
                   res = {atoms, bonds, ConstantArray["Single", Length[bonds]], pos};
                   If[validateChem[res], res, Message[MoleculeViewer::badfile]; $Failed]],

               (* normal mode *)
               Which[StringMatchQ[ext, "MOL", IgnoreCase -> True], If[TrueQ[$debug], Print["MOL"]]; molFromFile[fileName], (* special treatment for MOL files *)

                         StringMatchQ[ext, "HIN", IgnoreCase -> True],
                         If[TrueQ[$debug], Print["HIN"]];
                         res = Check[Import[fileName, {ext, {"VertexTypes", "EdgeRules", "EdgeTypes", "VertexCoordinates"}}],
                                            Message[MoleculeViewer::badfile]; Return[$Failed, Module]];
                         (* extra factor for Å *) res = Map[Function[r, MapAt[(100 $conversionFactor #) &, r, -1]], Transpose[res]];
                         If[validateChem[#], #, Message[MoleculeViewer::badfile]; $Failed] & /@ res,

                         True,
                         If[TrueQ[$debug], Print["SDF/MOL2"]];
                         If[$VersionNumber < 12.,
                             res = Check[Import[fileName, {ext, {"VertexTypes", "EdgeRules", "EdgeTypes", "VertexCoordinates"}}],
                                                Message[MoleculeViewer::badfile]; Return[$Failed, Module]];
                             If[StringMatchQ[ext, "MOL2" | "SDF", IgnoreCase -> True],
                                 res = Map[Function[r, MapAt[(# $conversionFactor) &, r, -1]], Transpose[res]]];
                             If[StringMatchQ[ext, "MOL2", IgnoreCase -> True], res = res /. "Amide" -> "Single"; (* special treatment for MOL2 unconventional bonds *)
                                 If[! FreeQ[res, "Dummy" | "Unknown" | "Not Connected"],
                                     res = With[{bp = Table[! MatchQ[bb, "Dummy" | "Unknown" | "Not Connected"], {bb, #[[3]]}]},
                                                      ReplacePart[#, Thread[{2, 3} -> {Pick[#[[2]], bp], Pick[#[[3]], bp]}]]] & /@ res]];
                             If[validateChem[#], #, Message[MoleculeViewer::badfile]; $Failed] & /@ res,
                             res = Check[Import[fileName, ext], Message[MoleculeViewer::badfile]; Return[$Failed, Module]];
                             processMoleculeObject /@ res]]]]

processChemicalURL[url_String, type : (_String | Automatic) : Automatic, fileName : (_String | Automatic) : Automatic] := Block[{fName = fileName, fnrx, tmp, urlrx},
           (* regular expression for parsing a URL, from https://tools.ietf.org/html/rfc3986#appendix-B *)
           urlrx = RegularExpression["^(([^:/?#]+):)?(//([^/?#]*))?([^?#]*)(\\?([^#]*))?(#(.*))?"];
           (* regular expression for invalid filename characters *)
           fnrx = RegularExpression["[\":|\\\\<*>/?]+"];

           If[fName === Automatic,
               fName = StringReplace[If[$VersionNumber >= 10.,
                                                      Last[Lookup[URLParse[url], "Path"]],
                                                      Last[StringSplit[First[StringCases[url, urlrx -> "$5"]], {"/", "\\"}]]], fnrx -> "_"]];
           fName = FileNameJoin[{$TemporaryDirectory, fName}]; If[TrueQ[$debug], Print[{url, fName}]];
           Internal`WithLocalSettings[tmp = Check[If[$VersionNumber >= 9., URLSave, Utilities`URLTools`FetchURL][url, fName], $Failed],
                                                   If[TrueQ[$debug], Print[{Round[FileByteCount[fName]/1024., 0.01], IntegerString[FileHash[fName, "MD5"], 16, 32]}]];
                                                   If[tmp =!= $Failed, Quiet[processChemicalFile[tmp, type], {MoleculeViewer::fnfnd}], $Failed],
                                                   If[tmp =!= $Failed, Quiet[DeleteFile[tmp], General::fdnfnd]]]]
                                                                                     
If[$VersionNumber >= 10., (* handling for Entity objects in version 10+ *)

processEntity[ent_Entity] := Module[{res},
           If[EntityTypeName[ent] =!= "Chemical", Return[$Failed, Module]];
           res = EntityValue[ent, EntityProperty["Chemical", #] & /@ {"VertexTypes", "EdgeRules", "EdgeTypes", "AtomPositions"}];
           res = MapAt[(# $conversionFactor) &, res, {4}];
           If[validateChem[res], res, $Failed]];

processEntity[__] := $Failed
]

If[$VersionNumber >= 12., (* handling for native Molecule objects in version 12+ *)

processMoleculeObject[mol_Molecule] := Block[{res},
           res = MoleculeValue[Molecule[mol, IncludeAromaticBonds -> False, IncludeHydrogens -> True],
                                         {"AtomList", "BondList", "AtomCoordinates"}];
           If[! (Head[Last[res]] === StructuredArray && StringMatchQ[Last[res]["Structure"], "QuantityArray"]) ||
               ! MatrixQ[Last[res]], Return[$Failed, Block]];
           res = Replace[res, {Atom[s_] :> s, Bond[i_, t_] :> {Rule @@ i, t}}, {2}];
           res = MapAt[QuantityMagnitude[UnitConvert[#, $DefaultLengthUnit]] &, res, {3}];
           res = MapAt[Composition[Apply[Sequence], Transpose], res, {2}];
           If[validateChem[res], res, $Failed]];

processMoleculeObject[__] := $Failed
]

Options[MakeMOL] = {Method -> Automatic};

MakeMOL[s_String, OptionsPattern[]] := Block[{babQ},
       If[! (testOpenBabel || PacletManager`$AllowInternet), Return[$Failed, Block]];
       babQ = Switch[System`Utilities`StringName[OptionValue[Method]],
                              "Automatic", testOpenBabel || (! PacletManager`$AllowInternet),
                              "OpenBabel", True,
                              "Beaker" | "ChEMBL", False,
                              _, Return[$Failed, Block]];
       If[babQ, obMakeString[s, "MOL"],
           Which[validSMILESQ[s], ChEMBLSMILESToMOL[s],
                     StringMatchQ[s, inchirx], ChEMBLInChIToMOL[s],
                     True, $Failed]]]

MakeMOL[mol_List?validateChem, opts___?OptionQ] := Block[{pts = Last[mol], raw},
       raw = Which[MatrixQ[pts, QuantityUnits`NumberQuantityQ],
                           Check[MapAt[QuantityMagnitude[UnitConvert[#, "Picometers"]] &, mol, {4}], Return[$Failed, Block], Quantity::compat],
                           MatrixQ[pts, NumberQ], MapAt[(#/$conversionFactor) &, mol, {4}],
                           True, Message[MakeMOL::ncnv]; Return[$Failed, Block]];
       ExportString[Thread[{"VertexTypes", "EdgeRules", "EdgeTypes", "VertexCoordinates"} -> raw], {"MOL", "Rules"}]]

MakeMOL[mol_List?validateChem, fn_String, opts___?OptionQ] := Block[{pts = Last[mol], raw},
       raw = Which[MatrixQ[pts, QuantityUnits`NumberQuantityQ],
                           Check[MapAt[QuantityMagnitude[UnitConvert[#, "Picometers"]] &, mol, {4}], Return[$Failed, Block], Quantity::compat],
                           MatrixQ[pts, NumberQ], MapAt[(#/$conversionFactor) &, mol, {4}],
                           True, Message[MakeMOL::ncnv]; Return[$Failed, Block]];
       Export[fn, Thread[{"VertexTypes", "EdgeRules", "EdgeTypes", "VertexCoordinates"} -> raw], {"MOL", "Rules"}]]

Options[MakeInChIKey] = {Method -> Automatic};

MakeInChIKey[s_String, OptionsPattern[]] := Block[{babQ, res},
       If[! (testOpenBabel || PacletManager`$AllowInternet), Return[$Failed, Block]];
       If[StringMatchQ[s, inkeyrx], Return[s, Block]];
       babQ = Switch[System`Utilities`StringName[OptionValue[Method]],
                              "Automatic", testOpenBabel || (! PacletManager`$AllowInternet),
                              "OpenBabel", True,
                              "Beaker" | "ChEMBL", False,
                              _, Return[$Failed, Block]];
       If[babQ, res = If[validSMILESQ[s] || StringMatchQ[s, inchirx], obMakeString[s, "InChIKey"], obLineNotation[s, "InChIKey"]];
           If[StringQ[res], StringTrim[res], res],
           ChEMBLMakeInChIKey[s]]]

MakeInChIKey[mol_List?validateChem, opts___] := MakeInChIKey[MakeMOL[mol, opts], opts]

(* OpenBabel interfacing routines *)

testOpenBabel := testOpenBabel = Quiet[Check[Import["!obabel -H", "Text"] =!= "", Message[RunOpenBabel::nobab]; False]]
$OpenBabelVersion := $OpenBabelVersion = Quiet[Check[StringReplace[Import["!obabel -V", "Text"], " -- " -> "\n"], Message[RunOpenBabel::nobab]; $Failed]]

(* defaults for Open Babel methods *)
$obConOpts = {"GeneratedConformers" -> 20, "Method" -> "Genetic"};
$obConGOpts = {"Children" -> Automatic, "Mutability" -> Automatic, "Converge" -> Automatic, "Score" -> Automatic};
$obMinOpts = {"CutOff" -> Automatic, "ForceField" -> "MMFF94s", "MaxIterations" -> 2000, "Newton" -> Automatic, "SteepestDescent" -> Automatic, "Tolerance" -> Automatic, "UpdateNonBondedFrequency" -> 10};

$obForceFields = {"GAFF", "Ghemical", "MMFF94", "MMFF94s", "UFF"};

makeStringRules[expr_List] := Block[{o, s}, Flatten[expr /. (o_ -> s_) :> (System`Utilities`StringName[o] -> s)]]

makeOBSwitches[met_, mopt_List] := Module[{mok, mon, res, tmp},
        If[! StringQ[met], Message[RunOpenBabel::bdmtd, met]; Return[$Failed, Module]];
        mon = makeStringRules[mopt]; res = {"--" <> ToLowerCase[met]};
        Switch[met,
                   "Generate", res = {},
                   "Conformer",
                   mok = FilterRules[Join[mon, $obConOpts], $obConOpts];
                   AppendTo[res, "--nconf " <> IntegerString[Round["GeneratedConformers" /. mok]]];
                   tmp = "Method" /. mok;
                   tmp = If[! ListQ[tmp], System`Utilities`StringName[tmp], 
                                 MapAt[System`Utilities`StringName, tmp, {1}]];
                   Switch[tmp,
                              "Automatic",
                              AppendTo[res, "--weighted"],
                              "Random" | "Systematic" | "Weighted",
                              AppendTo[res, "--" <> ToLowerCase[tmp]],
                              "ForceField",
                              AppendTo[res, "--ff MMFF94s"],
                              {"ForceField", _Rule},
                              tmp = makeStringRules[Rest[tmp]];
                              AppendTo[res, "--ff " <> ("Method" /. tmp /. Automatic | "Automatic" -> "MMFF94s")],
                              Alternatives @@ $obForceFields,
                              AppendTo[res, "--ff " <> tmp],
                              "Genetic" | {"Genetic", __Rule},
                              If[ListQ[tmp],
                                  tmp = makeStringRules[Rest[tmp]];
                                  tmp = FilterRules[Join[tmp, $obConGOpts], $obConGOpts],
                                  tmp = $obConGOpts];
                              AppendTo[res, "--score " <> ToLowerCase["Score" /. tmp /. Automatic -> "Energy"]];
                              tmp = DeleteCases[MapThread[If[MatchQ[#2, Automatic | "Automatic"], "", #1 <> IntegerString[Max[0, Round[#2]]]] &,
                                                                            {{"--children ", "--mutability ", "--converge "}, {"Children", "Mutability", "Converge"} /. tmp}], ""];
                              res = Join[res, tmp],
                              _, Message[RunOpenBabel::nsupm, tmp, "Conformer"]; Return[$Failed, Module]],
                   "Minimize",
                   mok = FilterRules[Join[mon, $obMinOpts], $obMinOpts];
                   tmp = "ForceField" /. mok /. Automatic | "Automatic" -> "MMFF94s";
                   If[MemberQ[$obForceFields, tmp], AppendTo[res, "--ff " <> tmp],
                       Message[RunOpenBabel::bdff, tmp]; Return[$Failed, Module]];
                   tmp = Round["MaxIterations" /. mok];
                   If[Internal`PositiveIntegerQ[tmp], AppendTo[res, "--steps " <> IntegerString[tmp]],
                       Message[RunOpenBabel::ioppm, "MaxIterations", tmp]; Return[$Failed, Module]];
                   tmp = "Tolerance" /. mok /. Automatic -> 1*^-7;
                   If[NumberQ[tmp] && Positive[tmp], AppendTo[res, "--crit " <> ToString[CForm[tmp], InputForm]],
                       Message[RunOpenBabel::tol2, tmp]; Return[$Failed, Module]];
                   tmp = Round["UpdateNonBondedFrequency" /. mok];
                   If[Internal`PositiveIntegerQ[tmp], AppendTo[res, "--freq " <> IntegerString[tmp]],
                       Message[RunOpenBabel::ioppm, "UpdateNonBondedFrequency", tmp]; Return[$Failed, Module]];
                   If[TrueQ["SteepestDescent" /. mok /. Automatic -> False], AppendTo[res, "--sd"]];
                   If[TrueQ["Newton" /. mok /. Automatic -> False], AppendTo[res, "--newton"]];
                   If[(tmp = "CutOff" /. mok /. Automatic -> False) =!= False,
                       Switch[tmp,
                                  True, AppendTo[res, "--cut"],
                                  {_?NumberQ} | {_?NumberQ, _?NumberQ}, 
                                  res = Join[res, Prepend[MapThread[StringJoin, {{"--rvdw ", "--rele "}, ToString /@ PadRight[tmp, 2, 10]}], "--cut"]],
                                  _, Message[RunOpenBabel::erropts, tmp, "CutOff"]; Return[$Failed, Module]]],
                   _, Message[RunOpenBabel::bdmtd, met, Automatic]; res = $Failed];
        res]

(* format conversion *)
openBabelConvert[fileName_String, Rule[ex_String, out_String], make3D : (True | False) : True, mc_: Automatic, tc_: Automatic, msgQ : (True | False) : True] :=
       Module[{ext = ToLowerCase[ex], fmt = ToLowerCase[out], pars = {"-c", "-h"}, args, ec, msg, proc, res},
                   If[! testOpenBabel, Return[$Failed, Module]];
                   If[make3D, PrependTo[pars, "--gen3d"], If[MatchQ[fmt, "smiles" | "inchi" | "inchikey"], pars = {}, PrependTo[pars, "--gen2d"]]];
                   constrainedEvaluate[If[$VersionNumber >= 10.,
                                                     args = {"obabel", "-i", ext, fileName, "-o" <> fmt} ~Join~ pars;
                                                     proc = Check[RunProcess[args], Return[$Failed, Module]];
                                                     res = proc["StandardOutput"]; ec = proc["ExitCode"]; msg = proc["StandardError"];
                                                     If[ec == 0 && res =!= "",
                                                         If[make3D, Check[molFromString[res], If[msgQ, Echo[StringTrim[msg]]]; $Failed], res],
                                                         msg = StringTrim[StringReplace[msg, {"0 molecules converted" -> "", "\r" | "\n" -> " "}]];
                                                         If[msgQ,
                                                             If[msg =!= "", Echo[msg],
                                                                 If[ec != 0, Echo["OpenBabel exit code: " <> IntegerString[ec]]]]]; $Failed],
                                                     args = StringJoin[Riffle[{"!obabel", "-i", ext, fileName, "-o" <> fmt} ~Join~ pars, " "]];
                                                     If[make3D, Check[molFromFile[args], If[msgQ, Print[StringTrim[msg]]]; $Failed], res]], RunOpenBabel,
                                                 mc /. Automatic -> 536870912 (* 512 MB *), tc /. Automatic -> 120 (* 2 min. *)]]

openBabelConvert[fileName_String, ex_String, rest___] := openBabelConvert[fileName, ex -> "MOL", rest]

(* generate data from MOL string *)
obConvertMOL[mol_String, out_String, rest___] := Block[{tmp = FileNameJoin[{$TemporaryDirectory, IntegerString[Hash[mol, "MD5"], 16, 32] <> ".mol"}]},
                                                                                     Internal`WithLocalSettings[tmp = Check[Export[tmp, mol, "Text"], $Failed],
                                                                                                                             If[tmp =!= $Failed, openBabelConvert[tmp, "MOL" -> out, rest], $Failed],
                                                                                                                             If[tmp =!= $Failed, Quiet[DeleteFile[tmp], General::fdnfnd]]]]

obConvertMOL[mol_List?validateChem, rest___] := obConvertMOL[MakeMOL[mol], rest]

(* generate 3D coordinates from MOL string *)
obMOLTo3D[mol_String, rest___] := obConvertMOL[mol, "MOL", True, rest]
obMOLTo3D[mol_List?validateChem, rest___] := obMOLTo3D[MakeMOL[mol], rest]

(* generate line notation from MOL string *)
obLineNotation[mol_String, type : ("SMILES" | "InChI" | "InChIKey") : "SMILES", rest___] := With[{res = obConvertMOL[mol, type, False, rest]}, If[StringQ[res], StringTrim[res], res]]
obLineNotation[mol_List?validateChem, rest__] := obLineNotation[MakeMOL[mol], rest]

(* generate InChIKey or 2D MOL string from line notation *)
obMakeString[chem_String, type : ("InChIKey" | "MOL") : "InChIKey", mc_: Automatic, tc_: Automatic] := Module[{args, ec, incq, msg, otyp, proc, res, sw},
    If[! testOpenBabel, Return[$Failed, Module]];
    incq = If[StringMatchQ[chem, inchirx], "-iinchi", Nothing]; otyp = ToLowerCase[type];
    sw = If[StringMatchQ[type, "MOL"], {"--gen2d", "-c", "-h"}, {}];
    constrainedEvaluate[If[$VersionNumber >= 10., 
                                      args = {"obabel", incq, "-:" <> chem, "-o" <> otyp} ~Join~ sw;
                                      proc = Check[RunProcess[args], Return[$Failed, Module]];
                                      res = proc["StandardOutput"]; ec = proc["ExitCode"]; msg = proc["StandardError"];
                                      If[ec == 0 && res =!= "",
                                          Check[ImportString[res, "Text"], Echo[StringTrim[msg]]; $Failed],
                                          msg = StringTrim[StringReplace[msg, {"0 molecules converted" -> "", "\r" | "\n" -> " "}]];
                                          If[msg =!= "", Echo[msg], 
                                              If[ec != 0, Echo["OpenBabel exit code: " <> IntegerString[ec]]]]; $Failed],
                                      args = StringJoin[Riffle[{"!obabel", incq, "-:" <> chem, "-o" <> otyp} ~Join~ sw, " "]];
                                      Check[Import[args, "Text"], Print[StringTrim[msg]]; $Failed]], RunOpenBabel,
                                  mc /. Automatic -> 67108864 (* 64 MB *), tc /. Automatic -> 30 (* 30 sec. *)]]

Options[RunOpenBabel] = {MemoryConstraint -> Automatic, Method -> Automatic, ProcessDirectory -> Inherited,
                                         ProcessEnvironment -> Inherited, RunProcess -> Automatic, TimeConstraint -> Automatic, Verbose -> True};

RunOpenBabel[mol_List?validateChem, opts : OptionsPattern[]] := obMOLTo3D[MakeMOL[mol], OptionValue[MemoryConstraint], OptionValue[TimeConstraint], OptionValue[Verbose]]

RunOpenBabel[list : {__String}, opts___] := RunOpenBabel[#, opts] & /@ list

RunOpenBabel[file_String /; Quiet[FindFile[file] =!= $Failed], opts : OptionsPattern[]] :=
     With[{fn = FindFile[file]}, openBabelConvert[fn, FileExtension[fn], True, OptionValue[MemoryConstraint], OptionValue[TimeConstraint], OptionValue[Verbose]]]

RunOpenBabel[chem_String, opts : OptionsPattern[]] := Module[{args, autoSet, ec, incq, met, mopt, msg, opo, pars, proc, res},
      If[! testOpenBabel, Return[$Failed, Module]];

      autoSet = {"Minimize", {"ForceField" -> "MMFF94s", "MaxIterations" -> 2000, "Tolerance" -> 1.*^-7}};
      met = OptionValue[Method];
      If[ListQ[met],
          {met, mopt} = {First[met], Flatten[Rest[met]]},
          mopt = {}];
      met = System`Utilities`StringName[met];
      If[met === "Automatic", {met, mopt} = autoSet];
      pars = makeOBSwitches[met, mopt];
      If[pars === $Failed, pars = makeOBSwitches @@ autoSet];

      incq = If[StringMatchQ[chem, inchirx], "-iinchi", Nothing];

      constrainedEvaluate[If[TrueQ[OptionValue[RunProcess] /. Automatic -> True] && $VersionNumber >= 10.,
                                        args = {"obabel", incq, "-:"<> chem, "-omol", "--gen3d", "-c", "-h"} ~Join~ pars;
                                        If[TrueQ[$debug], Print[Defer[RunProcess][args]]];
                                        proc = Check[RunProcess[args, FilterRules[Join[{opts}, Options[RunOpenBabel]], Options[RunProcess]]], 
                                                             Return[$Failed, Module]];
                                        res = proc["StandardOutput"]; ec = proc["ExitCode"]; msg = proc["StandardError"];
                                        If[ec == 0 && res =!= "", 
                                            opo = StringPosition[res, "OpenBabel"];
                                            If[opo =!= {},
                                                res = molFromString[StringDrop[res, opo[[1, 1]] - 3]];
                                                If[validateChem[res], res, $Failed],
                                                $Failed],
                                            If[TrueQ[OptionValue[Verbose] /. Automatic -> True],
                                                msg = StringTrim[StringReplace[msg, {"0 molecules converted" -> "", "\r" | "\n" -> " "}]];
                                                If[msg =!= "", Echo[msg], 
                                                    If[ec != 0, Echo["OpenBabel exit code: " <> IntegerString[ec]]]]]; $Failed],
                                        args = StringJoin[Riffle[{"!obabel", incq, "-:" <> chem, "-omol", "--gen3d", "-c", "-h"} ~Join~ pars, " "]];
                                        If[TrueQ[$debug], Print[Defer[Import][args, "Text"]]];
                                        res = Check[Import[args, "Text"], $Failed];
                                        If[res =!= $Failed,
                                            opo = StringPosition[res, "OpenBabel"];
                                            If[opo =!= {},
                                                res = molFromString[StringDrop[res, opo[[1, 1]] - 3]];
                                                If[validateChem[res], res, $Failed],
                                                $Failed],
                                            $Failed]], RunOpenBabel,
                                    OptionValue[MemoryConstraint] /. Automatic -> 536870912 (* 512 MB *),
                                    OptionValue[TimeConstraint] /. Automatic -> 120 (* 2 min. *)]]

(* Imago OCR interfacing routines *)

testImago := testImago = Quiet[Check[Import["!imago_console", "Text"] =!= "", Message[RunImago::noimg]; False]]

Options[RunImago] = {ImageResolution -> Automatic, ImageSize -> Automatic, MemoryConstraint -> Automatic, ProcessDirectory -> Inherited,
                                   ProcessEnvironment -> Inherited, RasterSize -> Automatic, RunProcess -> Automatic, TimeConstraint -> Automatic, Verbose -> Automatic};

RunImago[x_, opts___] := RunImago[Rasterize[x, "Image", FilterRules[Join[{ColorSpace -> "Grayscale", opts}, Options[RunImago]], Options[Rasterize]]], opts];

RunImago[img_Image, opts : OptionsPattern[{RunImago, Rasterize}]] := Module[{fn, fnIn, fnOut, msg, proc},
     If[! testImago, Return[$Failed]];
  
     fn = IntegerString[Hash[img, "MD5"], 16, 32];
     fnIn = FileNameJoin[{$TemporaryDirectory, fn <> ".png"}];
     fnOut = FileNameJoin[{$TemporaryDirectory, fn <> ".mol"}];
  
     Internal`WithLocalSettings[Export[fnIn, If[ImageColorSpace[img] =!= "Grayscale", ColorConvert[img, "Grayscale"], img], "PNG"],
                                             constrainedEvaluate[If[TrueQ[OptionValue[RunProcess] /. Automatic -> True] && $VersionNumber >= 10.,
                                                                               proc = Check[RunProcess[{"imago_console", "-images", "-o", fnOut, fnIn},
                                                                                                                      FilterRules[Join[{opts}, Options[RunImago]], Options[RunProcess]]],
                                                                                                    Return[$Failed, Module]];
                                                                               If[proc["ExitCode"] == 0, Check[Import[fnOut, "Text"], $Failed],
                                                                                   If[TrueQ[OptionValue[Verbose] /. Automatic -> True],
                                                                                       If[(msg = proc["StandardError"]) =!= "", Echo[msg], 
                                                                                           Echo[ToString[proc["StandardOutput"]] <> "\n" <>
                                                                                                   "Imago exit code: " <> IntegerString[proc["ExitCode"]]]]];
                                                                                   $Failed],
                                                                               proc = Check[Run["imago_console -images -o", fnOut, fnIn], Return[$Failed, Module]];
                                                                               If[proc == 0, Check[Import[fnOut, "Text"], $Failed], $Failed]], RunImago,
                                                                           OptionValue[MemoryConstraint] /. Automatic -> 536870912 (* 512 MB *),
                                                                           OptionValue[TimeConstraint] /. Automatic -> 120 (* 2 min. *)],
                                             Quiet[DeleteFile[{fnIn, fnOut}], General::fdnfnd]]]

(* CACTUS *)

SetAttributes[GetCACTUS, Listable];

GetCACTUS[arg_String] := Module[{in, url},
     If[$VersionNumber < 9., Message[GetCACTUS::unavail, GetCACTUS]; Return[$Failed, Module]];
     url = "https://cactus.nci.nih.gov/chemical/structure/``/file?format=sdf&get3d=true";
     in = If[StringMatchQ[arg, inchirx], arg, percentEncode[arg]];
     If[TrueQ[$debug], Print[Defer[Import][ToString[StringForm[url, in]], "MOL"]]];
     molFromFile[If[$VersionNumber >= 10., StringTemplate[url] @ in, ToString[StringForm[url, in]]]]]

SetAttributes[GetChemicalData, Listable];

Options[GetChemicalData] = {"MaxResults" -> 1};

GetChemicalData[str_String, out : (Automatic | "SMILES" | "InChI") : Automatic, OptionsPattern[]] := Module[{chk, name, res, tmp},

     name = StringJoin[StringReplacePart[#, ToUpperCase[StringTake[#, 1]], {1, 1}] & /@ StringSplit[str, " " | "-"]];
     If[TrueQ[$debug], Print[Defer[ChemicalData][name]]];
     If[! StringFreeQ[name, "*"],
         chk = Table[ChemicalData[c, "StandardName"], {c, ChemicalData[name]}];
         If[chk === {}, Return[$Failed, Module]],
         If[(chk = ChemicalData[name, "StandardName"]) === $Failed, Return[chk, Module]]];
     chk = Take[Flatten[{chk}], Max[1, Round[OptionValue["MaxResults"] /. All -> Infinity]] /. Infinity -> All];

     If[out === Automatic,
         res = Table[ChemicalData[c, prop], {c, chk}, {prop, {"VertexTypes", "EdgeRules", "EdgeTypes", "AtomPositions"}}];
         res = If[$VersionNumber >= 9.,
                      MapAt[If[FreeQ[#, _Missing], QuantityMagnitude[UnitConvert[#, $DefaultLengthUnit]], #] &, res, {All, -1}],
                      Function[r, MapAt[If[FreeQ[#, _Missing], # $conversionFactor, #] &, r, -1]] /@ res];
         res = Replace[res, x_ /; ! validateChem[x] :> $Failed, {1}],
         res = Table[If[out === "SMILES", 
                               tmp = Cases[ChemicalData[c, #] & /@ {"IsomericSMILES", "SMILES"}, _String]; 
                               If[Length[tmp] > 0, First[tmp], $Failed], 
                               ChemicalData[c, out] /. Except[_String] -> $Failed], {c, chk}]];
     If[Length[res] == 1, First[res], res]]

(* **ADVANCED USERS**:
     You can remove the code sandwiched by the row of dashes below and
     permanently assign your API token to $ChemSpiderAPIToken instead
     by editing this package file and uncommenting the line below:

     $ChemSpiderAPIToken = ;

 *)

(* ---------------- *)

$dir := DirectoryName[$InputFileName];
$ChemSpiderAPIToken := $ChemSpiderAPIToken = loadToken["RESTful"];
reloadToken["RESTful"] := ($ChemSpiderAPIToken = loadToken["RESTful"]);

loadToken[type_] := Block[{key = FileNameJoin[{$dir, "apikey-restful"}], res}, 

          If[FileExistsQ[key], Uncompress[Get[key]],
              res = DialogInput[DynamicModule[{api}, 
                                         Column[{Style["Please enter your ChemSpider API key:", Bold, 16],
                                                       InputField[Dynamic[api, (api = #) &], String, FieldSize -> 25], 
                                                       Item[Row[{DefaultButton[DialogReturn[api]], CancelButton[DialogReturn[$Failed]]}], 
                                                               Alignment -> Right]}, Left]]];
              If[StringQ[res], Put[Compress[res], key]; res, $Failed]]]

(* ---------------- *)

Options[GetChemSpider] = {"MaxResults" -> 1, Method -> Automatic};

GetChemSpider[arg_, out : (Automatic | "SMILES" | "InChI") : Automatic, OptionsPattern[]] := Block[{met, n},

     If[$VersionNumber < 9., Message[GetChemSpider::unavail, GetChemSpider]; Return[$Failed, Block]];
     met = OptionValue[Method];
     If[met === Automatic, met = "RESTful",
         met = met /. {"Legacy" -> "SOAP", _ -> "RESTful"}];

     n = Round[OptionValue["MaxResults"]];
     Switch[met,
                "RESTful", gCSRESTful[arg, out, "MaxResults" -> n],
                "SOAP", Message[MoleculeViewer::obsmet]; $Failed, (* for future removal *)
                _, $Failed]]

SetAttributes[csconv, Listable];
csconv[str_String, in_String, out_String] := Block[{p, q, r, url},
   url = "https://api.rsc.org/compounds/v1/tools/convert";
   p = {"Method" -> "POST", "Headers" -> {"apikey" -> $ChemSpiderAPIToken, "Content-Type" -> "application/json"}, 
           "Body" -> ExportString[{"input" -> str, "inputFormat" -> in, "outputFormat" -> out}, "JSON"]};

   If[$VersionNumber >= 11.,
       q = HTTPRequest[url, Association[p]];
       If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[q], "JSON"]]];
       r = URLExecute[q, "JSON"],
       If[TrueQ[$debug], Print[Defer[URLFetch][url, p]]];
       r = ImportString[URLFetch[url, p], "JSON"]];

   If[r =!= $Failed && Head[r] =!= Failure, r = "output" /. r;
       If[r =!= "output", r, Message[GetChemSpider::noconv, in, out]; $Failed], Message[GetChemSpider::noconv, in, out]; $Failed]]

gCSRESTful[csid_Integer?Positive, out : (Automatic | "SMILES" | "InChI") : Automatic, opts___] := Module[{pars, req, rsp, type, url},

      If[! ValueQ[$ChemSpiderAPIToken] || ! StringQ[$ChemSpiderAPIToken], 
          Message[GetChemSpider::token]; Return[$Failed, Module]];

      url = "https://api.rsc.org/compounds/v1/records/" <> IntegerString[csid] <> "/details";
      type = If[out === Automatic, "mol3D", "smiles"];
      pars = {"Method" -> "GET", "Headers" -> {"apikey" -> $ChemSpiderAPIToken, "Content-Type" -> "application/json"}};

      If[$VersionNumber >= 11.,
          AppendTo[pars, "Query" -> {"fields" -> type}];
          req = HTTPRequest[url, Association[pars]];
          If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req], "JSON"]]];
          rsp = URLExecute[req, "JSON"],
          AppendTo[pars, "Parameters" -> {"fields" -> type}];
          If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
          rsp = ImportString[URLFetch[url, pars], "JSON"]];

      If[rsp =!= $Failed && Head[rsp] =!= Failure, rsp = type /. rsp;
          If[! StringMatchQ[rsp, type],
              Which[out === Automatic, molFromString[rsp],
                        out === "InChI", csconv[rsp, "SMILES", "InChI"],
                        True, rsp],
              $Failed], $Failed]]

gCSRESTful[csidList_List /; VectorQ[csidList, Composition[Through, IntegerQ && Positive]], out : (Automatic | "SMILES" | "InChI") : Automatic, opts___] := 
      Module[{cnv, pars, req, rsp, type, url},

                  If[! ValueQ[$ChemSpiderAPIToken] || ! StringQ[$ChemSpiderAPIToken], 
                      Message[GetChemSpider::token]; Return[$Failed, Module]];

                  url = "https://api.rsc.org/compounds/v1/records/batch";
                  type = If[out === Automatic, "mol3D", "smiles"];
                  pars = {"Method" -> "POST", "Headers" -> {"apikey" -> $ChemSpiderAPIToken, "Content-Type" -> "application/json"}, 
                               "Body" -> ExportString[{"recordIds" -> csidList, "fields" -> {type}}, "JSON"]};

                  If[$VersionNumber >= 11.,
                      req = HTTPRequest[url, Association[pars]];
                      If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req], "JSON"]]];
                      rsp = URLExecute[req, "JSON"],
                      If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
                      rsp = ImportString[URLFetch[url, pars], "JSON"]];

                  If[rsp =!= $Failed && Head[rsp] =!= Failure, rsp = "records" /. rsp;
                      If[ListQ[rsp], rsp = type /. rsp;
                          If[VectorQ[rsp, ! StringMatchQ[#, type] &],
                              Which[out === Automatic, molFromString /@ rsp,
                                        out === "InChI", csconv[rsp, "SMILES", "InChI"],
                                        True, rsp],
                              $Failed], $Failed], $Failed]]

cscheck[id_String] := With[{url = "https://api.rsc.org/compounds/v1/filter/" <> id <> "/status", 
                                           pars = {"Method" -> "GET", "Headers" -> {"apikey" -> $ChemSpiderAPIToken, "Content-Type" -> "application/json"}}},
                                         If[TrueQ[$debug], Print["Checking ChemSpider status…"]];
                                         If[$VersionNumber >= 11.,
                                             URLExecute[HTTPRequest[url, Association[pars]], "JSON"],
                                             ImportString[URLFetch[url, pars], "JSON"]]]

gCSRESTful[str_String, out : (Automatic | "SMILES" | "InChI") : Automatic, OptionsPattern[{"MaxResults" -> 1}]] := 
      Module[{chk, nres, pars, qid, req, res, stat, type, url},

                  If[! ValueQ[$ChemSpiderAPIToken] || ! StringQ[$ChemSpiderAPIToken], 
                      Message[GetChemSpider::token]; Return[$Failed, Module]];

                  If[StringMatchQ[str, NumberString], 
                      Return[gCSRESTful[Round[FromDigits[str]], out], Module]];

                  type = Which[validSMILESQ[str], "smiles", StringMatchQ[str, inchirx], "inchi", StringMatchQ[str, inkeyrx], "inchikey", True, "name"];

                  url = "https://api.rsc.org/compounds/v1/filter/" <> type;
                  pars = {"Method" -> "POST", "Headers" -> {"apikey" -> $ChemSpiderAPIToken, "Content-Type" -> "application/json"}, 
                              "Body" -> ExportString[{type -> str}, "JSON"]};

                  If[$VersionNumber >= 11.,
                      req = HTTPRequest[url, Association[pars]];
                      If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req], "JSON"]]];
                      qid = URLExecute[req, "JSON"],
                      If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
                      qid = ImportString[URLFetch[url, pars], "JSON"]];

                  If[qid =!= $Failed && Head[qid] =!= Failure, 
                      qid = "queryId" /. qid, Return[$Failed, Module]];
                  If[qid === "queryId", Return[$Failed, Module]];

                  While[Pause[RandomReal[{4.5, 5.5}]]; chk = cscheck[qid];
                           chk =!= $Failed && Head[chk] =!= Failure && ("status" /. chk) === "Processing"];

                  If[chk =!= $Failed && Head[chk] =!= Failure,
                      {stat, nres} = {"status", "count"} /. chk /. "count" -> 1 (* temporary fix for ChemSpider API bug *);
                      If[stat === "Complete" && nres > 0,

                          nres = Min[Max[1, Round[OptionValue["MaxResults"] /. All -> Infinity]], nres];
                          url = "https://api.rsc.org/compounds/v1/filter/" <> qid <> "/results";
                          pars = {"Method" -> "GET", "Headers" -> {"apikey" -> $ChemSpiderAPIToken, "Content-Type" -> "application/json"}};

                          If[$VersionNumber >= 11.,
                              AppendTo[pars, "Query" -> {"start" -> 0, "count" -> nres}];
                              req = HTTPRequest[url, Association[pars]];
                              If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req], "JSON"]]];
                              res = URLExecute[req, "JSON"],
                              url = url <> "?start=0&count=" <> IntegerString[nres];
                              If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
                              res = ImportString[URLFetch[url, pars], "JSON"]];

                          If[res =!= $Failed && Head[res] =!= Failure,
                              If[(res = "results" /. res) === "results", Return[$Failed, Module]];
                              If[TrueQ[$debug], Print[res]];
                              If[res === {}, Return[res, Module]];
                              If[Length[res] == 1, res = First[res]];
                              gCSRESTful[res, out],
                              $Failed], $Failed], $Failed]]

gCSRESTful[l_List, rest___] := gCSRESTful[#, rest] & /@ l

SetAttributes[GetPubChem, Listable];

Options[GetPubChem] = {"MaxResults" -> 1};

GetPubChem[cid_Integer?Positive, out : (Automatic | "SMILES" | "InChI") : Automatic, opts___] := Module[{res, url},
     If[$VersionNumber < 9., Message[GetPubChem::unavail, GetPubChem]; Return[$Failed, Module]];
     If[out === Automatic,
         url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" <> IntegerString[cid] <> "/record/SDF?record_type=3d";
         If[TrueQ[$debug], Print[Defer[Import][url, "MOL"]]]; molFromFile[url],
         url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" <> IntegerString[cid] <> "/property/";
         res = Quiet[If[out === "SMILES",
                               Check[If[TrueQ[$debug], Print[Defer[Import][url <> "IsomericSMILES/CSV", "CSV"]]];
                                         Import[url <> "IsomericSMILES/CSV", "CSV"],
                                         Check[If[TrueQ[$debug], Print[Defer[Import][url <> "CanonicalSMILES/CSV", "CSV"]]];
                                                   Import[url <> "CanonicalSMILES/CSV", "CSV"],
                                                   $Failed]],
                               Check[If[TrueQ[$debug], Print[Defer[Import][url <> out <> "/CSV", "CSV"]]];
                                         Import[url <> out <> "/CSV", "CSV"], $Failed]]];
         If[res =!= $Failed, res[[2, 2]], $Failed]]]

GetPubChem[str_String, out : (Automatic | "SMILES" | "InChI") : Automatic, OptionsPattern[]] := Module[{res, id, pars, req, type, url},
     If[$VersionNumber < 9., Message[GetPubChem::unavail, GetPubChem]; Return[$Failed, Module]];
     If[StringMatchQ[str, NumberString], Return[GetPubChem[Round[FromDigits[res]], out], Module]];

     If[StringMatchQ[str, inchirx],
         url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON";
         pars = {Method -> "POST", "Body" -> {"inchi" -> str}};
         If[$VersionNumber >= 11., req = HTTPRequest[url, Association[pars]];
             If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req]]]];
             id = URLExecute[req],
             pars = ReplacePart[pars, 2 -> ("Body" -> "inchi=" <> str)];
             If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
             id = ImportString[URLFetch[url, pars], "JSON"]],
         type = Which[validSMILESQ[str], "smiles", StringMatchQ[str, inkeyrx], "inchikey", True, "name"];
         url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/" <> type <> "/" <> percentEncode[str] <> "/cids/JSON";
         If[TrueQ[$debug], Print[Defer[Import][url, "JSON"]]];
         id = Quiet[Check[Import[url, "JSON"], $Failed]]];

     If[id =!= $Failed,
         res = "CID" /. ("IdentifierList" /. id); If[TrueQ[$debug], Print[res]];
         res = Take[res, Max[1, Round[OptionValue["MaxResults"] /. All -> Infinity]] /. Infinity -> All];
         If[Length[res] == 1, res = First[res]]; GetPubChem[res, out],
         $Failed]]

(* common subroutine for processing GetPDB and GetProteinData results *)
extractPDBParts[{res_, het_}, type_String] := Switch[type,
                                                                                "All" | "Automatic", 
                                                                                If[! MatchQ[het, _Missing], Join[res, het, 2], res], 
                                                                                "Components" | "Parts", {res, het},
                                                                                "Main" | "Protein" | "Residue", res,
                                                                                "Additional" | "HetAtoms" | "Heterogen", het,
                                                                                _, $Failed]

SetAttributes[GetPDB, Listable];

Options[GetPDB] = {"Format" -> Automatic, "MaxResults" -> 1, "ShowPart" -> "Components"};

GetPDB[str_String, OptionsPattern[]] := Module[{het, pars, prts, query, req, res, type, u1, u2},
     If[$VersionNumber < 9., Message[GetPDB::unavail, GetPDB]; Return[$Failed, Module]];

     prts = System`Utilities`StringName[OptionValue["ShowPart"]];
     type = ToLowerCase[System`Utilities`StringName[OptionValue["Format"] /. {Automatic -> "PDB", "mmCIF" | "MMCIF" -> "CIF"}]];
     If[! MatchQ[type, "pdb" | "cif"], Return[$Failed, Module]];

     u1 = "https://files.rcsb.org/download/";
     If[StringMatchQ[str, RegularExpression["(\\w){4}"]],
         u1 = u1 <> str <> "." <> type;
         If[TrueQ[$debug], Print[Defer[Import][u1, ToUpperCase[type]]]];
         res = processChemicalURL[u1];
         Return[If[type === "pdb", extractPDBParts[res, prts], res], Module]];

     u2 = "http://www.rcsb.org/pdb/rest/search/?sortfield=rank%20Descending";
     query = XMLObject["Document"][{}, XMLElement["orgPdbQuery", {},
                                                                                {XMLElement["queryType", {}, {"org.pdb.query.simple.AdvancedKeywordQuery"}],
                                                                                  XMLElement["keywords", {}, {str}]}], {}];
     pars = {"Method" -> "POST", "Body" -> ExportString[query, "XML"]};

     If[$VersionNumber >= 11.,
         req = HTTPRequest[u2, Association[Append[pars, "ContentType" -> "application/x-www-form-urlencoded"]]];
         If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req]]]];
         res = URLExecute[req, "Text"],
         If[TrueQ[$debug], Print[Defer[URLFetch][u1, pars]]];
         res = URLFetch[u2, pars]];

     If[StringQ[res], res = ImportString[res, "List", "Numeric" -> False], Return[$Failed, Module]];
     If[MatchQ[res, {__String}] && res =!= {"null"},
         If[TrueQ[$debug], Print[res]];
         res = Take[DeleteDuplicates[res], Min[Max[1, Round[OptionValue["MaxResults"] /. All -> Infinity]], Length[res]] /. Infinity -> All];
         res = processChemicalURL[u1 <> # <> "." <> type] & /@ res;
         If[type === "pdb", res = extractPDBParts[#, prts] & /@ res];
         If[Length[res] == 1, res = First[res]],
         res = $Failed];
     res]

SetAttributes[GetProteinData, Listable];

Options[GetProteinData] = {"ShowPart" -> "Components"};

GetProteinData[str_String, OptionsPattern[]] := Module[{res, het},
     Check[res = Table[ProteinData[str, p], {p, {"AtomTypes", "AtomPositions"}}], Return[$Failed, Module]];
     If[! (VectorQ[res[[1]], StringQ] && MatrixQ[res[[2]]]), Return[$Failed, Module]];
     res = If[$VersionNumber >= 9.,
                  MapAt[QuantityMagnitude[UnitConvert[#, $DefaultLengthUnit]] &, res, 2],
                  MapAt[(# $conversionFactor) &, res, 2]];

     het = Table[ProteinData[str, p], {p, {"AdditionalAtomTypes", "AdditionalAtomPositions"}}];
     If[VectorQ[het[[1]], StringQ] && MatrixQ[het[[2]]],
         het = het /. Thread[ToUpperCase[$atoms] -> $atoms];
         het = If[$VersionNumber >= 9.,
                      MapAt[QuantityMagnitude[UnitConvert[#, $DefaultLengthUnit]] &, het, 2],
                      MapAt[(# $conversionFactor) &, het, 2]],
         het = Missing["NotAvailable"]];

     extractPDBParts[{res, het}, System`Utilities`StringName[OptionValue["ShowPart"]]]]

SetAttributes[GetZINC, Listable];

Options[GetZINC] = {"MaxResults" -> 1};

GetZINC[str_String, out : (Automatic | "SMILES" | "InChI") : Automatic, OptionsPattern[]] := Module[{cns, cnt, res},
     cnt = OptionValue["MaxResults"]; If[NumericQ[cnt], cnt = Round[cnt]];
     cns = If[IntegerQ[cnt], IntegerString[cnt], ToLowerCase[ToString[cnt]]];
     If[out === Automatic,
         res = processChemicalURL["http://zinc.docking.org/substances/search.mol2:mol3d?" <>
                                                 "count=" <> cns <> "&q=" <> percentEncode[str], "MOL2"];
         If[ListQ[res],
             cnt = Min[Max[1, cnt /. All -> Infinity], Length[res]];
             res = Take[res, cnt];
             If[Length[res] == 1, res = First[res]],
             res = $Failed],
         res = Import["http://zinc.docking.org/substances/search.json?" <>
                             "count=" <> cns <> "&q=" <> percentEncode[str], "JSON"];
         If[ListQ[res] && res =!= {},
             cnt = Min[Max[1, cnt /. All -> Infinity], Length[res]]; res = Take[res, cnt];
             If[out === "SMILES", res = "smiles" /. res,
                 res = Import["http://zinc.docking.org/apps/mol/convert?from=" <> # <> "&to=inchi&onfail=empty", "Text"] & /@ ("zinc_id" /. res)];
                 If[Length[res] == 1, res = First[res]],
             If[! ListQ[res], res = $Failed]]];
     res]

(* ChEMBL Beaker API functions *)

SetAttributes[ChEMBLSMILESTo3D, Listable];

Options[ChEMBLSMILESTo3D] = {Method -> Automatic};

ChEMBLSMILESTo3D[smiles_?validSMILESQ, OptionsPattern[]] := Block[{met, pars, raw, req, url},
            If[$VersionNumber < 9., Message[ChEMBLSMILESTo3D::unavail, ChEMBLSMILESTo3D]; Return[$Failed, Block]];
            url = "https://www.ebi.ac.uk/chembl/api/utils/";
            met = ToString[OptionValue[Method] /. Automatic -> ""];
            If[! MemberQ[{"ETKDG", "KDG", "MMFF", "UFF", ""}, met],
                Message[ChEMBLSMILESTo3D::bdmtd, met, Automatic]; met = ""];
            url = url <> (met /. "UFF" -> "") <> "smiles23D"; pars = {Method -> "POST", "Body" -> smiles};
            raw = If[$VersionNumber >= 11.,
                         req = HTTPRequest[url, Association[pars]];
                         If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req]]]];
                         URLExecute[req],
                         If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
                         URLFetch[url, pars]];
            If[StringFreeQ[raw, "error", IgnoreCase -> True] && raw =!= "" && raw =!= $Failed,
                Check[molFromString[raw], $Failed], $Failed]]

ChEMBLSMILESTo3D[__] := $Failed

Options[ChEMBLMOLTo3D] = {Method -> Automatic};

ChEMBLMOLTo3D[ms_String, OptionsPattern[]] := Block[{mol = ms, fn, met, pars, raw, req, url},
            If[$VersionNumber < 9., Message[ChEMBLMOLTo3D::unavail, ChEMBLMOLTo3D]; Return[$Failed, Block]];
            fn = If[FileExistsQ[mol], mol, FindFile[mol]]; If[fn =!= $Failed, mol = Import[fn, "Text"]];
            url = "https://www.ebi.ac.uk/chembl/api/utils/";
            met = ToString[OptionValue[Method] /. Automatic -> ""];
            If[! MemberQ[{"ETKDG", "KDG", "MMFF", "UFF", ""}, met],
                Message[ChEMBLMOLTo3D::bdmtd, met, Automatic]; met = ""];
            url = url <> (met /. "UFF" -> "") <> "ctab23D";
            pars = {Method -> "POST", "Headers" -> {"Expect" -> ""}, "Body" -> mol};
            raw = If[$VersionNumber >= 11.,
                         req = HTTPRequest[url, Association[pars]];
                         If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req]]]];
                         URLExecute[req],
                         If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
                         URLFetch[url, pars]];
            If[StringFreeQ[raw, "error", IgnoreCase -> True] && raw =!= "" && raw =!= $Failed,
                Check[molFromString[raw], $Failed], $Failed]]

ChEMBLMOLTo3D[mol_List?validateChem, opts___] := ChEMBLMOLTo3D[MakeMOL[mol], opts]

If[$VersionNumber >= 11.,
    ChEMBLMOLTo3D[f_File, opts___] := ChEMBLMOLTo3D[ExpandFileName[f], opts]
]

ChEMBLMOLTo3D[list : {(_String | _List?validateChem) ..}, opts___] := ChEMBLMOLTo3D[#, opts] & /@ list

ChEMBLMOLTo3D[__] := $Failed

SetAttributes[ChEMBLInChIToMOL, Listable];

ChEMBLInChIToMOL[ist_String /; StringMatchQ[ist, inchirx]] := Block[{pars, raw, req, url},
           If[$VersionNumber < 9., Message[ChEMBLInChIToMOL::unavail, ChEMBLInChIToMOL]; Return[$Failed, Block]];
           url = "https://www.ebi.ac.uk/chembl/api/utils/inchi2ctab";
           pars = {Method -> "POST", "Body" -> ist};
           raw = If[$VersionNumber >= 11.,
                        req = HTTPRequest[url, Association[pars]];
                        If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req]]]];
                        URLExecute[req],
                        If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
                        URLFetch[url, pars]];
           If[StringFreeQ[raw, "error", IgnoreCase -> True] && raw =!= "" && raw =!= $Failed, raw, $Failed]]

ChEMBLInChIToMOL[__] := $Failed

SetAttributes[ChEMBLSMILESToMOL, Listable];

ChEMBLSMILESToMOL[smiles_?validSMILESQ] := Block[{pars, raw, req, url},
           If[$VersionNumber < 9., Message[ChEMBLSMILESToMOL::unavail, ChEMBLSMILESToMOL]; Return[$Failed, Block]];
           url = "https://www.ebi.ac.uk/chembl/api/utils/smiles2ctab";
           pars = {Method -> "POST", "Body" -> smiles};
           raw = If[$VersionNumber >= 11.,
                        req = HTTPRequest[url, Association[pars]];
                        If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req]]]];
                        URLExecute[req],
                        If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
                        URLFetch[url, pars]];
           If[StringFreeQ[raw, "error", IgnoreCase -> True] && raw =!= "" && raw =!= $Failed, raw, $Failed]]

ChEMBLSMILESToMOL[__] := $Failed

SetAttributes[ChEMBLInChITo3D, Listable];

Options[ChEMBLInChITo3D] = {Method -> Automatic};

ChEMBLInChITo3D[ist_String /; StringMatchQ[ist, inchirx], opts___] := With[{mol = ChEMBLInChIToMOL[ist]}, If[mol =!= $Failed, ChEMBLMOLTo3D[mol, opts, Options[ChEMBLMOLTo3D]], $Failed]]

ChEMBLInChITo3D[__] := $Failed

ChEMBLAddHydrogens[ms_String] := Block[{mol = ms, fn, pars, raw, req, url},
            If[$VersionNumber < 9., Message[ChEMBLAddHydrogens::unavail, ChEMBLAddHydrogens]; Return[$Failed, Block]];
            fn = If[FileExistsQ[mol], mol, FindFile[mol]]; If[fn =!= $Failed, mol = Import[fn, "Text"]];
            url = "https://www.ebi.ac.uk/chembl/api/utils/addHs";
            pars = {Method -> "POST", "Headers" -> {"Expect" -> ""}, "Body" -> mol};
            raw = If[$VersionNumber >= 11.,
                         AppendTo[pars, "Query" -> {"addCoords" -> 1}];
                         req = HTTPRequest[url, Association[pars]];
                         If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req]]]];
                         URLExecute[req],
                         url = url <> "?addCoords=1";
                         If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
                         URLFetch[url, pars]];
            If[StringFreeQ[raw, "error", IgnoreCase -> True] && raw =!= "" && raw =!= $Failed, 
                Check[molFromString[raw], $Failed], $Failed]]

ChEMBLAddHydrogens[mol_List?validateChem] := ChEMBLAddHydrogens[MakeMOL[mol]]

If[$VersionNumber >= 11.,
    ChEMBLAddHydrogens[f_File] := ChEMBLAddHydrogens[ExpandFileName[f]]
]

ChEMBLAddHydrogens[list : {(_String | _List?validateChem) ..}] := ChEMBLAddHydrogens /@ list

ChEMBLAddHydrogens[__] := $Failed

ChEMBLSanitize[ms_String] := Block[{mol = ms, fn, pars, raw, req, url},
            If[$VersionNumber < 9., Message[ChEMBLSanitize::unavail, ChEMBLSanitize]; Return[$Failed, Block]];
            fn = If[FileExistsQ[mol], mol, FindFile[mol]]; If[fn =!= $Failed, mol = Import[fn, "Text"]];
            url = "https://www.ebi.ac.uk/chembl/api/utils/sanitize";
            pars = {Method -> "POST", "Headers" -> {"Expect" -> ""}, "Body" -> mol};
            raw = If[$VersionNumber >= 11., 
                         req = HTTPRequest[url, Association[pars]];
                         If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req]]]];
                         URLExecute[req],
                         If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
                         URLFetch[url, pars]];
            If[StringFreeQ[raw, "error", IgnoreCase -> True] && raw =!= "" && raw =!= $Failed, 
                Check[molFromString[raw], $Failed], $Failed]]

ChEMBLSanitize[mol_List?validateChem] := ChEMBLSanitize[MakeMOL[mol]]

If[$VersionNumber >= 11.,
    ChEMBLSanitize[f_File] := ChEMBLSanitize[ExpandFileName[f]]
]

ChEMBLSanitize[list : {(_String | _List?validateChem) ..}] := ChEMBLSanitize /@ list

ChEMBLSanitize[__] := $Failed

ChEMBLMakeLineNotation[ms_String, type : ("SMILES" | "InChI") : "SMILES"] := Block[{mol = ms, fn, pars, qry, raw, req, url},
           If[$VersionNumber < 9., Message[ChEMBLMakeLineNotation::unavail, ChEMBLMakeLineNotation]; Return[$Failed, Block]];
           fn = If[FileExistsQ[mol], mol, FindFile[mol]]; If[fn =!= $Failed, mol = Import[fn, "Text"]];
           url = "https://www.ebi.ac.uk/chembl/api/utils/ctab2" <> ToLowerCase[type];
           pars = {Method -> "POST", "Headers" -> {"Expect" -> ""}, "Body" -> mol};
           qry = {"isomericSmiles" -> 1, "kekuleSmiles" -> 0, "sanitize" -> 1, "removeHs" -> 0}; (* options for SMILES generation *)
           raw = If[$VersionNumber >= 11.,
                        If[type === "SMILES", AppendTo[pars, "Query" -> qry]];
                        req = HTTPRequest[url, Association[pars]];
                        If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req]]]];
                        URLExecute[req],
                        If[type === "SMILES",
                            url = url <> "?" <>
                                    StringJoin[Riffle[StringJoin[#1, "=", IntegerString[#2]] & @@@ qry, "&"]]];
                        If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
                        ImportString[URLFetch[url, pars], "Table"]];
           If[raw =!= "" && raw =!= $Failed,
               StringTrim[If[type === "SMILES", raw[[-1, 1]], raw]], $Failed]]

ChEMBLMakeLineNotation[mol_List?validateChem, type_] := ChEMBLMakeLineNotation[MakeMOL[mol], type]

If[$VersionNumber >= 11.,
    ChEMBLMakeLineNotation[f_File, type_] := ChEMBLMakeLineNotation[ExpandFileName[f], type]
]

ChEMBLMakeLineNotation[list : {(_String | _List?validateChem) ..}, type_] := ChEMBLMakeLineNotation[#, type] & /@ list

ChEMBLMakeLineNotation[__] := $Failed

ChEMBLMakeInChIKey[mol_String] := Block[{pars, raw, req, type, url},
            If[$VersionNumber < 9., Message[ChEMBLMakeInChIKey::unavail, ChEMBLMakeInChIKey]; Return[$Failed, Block]];
            type = Which[validSMILESQ[mol], "smiles", StringMatchQ[mol, inchirx], "inchi", True, "ctab"];
            url = "https://www.ebi.ac.uk/chembl/api/utils/" <> type <> "2inchiKey";
            pars = {Method -> "POST", "Headers" -> {"Expect" -> ""}, "Body" -> mol};
            raw = If[$VersionNumber >= 11., 
                         req = HTTPRequest[url, Association[pars]];
                         If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req]]]];
                         URLExecute[req],
                         If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
                         URLFetch[url, pars]];
            If[StringFreeQ[raw, "error", IgnoreCase -> True] && raw =!= "" && raw =!= $Failed, raw, $Failed]]

ChEMBLMakeInChIKey[mol_List?validateChem] := ChEMBLMakeInChIKey[MakeMOL[mol]]

ChEMBLMakeInChIKey[list : {(_String | _List?validateChem) ..}] := ChEMBLMakeInChIKey /@ list

ChEMBLMakeInChIKey[__] := $Failed

ChEMBLDescriptors[ms_String, type_String: "JSON"] := Block[{mol = ms, fn, pars, raw, req, url},
           If[$VersionNumber < 9., Message[ChEMBLDescriptors::unavail, ChEMBLDescriptors]; Return[$Failed, Block]];
           fn = If[FileExistsQ[mol], mol, FindFile[mol]]; If[fn =!= $Failed, mol = Import[fn, "Text"]];
           url = "https://www.ebi.ac.uk/chembl/api/utils/descriptors";
           pars = {Method -> "POST", "Headers" -> {"Expect" -> ""}, "Body" -> mol};
           raw = If[$VersionNumber >= 11.,
                        req = HTTPRequest[url, Association[pars]];
                        If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req]]]];
                        URLExecute[req],
                        If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
                        URLFetch[url, pars]];
           If[StringFreeQ[raw, "error", IgnoreCase -> True] && raw =!= "" && raw =!= $Failed,
               raw = ImportString[raw, type];
               If[ListQ[raw], raw = First[raw], Return[$Failed, Block]];
               Switch[raw, {__Rule}, SortBy[raw, First], _Association, KeySort[raw], _, $Failed], $Failed]]

ChEMBLDescriptors[mol_List?validateChem, type___] := ChEMBLDescriptors[MakeMOL[mol], type]

If[$VersionNumber >= 11., 
    ChEMBLDescriptors[f_File, type___] := ChEMBLDescriptors[ExpandFileName[f], type]]

ChEMBLDescriptors[list : {(_String | _List?validateChem) ..}, type___] := ChEMBLDescriptors[#, type] & /@ list

ChEMBLDescriptors[__] := $Failed

(* ChEMBL Beaker OCR methods *)

imgToB64[img_Image] := With[{raw = ExportString[If[ImageColorSpace[img] =!= "Grayscale", ColorConvert[img, "Grayscale"], img], {"Base64", "PNG"}]},
                                               StringJoin["data:image/png;base64,", (*StringReplace[*) raw (*, {"+" -> "-", "/" -> "_"}]*)]]

imgToMOL[img_Image, opts___] := Block[{pars, raw, req, url},
     If[$VersionNumber < 9., Return[$Failed, Block]];
     url = "https://www.ebi.ac.uk/chembl/api/utils/image2ctab";
     pars = {Method -> "POST", "Headers" -> {"Expect" -> ""}, "Body" -> imgToB64[img]};
     raw = If[$VersionNumber >= 11.,
                  req = HTTPRequest[url, Association[pars]];
                  If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req]]]];
                  URLExecute[req],
                  If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
                  URLFetch[url, pars]];
     If[StringFreeQ[raw, "error", IgnoreCase -> True] && raw =!= "" && raw =!= $Failed, raw, $Failed]]

imgToMOL[x_, opts___] := imgToMOL[Rasterize[x, "Image", ColorSpace -> "Grayscale", FilterRules[{opts}, Options[Rasterize]]]]

imgToSMILES[img_Image, opts___] := Block[{pars, raw, req, url},
     If[$VersionNumber < 9., Return[$Failed, Block]];
     url = "https://www.ebi.ac.uk/chembl/api/utils/image2smiles";
     pars = {Method -> "POST", "Headers" -> {"Expect" -> ""}, "Body" -> imgToB64[img]};
     raw = If[$VersionNumber >= 11.,
                  req = HTTPRequest[url, Association[pars]];
                  If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req]]]];
                  URLExecute[req, "Text"],
                  If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
                  URLFetch[url, pars]];
     If[StringFreeQ[raw, "error", IgnoreCase -> True] && raw =!= "" && raw =!= $Failed, Last[StringSplit[raw, "\n"]], $Failed]]

imgToSMILES[x_, opts___] := imgToSMILES[Rasterize[x, "Image", ColorSpace -> "Grayscale", FilterRules[{opts}, Options[Rasterize]]]]

(* OPSIN interface *)

OPSINNameToStructure[name_String, type : ("InChI" | "SMILES" | "CML" | "InChIKey") : "SMILES"] := Module[{ext},
         If[$VersionNumber < 9., Return[$Failed, Block]];
         ext = Switch[type, "SMILES", ".smi", "InChI", ".stdinchi", "CML", ".cml", "InChIKey", "stdinchikey", _, Return[$Failed, Module]];
         Import["https://opsin.ch.cam.ac.uk/opsin/" <> percentEncode[name] <> ext, "Text"]]

(* functions for Java-based chemical drawing *)

supportedApplets = {"JME", "JDraw", "JChemPaint"};
installedApplets := Pick[supportedApplets, Table[DirectoryQ[FileNameJoin[{$UserBaseDirectory, "Applications", ap}]], {ap, supportedApplets}]]
$DrawingApplets := Sort[installedApplets]

SetDrawingApplet[app_] := With[{ap = System`Utilities`StringName[app]},
     Check[If[MemberQ[supportedApplets, ap], SetOptions[{MoleculeDraw, GetMolecule, SetMolecule}, "Applet" -> ap];, $Failed], $Failed]]

Options[MoleculeDraw] = {"Applet" -> Automatic, "PreProcess" -> Automatic};

MoleculeDraw[applet_ /; MemberQ[supportedApplets, System`Utilities`StringName[applet]]][mol_, opts___] := MoleculeDraw[mol, "Applet" -> applet, opts]

MoleculeDraw[arg_String /; ! MemberQ[supportedApplets, arg], OptionsPattern[]] := Module[{mol = arg, app, procF},
            app = System`Utilities`StringName[OptionValue["Applet"] /. Automatic -> "JME"];
            If[! MemberQ[supportedApplets, app],
               Message[MoleculeDraw::noop, app, MoleculeDraw]; Return[$Failed, Module]];
            procF = OptionValue["PreProcess"] /. None -> Identity;
            If[MatchQ[procF, Automatic | "Automatic"],
                Which[validSMILESQ[mol],
                          procF = Which[testOpenBabel, obMakeString,
                                                PacletManager`$AllowInternet, ChEMBLSMILESToMOL,
                                                True, Message[MoleculeDraw::nopre]; Return[$Failed, Module]],
                          StringMatchQ[mol, inchirx],
                          procF = Which[testOpenBabel, obMakeString,
                                                PacletManager`$AllowInternet, ChEMBLInChIToMOL,
                                                True, Message[MoleculeDraw::nopre]; Return[$Failed, Module]],
                          True, procF = Identity]];
            openApplet[app, procF[mol]]]

MoleculeDraw[mol_List?validateChem, opts___] := MoleculeDraw[MakeMOL[mol], opts]

MoleculeDraw[(Rule | RuleDelayed)["Applet", applet_]] /; MemberQ[supportedApplets, System`Utilities`StringName[applet]] := openApplet[System`Utilities`StringName[applet]]

MoleculeDraw[] := openApplet[System`Utilities`StringName["Applet" /. Options[MoleculeDraw] /. Automatic -> "JME"]]

Options[GetMolecule] = {"Applet" -> Automatic, "Format" -> Automatic, "PostProcess" -> Automatic};

GetMolecule[OptionsPattern[]] := Module[{app, procF, type},
     app = System`Utilities`StringName[OptionValue["Applet"] /. Automatic -> "JME"];
     If[! MemberQ[supportedApplets, app],
         Message[GetMolecule::noop, app, GetMolecule]; Return[$Failed, Module]];
     procF = OptionValue["PostProcess"] /. None -> Identity;
     type = OptionValue["Format"] /. Automatic -> "MOL";
     If[MatchQ[procF, Automatic | "Automatic"],
         procF = Switch[type,
                                "MOL",
                                Which[# === $Failed, Return[$Failed, Module],
                                          testOpenBabel, obMOLTo3D[#],
                                          PacletManager`$AllowInternet, ChEMBLMOLTo3D[#, Method -> "ETKDG"],
                                          True, Message[GetMolecule::nopost]; #] &,
                                "SMILES",
                                Which[# === $Failed, Return[$Failed, Module],
                                          testOpenBabel, RunOpenBabel[#, Method -> "Generate"],
                                          PacletManager`$AllowInternet, ChEMBLSMILESTo3D[#, Method -> "ETKDG"],
                                          True, Message[GetMolecule::nopost]; #] &,
                                True, Identity]];
     Switch[type,
                "MOL", procF[getMOL[app]],
                "SMILES", procF[getSMILES[app]],
                _, Message[GetMolecule::nsupo, type, app]; $Failed]]

GetMolecule[applet_ /; MemberQ[supportedApplets, System`Utilities`StringName[applet]], opts___] := GetMolecule["Applet" -> applet, opts]

Options[SetMolecule] = {"Applet" -> Automatic};

SetMolecule[applet_ /; MemberQ[supportedApplets, System`Utilities`StringName[applet]]][mol_, opts___] := SetMolecule[mol, "Applet" -> applet, opts]

SetMolecule[mol_String /; ! MemberQ[supportedApplets, mol], OptionsPattern[]] := Module[{app, tmp},
     app = System`Utilities`StringName[OptionValue["Applet"] /. Automatic -> "JME"];
     If[! MemberQ[supportedApplets, app],
         Message[SetMolecule::noop, app, SetMolecule]; Return[False, Module]];
     Which[validSMILESQ[mol],
               setSMILES[app, mol] =!= $Failed ||
               With[{r = Which[testOpenBabel, obMakeString[mol, "MOL"],
                                         PacletManager`$AllowInternet, ChEMBLSMILESToMOL[mol],
                                         True, $Failed]},
                       r =!= $Failed && setMOL[app, r] =!= $Failed],
               StringMatchQ[mol, inchirx],
               Message[SetMolecule::nsupi, "InChI", app];
               tmp = Which[testOpenBabel, obMakeString[mol, "MOL"],
                                   PacletManager`$AllowInternet, ChEMBLInChIToMOL[mol],
                                   True, $Failed];
               tmp =!= $Failed && setMOL[app, tmp] =!= $Failed,
               True, setMOL[app, mol] =!= $Failed]]

SetMolecule[mol_List?validateChem, opts___] := SetMolecule[MakeMOL[mol], opts]

(* JME *)

installApplet["JME"] := Block[{dir = FileNameJoin[{$UserBaseDirectory, "Applications", "JME", "Java"}], applet, url},
        applet = FileNameJoin[{dir, "JME.jar"}];
        url = "http://www.molinspiration.com/jme/doc/JME.jar";
        If[FileExistsQ[applet], True,
            CreateDirectory[dir];
            Check[If[$VersionNumber >= 9., URLSave, Utilities`URLTools`FetchURL][url, applet]; True,
                      Message[MoleculeDraw::noapp, "JME"]; False]]]

openApplet["JME", args___] := Block[{dat = {args}},
                                                        If[dat =!= {} && StringQ[First[dat]], dat = {"mol=" <> StringReplace[First[dat], "\n" -> "|"]}, dat = {}];
                                                        JLink`JavaBlock[If[installApplet["JME"], JLink`InstallJava[];
                                                                                    Unprotect[$JME]; JLink`ReleaseJavaObject[$JME];
                                                                                    $JME = JLink`JavaNew["JME"]; Protect[$JME];
                                                                                    JLink`AppletViewer[$JME, {"width=360", "height=335",
                                                                                                                             "options=autoez,hidehs,noquery,norbutton,polarnitro"} ~Join~ dat];
                                                                                    $JME,
                                                                                    Message[MoleculeDraw::noapp, "JME"]; $Failed]]]

getMOL["JME"] := JLink`JavaBlock[If[ValueQ[$JME], $JME @ molFile[], Message[GetMolecule::noapp, "JME"]; $Failed]]

getSMILES["JME"] := JLink`JavaBlock[If[ValueQ[$JME], $JME @ smiles[], Message[GetMolecule::noapp, "JME"]; $Failed]]

setMOL["JME", mol_String] := JLink`JavaBlock[If[ValueQ[$JME], $JME @ readMolFile[StringReplace[mol, "\n" -> " |"]], Message[SetMolecule::noapp, "JME"]; $Failed]]

setSMILES["JME", mol_String] := (Message[SetMolecule::nsupi, "SMILES", "JME"]; $Failed)

(* for backward compatibility with earlier MoleculeViewer versions *)

InstallJME[] := installApplet["JME"]

LaunchJME[] := (Message[LaunchJME::obsfun, LaunchJME, MoleculeDraw]; openApplet["JME"])

GetJMEMOL[] := (Message[GetJMEMOL::obsfun, GetJMEMOL, GetMolecule]; getMOL["JME"])

GetJMESMILES[] := (Message[GetJMESMILES::obsfun, GetJMESMILES, GetMolecule]; getSMILES["JME"])

SetJMEMOL[mol_String] := (Message[SetJMEMOL::obsfun, SetJMEMOL, SetMolecule]; setMOL["JME", mol])

JME[] := If[ValueQ[$JME], getSMILES["JME"], openApplet["JME"]]

(* JChemPaint *)

installApplet["JChemPaint"] := Block[{applets, dir, tmp, url},
        dir = FileNameJoin[{$UserBaseDirectory, "Applications", "JChemPaint", "Java"}];
        applets = Uncompress["1:eJyVk8FOwzAMhnvgIWACwYQ4INExaTwF4sQeIDKZt6ZL4pB4XcXTU9IBlzRrL1GS//vtyLHnH/S+nRV\
                                           FES665U0FXt92m1pWaBwoyyU4p5FLSUbQF6Nd1ODXswHGY5SvkzJuFJOPwDwDlOQ4ZLK8goEoXyXlG\
                                           nrz/YDaQCsalAa4itxlkrOHKN4kRfI7sXKfkbgbJMBBd595aQdFNV3wnxA1tKeCLwcZcmiDVGglCrnZC0l2q/\
                                           rAL1NM7Elr9BOzKRLS6Gh6mmCKhsexhkgvRof/+9rn0RYm0n3bPIzxnOmN40qe6Y3W6N+ypWeFuEIfMu\
                                           3hMdDBS8wx3J018IlJz1yj8Ij/M/cNB6Nujg=="];
        url = "https://github.com/downloads/JChemPaint/jchempaint/jchempaint-applet-3.3-1210.zip";
        If[Apply[And, FileExistsQ[FileNameJoin[{dir, #}]] & /@ applets], True,
            Check[Internal`WithLocalSettings[CreateDirectory[dir]; tmp = FileNameJoin[{dir, "jchempaint.zip"}];
                                                              If[$VersionNumber >= 9., URLSave, Utilities`URLTools`FetchURL][url, tmp],
                                                              ExtractArchive[tmp, dir, "*.jar", CreateIntermediateDirectories -> False]; True,
                                                              DeleteFile[tmp]],
                      Message[MoleculeDraw::noapp, "JChemPaint"]; False]]]

openApplet["JChemPaint", args___] := Block[{dat = {args}, mol},
                                                                   JLink`JavaBlock[If[installApplet["JChemPaint"], JLink`InstallJava[];
                                                                                               Unprotect[$JChemPaint]; JLink`ReleaseJavaObject[$JChemPaint];
                                                                                               $JChemPaint = JLink`JavaNew["org.openscience.jchempaint.applet.JChemPaintEditorApplet"]; Protect[$JChemPaint];
                                                                                               JLink`AppletViewer[$JChemPaint, {"impliciths=true", "codebase_lookup=false", "scrollbars=true", "width=585", "height=455"}];
                                                                                               If[dat =!= {} && StringQ[First[dat]], mol = First[dat];
                                                                                                   If[validSMILESQ[mol], $JChemPaint @ setSmiles[mol], $JChemPaint @ setMolFile[mol]]]; $JChemPaint,
                                                                                               Message[MoleculeDraw::noapp, "JChemPaint"]; $Failed]]]

getMOL["JChemPaint"] := JLink`JavaBlock[If[ValueQ[$JChemPaint], $JChemPaint @ getMolFile[], Message[GetMolecule::noapp, "JChemPaint"]; $Failed]]

getSMILES["JChemPaint"] := JLink`JavaBlock[If[ValueQ[$JChemPaint], $JChemPaint @ getSmiles[], Message[GetMolecule::noapp, "JChemPaint"]; $Failed]]

setMOL["JChemPaint", mol_String] := JLink`JavaBlock[If[ValueQ[$JChemPaint], $JChemPaint @ setMolFile[mol], Message[SetMolecule::noapp, "JChemPaint"]; $Failed]]

setSMILES["JChemPaint", mol_?validSMILESQ] := JLink`JavaBlock[If[ValueQ[$JChemPaint], $JChemPaint @ setSmiles[mol], Message[SetMolecule::noapp, "JChemPaint"]; $Failed]]

JChemPaint[] := If[ValueQ[$JChemPaint], getSMILES["JChemPaint"], LaunchJChemPaint[]]

(* JDraw *)

installApplet["JDraw"] := Block[{dir = FileNameJoin[{$UserBaseDirectory, "Applications", "JDraw"}], applets, tmp, url},
        applets = Table[FileNameJoin[{dir, "Java", files}], {files, {"CsInline.jar", "jdrawapplet.jar", "jdrawcore.jar"}}];
        url = "http://download.accelrys.com/freeware/accelrys_jdraw/AccelrysJDraw-1.1.400.25.zip";
        If[Apply[And, FileExistsQ /@ applets], True,
            Check[Internal`WithLocalSettings[CreateDirectory[dir]; tmp = FileNameJoin[{dir, "jdraw.zip"}];
                                                              If[$VersionNumber >= 9., URLSave, Utilities`URLTools`FetchURL][url, tmp],
                                                              ExtractArchive[tmp, dir, "*.jar", CreateIntermediateDirectories -> False];
                                                              RenameDirectory[FileNameJoin[{dir, "AccelrysJDraw-1.1.400.25"}], FileNameJoin[{dir, "Java"}]]; True,
                                                              DeleteFile[tmp]],
                      Message[MoleculeDraw::noapp, "JDraw"]; False]]]

openApplet["JDraw", args___] := Block[{dat = {args}},
                                                           If[dat =!= {} && StringQ[First[dat]], dat = {"molString=" <> First[dat]}, dat = {}];
                                                           JLink`JavaBlock[If[installApplet["JDraw"], JLink`InstallJava[(* JLink`JVMArguments -> "-Xmx512m -Dsun.java2d.noddraw=true" *)];
                                                                                       Unprotect[$JDraw]; JLink`ReleaseJavaObject[$JDraw];
                                                                                       $JDraw = JLink`JavaNew["com.symyx.draw.JDrawEditor"]; Protect[$JDraw];
                                                                                       JLink`AppletViewer[$JDraw, {"hydrogenDisplayMode=1", "width=800", "height=500"} ~Join~ dat]; $JDraw,
                                                                                       Message[MoleculeDraw::noapp, "JDraw"]; $Failed]]]

getMOL["JDraw"] := JLink`JavaBlock[If[ValueQ[$JDraw], $JDraw @ getMolString[], Message[GetMolecule::noapp, "JDraw"]; $Failed]]

getSMILES["JDraw"] := (Message[GetMolecule::nsupo, "SMILES", "JDraw"]; $Failed)

setMOL["JDraw", mol_String] := JLink`JavaBlock[If[ValueQ[$JDraw], $JDraw @ setMolString[mol], Message[SetMolecule::noapp, "JDraw"]; $Failed]]

setSMILES["JDraw", mol_String] := (Message[SetMolecule::nsupi, "SMILES", "JDraw"]; $Failed)

JDraw[] := If[ValueQ[$JDraw], getMOL["JDraw"], openApplet["JDraw"]]

(* MoleculeRecognize *)

Options[MoleculeRecognize] = {Method -> Automatic, "PostProcess" -> Automatic};

MoleculeRecognize[img_Image, OptionsPattern[{MoleculeRecognize, Rasterize}]] := Module[{met, mopt, ocrF, procF, res},
            met = OptionValue[Method];
            If[ListQ[met],
               {met, mopt} = {First[met], Flatten[Rest[met]]},
               mopt = {}];
            ocrF = Switch[System`Utilities`StringName[met],
                                  "Automatic", Which[testImago, RunImago,
                                                               PacletManager`$AllowInternet, imgToMOL,
                                                               True, Message[MoleculeRecognize::noocr]; Return[$Failed, Module]],
                                  "Imago" | "RunImago", RunImago,
                                  "Beaker" | "ChEMBL", imgToMOL,
                                  _, Message[MoleculeRecognize::method, met]; Return[$Failed, Module]];
            res = ocrF[img, mopt];

            Switch[procF = OptionValue["PostProcess"],
                       Automatic | "Automatic",
                       procF = Which[testOpenBabel, obMOLTo3D,
                                             PacletManager`$AllowInternet, ChEMBLMOLTo3D[#, Method -> "ETKDG"] &,
                                             True, Message[MoleculeRecognize::badimg]; Return[$Failed, Module]],
                       None | "None", procF = Identity,
                       "OpenBabel", procF = obMOLTo3D,
                       "Beaker" | "ChEMBL", procF = ChEMBLMOLTo3D];
            procF[res]]

MoleculeRecognize[x_, opts___] := MoleculeRecognize[Rasterize[x, "Image", ColorSpace -> "Grayscale", FilterRules[{opts}, Options[Rasterize]]], opts]

(* MoleculeViewer *)

Options[MoleculeViewer] = SortBy[Join[{Background -> RGBColor[0., 43./255, 54./255], BaseStyle -> Automatic, Boxed -> False, ColorRules -> Automatic, Highlighted -> None, Lighting -> "Neutral",
                                                             PlotLegends -> None, PlotStyle -> Automatic, PreserveImageOptions -> True, SphericalRegion -> True, Tooltip -> None, ViewPoint -> Automatic},
                                                           Select[Options[Graphics3D], FreeQ[#, Background | BaseStyle | Boxed | ImageSizeRaw | Lighting | PreserveImageOptions | SphericalRegion | ViewPoint] &]], First];

SyntaxInformation[MoleculeViewer] = {"ArgumentsPattern" -> {_, OptionsPattern[MoleculeViewer]}};

MoleculeViewer[$Failed | {} | _Missing | {($Failed | {} | _Missing) ..}, opts___] := $Failed

MoleculeViewer[args___] := Block[{a, res},
            a = System`Private`Arguments[MoleculeViewer[args], 1];
            res /; TrueQ[a =!= {} &&
                               (res = iMoleculeViewer @@ Flatten[a, 1]) =!= $Failed &&
                               MatchQ[res, _Graphics3D | _Legended | _Row | {(_Graphics3D | _Legended | _Row) ..}]]]

iMoleculeViewer[dat : {__String}, opts___] := iMoleculeViewer[#, opts] & /@ dat
iMoleculeViewer[dat : {{_List, _?MatrixQ} ..}, opts___] := iMoleculeViewer[#, opts] & /@ dat
iMoleculeViewer[dat : {({_List, _List, _List, _?MatrixQ} | {_List, None, None, _?MatrixQ}) ..}, opts___] := iMoleculeViewer[#, opts] & /@ dat
iMoleculeViewer[dat : {{__Rule} ..}, opts___] := iMoleculeViewer[#, opts] & /@ dat

iMoleculeViewer[str_String, opts___] := Block[{res = $Failed},
 Quiet[Catch[Scan[Function[op, With[{r = op[str]}, If[r =!= $Failed, If[TrueQ[$debug], Print["Function used: ", op]]; Throw[Set[res, r /. _Missing -> Nothing], "$MoleculeViewer"]]]],
                            Which[FileExistsQ[str] || FindFile[str] =!= $Failed, {processChemicalFile},
                                      StringMatchQ[str, "http://*" | "ftp://*" | "https://*"], {processChemicalURL},
                                      StringMatchQ[str, inchirx] || validSMILESQ[str], {RunOpenBabel, If[validSMILESQ[str], ChEMBLSMILESTo3D, ChEMBLInChITo3D], GetCACTUS},
                                      True, {GetChemicalData, GetPubChem, GetCACTUS, GetChemSpider}]], "$MoleculeViewer"],
          {ChemicalData::notent, Import::fmterr, Utilities`URLTools`FetchURL::conopen, Utilities`URLTools`FetchURL::httperr, MoleculeViewer::badfile}];
 If[res =!= $Failed, iMoleculeViewer[res, opts], Message[MoleculeViewer::badstrng]; res]]

iMoleculeViewer[rl : {__Rule}, opts___] := iMoleculeViewer[DeleteCases[{"VertexTypes", "EdgeRules", "EdgeTypes", "VertexCoordinates"} /. rl, _String, {1}], opts];

If[$VersionNumber >= 10.,
    iMoleculeViewer[assoc_Association, opts___] :=
     With[{tmp = DeleteMissing[Lookup[assoc, {"VertexTypes", "EdgeRules", "EdgeTypes", "VertexCoordinates"}]]},
             If[MatchQ[Length[tmp], 2 | 4], iMoleculeViewer[tmp, opts], iMoleculeViewer[#, opts] & /@ assoc]];
    iMoleculeViewer[ent_Entity, opts___] := iMoleculeViewer[processEntity[ent], opts];
    iMoleculeViewer[dat : {__Association} | {__Entity}, opts___] := iMoleculeViewer[#, opts] & /@ dat
]

If[$VersionNumber >= 11.,
    iMoleculeViewer[f_File, opts___] := iMoleculeViewer[processChemicalFile[ExpandFileName[f]], opts];
    iMoleculeViewer[u_URL, opts___] := iMoleculeViewer[processChemicalURL[URLBuild[URLParse[u]]], opts];
    iMoleculeViewer[dat : {__File} | {__URL}, opts___] := iMoleculeViewer[#, opts] & /@ dat
]

If[$VersionNumber >= 12.,
    iMoleculeViewer[mol_Molecule, opts___] := iMoleculeViewer[processMoleculeObject[mol], opts];
    iMoleculeViewer[dat : {__Molecule}, opts___] := iMoleculeViewer[#, opts] & /@ dat
]

iMoleculeViewer[{vertexTypes_, atomPositions_}, opts___] := With[{bonds = inferBonds[vertexTypes, atomPositions]},
 iMoleculeViewer[{vertexTypes, bonds, ConstantArray["Single", Length[bonds]], atomPositions}, opts]]

iMoleculeViewer[{vertexTypes_, None | _Missing, None | _Missing, atomPositions_}, opts___] := iMoleculeViewer[{vertexTypes, {}, {}, atomPositions}, opts]

iMoleculeViewer[{vertexTypes_, edgeRules_, edgeTypes_, atomPositions_}, opts : OptionsPattern[{MoleculeViewer, Graphics3D}]] :=
 Block[{amlg, arQ, atomplot, atomQ, atoms, bondch, bondlist, bondplot, bondQ, bondtypes, bnwQ, cc, cent, cr, cRules, e, hiCol, hiSet, imat, ipos, lbl, leg, lchk,
            makebonds, mbp, mbstl, myLegendFunction, myOrient, mySwatch, nbn, nei, nmb, nrms, pos, pStyle, psOp, res, skel, shft, tips, tubeQ, v, view, vln, zpos},

          If[! validateChem[{vertexTypes, edgeRules, edgeTypes, atomPositions}], Message[MoleculeViewer::badchem]; Return[$Failed, Block]];
          Check[$conversionFactor, Return[$Failed, Block]];

          v = Length[vertexTypes]; e = Length[edgeTypes];

          bondQ = (edgeRules =!= {} && edgeTypes =!= {});

          pos = atomPositions;
          If[MatrixQ[pos, QuantityQ], pos = QuantityMagnitude[UnitConvert[pos, $DefaultLengthUnit]]];
          pos = SetPrecision[pos, MachinePrecision];

          pStyle = OptionValue[PlotStyle]; psOp = {};
          If[ListQ[pStyle],
              If[! MatchQ[pStyle, {(_Rule | _RuleDelayed) ..}],
                  {pStyle, psOp} = {First[pStyle], Flatten[Rest[pStyle]]};
                  pStyle = System`Utilities`StringName[pStyle],
                  psOp = pStyle],
              pStyle = System`Utilities`StringName[pStyle]];
          If[StringQ[pStyle], (* set up styles *)
              Switch[pStyle,
                         "Automatic",
                         atomQ = tubeQ = (Count[vertexTypes, Except["H"], {1}] < 500);,
                         "BallAndStick",
                         atomQ = tubeQ = True;,
                         "Tubes" | "Tube" | "Sticks" | "Stick",
                         If[bondQ, atomQ = False; tubeQ = True, Message[MoleculeViewer::nbnd, pStyle]; atomQ = True];,
                         "Wireframe" | "WireFrame",
                         If[bondQ, atomQ = tubeQ = False, Message[MoleculeViewer::nbnd, pStyle]; atomQ = True];,
                         "Spacefilling" | "SpaceFilling" | "Calotte" | "CPK",
                         atomQ = True; bondQ = False;,
                         _, Message[MoleculeViewer::nstl, pStyle];
                         atomQ = tubeQ = True],
              {atomQ, tubeQ} = Replace[{"Atoms", "RenderTubes"} /. psOp, {None -> False, Except[True | False] -> True}, {1}];
              If[bondQ,
                  bondQ = "Bonds" /. psOp /. {None -> False, Except[True | False] -> True},
                  Message[MoleculeViewer::nbnd, pStyle]; atomQ = True]];

          mbstl = 2; shft = False; arQ = False; bnwQ = False;
          If[psOp =!= {}, (* additional styles *)
              {bnwQ, mbstl, shft, arQ} = {"Monochrome", "MultipleBondStyle", "ShiftIsolatedAtoms", "ShowAromaticRings"} /. psOp;
              mbstl = System`Utilities`StringName[mbstl] /. {"Classic" | "Parabolic" -> 1, "Taut" | "Automatic" -> 2} /. Except[1 | 2] -> 2;
              shft = TrueQ[shft /. Automatic -> False]; arQ = TrueQ[arQ /. Automatic -> False]; bnwQ = TrueQ[bnwQ /. Automatic -> False]];

          myOrient = AffineTransform[RotationMatrix[-Pi/12, {1, 0, 0}].RotationMatrix[(1 - Sqrt[5]/2) Pi, {0, 0, 1}]];
          If[(view = OptionValue[ViewPoint]) === Automatic, (* reorient based on heavy atoms *)
              view = getNormal[If[Length[vertexTypes] - Count[vertexTypes, "H"] > 2, Pick[pos, Thread[! StringMatchQ[vertexTypes, "H"]]], pos]];
              view = 3.4 myOrient[PadRight[view, 3]]];

          If[bondQ && shft, pos = Last[ShiftIsolatedAtoms[{vertexTypes, edgeRules, edgeTypes, pos}]]];

          If[Last[Dimensions[pos]] == 2, (* flat/linear molecule *)
              Message[MoleculeViewer::is2d]; view = {0, -Infinity, 0}; pos = Insert[#, 0., 2] & /@ pos];

          cRules = Replace[Flatten[{OptionValue[ColorRules] /. Automatic -> {}}], (* generate color list *)
                                    HoldPattern[(s_String /;  MemberQ[$atoms, s]) -> Opacity[o_]] :> (s -> Opacity[o, s /. $colorRules]), {1}];
          cRules = Join[cRules, Append[$colorRules, _ -> RGBColor[211./255, 18./85, 26./51]] /.
                                          If[bnwQ, HoldPattern[a_ -> c_] :> (a -> ColorConvert[c, GrayLevel]), {}]];
          cr = Replace[vertexTypes, cRules, {1}];

          tips = OptionValue[Tooltip]; hiSet = OptionValue[Highlighted];

          If[hiSet =!= None, (* set up highlights *)
              If[! ListQ[hiSet], hiSet = {hiSet}]; hiCol = RGBColor[0.5, 1., 0., 0.6];
              If[VectorQ[hiSet, MatchQ[#, _Integer | _String] &], hiSet = Thread[hiSet -> hiCol]];
              hiSet = Function[l, MapAt[If[IntegerQ[#], {#}, #] &, l, {1}]] /@
                          Replace[hiSet, {Style[x_, c_] :> (x -> c), x : Except[_Rule | _RuleDelayed | _Style] :> (x -> hiCol)}, {1}];
              hiSet = DeleteCases[hiSet /. s_String :> Flatten[Position[vertexTypes, s]], HoldPattern[{} -> _] | HoldPattern[{} :> _]]];

          If[atomQ,
              atoms = MapThread[Sphere, {pos, (vertexTypes /. If[bondQ, $sizeRules, $vdWRules]) $conversionFactor}];
              atomplot = Transpose[{cr, atoms}];

              If[hiSet =!= None, (* highlights *)
                  atomplot = Fold[MapAt[Function[{l}, Join[l, {makeTranslucent[#2[[2]]], MapAt[If[bondQ, 1.2, 1.05] # &, l[[-1]], -1]}]],
                                                     #1, List /@ #2[[1]]] &, atomplot, hiSet]];

                  If[MatchQ[tips, All | "Atoms"], (* tooltips *)
                      atomplot = MapThread[Tooltip[#, Grid[{{Style[#2, Bold, Larger], ""}, {Style["Atom no.", Bold], #3},
                                                                                 {Style["Coordinates:", Bold], #4}}, Alignment -> {{Left, Left}}]] &,
                                                        {atomplot, vertexTypes, Range[v], atomPositions}]],
              (* no atoms *) atomplot = {}];

          cent = Mean[pos];

          If[bondQ, bondlist = List @@@ edgeRules;
              bondtypes = Replace[edgeTypes, {"Single" -> 1, "Double" -> 2, "Triple" -> 3, _ -> 1}, {1}];

              skel = Graph[Range[v], edgeRules, DirectedEdges -> False, GraphLayout -> None];
              vln = VertexDegree[skel];
              imat = Quiet[Check[IncidenceMatrix[skel], $Failed, IncidenceMatrix::ninc], IncidenceMatrix::ninc];
              If[! MatrixQ[imat], Message[MoleculeViewer::badchem]; Return[$Failed, Block]];

              ipos = GatherBy[imat["NonzeroPositions"], First];
              zpos = Position[Total[imat, {2}], 0]; (* find isolated atoms *)
              If[zpos =!= {}, zpos = List /@ PadRight[zpos, {Automatic, 2}]];
              ipos = #[[All, -1]] & /@ SortBy[Join[ipos, zpos], #[[1, 1]] &];

              (* set up normal directions for multiple bonds *)
              nrms = ConstantArray[Null, e]; mbp = Position[bondtypes, Except[1], 1, Heads -> False];
              If[mbp =!= {},
                  nei = Union[Flatten[bondlist[[Union @@ ipos[[#]]]]]] & /@ Extract[bondlist, mbp];
                  (* section of adjacency matrix of line graph *)
                  amlg = (Transpose[imat].imat - 2 IdentityMatrix[e, SparseArray])[[Flatten[mbp]]];
                  (* neighboring bonds *)
                  nbn = MapThread[Complement[#[[All, 2]], #2] &, {GatherBy[amlg["NonzeroPositions"], First], mbp}];

                  nmb = MapThread[With[{tst = getNormal[pos[[#]], True]},
                                                      If[#2 == 2 && Chop[tst[[1]], 0.1 $conversionFactor] == 0 &&
                                                          Complement[bondtypes[[#3]], {2}] === {}, {0., 0., 0.}, tst[[2]]]] &,
                                              {nei, Extract[bondtypes, mbp], nbn}];
                  If[Chop[Norm[nmb, Infinity]] == 0, nmb[[1]] = getNormal[pos[[nei[[1]]]]]];

                  (* reorient normals *)
                  nmb = NestWhile[Function[nset, MapThread[If[Norm[#, Infinity] != 0, #,
                                                                                        Normalize[Cross[Mean[#2 /. Thread[Flatten[mbp] -> nset]],
                                                                                                                 Subtract @@ pos[[Extract[bondlist, #3]]]]]] &,
                                                                                     {nset, nbn, mbp}]], nmb,
                                             ! FreeQ[Chop[Map[Norm[#, Infinity] &, #]], 0] &];

                  If[TrueQ[$adjustNormals], (* further adjustment of normals for adjacent double bonds *)
                      Block[{as, at, bn, bnt, bs, bt, ks, ms, mt, nid, ns, nt},
                               Do[If[(bt = Extract[bondtypes, mbp[[nn]]]) == 2,
                                        bnt = bondtypes[[nid = nbn[[nn]]]];
                                        If[Length[bnt] > 1 && Count[bnt, 2] == 1,
                                            bn = Extract[bondlist, mbp[[nn]]];
                                            bs = First[Pick[bondlist[[nid]], bnt, 2]];
                                            mt = Mean[pos[[bn]]]; ms = Mean[pos[[bs]]];
                                            bn = bn[[Ordering[vertexTypes[[bn]], 2, hillOrder[##] == 1 &]]];
                                            nt = nmb[[nn]] Norm[at = Subtract @@ pos[[bn]]]/2;
                                            ks = Position[mbp, Pick[nid, bnt, 2]][[1, 1]];
                                            bs = bs[[Ordering[vertexTypes[[bs]], 2, hillOrder[##] == 1 &]]];
                                            ns = nmb[[ks]] Norm[as = Subtract @@ pos[[bs]]]/2;
                                            If[Norm[Cross[nt, ns]] < 0.01 ||
                                                Abs[Det[PadRight[{mt, mt + nt, ms, ms + ns}, {4, 4}, 1]]] < 0.06 $conversionFactor^3,
                                                nmb[[nn]] = Normalize[RotationTransform[Pi/4, at] @ nmb[[nn]]];
                                                nmb[[ks]] = Normalize[RotationTransform[Pi/4, as] @ nmb[[ks]]]];]],
                                    {nn, Length[mbp]}]]];

                  nrms = ReplacePart[nrms, Thread[mbp -> nmb]]];

              (* bond constructor function *)
              makebonds[bid_, no_, mul_] := Block[{bo = no, en = pos[[bid]], d, del, fac, mids, mh, ml, pr, tm, tubeF},
                       tubeF = If[tubeQ, Tube[#, (12 - 2 mul) $conversionFactor] &, Identity];
                       del = Subtract @@ en; d = Norm[del]; del /= d;
                       mh = 1 + Boole[mul > 1] (1 + 4 Boole[mbstl == 2]); ml = 2 mh;
                       mids = Range[0, 1, 1/ml]; tm = Transpose[{1 - mids, mids}].en;

                       If[mul == 1, (* single bond *)
                           Transpose[{cr[[bid]], tubeF[Line[#]] & /@ Partition[tm, 2, 1]}],
                           (* multiple bond *)
                           If[Chop[Abs[pr = bo.del]] != 0, (* enforce orthogonality to bond *)
                               bo -= del pr; bo = Normalize[bo];
                               bo *= 2 UnitStep[Apply[Subtract, EuclideanDistance[cent, # bo] & /@ {1, -1}]] - 1];
                           fac = d If[mbstl == 1, 1/2, 1/6]/Sqrt[2]; bo *= fac; (* tetrahedral angle *)
                           If[mbstl == 1 && Max[vln[[bid]]] > 3, bo *= Sqrt[1/2]]; (* adjustment for multiple bonds on sulfur/phosphorus *)
                           Table[Transpose[{cr[[bid]], tubeF[BezierCurve[#, If[mbstl == 2, SplineDegree -> 6, Unevaluated[]]]] & /@ 
                                                      Partition[tm + ArrayPad[ConstantArray[RotationMatrix[2 Pi j/mul, del].bo, ml - 1], {{1}, {0}}], mh + 1, mh]}],
                                    {j, 0, mul - 1}]]];

              bondplot = MapThread[makebonds, {bondlist, nrms, bondtypes}];
              If[MatchQ[tips, All | "Bonds"], (* tooltips *)
                  bondch = Replace[#, Thread[{"Single", "Double", "Triple", _String} -> {"-", "=", "\[Congruent]", "-"}], {0, 1}] &;
                  bondplot = MapThread[Tooltip[#, Grid[{{Style[StringJoin[Insert[Sort[vertexTypes[[#4]], hillOrder[##] > 0 &], bondch[#2], 2]], Bold, Larger], ""},
                                                                             {Style["Bond no.", Bold], #3}, {Style["Bond type:", Bold], ToLowerCase[#2]},
                                                                             {Style["Bond length:", Bold], Row[{numDisplay[#5], unitLabel[$DefaultLengthUnit]}, " "]}}, Alignment -> {{Left, Left}}]] &,
                                                    {bondplot, edgeTypes, Range[e], bondlist, Apply[EuclideanDistance, pos[[#]]] & /@ bondlist}]];

              If[MemberQ[edgeTypes, "Aromatic"] && $VersionNumber >= 10. && arQ,
                  Block[{ars, cyc, rbt, rcol, ridx, rtu, swt}, (* aromatic ring depiction for explicit "Aromatic" bonds *)
                           cyc = FindCycle[skel, {5, 8}, All]; (* treat only five- to eight-membered rings *)
                           If[cyc =!= {},
                               swt = Fold[SetProperty[{#1, #2[[1]]}, EdgeWeight -> #2[[2]]] &,
                                                skel, Transpose[{UndirectedEdge @@@ edgeRules, edgeTypes}]];
                               rbt = Table[PropertyValue[{swt, ed}, EdgeWeight], {cc, cyc}, {ed, cc}];
                               ars = Pick[cyc, MatchQ[#, {"Aromatic" ..}] & /@ rbt];
                               If[ars =!= {}, ridx = VertexList /@ ars;
                                   rcol = First[Sort[Commonest[vertexTypes[[#]]], hillOrder[##] == 1 &]] & /@ ridx /. $colorRules;
                                   rtu = Composition[If[tubeQ, Tube[#, 8 $conversionFactor] &, Identity], incirc[pos[[#]], If[atomQ, 0.6, 0.7]] &] /@ ridx;
                                   bondplot = Join[bondplot, Transpose[{rcol, rtu}]]]]]];

              If[! atomQ, (* additional objects for atom-less styles *)
                  If[hiSet =!= None, (* atom highlights *)
                      bondplot = Join[If[tubeQ,
                                                  {makeTranslucent[#[[2]]], Sphere[pos[[#[[1]]]], 16 $conversionFactor]},
                                                  {Directive[AbsolutePointSize[12], #[[2]]], Point[pos[[#[[1]]]]]}] & /@ hiSet, bondplot]];

                  If[zpos =!= {}, zpos = zpos[[All, All, 1]]; (* force rendering of isolated atoms *)
                      bondplot = Join[MapThread[If[tubeQ, {#1, Sphere[#2, 14 $conversionFactor]}, {Directive[AbsolutePointSize[6], #1], Point[#2]}] &,
                                                               {Extract[cr, zpos], Extract[pos, zpos]}], bondplot]]],
              (* no bonds *) bondplot = {}];

          If[! FreeQ[(lbl = OptionValue[PlotLabel]), "Formula" | "MolarMass" | "MolecularWeight" | "ChemicalFormula" | "MolecularFormula"],
              lbl = lbl /. {"Formula" | "ChemicalFormula" | "MolecularFormula" :> ToString[buildFormula[vertexTypes], StandardForm, NumberMarks -> False],
                                "MolarMass" | "MolecularWeight" :> ToString[buildMolarMass[vertexTypes], TraditionalForm, NumberMarks -> False]}];

          res = Graphics3D[{atomplot, bondplot},
                                     FilterRules[Join[{PlotLabel -> lbl, ViewPoint -> view, opts}, DeleteCases[Options[MoleculeViewer], (PlotLabel -> _) | (ViewPoint -> _)]] /.
                                                     {(BaseStyle -> Automatic) -> (BaseStyle -> Apply[Directive, {If[atomQ && bondQ && tubeQ, CapForm["Butt"], Nothing],
                                                                                                                                                  If[bondQ && ! tubeQ, AbsoluteThickness[3], Nothing], Specularity[1, 50]}]),
                                                       (ImageSize -> Automatic) -> Unevaluated[],
                                                       (Method -> Automatic) -> (Method -> {"RotationControl" -> "Globe"})}, 
                                                     Options[Graphics3D]]];

          If[MatchQ[leg = OptionValue[PlotLegends], None | False], res,
              lchk = True; (* use special rules *)
              If[MatchQ[leg, Automatic | True] || ! FreeQ[leg, "Atoms" | Automatic],
                  If[$VersionNumber >= 11.,
                      atoms = Sort[DeleteDuplicates[vertexTypes], hillOrder],
                      atoms = Sort[DeleteDuplicates[vertexTypes], hillOrder[##] == 1 &]];
                  cr = Replace[atoms, cRules, {1}];
                  cc = Position[cr, (col_ /; translucentColorQ[col] && Last[col] == 0) | (Opacity[o_, r___] /; o == 0)];
                  atoms = Apply[Sequence, Delete[#, cc] & /@ {cr, atoms}],
                  If[(atoms = Cases[leg, l_List /; Complement[l, $atoms] === {}, {0, 2}]) =!= {},
                      atoms = DeleteDuplicates[First[atoms]];
                      atoms = Apply[Sequence, {Replace[atoms, cRules, {1}], atoms}],
                      lchk = False]];
 
              myLegendFunction = Framed[#, Background -> GrayLevel[0.863], FrameStyle -> None, RoundingRadius -> 5] &;
 
              If[$VersionNumber >= 9.,
                  mySwatch = SwatchLegend[##, LabelStyle -> {Bold, 12, GrayLevel[0.3]}, LegendFunction -> myLegendFunction, 
                                                           If[atomQ, LegendMarkers -> "SphereBubble", LegendMarkerSize -> 16]] &;
                  If[TrueQ[lchk],
                      If[FreeQ[leg, PointLegend | SwatchLegend | Placed],
                          Legended[res, mySwatch[atoms]],
                          Switch[leg,
                                     Placed[Except[_PointLegend | _SwatchLegend], _],
                                     Legended[res, ReplacePart[leg, 1 -> Append[mySwatch[atoms],
                                                                                                       LegendLayout -> If[MatchQ[Last[leg], Above | Below | Bottom | Top], "Row", "Column"]]]],
                                     PointLegend[__] | SwatchLegend[__], Legended[res, ReplacePart[leg, 1 -> atoms]],
                                     Placed[(PointLegend | SwatchLegend)[__], __], Legended[res, ReplacePart[leg, {1, 1} -> atoms]],
                                     _, Legended[res, leg]]],
                      Legended[res, leg]],
                  If[TrueQ[lchk],
                      Row[{Show[res, ImageSize -> Medium],
                               Deploy[myLegendFunction[Grid[MapThread[{Style["\[FilledSquare]", #1, 16],
                                                                                                  Style[#2, Bold, 12, GrayLevel[0.3], FontFamily -> "Helvetica"]} &, {atoms}]]]]},
                             Spacer[10], Editable -> False],
                      res]]]
          ]

End[ ]

SetAttributes[{MoleculeViewer, MoleculeRecognize, MakeMOL, MakeInChIKey,
                      GetCACTUS, GetChemicalData, GetChemSpider, GetPDB, GetProteinData, GetPubChem, GetZINC,
                      MoleculeDraw, GetMolecule, SetMolecule, JChemPaint, JDraw, JME, LaunchJME, GetJMEMOL, GetJMESMILES, SetJMEMOL,
                      ChEMBLMOLTo3D, ChEMBLSMILESTo3D, ChEMBLInChITo3D, ChEMBLInChIToMOL, ChEMBLSMILESToMOL,
                      ChEMBLAddHydrogens, ChEMBLSanitize, ChEMBLMakeLineNotation, ChEMBLMakeInChIKey, ChEMBLDescriptors,
                      RunImago, RunOpenBabel, OPSINNameToStructure},
                    ReadProtected];

Protect[MoleculeViewer, MoleculeRecognize, MakeMOL, MakeInChIKey,
            GetCACTUS, GetChemicalData, GetChemSpider, GetPDB, GetProteinData, GetPubChem, GetZINC,
            MoleculeDraw, GetMolecule, SetMolecule, JChemPaint, JDraw, JME,
            ChEMBLMOLTo3D, ChEMBLSMILESTo3D, ChEMBLInChITo3D, ChEMBLInChIToMOL, ChEMBLSMILESToMOL,
            ChEMBLAddHydrogens, ChEMBLSanitize, ChEMBLMakeLineNotation, ChEMBLMakeInChIKey, ChEMBLDescriptors,
            RunImago, RunOpenBabel];

EndPackage[ ]