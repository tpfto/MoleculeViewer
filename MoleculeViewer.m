
(* :Title: Molecule Viewer *)

(* :Author: J. M. *)

(* :Summary:

     This package implements an improved molecule viewer derived from Mathematica's default viewer,
     as well as providing a number of useful utility functions for obtaining and visualizing molecular models.

 *)

(* :Copyright:

     © 2017-2019 by J. M. (pleasureoffiguring(AT)gmail(DOT)com)
     This work is free. It comes without any warranty, to the extent permitted by applicable law.
     You can redistribute and/or modify it under the terms of the WTFPL (http://www.wtfpl.net/).

 *)

(* :Package Version: 2.1 *)

(* :Mathematica Version: 8.0 *)

(* :History:

     2.1 - removed support for deprecated ChemSpider SOAP API
     2.0 - improved handling of double bonds, added PubChem support, switched to new ChemSpider API, better handling of files and URLs
     1.1 - fixed a few reported bugs, added explicit input dialog for $ChemSpiderAPIToken, support suboptions for RunOpenBabel
     1.0 - initial release

*)

(* :Keywords:
     3D, chemistry, ChemSpider, graphics, JME, molecules, PubChem, Open Babel *)

(* :Limitations:
     currently limited handling of imported chemical file formats
     needs Open Babel to be installed for some functionality
*)

(* :Requirements:

     Requires Open Babel (http://openbabel.org/) and JME (http://www.molinspiration.com/jme/) to be installed.
     ChemSpider API key (https://developer.rsc.org/) needed for ChemSpider functionality.

*)   

BeginPackage["MoleculeViewer`"]
Unprotect[GetCACTUS, GetChemicalData, GetChemSpider, GetJMESMILES, GetPubChem, JME, LaunchJME, MoleculeViewer, RunOpenBabel];

(* usage messages *)

MoleculeViewer::usage = "MoleculeViewer[expr] displays a molecular model corresponding to expr. expr can be a file name of a chemical file format, a chemical name, a SMILES or an InChI string, or a list containing structural information."

GetChemicalData::usage = "GetChemicalData[expr] uses the built-in ChemicalData function to generate structural information suitable for MoleculeViewer. GetChemicalData[expr, \"SMILES\"] or GetChemicalData[expr, \"InChI\"] will generate the corresponding line notation."

GetChemSpider::usage = "GetChemSpider[expr] uses the ChemSpider API to generate structural information suitable for MoleculeViewer. GetChemSpider[expr, \"SMILES\"] or GetChemSpider[expr, \"InChI\"] will generate the corresponding line notation."

GetJMESMILES::usage = "GetJMESMILES[] extracts the SMILES of the chemical structure currently displayed in the JME applet.";

GetPubChem::usage = "GetPubChem[expr] uses the PubChem API to generate structural information suitable for MoleculeViewer. GetPubChem[expr, \"SMILES\"] or GetPubChem[expr, \"InChI\"] will generate the corresponding line notation."

LaunchJME::usage = "LaunchJME[] launches an instance of the JME applet.";

RunOpenBabel::usage = "RunOpenBabel[expr] is an interface to Open Babel, for generating structural information suitable for MoleculeViewer. expr can be a SMILES or an InChI string.";

$ChemSpiderAPIToken; $JME; GetCACTUS; JME;
(* declare symbols for older versions *)
Highlighted; MemoryConstraint; PlotLegends;
RunProcess; ProcessEnvironment; ProcessDirectory;

(* error/warning messages *)

MoleculeViewer::badchem = "Invalid structural information given to MoleculeViewer.";
MoleculeViewer::badstrng = "Valid structural information could not be generated from the input string.";
MoleculeViewer::ftconv = "The file `1` is not in a format supported by Import. Conversion with OpenBabel will be attempted.";
MoleculeViewer::is2d = "Two-dimensional coordinates detected; embedding to 3D.";

GetChemSpider::token = "ChemSpider API key in $ChemSpiderAPIToken not detected or invalid. Please obtain one from https://developer.rsc.org/.";
GetChemSpider::obsmet = "Method \[Rule] \"SOAP\" is now obsolete. Please obtain a new API key from https://developer.rsc.org/.";
JME::nojme = "Java Molecular Editor could not be loaded.";

RunOpenBabel::nobab = "Open Babel not installed; please download from http://openbabel.org/.";
MoleculeViewer::mem = RunOpenBabel::mem = "Memory allocated by Open Babel exceeded `1` bytes, and was then aborted. Increasing the value of the MemoryConstraint option may yield a result."; 
MoleculeViewer::time = RunOpenBabel::time = "Time spent by Open Babel exceeded `1` seconds, and was then aborted. Increasing the value of the TimeConstraint option may yield a result.";

Begin["`Private`"]

Needs["Utilities`URLTools`"]

$debug = False;
If[$VersionNumber < 10.2, Nothing := Unevaluated[]];

$atoms = Array[ElementData[#, "Symbol"] &, 118] ~Join~ {"D", "T"};

$colorRules = Thread[$atoms -> ((Apply[RGBColor, IntegerDigits[FromDigits[#, 16], 256, 3]/255.] & /@
                                (* colors from ElementData[elem, "IconColor"] *)
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

(* Pyykkö's covalent radii; https://doi.org/10.1002/chem.200800987 ; entries for deuterium and tritium were slightly enlarged from hydrogen *)
$elementCovalentRadii = {32, 46, 133, 102, 85, 75, 71, 63, 64, 67, 155, 139, 126, 116, 111, 
                                       103, 99, 96, 196, 171, 148, 136, 134, 122, 119, 116, 111, 110, 112, 
                                       118, 124, 121, 121, 116, 114, 117, 210, 185, 163, 154, 147, 138, 
                                       128, 125, 125, 120, 128, 136, 142, 140, 140, 136, 133, 131, 232, 
                                       196, 180, 163, 176, 174, 173, 172, 168, 169, 168, 167, 166, 165, 
                                       164, 170, 162, 152, 146, 137, 131, 129, 122, 123, 124, 133, 144, 
                                       144, 151, 145, 147, 142, 223, 201, 186, 175, 169, 170, 171, 172, 
                                       166, 166, 168, 168, 165, 167, 173, 176, 161, 157, 149, 143, 141, 
                                       134, 129, 128, 121, 122, 136, 143, 162, 175, 165, 157, (* D and T *) 34, 36};
$sizeRules = Thread[$atoms -> N[$elementCovalentRadii/2]];

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
                                                               _, False] && 0 <= Last[col] <1]
]

makeTranslucent[dir_] := If[translucentColorQ[dir] || Head[dir] === Opacity, dir, Opacity[0.6, dir]]

(* Hill system order for atoms *)
hillOrder["C", _] = 1; hillOrder[_, "C"] = -1;
hillOrder["H", _] = 1; hillOrder[_, "H"] = -1;
hillOrder["D", _] = 1; hillOrder[_, "D"] = -1;
hillOrder["T", _] = 1; hillOrder[_, "T"] = -1;
hillOrder[s1_, s2_] := Order[s1, s2];

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

percentEncode[s_String] := If[$VersionNumber >= 10., StringReplace[URLEncode[s], "+" -> "%20"], 
                                             StringReplace[s, RegularExpression["[^\\w/\\?=:._~]"] :> 
                                                                      StringJoin["%", ToUpperCase[IntegerString[First[ToCharacterCode["$0"]], 16]]]]]

getNormal[mat_?MatrixQ, sv : (True | False) : False] := With[{svd = SingularValueDecomposition[Standardize[mat, Mean, 1 &], {-1}, Tolerance -> 0]}, 
     If[TrueQ[sv], FlattenAt[Flatten[Rest[svd], {1, 3}], 1], Flatten[Last[svd]]]]

(* for XYZ format *)
inferBonds[vertexTypes_List, atomPositions_List, minDistance_: 40., bondTolerance_: 25., maxNeighbors_Integer: 8] /; Length[vertexTypes] == Length[atomPositions] :=
       Module[{atomP, curAtom, maxD, near, nearD, neighbors, nF, pos, rads},
                   pos = SetPrecision[atomPositions, MachinePrecision];
                   rads = N[vertexTypes /. Thread[$atoms -> $elementCovalentRadii]];
                   maxD = Max[rads] + bondTolerance;
  
                   nF = Nearest[pos -> Automatic];
                   neighbors = Function[idx,
                                                   curAtom = rads[[idx]]; atomP = pos[[idx]];
                                                   near = DeleteCases[nF[atomP, {maxNeighbors, curAtom + maxD}], idx];
                                                   nearD = SquaredEuclideanDistance[#, atomP] & /@ pos[[near]];
                                                   Pick[near, MapThread[minDistance < #1 < #2 &,
                                                                                    {nearD, (curAtom + rads[[near]] + bondTolerance)^2}]]];
  
                   Union[Sort /@ Flatten[Array[Thread[# -> neighbors[#]] &, Length[vertexTypes]]]]]

(* validate input format for viewer *)
validateChem[{vertexTypes_, edgeRules_, edgeTypes_, atomPositions_}] :=
           With[{dim = Dimensions[atomPositions], e = Length[edgeRules], v = Length[vertexTypes]},
                   (e == Length[edgeTypes] && v == First[dim]) &&
                   VectorQ[vertexTypes, StringQ] && Complement[vertexTypes, $atoms] === {} && 
                   VectorQ[edgeRules, MatchQ[#, _Integer?Positive -> _Integer?Positive] &] && 
                   VectorQ[edgeTypes, StringQ] && Complement[edgeTypes, {"Single", "Aromatic", "Double", "Triple"}] === {} &&
                   (MatrixQ[atomPositions, NumericQ] && MatchQ[Last[dim], 2 | 3] && Norm[atomPositions, Infinity] > 0)]

inchirx = RegularExpression["^(InChI=1)(S*)/(\\w+)/(.+)"];

smilesElems = With[{syms = Drop[$atoms, -2], pos = Transpose[{{1, 5, 6, 7, 8, 9, 15, 16, 19, 23, 39, 53, 74, 92}}]}, Delete[syms, pos] ~Join~ Extract[syms, pos]];
validSymbols = smilesElems ~Join~ {"c", "n", "o", "p", "s", "as", "se"} ~Join~ {"-", "=", "#", ":", "(", ")", "[", "]", "%", ".", "/", "\\", "@", "+", DigitCharacter};

validSMILESQ[s_] := StringQ[s] && s === StringJoin[StringCases[s, validSymbols]]

processChemicalFile[fName_String, type : (_String | Automatic) : Automatic] := 
           Module[{atoms, bonds, ext, fileName, pos, prop, res},

                       If[FileExistsQ[fName], fileName = fName,
                           fileName = FindFile[fName];
                           If[fileName === $Failed, Message[MoleculeViewer::fnfnd, fName]; Return[$Failed, Module]]];

                       ext = type;
                       If[ext === Automatic,
                           ext = FileFormat[fileName];
                           If[ext === "", Message[MoleculeViewer::fftype, fName]; Return[$Failed, Module]]];

                       prop = Quiet[Import[fileName, {ext, "Elements"}]];
                       If[prop === $Failed || 
                           FreeQ[prop, "VertexTypes" | "EdgeRules" | "EdgeTypes" | "VertexCoordinates"],
                           Message[MoleculeViewer::ftconv, fName]; Return[openBabelConvert[fileName, FileExtension[fileName]], Module]];

                       If[FreeQ[prop, "EdgeRules" | "EdgeTypes"],
                           (* xyz or PDB format *)
                           res = Import[fileName, {ext, {"VertexTypes", "VertexCoordinates"}}];
                           If[res =!= $Failed,
                              {atoms, pos} = res; bonds = inferBonds[atoms, pos];
                              {atoms, bonds, ConstantArray["Single", Length[bonds]], pos},
                               res],
                           (* normal mode *)
                           res = Import[fileName, {ext, {"VertexTypes", "EdgeRules", "EdgeTypes", "VertexCoordinates"}}];
                           If[res =!= $Failed && StringMatchQ[ext, "HIN" |"MOL2" | "SDF", IgnoreCase -> True], res = Transpose[res]];
                           If[StringMatchQ[ext, "HIN", IgnoreCase -> True], (* Å to pm *) res = Map[Function[mol, MapAt[100 # &, mol, {4}]], res]];
                           res]]

processChemicalURL[url_String, type : (_String | Automatic) : Automatic] := Block[{fileName, temp},
           fileName= FileNameJoin[{$TemporaryDirectory, If[$VersionNumber >= 10.,
                                                                                     Last[Lookup[URLParse[url], "Path"]],
                                                                                     Last[StringSplit[url, {"/", "\\"}]]]}];
           Internal`WithLocalSettings[temp = If[$VersionNumber >= 9., URLSave, Utilities`URLTools`FetchURL][url, fileName],
                                                   If[temp =!= $Failed, processChemicalFile[temp, type], $Failed],
                                                   DeleteFile[temp]]]

processEntity[ent_Entity] := Module[{res},
           If[EntityTypeName[ent] =!= "Chemical", Return[$Failed, Module]];
           res = EntityValue[ent, EntityProperty["Chemical", #] & /@ {"VertexTypes", "EdgeRules", "EdgeTypes", "AtomPositions"}];
           If[validateChem[res], res, $Failed]]

testOpenBabel := testOpenBabel = Quiet[Check[Import["!obabel -H", "Text"]; True, Message[RunOpenBabel::nobab]; False]]

(* defaults for Open Babel methods *)
$obConOpts = {"GeneratedConformers" -> 20, "Method" -> "Genetic"};
$obConGOpts = {"Children" -> Automatic, "Mutability" -> Automatic, "Converge" -> Automatic, "Score" -> Automatic};
$obMinOpts = {"CutOff" -> False, "ForceField" -> "MMFF94s", "MaxIterations" -> 2000, "Newton" -> False, "SteepestDescent" -> False, "Tolerance" -> Automatic, "UpdateNonBondedFrequency" -> 10};

makeStringRules[expr_List] := Block[{o, s}, Flatten[expr /. (o_ -> s_) :> (System`Utilities`StringName[o] -> s)]]

makeOBSwitches[met_, mopt_List] := Module[{mok, mon, res, tmp},
        If[! StringQ[met], Return[$Failed, Module]];
        mon = makeStringRules[mopt];
        res = {"--" <> ToLowerCase[met]};
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
                              AppendTo[res, "--ff " <> ("Method" /. tmp /. "Automatic" -> "MMFF94s")],
                              "GAFF" | "Ghemical" | "MMFF94" | "MMFF94s" | "UFF",
                              AppendTo[res, "--ff " <> tmp],
                              "Genetic" | {"Genetic", __Rule},
                              If[ListQ[tmp],
                                  tmp = makeStringRules[Rest[tmp]];
                                  tmp = FilterRules[Join[tmp, $obConGOpts], $obConGOpts],
                                  tmp = $obConGOpts];
                              AppendTo[res, "--score " <> ToLowerCase["Score" /. tmp /. Automatic -> "Energy"]];
                              tmp = DeleteCases[MapThread[If[MatchQ[#2, Automatic | "Automatic"], "", #1 <> IntegerString[Round[#2]]] &,
                                                                            {{"--children ", "--mutability ", "--converge "}, {"Children", "Mutability", "Converge"} /. tmp}], ""];
                              res = Join[res, tmp];];,
                   "Minimize",
                   mok = FilterRules[Join[mon, $obMinOpts], $obMinOpts];
                   AppendTo[res, "--ff " <> ("ForceField" /. mok /. Automatic | "Automatic" -> "MMFF94s")];
                   AppendTo[res, "--steps " <> IntegerString[Round["MaxIterations" /. mok]]];
                   AppendTo[res, "--crit " <> ToString[CForm["Tolerance" /. mok /. Automatic -> 1*^-7], InputForm]];
                   AppendTo[res, "--freq " <> IntegerString[Round["UpdateNonBondedFrequency" /. mok]]];
                   If[TrueQ["SteepestDescent" /. mok], AppendTo[res, "--sd"]];
                   If[TrueQ["Newton" /. mok], AppendTo[res, "--newton"]];
                   Switch[tmp = "CutOff" /. mok,
                              True, AppendTo[res, "--cut"],
                              _List, 
                              res = Join[res, Prepend[MapThread[StringJoin, {{"--rvdw ", "--rele "}, ToString /@ PadRight[tmp, 2, 10]}], "--cut"]]];,
                   _, res = $Failed];
        res]

openBabelConvert[fileName_String, ex_String] := 
       Module[{ext = ToLowerCase[ex], pars = {"-c", "-h"}, args, ec, msg, opo, out, proc, res},
                   If[! testOpenBabel, Return[$Failed, Module]];

                   out = {"MOL", {"VertexTypes", "EdgeRules", "EdgeTypes", "VertexCoordinates"}};

                   constrainedEvaluate[If[$VersionNumber >= 10.,
                                                     args = {"obabel", "-i", ext, fileName, "-omol", "--gen3d"} ~Join~ pars;
                                                     proc = Check[RunProcess[args], Return[$Failed, Module]];
                                                     res = proc["StandardOutput"]; ec = proc["ExitCode"]; msg = proc["StandardError"];
                                                     If[ec == 0 && res =!= "",
                                                         Check[ImportString[res, out], Echo[StringTrim[msg]]; $Failed],
                                                         msg = StringTrim[StringReplace[msg, {"0 molecules converted" -> "", "\r" | "\n" -> " "}]];
                                                         If[msg =!= "", Echo[msg],
                                                             If[ec != 0, Echo["OpenBabel exit code: " <> IntegerString[ec]]]]; $Failed],
                                                     args = StringJoin[Riffle[{"!obabel", "-i", ext, fileName, "-omol", "--gen3d"} ~Join~ pars, " "]];
                                                     Check[Import[args, out], Print[StringTrim[msg]]; $Failed]], MoleculeViewer,
                                                 536870912 (* 512 MB *), 120 (* 2 min. *)]]

Options[RunOpenBabel] = {MemoryConstraint -> Automatic, Method -> Automatic, ProcessDirectory -> Inherited,
                                         ProcessEnvironment -> Inherited, RunProcess -> Automatic, TimeConstraint -> Automatic};

RunOpenBabel[chems : {__String}, opts___] := RunOpenBabel[#, opts] & /@ chems

RunOpenBabel[chem_String, opts : OptionsPattern[]] := Module[{args, autoSet, ec, incq, met, mopt, msg, opo, out, pars, proc, res},
      If[! testOpenBabel, Return[$Failed, Module]];

      autoSet = {"Minimize", {"ForceField" -> "MMFF94s", "MaxIterations" -> 2000, "Tolerance" -> 1*^-7}};
      met = OptionValue[Method];
      If[ListQ[met],
          {met, mopt} = {First[met], Flatten[Rest[met]]},
          mopt = {}];
      met = System`Utilities`StringName[met];
      If[met === "Automatic", {met, mopt} = autoSet];
      pars = makeOBSwitches[met, mopt];
      If[pars === $Failed,
          Message[RunOpenBabel::mthd, met, Automatic];
          pars = makeOBSwitches @@ autoSet];

      incq = If[StringMatchQ[chem, inchirx], "-iinchi", Nothing];
      out = {"MOL", {"VertexTypes", "EdgeRules", "EdgeTypes", "VertexCoordinates"}};

      constrainedEvaluate[If[TrueQ[OptionValue[RunProcess] /. Automatic -> True] && $VersionNumber >= 10.,
                                        args = {"obabel", incq, "-:"<> chem, "-omol", "--gen3d", "-c", "-h"} ~Join~ pars;
                                        If[TrueQ[$debug], Print[Defer[RunProcess][args]]];
                                        proc = Check[RunProcess[args, FilterRules[Join[{opts}, Options[RunOpenBabel]], Options[RunProcess]]], 
                                                             Return[$Failed, Module]];
                                        res = proc["StandardOutput"]; ec = proc["ExitCode"]; msg = proc["StandardError"];
                                        If[ec == 0 && res =!= "", 
                                            opo = StringPosition[res, "OpenBabel"];
                                            If[opo =!= {},
                                                res = ImportString[StringDrop[res, opo[[1, 1]] - 3], out];
                                                If[validateChem[res], res, $Failed],
                                                $Failed],
                                            msg = StringTrim[StringReplace[msg, {"0 molecules converted" -> "", "\r" | "\n" -> " "}]];
                                            If[msg =!= "", Echo[msg], 
                                                If[ec != 0, Echo["OpenBabel exit code: " <> IntegerString[ec]]]]; $Failed],
                                        args = StringJoin[Riffle[{"!obabel", incq, "-:" <> chem,
                                                                             "-omol", "--gen3d", "-c", "-h"} ~Join~ pars, " "]];
                                        If[TrueQ[$debug], Print[Defer[Import][args, "Text"]]];
                                        res = Check[Import[args, "Text"], $Failed];
                                        If[res =!= $Failed,
                                            opo = StringPosition[res, "OpenBabel"];
                                            If[opo =!= {},
                                                res = ImportString[StringDrop[res, opo[[1, 1]] - 3], out];
                                                If[validateChem[res], res, $Failed],
                                                $Failed],
                                            $Failed]], RunOpenBabel,
                                    OptionValue[MemoryConstraint] /. Automatic -> 536870912 (* 512 MB *),
                                    OptionValue[TimeConstraint] /. Automatic -> 120 (* 2 min. *)]]

SetAttributes[GetCACTUS, Listable];

GetCACTUS[arg_String] := Module[{in, out, url},
     url = "https://cactus.nci.nih.gov/chemical/structure/``/file?format=sdf&get3d=true";
     in = If[StringMatchQ[arg, inchirx], arg, percentEncode[arg]];
     If[TrueQ[$debug], Print[Defer[Import][ToString[StringForm[url, in]]]]];
     out = {"MOL", {"VertexTypes", "EdgeRules", "EdgeTypes", "VertexCoordinates"}};
     If[$VersionNumber >= 10., Import[StringTemplate[url] @ in, out], Import[ToString[StringForm[url, in]], out]]]

SetAttributes[GetChemicalData, Listable];

GetChemicalData[str_String, out : (Automatic | "SMILES" | "InChI") : Automatic, OptionsPattern[{"MaxResults" -> 1}]] := Module[{chk, name, res, tmp},

     name = StringJoin[StringReplacePart[#, ToUpperCase[StringTake[#, 1]], {1, 1}] & /@ StringSplit[str, " " | "-"]];
     If[TrueQ[$debug], Print[Defer[ChemicalData][name]]];
     chk = ChemicalData[name, "StandardName"];
     If[chk === $Failed, Return[chk, Module]];
     chk = Take[Flatten[{chk}], Max[1, Round[OptionValue["MaxResults"] /. All -> Infinity]] /. Infinity -> All];

     If[out === Automatic,
         res = Table[ChemicalData[c, prop], {c, chk}, {prop, {"VertexTypes", "EdgeRules", "EdgeTypes", "AtomPositions"}}];
         If[$VersionNumber >= 9., res = MapAt[QuantityMagnitude, res, {All, -1}]];
         res = Replace[res, x_ /; ! validateChem[x] :> $Failed, {1}],
         res = Table[If[out === "SMILES", 
                               tmp = Cases[ChemicalData[c, #] & /@ {"IsomericSMILES", "SMILES"}, _String]; 
                               If[Length[tmp] > 0, First[tmp], $Failed], 
                               ChemicalData[c, out] /. Except[_String] -> $Failed], {c, chk}]];
     If[Length[res] == 1, First[res], res]]

(* **ADVANCED USERS**:
     You can remove the code sandwiched by the row of dashes below and
     permanently assign your API token to $ChemSpiderAPIToken instead
     by editing this package file:

     $ChemSpiderAPIToken = ;

 *)

(* ---------------- *)

$dir = DirectoryName[$InputFileName];
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

     met = OptionValue[Method];
     If[met === Automatic, met = "RESTful",
         met = met /. {"Legacy" -> "SOAP", _ -> "RESTful"}];

     n = Round[OptionValue["MaxResults"]];
     Switch[met,
                "RESTful",
                gCSRESTful[arg, out, "MaxResults" -> n],
                "SOAP",
	        Message[MoleculeViewer::obsmet]; $Failed,
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

   If[! (r === $Failed || Head[r] === Failure), r = "output" /. r;
       If[r =!= "output", r, $Failed], $Failed]]

gCSRESTful[csid_Integer?Positive, out : (Automatic | "SMILES" | "InChI") : Automatic, opts___] := Module[{pars, req, rsp, type, url},

      If[! ValueQ[$ChemSpiderAPIToken] || ! StringQ[$ChemSpiderAPIToken], 
          Message[GetChemSpider::token]; Return[$Failed, Module]];

      url = "https://api.rsc.org/compounds/v1/records/" <> IntegerString[csid] <> "/details";
      type = If[out === Automatic, "mol3D", "smiles"];
      pars = {"Method" -> "GET", "Headers" -> {"apikey" -> $ChemSpiderAPIToken, "Content-Type" -> "application/json"}};

      If[$VersionNumber >= 11.,
          pars = Association[Append[pars, "Query" -> {"fields" -> type}]];
          req = HTTPRequest[url, pars];
          If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req], "JSON"]]];
          rsp = URLExecute[req, "JSON"],
          pars = Append[pars, "Parameters" -> {"fields" -> type}];
          If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
          rsp = ImportString[URLFetch[url, pars], "JSON"]];

      If[! (rsp === $Failed || Head[rsp] === Failure), rsp = type /. rsp;
          If[! StringMatchQ[rsp, type],
              Which[out === Automatic,
                        ImportString[rsp, {"MOL", {"VertexTypes", "EdgeRules", "EdgeTypes", "VertexCoordinates"}}],
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

                  If[! (rsp === $Failed || Head[rsp] === Failure), rsp = "records" /. rsp;
                      If[ListQ[rsp], rsp = type /. rsp;
                          If[VectorQ[rsp, ! StringMatchQ[#, type] &],
                              Which[out === Automatic,
                                        ImportString[#, {"MOL", {"VertexTypes", "EdgeRules", "EdgeTypes", "VertexCoordinates"}}] & /@ rsp,
                                        out === "InChI", csconv[rsp, "SMILES", "InChI"],
                                        True, rsp],
                              $Failed], $Failed], $Failed]]

cscheck[id_String] := With[{url = "https://api.rsc.org/compounds/v1/filter/" <> id <> "/status", 
                                           pars = {"Method" -> "GET", "Headers" -> {"apikey" -> $ChemSpiderAPIToken, "Content-Type" -> "application/json"}}},
                                         If[$VersionNumber >= 11.,
                                             URLExecute[HTTPRequest[url, Association[pars]], "JSON"],
                                             ImportString[URLFetch[url, pars], "JSON"]]]

gCSRESTful[str_String, out : (Automatic | "SMILES" | "InChI") : Automatic, OptionsPattern[{"MaxResults" -> 1}]] := 
      Module[{chk, nres, pars, qid, req, res, stat, url},

                  If[! ValueQ[$ChemSpiderAPIToken] || ! StringQ[$ChemSpiderAPIToken], 
                      Message[GetChemSpider::token]; Return[$Failed, Module]];

                  If[StringMatchQ[str, NumberString], 
                      Return[gCSRESTful[Round[FromDigits[str]], out], Module]];

                  url = "https://api.rsc.org/compounds/v1/filter/name";
                  pars = {"Method" -> "POST", "Headers" -> {"apikey" -> $ChemSpiderAPIToken, "Content-Type" -> "application/json"}, 
                              "Body" -> ExportString[{"name" -> str}, "JSON"]};

                  If[$VersionNumber >= 11.,
                      req = HTTPRequest[url, Association[pars]];
                      If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req], "JSON"]]];
                      qid = URLExecute[req, "JSON"],
                      If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
                      qid = ImportString[URLFetch[url, pars], "JSON"]];

                  If[! (qid === $Failed || Head[qid] === Failure), 
                      qid = "queryId" /. qid, Return[$Failed, Module]];
                  If[qid === "queryId", Return[$Failed, Module]];

                  While[Pause[RandomReal[{1.5, 2.5}]]; chk = cscheck[qid];
                           ! (chk === $Failed || Head[chk] === Failure) && ("status" /. chk) === "Processing"];

                  If[! (chk === $Failed || Head[chk] === Failure),
                      {stat, nres} = {"status", "count"} /. chk;
                      If[stat === "Complete" && nres > 0,

                          nres = Min[Max[1, Round[OptionValue["MaxResults"] /. All -> Infinity]], nres];
                          url = "https://api.rsc.org/compounds/v1/filter/" <> qid <> "/results";
                          pars = {"Method" -> "GET", "Headers" -> {"apikey" -> $ChemSpiderAPIToken, "Content-Type" -> "application/json"}};

                          If[$VersionNumber >= 11.,
                              pars = Append[pars, "Query" -> {"start" -> 0, "count" -> nres}];
                              req = HTTPRequest[url, Association[pars]];
                              If[TrueQ[$debug], Print[Defer[URLExecute][InputForm[req], "JSON"]]];
                              res = URLExecute[req, "JSON"],
                              url = url <> "?start=0&count=" <> IntegerString[nres];
                              If[TrueQ[$debug], Print[Defer[URLFetch][url, pars]]];
                              res = ImportString[URLFetch[url, pars], "JSON"]];

                          If[! (res === $Failed || Head[res] === Failure),
                              res = "results" /. res;
                              If[TrueQ[$debug], Print[res]];
                              If[res === "results", Return[$Failed, Module]];
                              If[Length[res] == 1, res = First[res]];
                              gCSRESTful[res, out],
                              $Failed], $Failed], $Failed]]

gCSRESTful[l_List, rest___] := gCSRESTful[#, rest] & /@ l

SetAttributes[GetPubChem, Listable];

GetPubChem[cid_Integer?Positive, out : (Automatic | "SMILES" | "InChI") : Automatic, opts___] := Module[{res, url},
     If[out === Automatic,
         url = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" <> IntegerString[cid] <> "/record/SDF?record_type=3d";
         If[TrueQ[$debug], Print[Defer[Import][url, "MOL"]]];
         Import[url, {"MOL", {"VertexTypes", "EdgeRules", "EdgeTypes", "VertexCoordinates"}}],
         url = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" <> IntegerString[cid] <> "/property/";
         res = Quiet[If[out === "SMILES",
                               Check[If[TrueQ[$debug], Print[Defer[Import][url <> "IsomericSMILES/CSV", "CSV"]]];
                                         Import[url <> "IsomericSMILES/CSV", "CSV"],
                                         Check[If[TrueQ[$debug], Print[Defer[Import][url <> "CanonicalSMILES/CSV", "CSV"]]];
                                                   Import[url <> "CanonicalSMILES/CSV", "CSV"],
                                                   $Failed]],
                               Check[If[TrueQ[$debug], Print[Defer[Import][url <> out <> "/CSV", "CSV"]]];
                                         Import[url <> out <> "/CSV", "CSV"], $Failed]]];
         If[res =!= $Failed, res[[2, 2]], $Failed]]]

GetPubChem[str_String, out : (Automatic | "SMILES" | "InChI") : Automatic, OptionsPattern[{"MaxResults" -> 1}]] := Module[{res, id, url},

     If[StringMatchQ[str, NumberString],
         Return[GetPubChem[Round[FromDigits[res]], out], Module]];

     url = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" <> percentEncode[str] <> "/cids/JSON";
     If[TrueQ[$debug], Print[Defer[Import][url, "JSON"]]];
     id = Quiet[Check[Import[url, "JSON"], $Failed]];
     If[id =!= $Failed,
         res = Take["CID" /. ("IdentifierList" /. id), Max[1, Round[OptionValue["MaxResults"] /. All -> Infinity]] /. Infinity -> All];
         If[Length[res] == 1, res = First[res]]; GetPubChem[res, out],
         $Failed]]

InstallJME[] := Block[{dir = FileNameJoin[{$UserBaseDirectory, "Applications", "JME", "Java"}], applet, url},
         applet = FileNameJoin[{dir, "JME.jar"}];
         url = "http://www.molinspiration.com/jme/doc/JME.jar";
         If[FileExistsQ[applet], True,
             CreateDirectory[dir];
             Check[If[$VersionNumber >= 9., URLSave, Utilities`URLTools`FetchURL][url, applet]; True,
                       Message[JME::nojme]; False]]]

LaunchJME[] := JLink`JavaBlock[If[InstallJME[], JLink`InstallJava[]; Unprotect[$JME];
                                                     $JME = JLink`JavaNew["JME"]; Protect[$JME];
                                                     JLink`AppletViewer[$JME, {}]; $JME,
                                                     Message[JME::nojme]; $Failed]]

GetJMESMILES[] := JLink`JavaBlock[If[ValueQ[$JME], $JME @ smiles[], Message[JME::nojme]; $Failed]]

JME[] := If[ValueQ[$JME], GetJMESMILES[], LaunchJME[]];

(* MoleculeViewer *)

Options[MoleculeViewer] = SortBy[Join[{Background -> RGBColor[0., 43./255, 54./255], BaseStyle -> Automatic, Boxed -> False, ColorRules -> Automatic, 
                                                            Highlighted -> None, Lighting -> "Neutral", PlotLegends -> None, SphericalRegion -> True, Tooltip -> None, ViewPoint -> Automatic}, 
                                                           Select[Options[Graphics3D], FreeQ[#, Background | BaseStyle | Boxed | ImageSizeRaw | Lighting | SphericalRegion | ViewPoint] &]], First];

SyntaxInformation[MoleculeViewer] = {"ArgumentsPattern" -> {_, OptionsPattern[MoleculeViewer]}};

MoleculeViewer[$Failed, opts___] := $Failed

MoleculeViewer[args___] := Block[{a, res},
            a = System`Private`Arguments[MoleculeViewer[args], 1];
            res /; TrueQ[a =!= {} &&
                               (res = iMoleculeViewer @@ Flatten[a, 1]) =!= $Failed &&
                               MatchQ[Head[res], Graphics3D | Legended | List | Row]]]

iMoleculeViewer[dat : {__String} | {__Entity}, opts___] := iMoleculeViewer[#, opts] & /@ dat
iMoleculeViewer[dat : {{_List, _List} ..}, opts___] := iMoleculeViewer[#, opts] & /@ dat
iMoleculeViewer[dat : {{_List, _List, _List, _List} ..}, opts___] := iMoleculeViewer[#, opts] & /@ dat

iMoleculeViewer[str_String, opts___] := Block[{res = $Failed},
 Quiet[Catch[Scan[Function[op, With[{r = op[str]}, If[r =!= $Failed, If[TrueQ[$debug], Print["Function used: ", op]]; Throw[Set[res, r], "MoleculeViewer"]]]],
                            Which[FileExistsQ[str] || FindFile[str] =!= $Failed, {processChemicalFile},
                                      StringMatchQ[str, "http://*" | "ftp://*" | "https://*"], {processChemicalURL},
                                      StringMatchQ[str, inchirx] || validSMILESQ[str], {RunOpenBabel, GetCACTUS},
                                      True, {GetChemicalData, GetPubChem, GetCACTUS, GetChemSpider}]], "MoleculeViewer"],
          {ChemicalData::notent, Import::fmterr, Utilities`URLTools`FetchURL::conopen, Utilities`URLTools`FetchURL::httperr}];
 If[res =!= $Failed, iMoleculeViewer[res, opts], Message[MoleculeViewer::badstrng]; res]]

iMoleculeViewer[ent_Entity, opts___] := iMoleculeViewer[processEntity[ent], opts]
iMoleculeViewer[f_File, opts___] := iMoleculeViewer[processChemicalFile[ExpandFileName[f]], opts]
iMoleculeViewer[u_URL, opts___] := iMoleculeViewer[processChemicalURL[URLBuild[URLParse[u]]], opts]

iMoleculeViewer[{vertexTypes_, atomPositions_}, opts___] := With[{bonds = inferBonds[vertexTypes, atomPositions]},
 iMoleculeViewer[{vertexTypes, bonds, ConstantArray["Single", Length[bonds]], atomPositions}, opts]]

iMoleculeViewer[{vertexTypes_, edgeRules_, edgeTypes_, atomPositions_}, opts : OptionsPattern[MoleculeViewer]] :=
 Block[{amlg, atomplot, atoms, bondlist, bondplot, bondtypes, cc, cent, cr, cRules, e, hiCol, hiSet, imat, ipos, leg, makebonds, mbp,
            myLegendFunction, myOrient, nbn, nei, nmb, nrms, pos, res, skel, tips, v, view, vln, zpos},

          If[! validateChem[{vertexTypes, edgeRules, edgeTypes, atomPositions}], Message[MoleculeViewer::badchem]; Return[$Failed, Block]];

          v = Length[vertexTypes]; e = Length[edgeTypes]; pos = SetPrecision[atomPositions, MachinePrecision];
          myOrient = AffineTransform[RotationMatrix[Pi/20, {1, 0, 0}].RotationMatrix[-Pi/6, {0, 0, 1}]];
          view = OptionValue[ViewPoint] /. Automatic -> 3.4 myOrient[PadRight[getNormal[pos], 3]];
          If[Last[Dimensions[pos]] == 2,
              Message[MoleculeViewer::is2d]; view = {0, -Infinity, 0}; pos = Insert[#, 0., 2] & /@ pos];

          cRules = Join[Flatten[{OptionValue[ColorRules] /. Automatic -> {}}], $colorRules, {_ -> RGBColor[211./255, 18./85, 26./51]}];

          cr = Replace[vertexTypes, cRules, {1}];
          atoms = MapThread[Sphere[#2, #1 /. $sizeRules] &, {vertexTypes, pos}];
          atomplot = Transpose[{cr, atoms}];

          hiSet = OptionValue[Highlighted];
          If[hiSet =!= None,
              If[! ListQ[hiSet], hiSet = {hiSet}];
              hiCol = RGBColor[0.5, 1., 0., 0.6];
              If[VectorQ[hiSet, MatchQ[#, _Integer | _String] &], 
                  hiSet = Thread[hiSet -> hiCol]];
              hiSet = Function[l, MapAt[If[IntegerQ[#], {#}, #] &, l, {1}]] /@
                          Replace[hiSet, x : Except[_Rule | _RuleDelayed] :> (x -> hiCol), {1}];
              hiSet = DeleteCases[hiSet /. s_String :> Flatten[Position[vertexTypes, s]],
                                             HoldPattern[{} -> _] | HoldPattern[{} :> _]];
              atomplot = Fold[MapAt[Function[{l}, Join[l, {makeTranslucent[#2[[2]]],
                                                                                 MapAt[1.2 # &, l[[-1]], -1]}]],
                                                 #1, List /@ #2[[1]]] &, atomplot, hiSet]];

          tips = OptionValue[Tooltip];
          If[MatchQ[tips, All | "Atoms"],
              atomplot = MapThread[Tooltip[#, Grid[{{Style["Atom no.", Bold], #2}, {Style["Element:", Bold], #3}, {Style["Coordinates:", Bold], #4}}, 
                                                                       Alignment -> {{Left, Left}}]] &,
                                                {atomplot, Range[v], vertexTypes, atomPositions}]];

          cent = Mean[pos];

          If[edgeRules =!= {} && edgeTypes =!= {},
              bondlist = List @@@ edgeRules;
              bondtypes = Replace[edgeTypes, {"Single" -> 1, "Double" -> 2, "Triple" -> 3, _ -> 1}, {1}];

              skel = Graph[Range[v], edgeRules, DirectedEdges -> False, GraphLayout -> None];
              vln = VertexDegree[skel];
              imat = Quiet[Check[IncidenceMatrix[skel], $Failed, IncidenceMatrix::ninc], IncidenceMatrix::ninc];
              If[! MatrixQ[imat], Message[MoleculeViewer::badchem]; Return[$Failed, Block]];

              ipos = GatherBy[imat["NonzeroPositions"], First];
              zpos = Position[Total[imat, {2}], 0]; (* find isolated atoms *)
              If[zpos =!= {}, zpos = List /@ PadRight[zpos, {Automatic, 2}]];
              ipos = #[[All, -1]] & /@ SortBy[Join[ipos, zpos], #[[1, 1]] &];

              mbp = Position[bondtypes, Except[1], 1, Heads -> False];
              nrms = ConstantArray[Null, e];
              If[mbp =!= {},
                  nei = Union[Flatten[bondlist[[Union @@ ipos[[#]]]]]] & /@ Extract[bondlist, mbp];

                  (* section of adjacency matrix of line graph *)
                  amlg = (Transpose[imat].imat - 2 IdentityMatrix[e, SparseArray])[[Flatten[mbp]]];
                  (* neighboring bonds *)
                  nbn = MapThread[Complement[#[[All, 2]], #2] &, {GatherBy[amlg["NonzeroPositions"], First], mbp}];

                  nmb = MapThread[With[{tst = getNormal[pos[[#]], True]},
                                                      If[#2 == 2 && Chop[tst[[1]], 0.1] == 0 &&
                                                          Complement[bondtypes[[#3]], {2}] === {}, {0., 0., 0.}, tst[[2]]]] &,
                                              {nei, Extract[bondtypes, mbp], nbn}];
                  If[Chop[Norm[nmb, Infinity]] == 0, nmb[[1]] = getNormal[pos[[nei[[1]]]]]];

                  nmb = NestWhile[Function[nset, MapThread[If[Norm[#, Infinity] != 0, #,
                                                                                        Normalize[Cross[Mean[#2 /. Thread[Flatten[mbp] -> nset]],
                                                                                                                 Subtract @@ pos[[Extract[bondlist, #3]]]]]] &,
                                                                                     {nset, nbn, mbp}]], nmb,
                                             ! FreeQ[Chop[Map[Norm[#, Infinity] &, #]], 0] &];
                  nrms = ReplacePart[nrms, Thread[mbp -> nmb]]];

              makebonds[bid_, no_, mul_] := 
              Block[{bo = no, en = pos[[bid]], r = 12 - 2 mul, d, del, mids, pr, tm},
                       del = Subtract @@ en; d = Norm[del]; del /= d;
                       mids = Range[0, 1, 1/(2 + 2 Boole[mul > 1])];
                       tm = Transpose[{1 - mids, mids}].en;

                       If[mul == 1,
                           (* single bond *)
                           Transpose[{cr[[bid]], Tube[#, r] & /@ Partition[tm, 2, 1]}],
                           (* multiple bond *)
                           If[Chop[Abs[pr = bo.del]] != 0, (* enforce orthogonality to bond *)
                               bo -= del pr; bo = Normalize[bo];
                               bo *= 2 UnitStep[Apply[Subtract, EuclideanDistance[cent, # bo] & /@ {1, -1}]] - 1];
                           bo *= d/Sqrt[8]; (* tetrahedral angle *) If[Max[vln[[bid]]] > 3, bo *= 0.5]; (* adjustment for multiple bonds on sulfur/phosphorus *)
                           Table[Transpose[{cr[[bid]], Tube[BSplineCurve[#], r] & /@ 
                                                      Partition[tm + ArrayPad[ConstantArray[RotationMatrix[2 Pi j/mul, del].bo, 3], {{1, 1}, {0, 0}}], 3, 2]}],
                                    {j, 0, mul - 1}]]];

              bondplot = MapThread[makebonds, {bondlist, nrms, bondtypes}];
              If[MatchQ[tips, All | "Bonds"],
                  bondplot = MapThread[Tooltip[#, Grid[{{Style["Bond no.", Bold], #3}, {Style["Bond type:", Bold], ToLowerCase[#2]},
                                                                             {Style["Bond length:", Bold], Row[{#4, " pm"}]}}, Alignment -> {{Left, Left}}]] &,
                                                    {bondplot, edgeTypes, Range[e], Round[Apply[EuclideanDistance, pos[[#]]] & /@ bondlist, 0.01]}]],
              (* no bonds *) bondplot = {}];

          res = Graphics3D[{atomplot, bondplot},
                                     FilterRules[Join[{opts}, Options[MoleculeViewer]] /.
                                                     {(BaseStyle -> Automatic) -> (BaseStyle -> Directive[CapForm["Butt"], Specularity[1, 50]]),
                                                       (ViewPoint -> Automatic) -> (ViewPoint -> view)}, 
                                                     Options[Graphics3D]]];

          If[MatchQ[leg = OptionValue[PlotLegends], None | False], res,

              If[$VersionNumber >= 11.,
                  atoms = Sort[DeleteDuplicates[vertexTypes], hillOrder],
                  atoms = Sort[DeleteDuplicates[vertexTypes], hillOrder[##] == 1 &]];
              cr = Replace[atoms, cRules, {1}];
              cc = Position[cr, col_ /; translucentColorQ[col] && Last[col] == 0];
              atoms = Apply[Sequence, Delete[#, cc] & /@ {cr, atoms}];

              myLegendFunction = Framed[#, Background -> GrayLevel[0.863], FrameStyle -> None, RoundingRadius -> 5] &;

              If[$VersionNumber >= 9.,
                  Switch[leg,
                             Automatic | True | "Atoms",
                             Legended[res, SwatchLegend[atoms, LabelStyle -> {Bold, 12, GrayLevel[0.3]}, LegendFunction -> myLegendFunction, LegendMarkerSize -> 16]],
                             _PointLegend | _SwatchLegend, 
                             Legended[res, Replace[leg, Automatic | "Atoms" :> atoms, 1]],
                             Placed[_PointLegend | _SwatchLegend, _], 
                             Legended[res, Replace[leg, Automatic | "Atoms" :> atoms, 2]],
                             _, Legended[res, leg]],
                  If[MatchQ[leg, Automatic | True | "Atoms"],
                      Row[{Show[res, ImageSize -> Medium],
                                Deploy[myLegendFunction[Grid[MapThread[{Style["\[FilledSquare]", #1, 16],
                                                                                                   Style[#2, Bold, 12, GrayLevel[0.3], FontFamily -> "Helvetica"]} &, {atoms}]]]]},
                              Spacer[10], Editable -> False],
                      res]]]
          ]

End[ ]

SetAttributes[{GetCACTUS, GetChemicalData, GetChemSpider, GetJMESMILES, GetPubChem, JME, LaunchJME, MoleculeViewer, RunOpenBabel}, ReadProtected];
Protect[GetCACTUS, GetChemicalData, GetChemSpider, GetJMESMILES, GetPubChem, JME, LaunchJME, MoleculeViewer, RunOpenBabel];

EndPackage[ ]