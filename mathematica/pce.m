(* ::Package:: *)

(* ::Section:: *)
(*Component erasing maps package*)


(* ::Input::Initialization:: *)
(*Author: Jos\[EAcute] Alfredo de Le\[OAcute]n*)
(*Date: August 07, 2020*)
BeginPackage["pce`"]

Reshuffle::usage=
"Reshuffle[SqMatrix] reshuffles the matrix SqMatrix."
Pauli::usage=
"Pauli[Indices_List] gives the tensor product of Pauli Matrices with indices in Indices_List."
PCESuperoperator::usage=
"PCh[diagElements, qubitsNum] calculates the matrix representation of a map in the tensor product of Pauli matrices given the
diagonal elements of the matrix in computational basis."
CubePositions::usage=
"CubePositions[diagElPos] gives the positions of the painted cubes given the positions of the 1's in the diagonal of the quantum channel."
Cube3q::usage=
"Cube3q[pauliDiagonal] graphs the 3-qubit board given the indices of the painted cubes."
CPtest::usage=
"CPtest[points] returns True if CP, False if not."
PtestM::usage=
"Ptest[A] evaluates the positive-semidefiniteness of A using the principal minors criterion."
PositivityTest::usage=
"PositivityTest[\!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\)] evaluates the positive-semidefiniteness of A with its eigenvalues."
Dirac::usage=
"Dirac[vector] returns vector in Dirac notation in computational basis."
TwoQBoard::usage=
"TwoQBoard[diagonalPCE] returns a two qubits board."
Erase::usage=
"Erase[EigInfo,PCEfrom,invariantComponents] erases in all valid forms."
BlochSphTransformation::usage="Returns Bloch Ball transformation of a 1-qubit quantum chanenl given a
list with center point and the factor of x,y and z. 
Example: BlochSphTransformation[{{0,0,0.3},{1,1/2,1/2}}] returns 
a taco-like figure with center at (0,0,0.3).
BlochSphTransformation[coord_List]
"
PCEgenerators::usage=
"PCEgenerators[n] returns the diagonal superoperators of the generators
of PCE quantum channels of n qubits."
tensorPower::usage=
"."

Begin["`Private`"]
Reshuffle[SqMatrix_]:=ArrayFlatten[ArrayFlatten/@Partition[Partition[ArrayReshape[#,{Sqrt[Dimensions[SqMatrix][[1]]],Sqrt[Dimensions[SqMatrix][[1]]]}]&/@SqMatrix,Sqrt[Dimensions[SqMatrix][[1]]]],Sqrt[Dimensions[SqMatrix][[1]]]],1];

Pauli[0]=Pauli[{0}]={{1,0},{0,1}};
Pauli[1]=Pauli[{1}]={{0,1},{1,0}};
Pauli[2]=Pauli[{2}]={{0,-I},{I,0}};
Pauli[3]=Pauli[{3}]={{1,0},{0,-1}};
Pauli[Indices_List]:=KroneckerProduct@@(Pauli/@Indices);

PCESuperoperator[pauliDiagonal_List]:=Module[{indices,n,pauliToComp},
n=Log[4,Length[pauliDiagonal]];
indices=Tuples[Range[0,3],n];
pauliToComp=Transpose[Map[Flatten[Pauli[#]]&,indices]];
pauliToComp.DiagonalMatrix[pauliDiagonal].Inverse[pauliToComp]
]

CubePositions[diagElPos_]:=
Position[ArrayReshape[SparseArray[diagElPos->ConstantArray[1,Length[diagElPos]],{64}]//Normal,{4,4,4}],1]-1

Cube3q[pauliDiagonal_]:=Module[{cubeIndices},
cubeIndices=Position[ArrayReshape[pauliDiagonal,{4,4,4}],1]-1;
Graphics3D[{If[Count[#,0]==3,{Black,Cube[#]},
If[Count[#,0]==2,{RGBColor["#CC0000"],Cube[#]},
If[Count[#,0]==1,{RGBColor["#004C99"],Cube[#]},
If[Count[#,0]==0,{RGBColor["#99FF33"],Cube[#]}]]]]&/@cubeIndices,
{Thickness[0.012],Line[{{{-0.5,-0.5,-0.5},{-0.5,-0.5,3.5}},{{-0.5,-0.5,-0.5},{-0.5,3.5,-0.5}},{{-0.5,-0.5,-0.5},{3.5,-0.5,-0.5}},
{{3.5,-0.5,-0.5},{3.5,-0.5,3.5}},
{{-0.5,-0.5,3.5},{3.5,-0.5,3.5}},
{{-0.5,3.5,-0.5},{3.5,3.5,-0.5}},
{{3.5,3.5,-0.5},{3.5,3.5,3.5}},
{{3.5,3.5,3.5},{-0.5,3.5,3.5}},
{{-0.5,3.5,3.5},{-0.5,3.5,-0.5}},
{{-0.5,3.5,3.5},{-0.5,-0.5,3.5}},
{{3.5,3.5,3.5},{3.5,-0.5,3.5}},
{{3.5,3.5,-0.5},{3.5,-0.5,-0.5}}}]}},
Axes->False,AxesLabel->{"x","y","z"},LabelStyle->Directive[Bold,Medium,Black],PlotRange->{{-0.5,3.5},{-0.5,3.5},{-0.5,3.5}},AxesOrigin->{0.5,0.5,0.5},AxesStyle->Thickness[0.005],ImageSize->Medium,ImagePadding->45]
]

CPtest[points_]:=If[(PCE[SparseArray[points+1->ConstantArray[1,{points//Length}],{4,4,4}]//Normal//Flatten,3]//Reshuffle//Eigenvalues//Min)>=0,True,False]

PtestM[A_]:=AllTrue[(Diagonal[Map[Reverse,Minors[A,#],{0,1}]]&/@Range[Length[A]]),#>=0&,2]

PositivityTest[A_]:=Min[Eigenvalues[A]]>=0

Dirac[vector_List]:=(vector[[#]]Ket[IntegerString[(#-1),2,Log[2,Length[vector]]]])&/@Delete[Range[Length[vector]],Position[vector,0]]//Total

TwoQBoard[diagonalPCE_List]:=ArrayPlot[SparseArray[Position[ArrayReshape[diagonalPCE,{4,4}],1]->(If[#[[1]]==1\[And]#[[2]]==1,Black,If[#[[1]]==1\[Or]#[[2]]==1,RGBColor["#CC0000"],If[#[[1]]!=1\[And]#[[2]]!=1,RGBColor["#004C99"],Nothing]]]&/@Position[ArrayReshape[diagonalPCE,{4,4}],1]),{4,4}]]

Erase[EigInfo_,PCEfrom_,invariantComponents_]:=Module[{dimPCE},
dimPCE=Length[Dimensions[PCEfrom][[2;;]]];
If[Count[#//Flatten,1]==invariantComponents,#,Nothing]&/@
DeleteDuplicates[
Flatten[
Table[
DeleteDuplicates[ReplacePart[PCEfrom[[i]]+#-ConstantArray[1,ConstantArray[4,dimPCE]],#->0&/@Position[PCEfrom[[i]]+#-ConstantArray[1,ConstantArray[4,dimPCE]],-1]]&/@EigInfo]
,{i,Length[PCEfrom]}]
,1]
]
]

BlochSphTransformation[coord_]:=Module[{x0,y0,z0,a,b,c},
{x0,y0,z0}=coord[[1]];
{a,b,c}=If[#==0,0.02,#]&/@coord[[2]];
Style[Show[ContourPlot3D[x^2+y^2+z^2==1,{x,-1,1},{y,-1,1},{z,-1,1},
ContourStyle->{Yellow,Opacity[0.25]},Mesh->None],
ContourPlot3D[(x-x0)^2/a^2+(y-y0)^2/b^2+(z-z0)^2/c^2==1,{x,-1,1},{y,-1,1},{z,-1,1},
ContourStyle->{Dashed,Pink,Opacity[0.65]},Mesh->None],
Graphics3D[{
        Black, Arrow[Tube[{{0,0,0},1.3*Normalize[{1,0,0}]}],0.05],
        Black,Arrow[Tube[{{0,0,0},1.3*Normalize[{0,1,0}]}],0.05],
        Black,Arrow[Tube[{{0,0,0},1.3*Normalize[{0,0,1}]}],0.05],
Text["x",{1.3,0,0}],Text["y",{0,1.3,0}],Text["z",{0,0,1.3}] },
Boxed->False],Boxed->False,Axes->False,PlotRange->1.3],
RenderingOptions->{"3DRenderingMethod"->"HardwareDepthPeeling"}]
]

tensorPower[A_, n_] := Nest[KroneckerProduct[A, #] &, A, n - 1]

PCEgenerators[n_]:=Module[{a},
a={{1,1,1,1},{1,1,-1,-1},{1,-1,1,-1},{1,-1,-1,1}};
DiagonalMatrix[#]&/@ReplacePart[tensorPower[a,n],Position[tensorPower[a,n],-1]->0][[2;;]]
]
End[];
EndPackage[]


(* ::InheritFromParent:: *)
(*"PauliToComp[qbitsNum] calculates the change of basis matrix for a qubit-system map from computational to tensor product of Pauli matrices basis."*)
