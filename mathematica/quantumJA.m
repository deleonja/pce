(* ::Package:: *)

(* ::Section:: *)
(*Component erasing maps package*)


(* ::Input::Initialization:: *)
(*Author: Jos\[EAcute] Alfredo de Le\[OAcute]n*)
(*Date: August 07, 2020*)
BeginPackage["quantumJA`"]

Reshuffle::usage=
"Reshuffle[SqMatrix] reshuffles the matrix SqMatrix."
Pauli::usage=
"Pauli[Indices_List] gives the tensor product of Pauli Matrices with indices in Indices_List."
ChangeOfBasisMatrix::usage=
"ChangeOfBasisMatrix[qbitsNum] calculates the change of basis matrix for a qubit-system map from computational to tensor product of Pauli matrices basis."
QC::usage=
"QC[diagElements, qubitsNum] calculates the matrix representation of a map in the tensor product of Pauli matrices given the
diagonal elements of the matrix in computational basis."
CubePositions::usage=
"CubePositions[diagElPos] gives the positions of the painted cubes given the positions of the 1's in the diagonal of the quantum channel."
Cube3q::usage=
"Cube3q[indices] graphs the 3-qubit board given the indices of the painted cubes."

Begin["`Private`"]
Reshuffle[SqMatrix_]:=ArrayFlatten[ArrayFlatten/@Partition[Partition[ArrayReshape[#,{Sqrt[Dimensions[SqMatrix][[1]]],Sqrt[Dimensions[SqMatrix][[1]]]}]&/@SqMatrix,Sqrt[Dimensions[SqMatrix][[1]]]],Sqrt[Dimensions[SqMatrix][[1]]]],1];

Pauli[0]=Pauli[{0}]={{1,0},{0,1}};
Pauli[1]=Pauli[{1}]={{0,1},{1,0}};
Pauli[2]=Pauli[{2}]={{0,-I},{I,0}};
Pauli[3]=Pauli[{3}]={{1,0},{0,-1}};
Pauli[Indices_List]:=KroneckerProduct@@(Pauli/@Indices);

ChangeOfBasisMatrix[qbitsNum_]:=
Flatten/@(Pauli[#]&/@Tuples[Range[0,3],qbitsNum])//Transpose

QC[diagElements_, qubitsNum_]:=
ChangeOfBasisMatrix[qubitsNum].
DiagonalMatrix[diagElements].
Inverse[ChangeOfBasisMatrix[qubitsNum]]

CubePositions[diagElPos_]:=
Position[ArrayReshape[SparseArray[diagElPos->ConstantArray[1,Length[diagElPos]],{64}]//Normal,{4,4,4}],1]-1

Cube3q[indices_]:=Graphics3D[{Cube[#]&/@indices,{Opacity[0],Cube[{2,2,2},3]}},Axes->True,AxesLabel->{"x","y","z"},LabelStyle->Directive[Bold, Large,Black],PlotRange->{{-0.5,3.5},{-0.5,3.5},{-0.5,3.5}},AxesOrigin->{0.5,0.5,0.5},AxesStyle->Thickness[0.01],ImageSize->Medium]

End[];
EndPackage[]
