(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: One-Loop W boson Self-Energy in the SM.													*)

(*
	
	This software is covered by the GNU General Public License 3.
	Copyright (C) 2023 Edilson A. Reyes R.
        Copyright (C) 2023 Aldair A. Ramon B.  
	
*)

(* :Summary:  This is a Mathematica routine for the calculation of a W boson
              self-energy mediated by a Higgs boson in the SM at 
              one-loop level.                                               *)

(* ------------------------------------------------------------------------ *)


$LoadAddOns = {"FeynArts"};
Get["FeynCalc`"];

TopW = CreateTopologies[1,1->1,ExcludeTopologies->{Internal},Adjacencies->3 ];
Paint[TopW,ColumnsXRows->1,FieldNumbers->True];

excludefields={F[1, _],F[2,_], F[3, _], F[4, _],
               S[2], S[3],
               V[1], V[2], V[5], 
               U[1 | 2 | 3 | 4 | 5 ]};
               
FDW=InsertFields[TopW, V[3]->V[3], 
                         Model->"SM",
                         InsertionLevel -> {Particles},
                         ExcludeParticles->excludefields];
                         
Paint[FDW,ColumnsXRows->1];

AmpFAW=CreateFeynAmp[FDW,Truncated->True];

AmpFCW=FCFAConvert[AmpFAW,
                    ChangeDimension->D,
                    IncomingMomenta->{p},
                    OutgoingMomenta->{p},
                    LoopMomenta->{k},                              
                    UndoChiralSplittings->True] // Contract;
                    
Changes = {FCGV[b_] :> b, Lor1->\[Mu],Lor2->\[Nu]};
AmpFCW = AmpFCW /. Changes ;

res=OneLoop[k,AmpFCW[[1]],Dimension->D,FeynCalcExternal->True] /. D-> 4 // Simplify;

Quit[];
