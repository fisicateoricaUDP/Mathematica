#!/usr/bin/env wolframscript

(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: One-Loop Electron Self-Energy in QED.														*)

(*
	
	This software is covered by the GNU General Public License 3.
	Copyright (C) 2022 Edilson A. Reyes R.
        Copyright (C) 2022 Cristian A. Lobo B. 
	
*)

(* :Summary:  This is a Mathematica routine for the calculation of the electron
              self-energy in QED at one-loop level. *)

(* ------------------------------------------------------------------------ *)

PrependTo[$Path,ToFileName[{"/home","edilson","bin","Feyncalc931","FeynCalc"}]];


$LoadFeynArts=True;


Get["FeynCalc.m"]


TopSE=CreateTopologies[1,1->1,ExcludeTopologies->{Internal},Adjacencies->3 ]


Paint[TopSE,ColumnsXRows->1,FieldNumbers->True]


excludefields={F[1, _],F[2,{2}],F[2,{3}], F[3, _], F[4, _],
               S[1],S[2], S[3],
               V[2], V[3], V[5], 
               U[1 | 2 | 3 | 4 | 5 ]};


FDSE=InsertFields[TopSE, F[2,{1}]->F[2,{1}], 
                         Model->"SM",
                         InsertionLevel -> {Particles},
                         ExcludeParticles->excludefields]


Paint[FDSE,ColumnsXRows->1]


AmpFASE=CreateFeynAmp[FDSE,Truncated->True] 


AmpFCSE=FCFAConvert[AmpFASE,
                    ChangeDimension->D,
                    IncomingMomenta->{p},
                    OutgoingMomenta->{p},
                    LoopMomenta->{k},                              
                    UndoChiralSplittings->True]


Changes = {FCGV[a_]:>a,k->k+p}


Den=Flatten@@Cases[AmpFCSE,FeynAmpDenominator[__],Infinity] /. Changes // FCE 


Num=Block[{res},
          res=DeleteCases[AmpFCSE,FeynAmpDenominator[__],Infinity];
          res=res /. Changes // FCE ]


Amp=Contract[Flatten@@Times[Num,Den]] // FCE


AmpOL=OneLoopSimplify[Amp,k,FeynCalcExternal->True]


res=OneLoop[k,Amp,Dimension->D,FeynCalcExternal->True] /. D-> 4 // Simplify 


Changes2 = {"ME" -> m,"EL"^2->4*Pi*\[Alpha],SP[a_, a_]->a^2}


resmod = res /. Changes2; 
resmod = Collect[resmod,GS[p]]
(** The final result is written in terms of Passarino-Veltmann functions. **)

Print["The End !"];

Quit[];
