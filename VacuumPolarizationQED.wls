#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: QED Vacuum Polarization															*)

(*
	
	This software is covered by the GNU General Public License 3.
	Copyright (C) 2022 Edilson A. Reyes R.
  Copyright (C) 2022 Fabian C. Barajas R. 
	
*)

(* :Summary:	This is a Mathematica routine for the calculation of the vacuum
              polarization in QED at one-loop level. *)

(* ------------------------------------------------------------------------ *)


PrependTo[$Path,ToFileName[{"/Path","to","FeynCalc"}]];


$LoadFeynArts=True;
Get["FeynCalc.m"];


TopVP = CreateTopologies[ 1, 1->1,
                             ExcludeTopologies -> {Internal},
                             Adjacencies -> 3 ];

(** If you want to visualize the topology, you must enable this command **)
(** Paint[TopVP,ColumnsXRows->1,FieldNumbers->True] **)


excludefields = {F[1, _],F[2,{2}],F[2,{3}], F[3, _], F[4, _],
                 S[1],S[2], S[3],
                 V[2], V[3], V[5], 
                 U[1 | 2 | 3 | 4 | 5 ]};


FDVP = InsertFields[TopVP, V[1]->V[1], 
                           Model->"SM",
                           InsertionLevel -> {Particles},
                           ExcludeParticles->excludefields]

(** If you want to visualize the diagram, you must enable this command **)
(** Paint[FDVP,ColumnsXRows->1] **)


AmpFAVP = CreateFeynAmp[FDVP, Truncated->True];


AmpFCVP = FCFAConvert[AmpFAVP,
                      ChangeDimension->D,
                      IncomingMomenta->{p},
                      OutgoingMomenta->{p},
                      LoopMomenta->{k},                              
                      UndoChiralSplittings->True];


Changes = {FCGV[a_]:>a};

(** Denominator of the VP amplitude **)
Den = Flatten@@Cases[AmpFCVP,
                     FeynAmpDenominator[__],Infinity] /. Changes // FCE ;

(** Numerator of the VP amplitude**)
Num = Block[{res},
            SetOptions[DiracTrace, DiracTraceEvaluate->True,
                                   FeynCalcExternal->True];
            res = DeleteCases[AmpFCVP,FeynAmpDenominator[__],Infinity];
            res = res /. Changes];

(** VP amplitude at one-loop level in QED **)
Amp = Flatten@@Times[Num,Den];

AmpOL = OneLoopSimplify[Amp,k,FeynCalcExternal->True];

res = OneLoop[k,Amp,Dimension->D,
                    FeynCalcExternal->True] /. D-> 4 // Simplify ;


Changes2 = {FV[a_,b_]:> a[b],
            Lor1->u,Lor2->v,
            SP[a_,a_]:>a^2,
            "ME"->m,"EL"^2->4*Pi*\[Alpha],
            MT[a_,b_]:>g[a,b]};


resmod = res //. Changes2;
(** The final result is written in terms of Passarino-Veltmann functions. **)

Quit[];
