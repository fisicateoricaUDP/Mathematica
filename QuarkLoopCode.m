#!/home/danielm/Programas -script

(* ::Package:: *)

(* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: Short distance constraints from HLbL contribution to the muon 
           anomalous magnetic moment.                                      *)

(*
	
	This software is covered by the GNU General Public License 3.
	Copyright (C) 2023 Daniel G. Melo P. 
	Copyright (C) 2023 Edilson A. Reyes R.
	Copyright (C) 2023 Angelo R. Fazio 
        
	
*)

(* :Summary:  This code was writted by Daniel Melo as part of his Master dissertation 
              under the supervision of Prof. Raffaele Fazio and Prof. Edilson Reyes. *)

(* ----------------------------------------------------------------------------------- *)


$FeynCalcStartupMessages=False;
Get["FeynCalc.m",Path -> "/Path/to/FeynCalc"];

Off[ParrallelCombine::nopar1];

SetOptions[DiracSimplify, DiracSubstitute67 -> True, Expanding -> True, Factoring -> False];
SetOptions[DiracTrace, DiracTraceEvaluate -> False];
SetOptions[Simplify, TimeConstraint -> 300];

LaunchKernels[50];

ParallelEvaluate[
 		 $FeynCalcStartupMessages=False;
 		 Get["FeynCalc.m",Path -> "/Path/to/FeynCalc"];
 		 Off[ParrallelCombine::nopar1];
 		 SetOptions[DiracSimplify, DiracSubstitute67 -> True, Expanding -> True, Factoring -> False];
 		 SetOptions[DiracTrace, DiracTraceEvaluate -> False];
 		 SetOptions[Simplify, TimeConstraint -> 300];
 		 ];

a = {a1, a2, a3, a4, a5};

q = {q1, q2, q3, q4};

Q = {Q1, Q2, Q3};

c = {c1, c2, c3, c4};

Den = {D1, D2, D3};

ChangeS = {S[x_, m_] :> (DiracGamma[Momentum[x, D], D] + m)*
	   FeynAmpDenominator[PropagatorDenominator[Momentum[x, D], m]]};

H[i1_, i2_, i3_, i4_, m_] :=(
  DiracTrace[
	     DiracMatrix[a[[i3]], Dimension -> D] .
	     S[p + q[[i1]] + q[[i2]] + q[[i4]], m] .
	     DiracMatrix[a[[i4]], Dimension -> D] .
	     S[p + q[[i1]] + q[[i2]], m] .
	     DiracMatrix[a[[i1]], Dimension -> D] . S[p + q[[i2]], m] .
	     DiracMatrix[a[[i2]], Dimension -> D] . S[p, m]]);

Integrando[x_] :=Flatten @@ Cases[{x}, DiracTrace[int__] :> int, Infinity];

Integrando[H[1, 2, 3, 4, m]];

DFP[exp_, x_, li_] :=  D[exp, x] /. {Derivative[1, 0][a_][b_,m_] :> -a[b, m] . DiracMatrix[li, Dimension -> D] . a[b, m]};
  
fun1[n_] := Pair[Momentum[p, D], LorentzIndex[a[[n]], D]];
fun2[x_, y_] :=Pair[Momentum[p, D], LorentzIndex[a[[x]], D]] Pair[Momentum[p, D],LorentzIndex[a[[y]], D]];
fun3[x_, y_, z_] :=Pair[Momentum[p, D], LorentzIndex[a[[x]], D]] Pair[Momentum[p, D],LorentzIndex[a[[y]], D]] Pair[Momentum[p, D],LorentzIndex[a[[z]], D]];
fun4[x_, y_, z_, w_] :=(Pair[Momentum[p, D], LorentzIndex[a[[x]], D]] Pair[Momentum[p, D],LorentzIndex[a[[y]], D]] Pair[Momentum[p, D],LorentzIndex[a[[z]], D]]*
			Pair[Momentum[p, D],LorentzIndex[a[[w]], D]]);
fun5[x_, y_, z_, w_, v_] :=(Pair[Momentum[p, D], LorentzIndex[a[[x]], D]] Pair[Momentum[p, D],LorentzIndex[a[[y]], D]] Pair[Momentum[p, D],LorentzIndex[a[[z]], D]]*
			    Pair[Momentum[p, D],LorentzIndex[a[[w]], D]] Pair[Momentum[p, D],LorentzIndex[a[[v]], D]]);

IntegralList =( Flatten[Table[Int2[{b1, b2, b3, b4, b5}, {e1, e2, e3, e4}, {d1, d2, d3}] , {b1,0, 1}, {b2, 0, 1}, {b3, 0, 1}, {b4, 0, 1}, {b5, 0, 1}, {e1, 0,1},
    {e2, 0, 1}, {e3, 0, 1}, {e4, 0, 1}, {d1, 0, 3}, {d2, 0,3}, {d3, 0,3}]]);

(*Contiene todas las posibles integrales con cinco índices libres, cuatro índices contraídos y 3 denominadores de propagador de 
 cada uno. Los primeros 5 cinco números indican los posibles cinco índices libres de Lorentz con que puede aparecer el momento de loop. 
 Los siguientes 4 números se refieren a índices contraídos. Los últimos tres números indican la potencia con que aparece cada	
 propagador con su respectivo momento externo.*)
  
Int[{b1_, b2_, b3_, b4_, b5_}, {e1_, e2_, e3_, e4_}, {d1_, d2_,d3_}] := (
									 Pair[Momentum[p, D], LorentzIndex[a[[1]], D]]^b1*
									 Pair[Momentum[p, D], LorentzIndex[a[[2]], D]]^b2*
									 Pair[Momentum[p, D], LorentzIndex[a[[3]], D]]^b3*
									 Pair[Momentum[p, D], LorentzIndex[a[[4]], D]]^b4*
									 Pair[Momentum[p, D], LorentzIndex[a[[5]], D]]^b5*
									 Pair[Momentum[p, D], LorentzIndex[c[[1]], D]]^e1*
									 Pair[Momentum[p, D], LorentzIndex[c[[2]], D]]^e2*
									 Pair[Momentum[p, D], LorentzIndex[c[[3]], D]]^e3*
									 Pair[Momentum[p, D], LorentzIndex[c[[4]], D]]^e4*
									 (1/Den[[1]]^d1) (1/Den[[2]]^d2)*(1/Den[[3]]^d3) );
  
(*Da como resultado el integrando de integrales de loop donde la variable es el momento p. El significado de las variables de entrada se encuentra en el objeto
 IntegralList*)

EstructuraTensorial[{b1_, b2_, b3_, b4_, b5_}, {e1_, e2_, e3_,e4_}, {k_, k1_, k2_, k3_}] :=
  (M = b1 + b2 + b3 + b4 + b5 + e1 + e2 + e3 + e4;
   If[M == 0, Return[1]];
   ind = {a[[1]]*b1, a[[2]]*b2, a[[3]]*b3, a[[4]]*b4, a[[5]]*b5,c[[1]]*e1, c[[2]]*e2, c[[3]]*e3,
	  c[[4]]*e4};
   ind = DeleteCases[ind, 0];
   vector =Table[Pair[Momentum[Q[[i]], D], LorentzIndex[ind[[j]], D]], {i, 3}, {j, M}];
   g[i1_, i2_] :=Pair[LorentzIndex[ind[[i1]], D], LorentzIndex[ind[[i2]], D]];
   tensor =Product[vector[[1]][[i]], {i, k1}] Product[vector[[2]][[i]], {i, k1 + 1, k1 + k2}] Product[vector[[3]][[i]], {i, k1 + k2 + 1, k1 + k2 + k3}]
   Product[g[i, i + 1], {i, k1 + k2 + k3 + 1, k1 + k2 + k3 + 2*k-1, 2}];
   Return[FCSymmetrize[tensor, ind]]
   );
(*Arroja la estructura tensorial que puede tener una integral de loop de acuerdo con la descomposición de Davydychev. El signficado de las primeras nueve variables
 de entrada es el mismo descrito en IntegralList. Las últimas cuatro variables indican el número de veces que está permitido que salgan cada uno de los tres momentos
 externos y la métrica, respectivamente.*)

Int3[{b1_, b2_, b3_, b4_, b5_}, {e1_, e2_, e3_, e4_}, {d1_, d2_,d3_}] :=
  (M = b1 + b2 + b3 + b4 + b5 + e1 + e2 + e3 + e4;
   temp = Flatten[Table[If[M - 2 k - k1 - k2 >= 0, (-1/2)^k *Pochhammer[d1, k1]*Pochhammer[d2, k2]*Pochhammer[d3, M - 2 k - k1 - k2]*(Pi)^(k-M)*
			   EstructuraTensorial[{b1, b2, b3, b4, b5}, {e1, e2, e3, e4}, {k, k1, k2, M - 2 k - k1 - k2}]*
			   IntSca[3 - Length[Cases[{d1, d2, d3}, 0]], 4 + 2 (M - k), d1 + k1, d2 + k2, d3 + M - 2 k - k1 - k2], 0],
      {k, 0, Floor[M/2]}, {k1, 0, M}, {k2, 0, M}]];
   Return[Sum[temp[[i]], {i, Length[temp]}]]
   );

(*Realiza la descomposición de una integral de loop de acuerdo con	\
 Davydychev. Incluye descomposición tensorial y shift dimensional.	\
 Deja las integrales escalares expresadas en términos de una función	\
 auxiliar IntSca. Las variables de entrada de IntSca funcionan de la	\
 siguiente manera: la primera es la dimensión de la integral y las	\
 demás son la potencia con la que aparece cada tipo de denominador de	\
 propagador escalar.*)

MomExt[amp_] := DeleteDuplicates[Cases[amp /. q[[4]] -> 0 /. p -> p - zero, S[p + x_, m] :>  x, 4] /. zero -> 0];

(*Arroja una lista con los momentos externos (en la notación de Davydychev) de una amplitud en el límite q4-> 0. Identifica incluso cuando el alguno de esos momentos
es cero. SE DEBE INGRESAR LA AMPLITUD ANTES DE REEMPLAZAR LA FUNCIÓN AUXILIAR S[]*)
    
CondEntrada[amp_] := Table[MomExt[amp][[i]] -> -Q[[i]], {i, 3}];

(*Arroja una lista de condiciones que se deben aplicar a la amplitud de entrada para convertirla a la notación de Davydychev.*)

CondSalida[amp_] :=Table[Q[[i]] -> -MomExt[amp][[i]], {i, 3}];

(*Arroja una lista de condiciones que se deben aplicar sobre la amplitud de entrada para devolverla a ella o sus derivados a la
 notación original, cuando se ha aplicado CondEntrada previamente.*)


(*Realiza la traza y el álgebra necesarios para identificar qué integrales de loop son realmente necesarias.*)

TrazaYAlgebra[amp_] :=(
  Block[{AmpNum, AmpDen, AmpDS, AmpDT, AmpDT1, AmpDT2, AmpDT3, AmpDT4, Den, TempDen, c, Externos, Entrada, Salida},
	Externos = MomExt[amp];
	AmpDS = DiracSimplify[amp /. ChangeS];
	AmpDT = DiracTrace[AmpDS, DiracTraceEvaluate -> True];
	Den = {D1, D2, D3};
	TempDen = {TempD1, TempD2, TempD3};
	c = {c1, c2, c3, c4};
	
	ChangeDeno = Union[Table[FeynAmpDenominator[PropagatorDenominator[Momentum[p - Q[[i]], D], m]] ->1/Den[[i]], {i, 3}],
			   Table[FeynAmpDenominator[PropagatorDenominator[Momentum[p - Externos[[i]], D], m]] ->1/Den[[i]], {i, 3}]];
	TempChangeDeno = Flatten[Table[1/Den[[i]]^j -> 1/TempDen[[i]]^j, {i, Length[Den]}, {j, 5}]];
	TempUnchangeDeno = Flatten[Table[1/TempDen[[i]]^j -> 1/Den[[i]]^j, {i, Length[Den]}, {j, 5}]];
	ChangeNum = Union[Table[Pair[Momentum[p - Q[[i]], D], Momentum[p - Q[[i]], D]] -> Den[[i]] + m^2, {i, 3}],
			  Table[Pair[Momentum[p - Externos[[i]], D], Momentum[p - Externos[[i]], D]] -> Den[[i]] + m^2, {i, 3}]];
	UnchangeNum = DeleteCases[Table[If[Externos[[i]] == 0, Den[[i]] -> Pair[Momentum[p, D], Momentum[p, D]] - m^2, 0, 0], {i, Length[Externos]}], 0];
	ExpandProduct = Table[Pair[Momentum[p, D], Momentum[Q[[i]], D]] -> (1/2)*(Pair[Momentum[p - Q[[i]], D], Momentum[p - Q[[i]], D]] -
										  Pair[Momentum[p, D], Momentum[p, D]] -
										  Pair[Momentum[Q[[i]], D], Momentum[Q[[i]], D]]), {i, 3}];
	ExpandContraction = {Pair[Momentum[p, D], Momentum[p, D]]^2 -> Pair[Momentum[p, D], LorentzIndex[c[[1]], D]] Pair[Momentum[p, D], LorentzIndex[c[[2]], D]]*
			     Pair[Momentum[p, D], LorentzIndex[c[[3]], D]] Pair[Momentum[p, D], LorentzIndex[c[[4]], D]], Pair[Momentum[p, D], Momentum[p, D]] ->
			     Pair[Momentum[p, D], LorentzIndex[c[[1]], D]] Pair[Momentum[p, D], LorentzIndex[c[[2]], D]]};
	AmpDT1 =Collect[AmpDT /. ExpandProduct, {fun5[1, 2, 3, 4, 5], fun4[1, 2, 3, 4], fun4[1, 2, 3, 5], fun4[1, 2, 5, 4], fun4[1, 5, 3, 4], fun4[5, 2, 3, 4],
						 fun3[1, 2, 3], fun3[1, 2, 4], fun3[1, 4, 3], fun3[4, 2, 3], fun3[1, 2, 5], fun3[1, 5, 3], fun3[5, 2, 3],
						 fun3[1, 4, 5], fun3[2, 4, 5], fun3[3, 4, 5], fun2[1, 2], fun2[1, 3], fun2[1, 4], fun2[1, 5], fun2[2, 4], fun2[2, 3],
						 fun2[2, 5], fun2[3, 4], fun2[3, 5], fun1[1], fun1[2], fun1[3], fun1[4], fun1[5]}];
	AmpDT2 = AmpDT1 /. ChangeNum;
	AmpDT3 = Expand[Expand[AmpDT2 (*]] Este símbolo es para que Mathematica no se confunda*) //. ChangeDeno] //. TempChangeDeno /. UnchangeNum /. TempUnchangeDeno] /. ExpandContraction;
	Return[AmpDT3];
	]);

(*amp = Integrando[H[1, 4, 3, 2, m]];
 amp2 = DFP[amp, q[[4]], a[[5]]] /. q[[4]] -> 0 /. CondEntrada[amp];
 res = TrazaYAlgebra[amp2];
 Coefficient[res, Pair[LorentzIndex[a1, D], Momentum[p, D]]]; Esto hace algunos checks básicos sobre las funciones previamente definidas.*)

ReplaceIntegrals = Reverse[ Flatten[ Table[ Int[{b1, b2, b3, b4, b5}, {e1, e2, e3, e4}, {d1, d2, d3}] -> Int2[{b1, b2, b3, b4, b5}, {e1, e2, e3, e4}, {d1, d2, d3}],
  {b1, 0, 1}, {b2, 0, 1}, {b3, 0, 1}, {b4, 0, 1}, {b5, 0, 1}, {e1, 0, 1}, {e2, 0, 1}, {e3, 0, 1}, {e4, 0, 1}, {d1, 0, 3}, {d2, 0, 3}, {d3, 0, 3}]]];

ActualIntegralList[amp_] := Union@Cases[amp, Int2[__], Infinity];
ScalarIntegralList[amp_] := DeleteDuplicates[Cases[ActualIntegralList[amp] /. Int2[x__] -> Int3[x], IntSca[__], Infinity]];

Combinada[mu1_, mu2_, mu3_, mu4_, nu4_] :=
  (
   amp = Integrando[H[mu1, mu2, mu3, mu4, m]];
   (*Genera la amplitud sin derivar. Los propagadores fermiónicos quedan expresados en una función auxiliar S[]*)
   
   condicionesEntrada = CondEntrada[amp];
   condicionesSalida = CondSalida[amp];
   (*Guarda las condiciones de entrada y de salida*)
   
   amp = DFP[amp, q[[4]], a[[nu4]]] /. q[[4]] -> 0 /. condicionesEntrada;
   (*Deriva la amplitud con respecto al momento q4, toma el límite q4->0, traduce a notación de Davydychev.*)
   
   amp = FCAntiSymmetrize[TrazaYAlgebra[amp], {a[[mu4]], a[[nu4]]}] /. ReplaceIntegrals;
   (*Realiza la traza y el álgebra necesaria para poder identificar las integrales de loop realmente necesarias. Luego antisimetrizo. Luego identifico las integrales y
    las reemplazo con una función auxiliar.*)
   
   amp = amp /. Int2[x__] -> Int3[x] /. condicionesSalida; (*Al final aplico la descmposición tensorial de Davydychev y traduzco de nuevo a la notación de Bijnens*)
   amp = Contract[amp /. c2 -> c1 /. c4 -> c3]; (*Contraigo los índices*)
   Return[amp];
   );


(*resultado = Combinada[2, 4, 3, 1, 5];
 Cases[resultado /. Int2[x__] -> Int3[x], Pair[LorentzIndex[_, D], Momentum[p, D]], 5] Una pruebita de la función Combinada*)

(*Cases[TodasLasAmplitudes /. Int2[x__] -> Int3[x], Den, 6]; (*Pequeña prueba de que todos los denominadores están en integrales escalares*)*)
   
(*TodasLasEscalares = DeleteDuplicates[Flatten[Table[ScalarIntegralList[TodasLasAmplitudes[[i]]], {i, Length[TodasLasAmplitudes]}]]];
 (*Antes buscaba los distintos tipos de integrales escalares que aparecían en todas las permutaciones. Ya no funciona pq en la función
  Combinada se está haciendo ya la descomposición de Davydychev*)*)


ProEsc[pot1_, pot2_, pot3_] := (Return[Pair[Momentum[q[[1]], D], Momentum[q[[1]], D]]^(pot1/2) Pair[Momentum[q[[2]], D], Momentum[q[[2]], D]]^(pot2/2)*
				       Pair[Momentum[q[[3]], D], Momentum[q[[3]], D]]^(pot3/2)]);

(*Devuelve un producto de las tres normas de q1, q2 y q3 elevadas a la potencia poti/2*)


ProTen[mom_, ind_] := (Return[Pair[Momentum[q[[mom]], D], LorentzIndex[a[[ind]], D]] ]);

(*Devuelve el cuadrimomento mom con índice ind*)

SP[mom1_, mom2_] := (Return[Pair[Momentum[q[[mom1]], D], Momentum[q[[mom2]], D]]]);

(*Una función que devuelve el producto escalar de dos momentos externos*)

g[ind1_, ind2_] := (Return[Pair[LorentzIndex[a[[ind1]], D], LorentzIndex[a[[ind2]], D]]]);

(*Devuelve la métrica con índices ind1 e ind2*)

lambda = ProEsc[0, 0, 4] + ProEsc[0, 4, 0] + ProEsc[4, 0, 0] -2 ProEsc[0, 2, 2] - 2ProEsc[2, 0, 2] -2 ProEsc[2, 2, 0];

(*Definimos lambda de Källen*)

Crossing[x_, ind1_, ind2_] := (
			       res = x /. {q[[ind1]] -> q[[ind2]], q[[ind2]] -> q[[ind1]], a[[ind1]] -> a[[ind2]], a[[ind2]] -> a[[ind1]]};
			       Return[res];
			       );

(*Implementa intercambio q_i<->q_j e ind_i<->ind_j*)

Crossing[ProTen[1, 1] ProTen[2, 1], 1, 2]

(*A continuación se presentan los proyectores que propone Bijnens.*)

(*Simétrico bajo crossing 1-2*)
Proyector1 =( -8/lambda^2*ProTen[1, 2] ProTen[1, 5] ProTen[2, 3] ProTen[2, 4] ProTen[3, 1]
	      +2/lambda*g[1, 5] g[2, 4] ProTen[2, 3]
	      -8/lambda^2*ProEsc[0, 2, 0] g[1, 5] ProTen[1, 2] ProTen[2, 3] ProTen[3, 4]
	      -4/lambda^2*(ProEsc[0, 0, 2] + ProEsc[0, 2, 0] - ProEsc[2, 0, 0]) g[1, 5] ProTen[1, 3] ProTen[2, 4] ProTen[3, 2]
	      -8/lambda^2*ProEsc[2, 0, 0] g[2, 5] ProTen[1, 3] ProTen[2, 1] ProTen[3, 4]
	      -4/lambda^2*(ProEsc[0, 0, 2] - ProEsc[0, 2, 0] + ProEsc[2, 0, 0]) g[2, 5] ProTen[1, 4] ProTen[2, 3] ProTen[3,1]);

(*Simétrico bajo crossing 1-2*)

(*Simétrico bajo crossing 1-2*)
Proyector4 =
  ( 8/lambda^4 (6 ProEsc[0, 0, 8] + 11 ProEsc[0, 2, 6] -
		29 ProEsc[0, 4, 4] + ProEsc[0, 6, 2] + 11 ProEsc[0, 8, 0] +
		11 ProEsc[2, 0, 6] + 14 ProEsc[2, 2, 4] - ProEsc[2, 4, 2] -
		44 ProEsc[2, 6, 0] - 29 ProEsc[4, 0, 4] - ProEsc[4, 2, 2] +
		66 ProEsc[4, 4, 0] + ProEsc[6, 0, 2] - 44 ProEsc[6, 2, 0] +
		11 ProEsc[8, 0, 0]) ProTen[1, 2] ProTen[1, 5] ProTen[2, 3] ProTen[2, 4] ProTen[3, 1]
    +1/lambda^2 (ProEsc[0, 0, 4] - 6 ProEsc[0, 2, 2] -ProEsc[0, 4, 0] + 2 ProEsc[2, 2, 0] - ProEsc[4, 0, 0])*
    g[1, 2] g[3, 5] ProTen[1, 4]
    +1/lambda^2 (ProEsc[0, 0, 4] - ProEsc[0, 4, 0] -6 ProEsc[2, 0, 2] + 2 ProEsc[2, 2, 0] - ProEsc[4, 0, 0])*
    g[1, 2] g[3, 5] ProTen[2, 4]
    -4/lambda^3 (ProEsc[0, 0, 6] + 3 ProEsc[0, 2, 4] -4 ProEsc[0, 4, 2] + 3 ProEsc[2, 0, 4] + 8 ProEsc[2, 2, 2]
		 -4 ProEsc[4, 0, 2])*g[1, 2] ProTen[1, 3] ProTen[1, 5] ProTen[2, 4]
    +2/lambda^2 ProEsc[0, 0, 2] (ProEsc[0, 2, 0] - ProEsc[2, 0, 0]) g[1, 3] g[2, 5] ProTen[1, 4]
    + 1/lambda^2 (ProEsc[0, 0, 4] - ProEsc[0, 4, 0] -2 ProEsc[2, 0, 2] + 2 ProEsc[2, 2, 0] - ProEsc[4, 0, 0])*
    g[1, 3] g[2, 5] ProTen[3, 4]
    +2/lambda^3*(3 ProEsc[0, 0, 6] + ProEsc[0, 2, 4] -ProEsc[0, 4, 2] - 3 ProEsc[0, 6, 0] - 3 ProEsc[2, 0, 4]
	       +4 ProEsc[2, 2, 2] + 9 ProEsc[2, 4, 0] - 3 ProEsc[4, 0, 2] -9 ProEsc[4, 2, 0] + 3 ProEsc[6, 0, 0])*
    g[1, 3] ProTen[1, 4] ProTen[3, 2] ProTen[3, 5]
    +2/lambda^2*(ProEsc[0, 0, 4] - 4 ProEsc[0, 2, 2] -ProEsc[0, 4, 0] + 2 ProEsc[2, 2, 0] - ProEsc[4, 0, 0])*
    g[1, 4] g[3, 5] ProTen[1, 2]
    -2/lambda^2 ProEsc[0, 0, 2] (ProEsc[0, 2, 0] - ProEsc[2, 0, 0]) g[1, 5] g[2, 3] ProTen[2, 4]
    +1/lambda^2 (ProEsc[0, 0, 4] - 2 ProEsc[0, 2, 2] -ProEsc[0, 4, 0] + 2 ProEsc[2, 2, 0] - ProEsc[4, 0, 0])*
    g[1, 5] g[2, 3] ProTen[3, 4]
    -10/lambda^3 (ProEsc[0, 0, 6] - ProEsc[0, 2, 4] -3 ProEsc[0, 4, 2] - ProEsc[0, 6, 0] - ProEsc[2, 0, 4] +
		  4 ProEsc[2, 2, 2] + 3 ProEsc[2, 4, 0] - ProEsc[4, 0, 2] -3 ProEsc[4, 2, 0] + ProEsc[6, 0, 0])*
    g[1, 5] ProTen[1, 2] ProTen[2, 3] ProTen[3, 4]
    -4/lambda^3 (2 ProEsc[0, 0, 6] - 9 ProEsc[0, 2, 4] -3 ProEsc[0, 4, 2] + ProEsc[2, 0, 4] + 6 ProEsc[2, 2, 2]
		 -3 ProEsc[4, 0, 2]) g[1, 5] ProTen[1, 3] ProTen[2, 4] ProTen[3, 2]
    +2/lambda^3 (3 ProEsc[0, 0, 6] - 3 ProEsc[0, 2, 4] -3 ProEsc[0, 4, 2] + 3 ProEsc[0, 6, 0] + ProEsc[2, 0, 4] +
		 +4 ProEsc[2, 2, 2] - 9 ProEsc[2, 4, 0] - ProEsc[4, 0, 2] + 9 ProEsc[4, 2, 0] - 3 ProEsc[6, 0, 0])*
    g[2, 3] ProTen[2, 1] ProTen[2, 5] ProTen[3, 4]
    +2/lambda^2 (ProEsc[0, 0, 4] - ProEsc[0, 4, 0] -4 ProEsc[2, 0, 2] + 2 ProEsc[2, 2, 0] - ProEsc[4, 0, 0])*
    g[2, 5] g[3, 4] ProTen[3, 1] -
    10/lambda^3 (ProEsc[0, 0, 6] - ProEsc[0, 2, 4] - ProEsc[0, 4, 2] +ProEsc[0, 6, 0] - ProEsc[2, 0, 4]
	       + 4 ProEsc[2, 2, 2] -3 ProEsc[2, 4, 0] - 3 ProEsc[4, 0, 2] + 3 ProEsc[4, 2, 0] -ProEsc[6, 0, 0])*
    g[2, 5] ProTen[1, 3] ProTen[2, 1] ProTen[3, 4]
    - 4/lambda^3 (2 ProEsc[0, 0, 6] + ProEsc[0, 2, 4] -3 ProEsc[0, 4, 2] - 9 ProEsc[2, 0, 4] + 6 ProEsc[2, 2, 2]
		  -3 ProEsc[4, 0, 2]) g[2, 5] ProTen[1, 4] ProTen[2, 3] ProTen[3, 1]
    + 1/lambda^3 (6 ProEsc[0, 0, 6] - 6 ProEsc[0, 2, 4] -6 ProEsc[0, 4, 2] + 6 ProEsc[0, 6, 0] - 50 ProEsc[2, 0, 4]
		+72 ProEsc[2, 2, 2] + 10 ProEsc[2, 4, 0] + 22 ProEsc[4, 0, 2] -38 ProEsc[4, 2, 0]
		+ 22 ProEsc[6, 0, 0]) g[3, 5] ProTen[1, 2] ProTen[2, 4] ProTen[3, 1]
    +1/lambda^3 (6 ProEsc[0, 0, 6] - 50 ProEsc[0, 2, 4] +22 ProEsc[0, 4, 2] + 22 ProEsc[0, 6, 0] - 6 ProEsc[2, 0, 4]
	       +72 ProEsc[2, 2, 2] - 38 ProEsc[2, 4, 0] - 6 ProEsc[4, 0, 2] +10 ProEsc[4, 2, 0] + 6 ProEsc[6, 0, 0])*
    g[3, 5] ProTen[1, 4] ProTen[2, 1] ProTen[3, 2] );

(*Este proyector sí corresponde a uno de los 19 \hat{T}*)

(*Sin simetría de crossing*)
Proyector7 =
  ( 80/lambda^4 (-2 ProEsc[0, 0, 6] - ProEsc[0, 2, 4] +ProEsc[0, 4, 2] + 2 ProEsc[0, 6, 0] + 2 ProEsc[2, 0, 4]
		 -3 ProEsc[2, 2, 2] - 6 ProEsc[2, 4, 0] + 2 ProEsc[4, 0, 2] +6 ProEsc[4, 2, 0]
		 - 2 ProEsc[6, 0, 0]) ProTen[1, 2] ProTen[1, 5] ProTen[2, 3] ProTen[2, 4] ProTen[3, 1]
    -2/lambda^2 (ProEsc[0, 0, 2] + ProEsc[0, 2, 0] - ProEsc[2, 0, 0]) g[1, 2] g[3, 5] ProTen[1, 4]
    +2/lambda^2 (2 ProEsc[0, 0, 2] - ProEsc[0, 2, 0] +ProEsc[2, 0, 0]) g[1, 2] g[3, 5] ProTen[2, 4]
    +20/lambda^3 (ProEsc[0, 0, 4] + ProEsc[0, 2, 2] -ProEsc[2, 0, 2]) g[1, 2] ProTen[1, 3] ProTen[1, 5] ProTen[2, 4]
    +2/lambda^2 ProEsc[0, 0, 2] g[1, 3] g[2, 5] ProTen[1, 4]
    -2/lambda^2 (ProEsc[0, 2, 0] - ProEsc[2, 0, 0]) g[1, 3] g[2, 5] ProTen[3, 4]
    -4/lambda^3 (3 ProEsc[0, 0, 4] + 4 ProEsc[0, 2, 2] +3 ProEsc[0, 4, 0] - 6 ProEsc[2, 0, 2] - 6 ProEsc[2, 2, 0]
		 +3 ProEsc[4, 0, 0]) g[1, 3] ProTen[1, 4] ProTen[3, 2] ProTen[3, 5]
    -4/lambda^2 (ProEsc[0, 0, 2] + ProEsc[0, 2, 0] -ProEsc[2, 0, 0]) g[1, 4] g[3, 5] ProTen[1, 2]
    -2/lambda^2 ProEsc[0, 0, 2] g[1, 5] g[2, 3] ProTen[2, 4]
    -2/lambda^2 (ProEsc[0, 0, 2] + ProEsc[0, 2, 0] - ProEsc[2, 0, 0]) g[1, 5] g[2, 3] ProTen[3, 4]
    +20/lambda^3 (ProEsc[0, 0, 4] + 2 ProEsc[0, 2, 2] + ProEsc[0, 4, 0] - 2 ProEsc[2, 0, 2] - 2 ProEsc[2, 2, 0]
		  + ProEsc[4, 0, 0]) g[1, 5] ProTen[1, 2] ProTen[2, 3] ProTen[3, 4]
    + 20/lambda^3 (ProEsc[0, 0, 4] + ProEsc[0, 2, 2] -ProEsc[2, 0, 2]) g[1, 5] ProTen[1, 3] ProTen[2, 4] ProTen[3, 2]
    - 4/lambda^3 (2 ProEsc[0, 0, 4] + ProEsc[0, 2, 2] -3 ProEsc[0, 4, 0] + ProEsc[2, 0, 2] + 6 ProEsc[2, 2, 0]
		  -3 ProEsc[4, 0, 0]) g[2, 3] ProTen[2, 1] ProTen[2, 5] ProTen[3, 4]
    + 4/lambda^2 (ProEsc[0, 0, 2] - ProEsc[0, 2, 0] +ProEsc[2, 0, 0]) g[2, 5] g[3, 4] ProTen[3, 1]
    +20/lambda^3 (ProEsc[0, 2, 2] - ProEsc[0, 4, 0] - ProEsc[2, 0, 2] +2 ProEsc[2, 2, 0] - ProEsc[4, 0, 0])*
    g[2, 5] ProTen[1, 3] ProTen[2, 1] ProTen[3, 4]
    -20/lambda^3 (ProEsc[0, 0, 4] - ProEsc[0, 2, 2] +ProEsc[2, 0, 2]) g[2, 5] ProTen[1, 4] ProTen[2, 3] ProTen[3, 1]
    + 4/lambda^3 (6 ProEsc[0, 0, 4] - 7 ProEsc[0, 2, 2] +ProEsc[0, 4, 0] + 3 ProEsc[2, 0, 2] + 8 ProEsc[2, 2, 0]
		-9 ProEsc[4, 0, 0]) g[3, 5] ProTen[1, 2] ProTen[2, 4] ProTen[3, 1]
    - 4/lambda^3 (3 ProEsc[0, 0, 4] + 4 ProEsc[0, 2, 2] -7 ProEsc[0, 4, 0] - 6 ProEsc[2, 0, 2] + 4 ProEsc[2, 2, 0]
		  +3 ProEsc[4, 0, 0]) g[3, 5] ProTen[1, 4] ProTen[2, 1] ProTen[3,2]);

(*Este proyector sí corresponde a uno de los 19 \hat{T}*)

(*Simétrico bajo crossing 1-2*)

Proyector17 =
  ( 80/lambda^3 (ProEsc[0, 0, 2] - ProEsc[0, 2, 0] -ProEsc[2, 0, 0])
    ProTen[1, 2] ProTen[1, 5] ProTen[2, 3] ProTen[2, 4] ProTen[3, 1]
    -2/lambda^2 (ProEsc[0, 0, 2] - 3 ProEsc[0, 2, 0] -ProEsc[2, 0, 0]) g[1, 2] g[3, 5] ProTen[1, 4]
    -2/lambda^2 (ProEsc[0, 0, 2] - ProEsc[0, 2, 0] -3 ProEsc[2, 0, 0]) g[1, 2] g[3, 5] ProTen[2, 4]
    -8/lambda^2 g[1, 2] ProTen[1, 3] ProTen[1, 5] ProTen[2, 4]
    -2/lambda^2 (ProEsc[0, 0, 2] + ProEsc[0, 2, 0] - ProEsc[2, 0, 0]) g[1, 3] g[2, 5] ProTen[1, 4]
    -2/lambda^2 (ProEsc[0, 0, 2] - ProEsc[0, 2, 0] - ProEsc[2, 0, 0]) g[1, 3] g[2, 5] ProTen[3, 4]
    -4/lambda^2 (ProEsc[0, 0, 2] - 2 ProEsc[0, 2, 0] -ProEsc[2, 0, 0]) g[1, 4] g[3, 5] ProTen[1, 2]
    -2/lambda^2 (ProEsc[0, 0, 2] - ProEsc[0, 2, 0] + ProEsc[2, 0, 0]) g[1, 5] g[2, 3] ProTen[2, 4]
    -2/lambda^2 (ProEsc[0, 0, 2] - ProEsc[0, 2, 0] - ProEsc[2, 0, 0]) g[1, 5] g[2, 3] ProTen[3, 4]
    -4/lambda^2 ProEsc[0, 0, 2] g[1, 5] g[2, 4] ProTen[2, 3]
    +8/lambda^3 (2 ProEsc[0, 0, 4] + ProEsc[0, 2, 2] -3 ProEsc[0, 4, 0] - 4 ProEsc[2, 0, 2] + ProEsc[2, 2, 0]
		 +2 ProEsc[4, 0, 0]) g[1, 5] ProTen[1, 2] ProTen[2, 3] ProTen[3, 4]
    +8/lambda^3 (4 ProEsc[0, 0, 4] - 3 ProEsc[0, 2, 2] -ProEsc[0, 4, 0] - 3 ProEsc[2, 0, 2] + 2 ProEsc[2, 2, 0]
		 -ProEsc[4, 0, 0])*g[1, 5] ProTen[1, 3] ProTen[2, 4] ProTen[3, 2]
    -4/lambda^2 (ProEsc[0, 0, 2] - ProEsc[0, 2, 0] -2 ProEsc[2, 0, 0]) g[2, 5] g[3, 4] ProTen[3, 1]
    +8/lambda^3 (2 ProEsc[0, 0, 4] - 4 ProEsc[0, 2, 2] +2 ProEsc[0, 4, 0] + ProEsc[2, 0, 2] + ProEsc[2, 2, 0]
		 -3 ProEsc[4, 0, 0]) g[2, 5] ProTen[1, 3] ProTen[2, 1] ProTen[3, 4]
    + 8/lambda^3 (4 ProEsc[0, 0, 4] - 3 ProEsc[0, 2, 2] -ProEsc[0, 4, 0] - 3 ProEsc[2, 0, 2] + 2 ProEsc[2, 2, 0]
		  -ProEsc[4, 0, 0])*
    g[2, 5] ProTen[1, 4] ProTen[2, 3] ProTen[3, 1]
    - 8/lambda^3 (2 ProEsc[0, 0, 4] - 4 ProEsc[0, 2, 2] +2 ProEsc[0, 4, 0] - 9 ProEsc[2, 0, 2] + 11 ProEsc[2, 2, 0]
		  +7 ProEsc[4, 0, 0]) g[3, 5] ProTen[1, 2] ProTen[2, 4] ProTen[3, 1]
    - 8/lambda^3 (2 ProEsc[0, 0, 4] - 9 ProEsc[0, 2, 2] +7 ProEsc[0, 4, 0] - 4 ProEsc[2, 0, 2] + 11 ProEsc[2, 2, 0]
		  +2 ProEsc[4, 0, 0]) g[3, 5] ProTen[1, 4] ProTen[2, 1] ProTen[3, 2] );


(*Simétrico bajo crossing 1-2, 1-3, 2-3*)

Proyector39 =
  -( 160/lambda^4 (ProEsc[0, 0, 6] - ProEsc[0, 2, 4] - ProEsc[0, 4, 2] +ProEsc[0, 6, 0] - ProEsc[2, 0, 4]
		   + ProEsc[2, 2, 2] -ProEsc[2, 4, 0] - ProEsc[4, 0, 2] - ProEsc[4, 2, 0]
		   +ProEsc[6, 0, 0])*
     ProTen[1, 2] ProTen[1, 5] ProTen[2, 3] ProTen[2, 4] ProTen[3, 1]
     -2/lambda^2 ProEsc[0, 2, 0] g[1, 2] g[3, 5] ProTen[1, 4]
     -2/lambda^2 ProEsc[2, 0, 0] g[1, 2] g[3, 5] ProTen[2, 4]
     -4/lambda^3 (2 ProEsc[0, 0, 4] + ProEsc[0, 2, 2] -3 ProEsc[0, 4, 0] + ProEsc[2, 0, 2] + 6 ProEsc[2, 2, 0]
		  -3 ProEsc[4, 0, 0]) g[1, 2] ProTen[1, 3] ProTen[1, 5] ProTen[2, 4]
     - 2/lambda^2 ProEsc[0, 0, 2] g[1, 3] g[2, 5] ProTen[1, 4]
     -2/lambda^2 ProEsc[2, 0, 0] g[1, 3] g[2, 5] ProTen[3, 4]
     +4/lambda^3 (3 ProEsc[0, 0, 4] - ProEsc[0, 2, 2] -2 ProEsc[0, 4, 0] - 6 ProEsc[2, 0, 2] - ProEsc[2, 2, 0]
	       +3 ProEsc[4, 0, 0]) g[1, 3] ProTen[1, 4] ProTen[3, 2] ProTen[3, 5]
     - 4/lambda^2 ProEsc[0, 2, 0] g[1, 4] g[3, 5] ProTen[1, 2]
     -2/lambda^2 ProEsc[0, 0, 2] g[1, 5] g[2, 3] ProTen[2, 4]
     -2/lambda^2 ProEsc[0, 2, 0] g[1, 5] g[2, 3] ProTen[3, 4]
     -4/lambda^2 ProEsc[0, 0, 2] g[1, 5] g[2, 4] ProTen[2, 3]
     -4/lambda^3 (ProEsc[0, 0, 4] - 7 ProEsc[0, 2, 2] -4 ProEsc[0, 4, 0] - 2 ProEsc[2, 0, 2] + 3 ProEsc[2, 2, 0]
		  +ProEsc[4, 0, 0])*
     g[1, 5] ProTen[1, 2] ProTen[2, 3] ProTen[3, 4]
     + 4/lambda^3 (4 ProEsc[0, 0, 4] + 7 ProEsc[0, 2, 2] -ProEsc[0, 4, 0] - 3 ProEsc[2, 0, 2] + 2 ProEsc[2, 2, 0]
		   -ProEsc[4, 0, 0])*
     g[1, 5] ProTen[1, 3] ProTen[2, 4] ProTen[3, 2]
     + 4/lambda^3 (3 ProEsc[0, 0, 4] - 6 ProEsc[0, 2, 2] +3 ProEsc[0, 4, 0] - ProEsc[2, 0, 2] - ProEsc[2, 2, 0]
		   -2 ProEsc[4, 0, 0]) g[2, 3] ProTen[2, 1] ProTen[2, 5] ProTen[3, 4]
     - 4/lambda^2 ProEsc[2, 0, 0] g[2, 5] g[3, 4] ProTen[3, 1]
     -4/lambda^3 (ProEsc[0, 0, 4] - 2 ProEsc[0, 2, 2] + ProEsc[0, 4, 0] -7 ProEsc[2, 0, 2] + 3 ProEsc[2, 2, 0]
		  - 4 ProEsc[4, 0, 0]) g[2, 5] ProTen[1, 3] ProTen[2, 1] ProTen[3, 4]
     +4/lambda^3 (4 ProEsc[0, 0, 4] - 3 ProEsc[0, 2, 2] -ProEsc[0, 4, 0] + 7 ProEsc[2, 0, 2] + 2 ProEsc[2, 2, 0]
		  -ProEsc[4, 0, 0])*
     g[2, 5] ProTen[1, 4] ProTen[2, 3] ProTen[3, 1]
     - 4/lambda^3 (ProEsc[0, 0, 4] - 2 ProEsc[0, 2, 2] +ProEsc[0, 4, 0] + 3 ProEsc[2, 0, 2] - 7 ProEsc[2, 2, 0]
		   -4 ProEsc[4, 0, 0]) g[3, 5] ProTen[1, 2] ProTen[2, 4] ProTen[3,1]
     - 4/lambda^3 (ProEsc[0, 0, 4] + 3 ProEsc[0, 2, 2] -4 ProEsc[0, 4, 0] - 2 ProEsc[2, 0, 2] - 7 ProEsc[2, 2, 0]
		   +ProEsc[4, 0, 0]) g[3, 5] ProTen[1, 4] ProTen[2, 1] ProTen[3, 2] );


(*Simétrico bajo crossing 1-2*) (*Añadí un signo menos global para que proyectara sobre T54 y no sobre -T54*)
Proyector54 =
  -( -40/lambda^3 (ProEsc[0, 2, 0] - ProEsc[2, 0, 0])*
     ProTen[1, 2] ProTen[1, 5] ProTen[2, 3] ProTen[2, 4] ProTen[3, 1]
     +1/lambda^2 (ProEsc[0, 0, 2] + ProEsc[0, 2, 0] - ProEsc[2, 0, 0]) g[1, 2] g[3, 5] ProTen[1, 4]
     -1/lambda^2 (ProEsc[0, 0, 2] - ProEsc[0, 2, 0] + ProEsc[2, 0, 0]) g[1, 2] g[3, 5] ProTen[2, 4]
     +2/lambda^2 (ProEsc[0, 2, 0] - ProEsc[2, 0, 0]) g[1, 3] g[2, 5] ProTen[1, 4]
     +1/lambda^2 (ProEsc[0, 0, 2] - ProEsc[0, 2, 0] -3 ProEsc[2, 0, 0]) g[1, 3] g[2, 5] ProTen[3, 4]
     +2/lambda^2 g[1, 3] ProTen[1, 4] ProTen[3, 2] ProTen[3, 5]
     +4/lambda^2 ProEsc[0, 2, 0] g[1, 4] g[3, 5] ProTen[1, 2]
     +2/lambda^2 (ProEsc[0, 2, 0] - ProEsc[2, 0, 0]) g[1, 5] g[2, 3] ProTen[2, 4]
     -1/lambda^2 (ProEsc[0, 0, 2] - 3 ProEsc[0, 2, 0] -ProEsc[2, 0, 0]) g[1, 5] g[2, 3] ProTen[3, 4]
     +4/lambda^2 (ProEsc[0, 2, 0] - ProEsc[2, 0, 0]) g[1, 5] g[2, 4] ProTen[2, 3]
     +2/lambda^3 (3 ProEsc[0, 0, 4] - 6 ProEsc[0, 2, 2] -17 ProEsc[0, 4, 0] - 6 ProEsc[2, 0, 2] + 14 ProEsc[2, 2, 0]
	       +3 ProEsc[4, 0, 0]) g[1, 5] ProTen[1, 2] ProTen[2, 3] ProTen[3, 4]
     + 2/lambda^3 (2 ProEsc[0, 0, 4] - 14 ProEsc[0, 2, 2] -8 ProEsc[0, 4, 0] + 6 ProEsc[2, 0, 2]
		   + 16 ProEsc[2, 2, 0] -8 ProEsc[4, 0, 0]) g[1, 5] ProTen[1, 3] ProTen[2, 4] ProTen[3, 2]
     - 2/lambda^2 g[2, 3] ProTen[2, 1] ProTen[2, 5] ProTen[3, 4]
     -4/lambda^2 ProEsc[2, 0, 0] g[2, 5] g[3, 4] ProTen[3, 1]
     -2/lambda^3 (3 ProEsc[0, 0, 4] - 6 ProEsc[0, 2, 2] +3 ProEsc[0, 4, 0] - 6 ProEsc[2, 0, 2] + 14 ProEsc[2, 2, 0]
	       -17 ProEsc[4, 0, 0]) g[2, 5] ProTen[1, 3] ProTen[2, 1] ProTen[3, 4]
     - 2/lambda^3 (2 ProEsc[0, 0, 4] + 6 ProEsc[0, 2, 2] -8 ProEsc[0, 4, 0] - 14 ProEsc[2, 0, 2] + 16 ProEsc[2, 2, 0]
		   -8 ProEsc[4, 0, 0]) g[2, 5] ProTen[1, 4] ProTen[2, 3] ProTen[3, 1]
     - 2/lambda^3 (3 ProEsc[0, 0, 4] - 6 ProEsc[0, 2, 2] +3 ProEsc[0, 4, 0] + 4 ProEsc[2, 0, 2] + 4 ProEsc[2, 2, 0]
		   -7 ProEsc[4, 0, 0]) g[3, 5] ProTen[1, 2] ProTen[2, 4] ProTen[3, 1]
     + 2/lambda^3 (3 ProEsc[0, 0, 4] + 4 ProEsc[0, 2, 2] -7 ProEsc[0, 4, 0] - 6 ProEsc[2, 0, 2] + 4 ProEsc[2, 2, 0]
		   +3 ProEsc[4, 0, 0]) g[3, 5] ProTen[1, 4] ProTen[2, 1] ProTen[3, 2]);


(*Nótese que algunos de los proyectores anteriores no corresponden a los \hat{T}_{i} de Colangelo.
 A continuación mostraremos los proyectores correspondientes a las demás estrctures \hat{T}_{i}
 que usa Bijnens basados en las relaciones de crossing entre ellas: *)

Proyector2 = Crossing[Proyector1, 2, 3];
Proyector3 = Crossing[Proyector1, 1, 3];
Proyector5 = Crossing[Proyector4, 2, 3];
Proyector6 = Crossing[Proyector4, 1, 3];
Proyector8 = Crossing[Proyector7, 1, 2];
Proyector9 = Crossing[Crossing[Proyector7, 2, 3], 1, 3];
Proyector10 = Crossing[Proyector7, 2, 3];
Proyector11 = Crossing[Proyector17, 1, 3];
Proyector13 = Crossing[Proyector7, 1, 3];
Proyector14 = Crossing[Crossing[Proyector7, 1, 3], 2, 3];
Proyector16 = Crossing[Proyector17, 2, 3];
Proyector50 = -Crossing[Proyector54, 2, 3];
Proyector51 = Crossing[Proyector50, 1, 2];

T1 = Contract[LCD[a[[1]], a[[2]], cc1, cc2] FVD[q[[1]], cc1] FVD[q[[2]],cc2]] *Contract[LCD[a[[3]], a[[4]], cc3, cc4] FVD[q[[3]], cc3] FVD[q[[4]],cc4]];

T4 =( (ProTen[2, 1] ProTen[1, 2] -g[1, 2] SP[1, 2]) (ProTen[4, 3] ProTen[3, 4] -g[3, 4] SP[3, 4]) );

T7 =( (ProTen[2, 1] ProTen[1, 2] -g[1, 2] SP[1, 2])*
      (SP[1,4] (ProTen[1, 3] ProTen[3, 4] - g[3, 4] SP[1, 3])
       +SP[1, 3] ProTen[4, 3] ProTen[1, 4] -ProTen[1, 3] ProTen[1, 4] SP[3, 4]) );


T19 =( (ProTen[2, 1] ProTen[1, 2] -g[1, 2] SP[1, 2])*
       (SP[2, 4] (ProTen[1, 3] ProTen[3, 4] - g[3, 4] SP[1, 3])
	+ SP[1, 3] ProTen[4, 3] ProTen[2, 4]-ProTen[1, 3] ProTen[2, 4] SP[3, 4]) );

T31 =( (ProTen[2, 1] ProTen[1, 2] -g[1, 2] SP[1, 2])*
       (ProTen[2, 3] SP[1, 3] -ProTen[1, 3] SP[2, 3]) (ProTen[2, 4] SP[1, 4] -ProTen[1, 4] SP[2, 4]) );

T37 =( (ProTen[3, 1] SP[1, 4] -ProTen[4, 1] SP[1, 3])*
       (ProTen[3, 2] ProTen[4, 3] ProTen[2, 4] -ProTen[4, 2] ProTen[2, 3] ProTen[3, 4]
	+g[3, 4] (ProTen[4, 2] SP[2, 3] - ProTen[3, 2] SP[2, 4])
	+g[2, 4] (ProTen[2, 3] SP[3, 4] - ProTen[4, 3] SP[2, 3])
	+g[2, 3] (ProTen[3, 4] SP[2, 4] - ProTen[2, 4] SP[3, 4])) );

T49 =( ProTen[3,4] (SP[1, 3] SP[2, 4] ProTen[4, 1] g[3, 2] - SP[2, 3] SP[1, 4] ProTen[4, 2] g[1, 3]
		    +ProTen[4, 1] ProTen[4,2]*(ProTen[1, 3] SP[2, 3] - ProTen[2, 3] SP[1, 3])
		    +SP[1, 4] ProTen[3, 1] ProTen[4, 2] ProTen[2, 3] -SP[2, 4] ProTen[4, 1] ProTen[3, 2] ProTen[1, 3]
		    +SP[1, 4] SP[2,4] (ProTen[3, 2] g[1, 3] - ProTen[3, 1] g[2, 3]))
       -ProTen[4,3] (SP[1, 4] SP[2, 3] ProTen[3, 1] g[2, 4] -SP[2, 4] SP[1, 3] ProTen[3, 2] g[1, 4]
		     +ProTen[3, 1] ProTen[3, 2] (ProTen[1, 4] SP[2, 4] - ProTen[2, 4] SP[1, 4])
		     +SP[1, 3] ProTen[4, 1] ProTen[3, 2] ProTen[2, 4] -SP[2, 3] ProTen[3, 1] ProTen[4, 2] ProTen[1, 4] +
		     SP[1, 3] SP[2, 3] (ProTen[4, 2] g[1, 4] - ProTen[4, 1] g[2, 4]))
       +SP[3,4] ((ProTen[1, 3] ProTen[4, 1] -SP[1, 4] g[1, 3]) (ProTen[3, 2] ProTen[2, 4] -SP[2, 3] g[2, 4])
		 - (ProTen[2, 3] ProTen[4, 2] -SP[2, 4] g[2, 3]) (ProTen[3, 1] ProTen[1, 4] -SP[1, 3] g[1, 4])) );

(*Aquí voy a poner las \hat{T}_{i} de Colangelo que contribuyen
 después de derivar con respecto a q4 y del límite estático*)

T2 = Crossing[T1, 1, 4];
T3 = Crossing[T1, 1, 3];
T5 = Crossing[T4, 1, 4];
T6 = Crossing[T4, 1, 3];
T8 = Crossing[T7, 1, 2];
T9 = Crossing[Crossing[T7, 2, 3], 1, 3];
T10 = Crossing[T7, 2, 3];
T11 = Crossing[T7, 2, 4];
T13 = Crossing[T7, 1, 3];
T14 = Crossing[Crossing[T7, 1, 3], 2, 3];
T16 = Crossing[Crossing[T7, 1, 4], 2, 4];
T17 = Crossing[Crossing[T7, 1, 3], 2, 4];
T39 = Crossing[T37, 1, 4];
T40 = Crossing[Crossing[T37, 1, 4], 1, 2];
T42 = Crossing[Crossing[T37, 2, 4], 1, 2];
T43 = Crossing[T37, 2, 4];
T46 = Crossing[Crossing[T37, 2, 3], 1, 4];
T50 = Crossing[Crossing[T49, 2, 4], 1, 2];
T51 = Crossing[T49, 2, 4];

(*Aquí van las derivadas de las \hat{T}_{i} de Colangelo*)

dT1 = FourDivergence[T1, ProTen[4, 5]] /. q[[4]] -> 0;
dT2 = FourDivergence[T2, ProTen[4, 5]] /. q[[4]] -> 0;
dT3 = FourDivergence[T3, ProTen[4, 5]] /. q[[4]] -> 0;
dT4 = FourDivergence[T4, ProTen[4, 5]] /. q[[4]] -> 0;
dT5 = FourDivergence[T5, ProTen[4, 5]] /. q[[4]] -> 0;
dT6 = FourDivergence[T6, ProTen[4, 5]] /. q[[4]] -> 0;
dT7 = FourDivergence[T7, ProTen[4, 5]] /. q[[4]] -> 0;
dT8 = FourDivergence[T8, ProTen[4, 5]] /. q[[4]] -> 0;
dT9 = FourDivergence[T9, ProTen[4, 5]] /. q[[4]] -> 0;
dT10 = FourDivergence[T10, ProTen[4, 5]] /. q[[4]] -> 0;
dT11 = FourDivergence[T11, ProTen[4, 5]] /. q[[4]] -> 0;
dT13 = FourDivergence[T13, ProTen[4, 5]] /. q[[4]] -> 0;
dT14 = FourDivergence[T14, ProTen[4, 5]] /. q[[4]] -> 0;
dT16 = FourDivergence[T16, ProTen[4, 5]] /. q[[4]] -> 0;
dT17 = FourDivergence[T17, ProTen[4, 5]] /. q[[4]] -> 0;
dT39 = FourDivergence[T39, ProTen[4, 5]] /. q[[4]] -> 0;
dT40 = FourDivergence[T40, ProTen[4, 5]] /. q[[4]] -> 0;
dT42 = FourDivergence[T42, ProTen[4, 5]] /. q[[4]] -> 0;
dT43 = FourDivergence[T43, ProTen[4, 5]] /. q[[4]] -> 0;
dT46 = FourDivergence[T46, ProTen[4, 5]] /. q[[4]] -> 0;
dT50 = FourDivergence[T50, ProTen[4, 5]] /. q[[4]] -> 0;
dT51 = FourDivergence[T51, ProTen[4, 5]] /. q[[4]] -> 0;
dTColangelo = {dT1, dT2, dT3, dT4, dT5, dT6, dT7, dT8, dT9, dT10,
	       dT11, dT13, dT14, dT16, dT17, dT39, dT40, dT42, dT43, dT50,
	       dT51};
(*Nótese que separé T39+T40 y T50-T51 por facilidad*)



(*Estas son las \hat{T}_{i} adicionales de Bijnens y sus derivadas en	\
 el límite estático.*)(*En lugar de T39+T40, T42, T43, T50-T51,		\
		       Bijnens utiliza: T39+T40+T46, T50, T51, T54*)
T54 = Crossing[Crossing[T49, 1, 4], 2, 3];
dT54 = FourDivergence[T54, ProTen[4, 5]] /. q[[4]] -> 0;
T50 = -Crossing[T54, 2, 3];
dT50 = FourDivergence[T50, ProTen[4, 5]] /. q[[4]] -> 0;
T51 = Crossing[T50, 1, 2];
dT51 = FourDivergence[T51, ProTen[4, 5]] /. q[[4]] -> 0;
dTBijnens = {dT1, dT2, dT3, dT4, dT5, dT6, dT7, dT8, dT9, dT10, dT11,
	     dT13, dT14, dT16, dT17, (1/3)*(dT39+dT40+dT46), dT50, dT51, dT54};

(*Estos son los proyectores sobre las derivadas de las \hat{T}_{i} de	\
 Bijnens.*)

ProyectoresBijnens = {Proyector1, Proyector2, Proyector3, Proyector4,
		      Proyector5, Proyector6, Proyector7, Proyector8, Proyector9,
		      Proyector10, Proyector11, Proyector13, Proyector14, Proyector16,
		      Proyector17, Proyector39, Proyector50, Proyector51, Proyector54};

(*PruebaProyectores2 = Table[Simplify[ExpandScalarProduct[Contract[Expand[ProyectoresBijnens[[1]]*dTBijnens[[i]]]]
							  /.D -> 4 /. q[[3]] -> -q[[1]] - q[[2]]]], {i, Length[dTBijnens]}];*)
 (*Prueba un proyector de Bijnens vs todas las estructuras de Bijnens*)

(*"El resultado de probar todos los proyectores de Bijnens contras todas las estructuras de Bijnens es el
siguiente:\n\n" // OutputForm >>> res.txt;

PruebaProyectores =(
  Parallelize[Table[Simplify[ExpandScalarProduct[Contract[Expand[ProyectoresBijnens[[j]]*dTBijnens[[i]]]]
						 /. D -> 4 /. q[[3]] -> -q[[1]] - q[[2]]]], {i, Length[dTBijnens]}, {j, Length[ProyectoresBijnens]}]]) >>> res.txt;
 (*Esta prueba todas los proyectores de Bijnens contra todas las estructuras de Bijnens*)*)


Permu = Permutations[{1, 2, 4}];

TodasLasAmplitudes = Table[Parallelize[Combinada[Permu[[i]][[1]], Permu[[i]][[2]], 3, Permu[[i]][[3]], 5]], {i, Length[Permu]}];
(*Calcula todas las amplitudes del quark loop con cada una de las permutaciones posibles.*)

AmplitudCompleta = Sum[TodasLasAmplitudes[[i]],{i,Length[TodasLasAmplitudes]}];
"\n\n\n\n La amplitud completa del quark loop es igual a:\n\n" // OutputForm >> ampCompleta.txt; AmplitudCompleta >>> ampCompleta.txt;

(*hjk = Parallelize[Expand[TodasLasAmplitudes[[6]]]];
 hjk = Parallelize[Table[hjk[[i]], {i, Length[hjk]}]];
 

 "Los índices de las partes que tienen más o menos de 1 índice de tipo a[[5]] son: \n\n" // OutputForm >> errores6.txt;
 errores = DeleteDuplicates[DeleteCases[Flatten[Parallelize[Table[If[Length[Cases[hjk[[i]], Pair[LorentzIndex[a[[j]], D], _],4]] !=1, i, 0],
   {i, Length[hjk]}, {j, Length[a]}]]], 0]] >>> errores6.txt;
 "\n\n\n\n Las partes problemáticas correspondientes a los índices anteriores son: \n\n" // OutputForm >>> errores6.txt;
 hjk[[errores]] >>> errores6.txt;
 (*Prueba que me da los elementos de la amplitud completa  donde el índice a[[5]] aparece más o menos de una vez*)
 
 hjk >> errores.txt;*)

(*Suma todas las permutaciones del quark loop*)

Proyector8 = Crossing[Proyector7, 1, 2];
Proyector9 = Crossing[Crossing[Proyector7, 2, 3], 1, 3];
Proyector10 = Crossing[Proyector7, 2, 3];
Proyector11 = Crossing[Proyector17, 1, 3];
Proyector13 = Crossing[Proyector7, 1, 3];
Proyector14 = Crossing[Crossing[Proyector7, 1, 3], 2, 3];
Proyector16 = Crossing[Proyector17, 2, 3];
Proyector50 = -Crossing[Proyector54, 2, 3];
Proyector51 = Crossing[Proyector50, 1, 2];



hatPi = Range[19]; (*Aquí voy a guardar los 19 coeficientes escalare que voy a evaluar a continuación.*)
(*
hatPi[[1]] = ExpandScalarProduct[Contract[Parallelize[Expand[ProyectoresBijnens[[1]]*AmplitudCompleta]]
					  /. D -> 4 /. q[[3]] -> -q[[1]] - q[[2]]]];
"\n\n\n\n La Pi 1 es igual a:\n\n" // OutputForm >> hatPi1.txt; hatPi[[1]] >>> hatPi1.txt;

hatPi[[2]] = Crossing[hatPi[[1]], 2, 3];

"\n\n\n\n La Pi 2 es igual a:\n\n" // OutputForm >> hatPi2.txt; hatPi[[2]] >>> hatPi2.txt;

hatPi[[3]] = Crossing[hatPi[[1]], 1, 3];

"\n\n\n\n La Pi 3 es igual a:\n\n" // OutputForm >> hatPi3.txt; hatPi[[3]] >>> hatPi3.txt;

hatPi[[4]] = ExpandScalarProduct[Contract[Parallelize[Expand[ProyectoresBijnens[[4]]*AmplitudCompleta]]
                                                   /. D -> 4 /. q[[3]] -> -q[[1]] - q[[2]]]];
"\n\n\n\n La Pi 4 es igual a:\n\n" // OutputForm >> hatPi4.txt; hatPi[[4]] >>> hatPi4.txt;

hatPi[[5]] = Crossing[hatPi[[4]], 2, 3];
  
"\n\n\n\n La Pi 5 es igual a:\n\n" // OutputForm >> hatPi5.txt; hatPi[[5]] >>> hatPi5.txt;

hatPi[[6]] = Crossing[hatPi[[4]], 1, 3]:
  
"\n\n\n\n La Pi 6 es igual a:\n\n" // OutputForm >> hatPi6.txt; hatPi[[6]] >>> hatPi6.txt;
  
hatPi[[7]] = ExpandScalarProduct[Contract[Parallelize[Expand[ProyectoresBijnens[[7]]*AmplitudCompleta]]
                                                   /. D -> 4 /. q[[3]] -> -q[[1]] - q[[2]]]];
"\n\n\n\n La Pi 7 es igual a:\n\n" // OutputForm >> hatPi7.txt; hatPi[[7]] >>> hatPi7.txt;

hatPi[[8]] = Crossing[hatPi[[7]], 1, 2];
                                        
"\n\n\n\n La Pi 8 es igual a:\n\n" // OutputForm >> hatPi8.txt; hatPi[[8]] >>> hatPi8.txt;

hatPi[[9]] = Crossing[Crossing[hatPi[[7]], 2, 3], 1, 3];
  
"\n\n\n\n La Pi 9 es igual a:\n\n" // OutputForm >> hatPi9.txt; hatPi[[9]] >>> hatPi9.txt;

hatPi[[10]] = Crossing[hatPi[[7]], 2, 3];
  
"\n\n\n\n La Pi 10 es igual a:\n\n" // OutputForm >> hatPi10.txt; hatPi[[10]] >>> hatPi10.txt;

hatPi[[15]] = Simplify[ExpandScalarProduct[Contract[Parallelize[Expand[ProyectoresBijnens[[15]]*AmplitudCompleta]] (*17*)
						    /. D -> 4 /. q[[3]] -> -q[[1]] - q[[2]]]]];
"\n\n\n\n La Pi 15 es igual a:\n\n" // OutputForm >> hatPi15.txt; hatPi[[15]] >>> hatPi15.txt;

hatPi[[11]] = Crossing[hatPi[[15]], 1, 3];
  
"\n\n\n\n La Pi 11 es igual a:\n\n" // OutputForm >> hatPi11.txt; hatPi[[11]] >>> hatPi11.txt;

hatPi[[12]] = Crossing[hatPi[[7]], 1, 3]; (*13*)
  
"\n\n\n\n La Pi 12 es igual a:\n\n" // OutputForm >> hatPi12.txt; hatPi[[12]] >>> hatPi12.txt;

hatPi[[13]] = Crossing[Crossing[hatPi[[7]], 1, 3], 2, 3]; (*14*)
  
"\n\n\n\n La Pi 13 es igual a:\n\n" // OutputForm >> hatPi13.txt; hatPi[[13]] >>> hatPi13.txt;

hatPi[[14]] = Crossing[hatPi[[15]], 2, 3]; (*16*)
							    
"\n\n\n\n La Pi 14 es igual a:\n\n" // OutputForm >> hatPi14.txt; hatPi[[14]] >>> hatPi14.txt;*)
 Print["Hola"];
hatPi[[16]] = ExpandScalarProduct[Contract[Parallelize[Expand[ProyectoresBijnens[[16]]*AmplitudCompleta]] (*39*)
                                                   /. D -> 4 /. q[[3]] -> -q[[1]] - q[[2]]]];
"\n\n\n\n La Pi 16 es igual a:\n\n" // OutputForm >> hatPi16.txt; hatPi[[16]] >>> hatPi16.txt;

hatPi[[19]] = ExpandScalarProduct[Contract[Parallelize[Expand[ProyectoresBijnens[[19]]*AmplitudCompleta]]
                                                   /. D -> 4 /. q[[3]] -> -q[[1]] - q[[2]]]];
"\n\n\n\n La Pi 19 es igual a:\n\n" // OutputForm >> hatPi19.txt; hatPi[[19]] >>> hatPi19.txt;

hatPi[[17]] = -Crossing[hatPi[[19]], 2, 3]; (*50*)
  
"\n\n\n\n La Pi 17 es igual a:\n\n" // OutputForm >> hatPi17.txt; hatPi[[17]] >>> hatPi17.txt;

hatPi[[18]] = Crossing[hatPi[[17]], 1, 2]; (*51*)
  
"\n\n\n\n La Pi 18 es igual a:\n\n" // OutputForm >> hatPi18.txt; hatPi[[18]] >>> hatPi18.txt;


 (*adicionales = Sum[Parallelize[Expand[AmplitudCompleta -hatPi[[i]]*dTBijnens[[i]],{i,Length[hatPi]}]]];

"\n\n\n\n\n\n Si se restan las contribuciones de los proyectores a la amplitud se obtiene:"] // OutputForm >> residuos.txt; adicionales >>> residuos.txt;*)

(*
Clear[IntSca2]

IntSca2[N_, n_, nus_] := Block[{k2, L, s},
			       L = N*(N - 1)/2;
			       k2 = Table[
					  If[i == j || i > j, 0,
					     ToExpression[ToString[k] <> ToString[i] <> ToString[j]]], {i,
													N}, {j, N}];
			       (*k2=Table[If[i\[Equal]j||i>j,0,Pair[Momentum[Q[[i]]-Q[[j]],D],
								    Momentum[Q[[i]]-Q[[j]],D]]],{i,3},{j,3}];*)
  s = Table[
	    If[i == j || i > j, 0,
	       ToExpression[ToString[s] <> ToString[i] <> ToString[j]]], {i,
									  N}, {j, N}];(*Estas son las s de Davydychev y las z de Friot*)
									  sFlat = DeleteCases[Flatten[s], 0];
x = DeleteCases[Flatten[Table[-k2[[j]][[l]], {j, N}, {l, N}]],
		0];(*Las x de Friot*)
		numerator1 = -sFlat;
numerator2 =
  Sum[nus[[i]], {i, N}] - n/2 + Sum[sFlat[[j]], {j, Length[sFlat]}];
numerator3 =
  Table[nus[[i]] + Sum[s[[j]][[i]], {j, Length[s]}] +
	  Sum[If[l > i, s[[i]][[l]], 0], {l, Length[s[[i]]]}], {i, N}];
numerator = Flatten[{numerator1, numerator2, numerator3}];
denominator = {Sum[nus[[i]], {i, N}] +
	       2 Sum[sFlat[[j]], {j, N*(N - 1)/2}]};
  const =
    Pi^(n/2)*
    I^(1 - n)*(-m)^(n/2 - Sum[nus[[i]], {i, 3}])*(1/
						  Product[Gamma[nus[[i]]], {i, 3}])(*(1/(2*Pi*I)^L)*);
MB = MBRep[const,
	   sFlat, -DeleteCases[Flatten[x], 0], {numerator, denominator}];
Resol = ResolveMB[MB, 16];
Eval = EvaluateSeries[Resol, {}, 2];
  Final =
    SumAllSeries[Eval,
		     Flatten[Table[
				   k2[[i]][[j]] -> -1000, {i, Length[k2]}, {j, Length[k2[[1]]]}]],
		 5];(*Eval2=EvaluateSeries[Resol,{},2];
		     Final2=SumAllSeries[Eval,Flatten[Table[k2[[i]][[j]]\[Rule] -1000,{i,
										       Length[k2]},{j,Length[k2[[1]]]}]],5];Eval3=EvaluateSeries[Resol,{},
																		 3];
Final4=SumAllSeries[Eval,Flatten[Table[k2[[i]][[j]]\[Rule] -1000,{i,
								  Length[k2]},{j,Length[k2[[1]]]}]],5];Eval4=EvaluateSeries[Resol,{},
															    4];
Final4=SumAllSeries[Eval,Flatten[Table[k2[[i]][[j]]\[Rule] -1000,{i,
								  Length[k2]},{j,Length[k2[[1]]]}]],5];Eval5=EvaluateSeries[Resol,{},
															    5];
Final5=SumAllSeries[Eval,Flatten[Table[k2[[i]][[j]]\[Rule] -1000,{i,
								  Length[k2]},{j,Length[k2[[1]]]}]],5];Eval6=EvaluateSeries[Resol,{},
															    6];
Final6=SumAllSeries[Eval,Flatten[Table[k2[[i]][[j]]\[Rule] -1000,{i,
								  Length[k2]},{j,Length[k2[[1]]]}]],5];Eval7=EvaluateSeries[Resol,{},
															    7];
Final7=SumAllSeries[Eval,Flatten[Table[k2[[i]][[j]]\[Rule] -1000,{i,
								  Length[k2]},{j,Length[k2[[1]]]}]],5];Eval8=EvaluateSeries[Resol,{},
															    8];
Final8=SumAllSeries[Eval,Flatten[Table[k2[[i]][[j]]\[Rule] -1000,{i,
								  Length[k2]},{j,Length[k2[[1]]]}]],5];Eval9=EvaluateSeries[Resol,{},
															    9];
Final9=SumAllSeries[Eval,Flatten[Table[k2[[i]][[j]]\[Rule] -1000,{i,
								  Length[k2]},{j,Length[k2[[1]]]}]],5];Eval10=EvaluateSeries[Resol,{},
															     10];
Final10=SumAllSeries[Eval,Flatten[Table[k2[[i]][[
						 j]]\[Rule] -1000,{i,Length[k2]},{j,Length[k2[[1]]]}]],5];Eval11=
													    EvaluateSeries[Resol,{},11];
Final11=SumAllSeries[Eval,Flatten[Table[k2[[i]][[
						 j]]\[Rule] -1000,{i,Length[k2]},{j,Length[k2[[1]]]}]],5];Eval12=
													    EvaluateSeries[Resol,{},12];
Final12=SumAllSeries[Eval,Flatten[Table[k2[[i]][[
						 j]]\[Rule] -1000,{i,Length[k2]},{j,Length[k2[[1]]]}]],5];Eval13=
													    EvaluateSeries[Resol,{},13];
Final13=SumAllSeries[Eval,Flatten[Table[k2[[i]][[
						 j]]\[Rule] -1000,{i,Length[k2]},{j,Length[k2[[1]]]}]],5];Eval14=
													    EvaluateSeries[Resol,{},14];
Final14=SumAllSeries[Eval,Flatten[Table[k2[[i]][[
						 j]]\[Rule] -1000,{i,Length[k2]},{j,Length[k2[[1]]]}]],5];Eval15=
													    EvaluateSeries[Resol,{},15];
Final15=SumAllSeries[Eval,Flatten[Table[k2[[i]][[
						 j]]\[Rule] -1000,{i,Length[k2]},{j,Length[k2[[1]]]}]],5];Eval16=
													    EvaluateSeries[Resol,{},16];
Final16=SumAllSeries[Eval,Flatten[Table[k2[[i]][[
						 j]]\[Rule] -1000,{i,Length[k2]},{j,Length[k2[[1]]]}]],5];*)
						 Return[Final];
(*Return[{const,s,k2,{numerator,denominator},sFlat}];*)
  ]

  Length[TodasLasEscalares[[1]]]


  (* IntSca2[3, 4, {1, 1, 1}] // InputForm

  MB = MBRep[1, {s12, s13, s23}, {-k12, -k13, -k23}, {{-s12, -s13, -s23}, {s12}}]

  Resol = ResolveMB[MB]

  A1 = MBRep[1, {z1, z2}, {-u1, -u2}, {{-z1, -z2, 2 + z1 + z2, 1 + z1, 1 + z2}, {0.5 + z1 + z2}}]

  (*SumAllSeries[C1,{u1->-0.3,u2->-10.1},5]*)

  B1 = ResolveMB[A1, 5]

  Flatten[Cases[B1, StringForm[__], 9]] /. StringForm[x_] -> x

  C1 = EvaluateSeries[B1, {}, 2] // OutputForm *)
   *)
