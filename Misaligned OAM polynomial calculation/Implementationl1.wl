(* ::Package:: *)

(* ::Title:: *)
(*Misaligned OAM-mode calculation  for l0 = 1*)


(* ::Text:: *)
(*Project: META-LiFi*)
(*Date: 24/07/2025*)
(*Version 1*)


(* ::Input:: *)
(*ClearAll["Global`*"]*)


(* ::Input:: *)
(*signfunction[k_]:= If[k<0,(-1)^k,1]*)
(*Signs[n_,Dl_,l0_]:= signfunction[n]*signfunction[Dl-n]*signfunction[l0+Dl-n]*)
(*Sfunc[n_,Dl_,l0_,j_,k_]:= Signs[n,Dl,l0]*(-1)^(j+k)*)
(*s[l_,l0_,n_,k_]:=1/2 (Abs[l-n]-Abs[l-l0-n]-Abs[l0]-2k)*)
(*B[l_,l0_,n_,j_,k_]:=2^Abs[l0] 2^(-1/2 (Abs[l-n]+ Abs[n]+2j)) (Factorial[1/2 (Abs[l0]+Abs[l-n]+Abs[l-l0-n]+2k)]/(Factorial[Abs[l0]]Factorial[Abs[l-l0-n]]Factorial[Abs[n]]Factorial[Abs[l-n]]Factorial[j]Factorial[ Abs[n]+j]Factorial[k]Factorial[Abs[l-l0-n]+k]))*)
(*Fullterm[l_,l0_,n_,m_,j_,jp_,k_,kp_]:=Sfunc[n,l-l0,l0,j,k]Sfunc[m,l-l0,l0,jp,kp]B[l,l0,n,j,k]B[l,l0,m,jp,kp] *)
(*Ifunc[l_,l0_,n_,m_,j_,jp_,k_,kp_,x_]:= x^(Abs[l-n]+ Abs[n] + Abs[l-m]+ Abs[m]+2*(j+jp)+1)*Exp[-x^2]*Hypergeometric1F1[s[l,l0,n,k],Abs[l-n]+1,x^2/2]*Hypergeometric1F1[s[l,l0,m,kp],Abs[l-m]+1,x^2/2]*)
(*Polynomial[l_,l0_,n_,m_,j_,jp_,k_,kp_] := Fullterm[l,l0,n,m,j,jp,k,kp]*Integrate[Ifunc[l,l0,n,m,j,jp,k,kp,x],{x,0,\[Infinity]}]*Cos[(n-m)*(\[Phi]+\[Pi]/2)]*X0^(Abs[n]+Abs[m]+2(j+jp) ) T0^(Abs[l-l0-n]+Abs[l-l0-m]+2(k+kp))*)
(**)


(* ::Input:: *)
(*l0test = 1;*)
(*l1  = 2;*)
(*l2 = 0;*)


(* ::Subchapter:: *)
(*Lets first check l = l0 = 1*)


(* ::Input:: *)
(*ltest= l0test;*)


(* ::Text:: *)
(*Lets start with term c = 1, m=n=0, all other terms terms also 0*)


(* ::Input:: *)
(*ntest = 0;*)
(*mtest = 0;*)
(*jtest = 0;*)
(*jptest = 0;*)
(*ktest = 0;*)
(*kptest = 0;*)


(* ::Input:: *)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*Pol1 = Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)


(* ::Text:: *)
(*Next one is c = 0: First n = +-1. This contribution is twice (takes into account m = +-1, n= 0)*)


(* ::Input:: *)
(*ntest = 1;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*Pol2 = 2*Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)


(* ::Input:: *)
(*ntest = -1;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*Pol3 = 2*Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(**)


(* ::Text:: *)
(*Next one is c = 0, and j = 1. we can take it twice to take  jp = 1 into account*)


(* ::Input:: *)
(*ntest  = 0;*)
(*jtest = 1;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*Pol4 = 2*Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)


(* ::Text:: *)
(*We do the same with k = 1, and twice because kp = 1 is same*)


(* ::Input:: *)
(*jtest = 0;*)
(*ktest = 1;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*Pol5= 2*Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*ktest = 0*)


(* ::Input:: *)
(*FinalPol = Pol1+Pol2+Pol3+Pol4+Pol5*)


(* ::Text:: *)
(*We  see the final input is something logical. If phi = 0, then the last term is 0, and does not contribute to the error. If Phi = pi/2, then the last term is 1, and adds to error. If phi= -pi/2, it reduces error*)


(* ::Subchapter:: *)
(*Then check l = 2*)


(* ::Input:: *)
(*ltest= l1;*)


(* ::Input:: *)
(*ntest = 0;*)
(*mtest = 0;*)
(*jtest = 0;*)
(*jptest = 0;*)
(*ktest = 0;*)
(*kptest = 0;*)


(* ::Text:: *)
(*In this case, we \[CapitalDelta]l = 1  = N, so only c = 0, and bessel expansion term is zero. So now, we have to consider four (actually three) terms: n = 0, m = 0 then  n = 1, m=0, (twice), as m=1 n= 0 same, and then n=1,m=1*)
(*We start with n = 0, m = 0*)


(* ::Input:: *)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*l2Pol1 = Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)


(* ::Text:: *)
(*Next is twice n = 1, m = 0*)


(* ::Input:: *)
(*ntest = 1;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*l2Pol2 = 2*Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)


(* ::Text:: *)
(*Finally,n =1, m = 1*)


(* ::Input:: *)
(*ntest = 1;*)
(*mtest = 1;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*l2Pol3= Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)


(* ::Input:: *)
(*Finall2Pol = l2Pol1 + l2Pol2 + l2Pol3*)


(* ::Subchapter:: *)
(*Finally check l = 0*)


(* ::Input:: *)
(*ltest= 0;*)


(* ::Input:: *)
(*ntest = 0;*)
(*mtest = 0;*)
(*jtest = 0;*)
(*jptest = 0;*)
(*ktest = 0;*)
(*kptest = 0;*)


(* ::Text:: *)
(*In this case, we \[CapitalDelta]l = -1 = N, so only c = 0, and bessel expansion term is zero. So now, we have to consider four (actually three) terms: n = 0, m = 0 then  n = -1, m=0, (twice), as m=-1 n= 0 same, and then n=-1,m=-1*)
(*We start with n = 0, m = 0*)


(* ::Input:: *)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*l0Pol1 = Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)


(* ::Text:: *)
(*Next is twice n =- 1, m = 0*)


(* ::Input:: *)
(*ntest = -1;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*l0Pol2 = 2*Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)


(* ::Text:: *)
(*Finally, n=-1, m = -1*)


(* ::Input:: *)
(*ntest = -1;*)
(*mtest = -1;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*l0Pol3= Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)


(* ::Input:: *)
(*Finall0Pol = l0Pol1 + l0Pol2 + l0Pol3*)


(* ::Subchapter:: *)
(*Lets check l = -1*)


(* ::Input:: *)
(*ltest= -1;*)


(* ::Input:: *)
(*ntest = 0;*)
(*mtest = 0;*)
(*jtest = 0;*)
(*jptest = 0;*)
(*ktest = 0;*)
(*kptest = 0;*)


(* ::Text:: *)
(*In this case, we \[CapitalDelta]l = -2 = N, so only c = 0, and bessel expansion term is zero. We now have much more terms to consider: n = 0 and m = 0, then n = -1, m= 0 (twice), n = -2, m = 0 (twice), n=-1,m=-1, n=-2,m=-1 (twice), and finally n=-2,m=-2. So we have 6 terms*)
(*First lets check all 0 term*)


(* ::Input:: *)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*lneg1Pol1 = Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)


(* ::Text:: *)
(*Next is twice n = -1, m = 0 twice*)


(* ::Input:: *)
(*ntest = -1;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*lneg1Pol2 = 2*Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)


(* ::Text:: *)
(*Then, n = -2, m = 0 twice*)


(* ::Input:: *)
(*ntest = -2;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*lneg1Pol3= FullSimplify[2*Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]]*)


(* ::Text:: *)
(*Then, n=-1,m=-1*)


(* ::Input:: *)
(*ntest = -1;*)
(*mtest  =-1;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*lneg1Pol4= FullSimplify[Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]]*)


(* ::Text:: *)
(*Then, n = -2,m=-1 twice*)


(* ::Input:: *)
(*ntest = -2;*)
(*mtest  =-1;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*lneg1Pol5= FullSimplify[2*Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]]*)


(* ::Text:: *)
(*Finally, n = -2, m = -2*)


(* ::Input:: *)
(*ntest = -2;*)
(*mtest  =-2;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*lneg1Pol6= FullSimplify[Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]]*)


(* ::Input:: *)
(*Finall1negPol  = lneg1Pol1 + lneg1Pol2 + lneg1Pol3 + lneg1Pol4 + lneg1Pol5 + lneg1Pol6*)


(* ::Subchapter:: *)
(*For completeness l = 3*)


(* ::Input:: *)
(*ltest= 3;*)


(* ::Input:: *)
(*ntest = 0;*)
(*mtest = 0;*)
(*jtest = 0;*)
(*jptest = 0;*)
(*ktest = 0;*)
(*kptest = 0;*)


(* ::Text:: *)
(*In this case, we \[CapitalDelta]l = 2 = N, so only c = 0, and bessel expansion term is zero. We now have much more terms to consider: n = 0 and m = 0, then n = 1, m= 0 (twice), n = 2, m = 0 (twice), n=1,m=1, n=2,m=1 (twice), and finally n=2,m=2. So we have 6 terms*)
(*First lets check all 0 term*)


(* ::Input:: *)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*l3Pol1 = Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)


(* ::Text:: *)
(*Next is twice n = 1, m = 0 twice*)


(* ::Input:: *)
(*ntest = 1;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*l3Pol2 = 2*Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)


(* ::Text:: *)
(*Then, n = 2, m = 0 twice*)


(* ::Input:: *)
(*ntest = 2;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*l3Pol3= FullSimplify[2*Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]]*)


(* ::Text:: *)
(*Then, n=1,m=1*)


(* ::Input:: *)
(*ntest = 1;*)
(*mtest  =1;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*l3Pol4= FullSimplify[Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]]*)


(* ::Text:: *)
(*Then, n = 2,m=1 twice*)


(* ::Input:: *)
(*ntest = 2;*)
(*mtest  =1;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*l3Pol5= FullSimplify[2*Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]]*)


(* ::Text:: *)
(*Finally, n = 2, m = 2*)


(* ::Input:: *)
(*ntest = 2;*)
(*mtest  =2;*)
(*term1 =  Fullterm[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]*)
(*Ifuncval=Ifunc[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest,x]*)
(*Ifuncintval = Integrate[Ifuncval,{x,0,\[Infinity]}]*)
(*l3Pol6= FullSimplify[Polynomial[ltest,l0test,ntest,mtest,jtest,jptest,ktest,kptest]]*)


(* ::Input:: *)
(*Finall3Pol  = l3Pol1 + l3Pol2 + l3Pol3 + l3Pol4 + l3Pol5 + l3Pol6*)
