\\ \r lfl3\lfl-utils-alpha3Variable.gp

read("lfl3\\step2-utils.gp");
read("lfl3\\step3-utils.gp");

\\ assume that \alpha_3 is the variable \alpha_i (i.e., \alpha_1 and \alpha_2 are fixed numbers)

\\ uses:
\\ be used for mignotte-eg1 and the Fibonacci perfect powers linear form

\\ 11 Dec 2021
alpha3_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,dbg=0)={
	local(d1,d2);
	
	d1=1;
	d2=1;
	return(alpha3_check_params_with_d1d2(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,d1,d2,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,dbg));
}

\\ 29 June 2022
alpha3_check_params_with_d1d2(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,d1,d2,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,dbg=0)={
	my(a,aP,bigK,bigR,bigR1,bigR2,bigR3,bigS,bigS1,bigS2,bigS3,bigT,bigT1,bigT2,bigT3,c1,c2,c3,chiV,chiVSqr,cM,eqn42,eqn42LHS,eqn42RHS,g,isComplex,logBUB,logLambdaLB,rt,two13,vSqr);

	if(type(al3)!="t_POL",
		print("ERROR: type(al3)=",type(al3)," must be t_POL");
		return([]);
	);
	if(type(al1)!="t_INT" && type(al1)!="t_FRAC" && type(al1)!="t_REAL" && type(al1)!="t_COMPLEX",
		print("ERROR: type(al1)=",type(al1)," must be t_INT, t_FRAC, t_COMPLEX or t_REAL");
		return([]);
	);
	if(type(al2)!="t_INT" && type(al2)!="t_FRAC" && type(al2)!="t_REAL" && type(al2)!="t_COMPLEX",
		print("ERROR: type(al2)=",type(al2)," must be t_INT, t_FRAC, t_COMPLEX or t_REAL");
		return([]);
	);
	if(nLB<1,
		print("ERROR: nLB=",nLB," must be a positive integer");
		return([]);
	);

	\\ al3 not used again, al1 and al2 only used for c2
	a=min(a1,a2);
	\\print("a1=",a1,", a3=",a3,", a=",a);
	aP=max(a1,a2); \\ aPrime=second largest
	
	two13=1.2599210498948731647672106072782283506; \\ 2^(1/3) to save recalculation every time
	c1=two13;
	c1=max(c1,(chi*m*bigL)^(2/3));
	c1=max(c1,sqrt(2*m*bigL/a));
	if(type(al1)=="t_REAL" && type(al2)=="t_REAL",
		c2=two13*(m*bigL)^(2/3);
	);
	if(type(al1)!="t_REAL" || type(al2)!="t_REAL",
		c2=max(two13*(m*bigL)^(2/3),sqrt(m/a)*bigL);
	);
	c3=(6*m*m)^(1/3)*bigL;
	if(dbg!=0,
		printf("in alpha3_check_params_with_d1d2(): c1=%10.6f\n",c1);
		printf("in alpha3_check_params_with_d1d2(): c2=%10.6f\n",c2);
		printf("in alpha3_check_params_with_d1d2(): c3=%10.6f\n",c3);
	);
	
	bigR1=c1*a2*a3;
	bigR1=(polcoef(bigR1,1)+polcoef(bigR1,0)/logXLB)*logX;
	bigR2=c2*a2*a3;
	bigR2=(polcoef(bigR2,1)+polcoef(bigR2,0)/logXLB)*logX;
	bigR3=c3*a2*a3;
	bigR3=(polcoef(bigR3,1)+polcoef(bigR3,0)/logXLB)*logX;
	bigR=bigR1+bigR2+bigR3+1;
	bigR=(polcoef(bigR,1)+polcoef(bigR,0)/logXLB)*logX;

	bigS1=c1*a1*a3;
	bigS1=(polcoef(bigS1,1)+polcoef(bigS1,0)/logXLB)*logX;
	bigS2=c2*a1*a3;
	bigS2=(polcoef(bigS2,1)+polcoef(bigS2,0)/logXLB)*logX;
	bigS3=c3*a1*a3;
	bigS3=(polcoef(bigS3,1)+polcoef(bigS3,0)/logXLB)*logX;
	bigS=bigS1+bigS2+bigS3+1;
	bigS=(polcoef(bigS,1)+polcoef(bigS,0)/logXLB)*logX;

	bigT1=c1*a1*a2;
	bigT1=floor(bigT1);
	bigT2=c2*a1*a2;
	bigT2=floor(bigT2);
	bigT3=c3*a1*a2;
	bigT3=floor(bigT3);
	bigT=bigT1+bigT2+bigT3+1;

	bigK=bigL*m*a1*a2*a3;
	bigK=(polcoef(bigK,1)+polcoef(bigK,0)/logXLB)*logX;
	logLambdaLB=-bigK*bigL*log(rho); \\-log(bigK*bigL);
	if(dbg!=0, \\ bigL==24 && m==89 && rho==100, \\
		printf("in alpha3_check_params_with_d1d2(): nUB=%8.6e\n",nUB);
		printf("bigR1=%11.6f*logX\n",polcoef(bigR1,1));
		printf("bigR2=%11.6f*logX\n",polcoef(bigR2,1));
		printf("bigR3=%11.6f*logX\n",polcoef(bigR3,1));

		printf("bigS1=%11.6f*logX\n",polcoef(bigS1,1));
		printf("bigS2=%11.6f*logX\n",polcoef(bigS2,1));
		printf("bigS3=%11.6f*logX\n",polcoef(bigS3,1));

		printf("bigT1=%8d\n",bigT1);
		printf("bigT2=%8d\n",bigT2);
		printf("bigT3=%8d\n",bigT3);

		printf("rho =%10.6f\n",rho);
		printf("bigK=%10.6f*logX\n",polcoef(bigK,1));
		printf("bigL=%6d\n",bigL);

		printf("bigR=%10.6f*logX\n",polcoef(bigR,1));
		printf("bigS=%10.6f*logX\n",polcoef(bigS,1));
		printf("bigT=%10.6f\n",bigT);
		printf("log |Lambda|>%9.6e*logX\n\n",polcoef(logLambdaLB,1));
	);

	logBUB=get_logBPrime_UB(bigK,bigR,bigS,bigT,b1,b2,b3,d1,d2,logXLB,nUB,dbg);
	if(dbg!=0,
		printf("log(b')<%9.6f\n",logBUB);
	);

	eqn42=get_eqn42(a1,a2,a3,bigK,bigL,bigR,bigS,bigT,d,rho,logBUB,nUB,logXLB,dbg);
	rt=polrootsreal(eqn42)[1];
	if(dbg!=0,
		printf("eqn42LHS-eqn42RHS=%s\n",eqn42);
		printf("eqn42(x)=0 at x=%9.6f\n",rt);
		printf("subst(eqn42,logX,logXLB)=%9.6f\n",subst(eqn42,logX,logXLB));
	);

	\\ not sure why first condition is here??
	\\ removed without harm
	\\if(polcoef(xvLHS,1)>polcoef(xvRHS,1) && rt<logXLB && subst(f,logX,logXLB)>0,
	if(rt<logXLB && subst(eqn42,logX,logXLB)>0,
		\\ we only handle lower bounds for lambda that are linear in logX at the moment (May 2022)
		if(poldegree(logLambdaLB,logX)==1 && polcoef(logLambdaLB,0)==0,
			newNonDegenNUB=logLambdaLB/lamUB1;
			newNonDegenNUB=polcoef(newNonDegenNUB,0); \\ need this to make bR real, not a poly
			if(dbg!=0,
				printf("newNonDegenNUB=%9.6e\n",newNonDegenNUB);
			);
			isComplex = type(al1)=="t_COMPLEX" || type(al2)=="t_COMPLEX" || type(al3)=="t_COMPLEX";
			newDegenNUB=alpha3_do_degenerate_case(d,hgtA1,absLogA1,hgtA2,absLogA2,hgtA3,absLogA3,bigR1,bigS1,bigT1,bigT2,chi,logXLB,nLB,lamUB1,lamUB0,isComplex,dbg);
			if(dbg!=0,
				printf("newDegenNUBArray=%12.6e\n",newDegenNUB);
			);
			if(length(newDegenNUB)>0,
				return([polcoef(bigK,1),newNonDegenNUB,newDegenNUB]);
			);
		);
		if(poldegree(logLambdaLB,logX)!=1 || polcoef(logLambdaLB,0)!=0,
			print("BAD");
		);
	);
	if(dbg!=0,
		if(rt>logXLB || subst(eqn42,logX,logXLB)<0,
			printf("alpha3_check_params_with_d1d2(): rt=%9.6e, logXLB=%9.6e, subst(eqn42,logX,logXLB)=%9.6e\n",rt,logXLB,subst(eqn42,logX,logXLB));
		);
	);
	return([]);
}

\\ assume the S_i's are constant and R_i's and T_i's are both linear in logZ
\\ also needs to have the alpha_i's
\\ 21 Nov 2021
alpha3_do_degenerate_case(d,hgtA1,absLogA1,hgtA2,absLogA2,hgtA3,absLogA3,bigR1,bigS1,bigT1,bigT2,chi,logXLB,nLB,lamUB1,lamUB0,isComplex,dbg=0)={
	my(a,aiArray1,aiArray2,aiArray3,b,b1T,bPCnst,bPLB,chiV,chiVSqr,cM,expB1T,h,lflCnst,log2,log5,logA1,logA2,nUB1,nUB1a,nUB1b,nUB2,nUB2a,nUB2b,nUB3,nUB3a,nUB3b,u1UB,u2UB,u3UB,vSqr);

	if(poldegree(bigR1)!=1,
		print("bad R1 in alpha3_do_degenerate_case(): R1=",bigR1);
		return([]);
	);
	if(poldegree(bigS1)!=1,
		print("bad S1 in alpha3_do_degenerate_case(): S1=",bigS1);
		return([]);
	);

	if(poldegree(bigT1)!=0,
		print("bad T1 in alpha3_do_degenerate_case(): T1=",bigT1);
		return([]);
	);
	
	if(dbg!=0,
		print("\nalpha3_do_degenerate_case(): START");
		printf("(C1) and (C2) bounds: n<max(T1,T2)=%7d\n",max(bigT1,bigT2));
	);

	\\ these look like degenerate case values	
	vSqr=(bigR1+1)*(bigS1+1)*(bigT1+1);
	chiVSqr=chi*chi*vSqr;
	\\ chiV is an upper bound for chi*\cV, written as number*log(x)
	chiV=sqrt(polcoef(chiVSqr,2)+polcoef(chiVSqr,1)/logXLB+polcoef(chiVSqr,0)/logXLB^2)*logX;

	\\ cM is the other terms in the definition of \cM
	\\ looks like at end we impose the condition that chiV is larger and
	\\ use that as our value of \cM
	\\ we do that because that is what Mignotte does in the degenerate case for eg 1
	cM=polcoef(bigR1+bigT1+1,1)+polcoef(bigR1+bigT1+1,0)/logXLB;
	\\printf("cM=%s\n",cM);
	cM=max(cM,polcoef(bigS1+bigT1+1,1)+polcoef(bigS1+bigT1+1,0)/logXLB);
	\\printf("cM=%s\n",cM);
	cM=max(cM,polcoef(bigR1+bigS1+1,1)+polcoef(bigR1+bigS1+1,0)/logXLB)*logX;
	if(dbg!=0,
		printf("alpha3_do_degenerate_case(): chi*cV=%9.6f*log(x)\n",polcoef(chiV,1));
		printf("alpha3_do_degenerate_case(): initial cM=%9.6f*log(x)\n",polcoef(cM,1));
	);
	if(polcoef(chiV,1)>polcoef(cM,1),
		cM=chiV;
	);
	if(dbg!=0,
		printf("alpha3_do_degenerate_case(): cM=%9.6f*log(x)\n",polcoef(cM,1));
	);

	\\ u1UB=bR=B_R in Step 3
	\\ we assume that S_1>T_1 there
	if(subst(bigS1,logX,logXLB)<bigT1,
		return([]);
	);
	u1UB=(bigS1+1)*(bigT1+1);
	u1UB=subst(u1UB,logX,logXLB)/logXLB*logX;
	u1UB=polcoef(u1UB/(cM-bigS1),0); \\ need this to make bR real, not a poly

	\\ doing u2UB (bT=B_T) calculation: we assume that R_1>T_1 there
	if(subst(bigR1,logX,logXLB)<bigT1,
		return([]);
	);
	u2UB=(bigR1+1)*(bigT1+1);
	u2UB=subst(u2UB,logX,logXLB)/logXLB*logX;
	u2UB=polcoef(u2UB/(cM-bigR1),0); \\ need this to make bT real, not a poly (otherwise, floor() and other functions do not work)

	\\ u3UB (bS=B_S) calc:
	u3UB=(bigR1+1)*(bigS1+1);
	u3UB=subst(u3UB,logX,logXLB)/logXLB^2*logX^2;
	if(subst(bigR1,logX,logXLB)<subst(bigS1,logX,logXLB),
		u3UB=u3UB/(cM-bigS1);
	);
	if(subst(bigS1,logX,logXLB)<subst(bigR1,logX,logXLB),
		u3UB=u3UB/(cM-bigR1);
	);

	if(dbg!=0,
		printf("alpha3_do_degenerate_case(): before rounding down, u1UB=%10.6f\n",u1UB); \\ was B_R, bR
		printf("alpha3_do_degenerate_case(): before rounding down, u2UB=%10.6f\n",u2UB); \\ was B_S, bS
	);

	u1UB=floor(u1UB);
	u2UB=floor(u2UB);
	if(dbg!=0,
		printf("alpha3_do_degenerate_case(): u1UB=%10.6f\n",u1UB); \\ was B_R, bR
		printf("alpha3_do_degenerate_case(): u2UB=%10.6f\n",u2UB); \\ was B_S, bS
		printf("alpha3_do_degenerate_case(): u3UB=%10.6f*logX\n",polcoef(u3UB,1)); \\ was B_S, bS
	);
	
	nUB1=0;
	if(type(u1UB)!="t_POL",
		if(dbg!=0,
			print("\nalpha3_do_degenerate_case(): about to check eliminating b1");
		);
		aiArray1=get_A1A2_from_b1(d,absLogA1,hgtA1,absLogA2,hgtA2,absLogA3,hgtA3,logXLB,u1UB,u2UB,u3UB,dbg);
		nUB1a=get_degen_nUB([aiArray1[1],aiArray1[2]],d,nLB,lamUB0,lamUB1,u1UB,logXLB,isComplex,dbg);
		\\ nUB1b is the value by eliminating b1 when u1=0
		nUB1b=get_degen_nUB([aiArray1[3],aiArray1[4]],d,nLB,lamUB0,lamUB1,u1UB,logXLB,isComplex,dbg);
		\\ value from trying to eliminate b1
		nUB1=max(nUB1a,nUB1b);
		if(dbg!=0,
			printf("alpha3_do_degenerate_case(): nUB1=%9.6e, nUB1a=%9.6e, nUB1b=%9.6e\n",nUB1,nUB1a,nUB1b);
		);
	);

	nUB2=0;
	if(type(u2UB)!="t_POL",
		if(dbg!=0,
			print("\nalpha3_do_degenerate_case(): about to check eliminating b2");
		);
		aiArray2=get_A1A2_from_b2(d,absLogA1,hgtA1,absLogA2,hgtA2,absLogA3,hgtA3,logXLB,u1UB,u2UB,u3UB,dbg);
		nUB2a=get_degen_nUB([aiArray2[1],aiArray2[2]],d,nLB,lamUB0,lamUB1,u2UB,logXLB,isComplex,dbg);
		\\ nUB2b is the value by eliminating b2 when u2=0
		nUB2b=get_degen_nUB([aiArray2[3],aiArray2[4]],d,nLB,lamUB0,lamUB1,u2UB,logXLB,isComplex,dbg);
		\\ value from trying to eliminate b2
		nUB2=max(nUB2a,nUB2b);
		if(dbg!=0,
			printf("alpha3_do_degenerate_case(): nUB2=%9.6e, nUB2a=%9.6e, nUB2b=%9.6e\n",nUB2,nUB2a,nUB2b);
		);
	);
	
	nUB3=0;
	if(type(u3UB)!="t_POL",
		if(dbg!=0,
			print("\nalpha3_do_degenerate_case(): about to check eliminating b3");
			print("u3UB=",u3UB,", type(u3UB)=",type(u3UB));
		);
		aiArray3=get_A1A2_from_b3(d,absLogA1,hgtA1,absLogA2,hgtA2,absLogA3,hgtA3,logXLB,u1UB,u2UB,u3UB,dbg);
		nUB3a=get_degen_nUB([aiArray3[1],aiArray3[2]],d,nLB,lamUB0,lamUB1,u3UB,logXLB,isComplex,dbg);
		\\ nUB3b is the value by eliminating b3 when u3=0
		nUB3b=get_degen_nUB([aiArray3[3],aiArray3[4]],d,nLB,lamUB0,lamUB1,u3UB,logXLB,isComplex,dbg);
		\\ value from trying to eliminate b3
		nUB3=max(nUB3a,nUB3b);
		if(dbg!=0,
			printf("alpha3_do_degenerate_case(): nUB3=%9.6e, nUB3a=%9.6e, nUB3b=%9.6e\n",nUB3,nUB3a,nUB3b);
		);
	);

	\\ check our guess that using the b_i for which u_iUB is smallest is best
	if(nUB1*nUB3!=0 && u1UB<u3UB && nUB3<nUB1,
		printf("BAD: in alpha3_do_degenerate_case(), u1UB=%10.6f<u3UB=%10.6f, but nUB1=%10.6f>nUB3=%10.6f\n",u1UB,u3UB,nUB1,nUB3);
		1/0;
	);
	if(nUB1*nUB3!=0 && u3UB<u1UB && nUB1<nUB3,
		printf("BAD: in alpha3_do_degenerate_case(), u3UB=%10.6f<u1UB=%10.6f, but nUB3=%10.6f>nUB1=%10.6f\n",u3UB,u1UB,nUB3,nUB1);
		1/0;
	);
	if(dbg!=0,
		printf("nUB1=%9.6e\n",nUB1);
		printf("nUB2=%9.6e\n",nUB2);
		printf("nUB3=%9.6e\n",nUB3);
		print("alpha3_do_degenerate_case(): END\n");
	);
	\\return(min(nUB1,nUB3));
	return([nUB1,nUB2,nUB3]);
}
