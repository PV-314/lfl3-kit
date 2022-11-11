\\ \r lfl3\lfl-utils-alpha1Variable.gp

read("lfl3\\lfl-utils-general.gp");

\\ assume that \alpha_1 is the variable \alpha_i (i.e., \alpha_2 and \alpha_3 are fixed numbers)

\\ uses:
\\ can be used for mignotte-eg2

\\ 4 Jan 2022
alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,dbg=0)={
	local(d1,d2);
	
	d1=1;
	d2=1;
	alpha1_check_params_with_d1d2(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,d1,d2,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,dbg);
}

\\ 29 June 2022
alpha1_check_params_with_d1d2(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,d1,d2,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,dbg=0)={
	my(a,aP,bigK,bigR,bigR1,bigR2,bigR3,bigS,bigS1,bigS2,bigS3,bigT,bigT1,bigT2,bigT3,c1,c2,c3,chiV,chiVSqr,cM,eqn42LHS,eqn42RHS,f,g,isComplex,logBUB,logLambdaLB,rt,two13,vSqr);

	if(type(al1)!="t_POL",
		print("ERROR: type(al1)=",type(al1)," must be t_POL");
		return([]);
	);
	if(type(al2)!="t_INT" && type(al2)!="t_FRAC" && type(al2)!="t_REAL" && type(al2)!="t_COMPLEX",
		print("ERROR: type(al2)=",type(al2)," must be t_INT, t_FRAC, t_COMPLEX or t_REAL");
		return([]);
	);
	if(type(al3)!="t_INT" && type(al3)!="t_FRAC" && type(al3)!="t_REAL" && type(al3)!="t_COMPLEX",
		print("ERROR: type(al3)=",type(al3)," must be t_INT, t_FRAC, t_COMPLEX or t_REAL");
		return([]);
	);
	
	\\ al1 not used again, al2 and al3 only used for c2
	a=min(a2,a3); \\ the assumption here is that a1 will be c*logX and logX large, so a2 or a3, the constant alphas, must be smaller
	\\print("a1=",a1,", a3=",a3,", a=",a);
	aP=max(a2,a3); \\ aPrime=second largest
	if(dbg!=0,
		print("a=",a);
		print("aPrime=",aP);
	);
	
	two13=1.2599210498948731647672106072782283506; \\ 2^(1/3) to save recalculation every time
	c1=two13;
	c1=max(c1,(chi*m*bigL)^(2/3));
	c1=max(c1,sqrt(2*m*bigL/a));
	if(type(al2)=="t_REAL" && type(al3)=="t_REAL",
		c2=two13*(m*bigL)^(2/3);
	);
	if(type(al2)!="t_REAL" || type(al3)!="t_REAL",
		c2=max(two13*(m*bigL)^(2/3),sqrt(m/a)*bigL);
	);
	c3=(6*m*m)^(1/3)*bigL;
	if(dbg!=0,
		printf("in alpha1_check_params(): c1=%10.6f\n",c1);
		printf("in alpha1_check_params(): c2=%10.6f\n",c2);
		printf("in alpha1_check_params(): c3=%10.6f\n",c3);
	);
	
	bigR1=c1*a2*a3;
	bigR1=floor(bigR1);
	bigR2=c2*a2*a3;
	bigR2=floor(bigR2);
	bigR3=c3*a2*a3;
	bigR3=floor(bigR3);
	bigR=bigR1+bigR2+bigR3+1;

	bigS1=c1*a1*a3;
	bigS1=(polcoef(bigS1,1)+polcoef(bigS1,0)/logXLB)*logX;
	bigS2=c2*a1*a3;
	bigS2=(polcoef(bigS2,1)+polcoef(bigS2,0)/logXLB)*logX;
	bigS3=c3*a1*a3;
	bigS3=(polcoef(bigS3,1)+polcoef(bigS3,0)/logXLB)*logX;
	bigS=bigS1+bigS2+bigS3+1;
	bigS=(polcoef(bigS,1)+polcoef(bigS,0)/logXLB)*logX;

	bigT1=c1*a1*a2;
	bigT1=(polcoef(bigT1,1)+polcoef(bigT1,0)/logXLB)*logX;
	bigT2=c2*a1*a2;
	bigT2=(polcoef(bigT2,1)+polcoef(bigT2,0)/logXLB)*logX;
	bigT3=c3*a1*a2;
	bigT3=(polcoef(bigT3,1)+polcoef(bigT3,0)/logXLB)*logX;
	bigT=bigT1+bigT2+bigT3+1;
	bigT=(polcoef(bigT,1)+polcoef(bigT,0)/logXLB)*logX;

	bigK=bigL*m*a1*a2*a3;
	bigK=(polcoef(bigK,1)+polcoef(bigK,0)/logXLB)*logX;
	logLambdaLB=-bigK*bigL*log(rho); \\-log(bigK*bigL);
	if(dbg!=0,
		printf("in alpha1_check_params(): nUB=%8.6e\n",nUB);
		printf("bigR1=%8d\n",bigR1);
		printf("bigR2=%8d\n",bigR2);
		printf("bigR3=%8d\n",bigR3);

		printf("bigS1=%11.6f*logX\n",polcoef(bigS1,1));
		printf("bigS2=%11.6f*logX\n",polcoef(bigS2,1));
		printf("bigS3=%11.6f*logX\n",polcoef(bigS3,1));
		
		printf("bigT1=%11.6f*logX\n",polcoef(bigT1,1));
		printf("bigT2=%11.6f*logX\n",polcoef(bigT2,1));
		printf("bigT3=%11.6f*logX\n",polcoef(bigT3,1));
		printf("rho=%10.6f\n",rho);
		printf("bigK=%s\n",bigK);
		printf("bigL=%s\n",bigL);

		printf("bigR=%10.6f\n",bigR);
		printf("bigS=%10.6f*logX\n",polcoef(bigS,1));
		printf("bigT=%10.6f*logX\n",polcoef(bigT,1));
		printf("log |Lambda|>%9.6e*logX\n\n",polcoef(logLambdaLB,1));
	);

	logBUB=get_logBPrime_UB(bigK,bigR,bigS,bigT,b1,b2,b3,d1,d2,logXLB,nUB,dbg);
	if(dbg!=0,
		printf("log(b')<%9.6f\n",logBUB);
	);

	\\ with K_0=2(K-1):
	\\eqn42LHS=(bigK*bigL/2+bigL/4-1-2*bigK/bigL)*log(rho);
	\\ with K_0=K-1, rather than K_0=2(K-1)
	eqn42LHS=(bigK*bigL-bigK-bigK/3/bigL)*log(rho);
	\\eqn42LHS=(polcoef(eqn42LHS,1)+polcoef(eqn42LHS,0)/logXLB)*logX;

	gDenom=12*bigR*bigS*bigT/logX/logX;
	gDenom=subst(gDenom,logX,logXLB)*logX*logX;
	gNumer=bigK*bigK*bigL;
	gNumer=polcoef(gNumer,2,logX)*logX*logX;
	g=1/4-gNumer/gDenom;
	if(dbg!=0,
		printf("g=%8.6f\n",g);
	);
	
	\\ first the (cD+1)*log(N) term
	eqn42RHS=(d+1)*log(polcoef(bigK,1)*polcoef(bigK,1)*bigL)+(d+1)*log(logXLB)/logXLB*logX;
	\\ next the gL(a1*R+a2*S+a3*T) term
	eqn42RHS=eqn42RHS+g*bigL*(a1*bigR+a2*bigS+a3*bigT);
	\\ lastly the D(K-1)*log(b) term
	eqn42RHS=eqn42RHS+d*(bigK-1)*logBUB;
	eqn42RHS=subst(eqn42RHS,logN,log(nUB));
	eqn42=eqn42LHS-eqn42RHS;
	eqn42=pollead(polcoef(eqn42,1))*logX+pollead(polcoef(eqn42,0));
	rt=polrootsreal(eqn42)[1];
	if(dbg!=0,
		printf("eqn42LHS=%s\n",eqn42LHS);
		printf("eqn42RHS D(K-1)*log(b) term=%s\n",d*(bigK-1)*logBUB);
		printf("eqn42RHS gL(a1*R+a2*S+a3*T) term=%s\n",g*bigL*(a1*bigR+a2*bigS+a3*bigT));
		printf("eqn42RHS (cD+1)*log(N) term=%s\n",(d+1)*log(polcoef(bigK,1)*polcoef(bigK,1)*bigL)+(d+1)*log(logXLB)/logXLB*logX);
		printf("eqn42RHS=%s\n",eqn42RHS);
		printf("eqn42LHS-eqn42RHS=%s\n",eqn42);
		printf("eqn42(x)=0 at x=%9.6f\n",rt);
		printf("subst(eqn42,logX,logXLB)=%9.6f\n",subst(eqn42,logX,logXLB));
	);

	if(rt<logXLB && subst(eqn42,logX,logXLB)>0,
		\\ we only handle lower bounds that are linear in logX at the moment (May 2022)
		if(poldegree(logLambdaLB,logX)==1 && polcoef(logLambdaLB,0)==0,
			newNonDegenNUB=logLambdaLB/lamUB1;
			newNonDegenNUB=polcoef(newNonDegenNUB,0); \\ need this to make newNonDegenNUB real, not a poly
			if(dbg!=0,
				printf("newNonDegenNUB=%9.6e\n",newNonDegenNUB);
			);
			isComplex = type(al1)=="t_COMPLEX" || type(al2)=="t_COMPLEX" || type(al3)=="t_COMPLEX";
			newDegenNUB=alpha1_do_degenerate_case(d,hgtA1,absLogA1,hgtA2,absLogA2,hgtA3,absLogA3,bigR1,bigR2,bigS1,bigT1,chi,logXLB,nLB,lamUB1,lamUB0,isComplex,dbg);
			if(dbg!=0,
				printf("newDegenNUB=%9.6e\n",newDegenNUB);
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
		print("eqn42 does not hold. Stopping here.");
	);
	return([]);
}

\\ assume the R_i's are constant and S_i's and T_i's are both linear in logZ
\\ also needs to have the alpha_i's
\\ 21 Nov 2021
alpha1_do_degenerate_case(d,hgtA1,absLogA1,hgtA2,absLogA2,hgtA3,absLogA3,bigR1,bigR2,bigS1,bigT1,chi,logXLB,nLB,lamUB1,lamUB0,isComplex,dbg=0)={
	my(a,aiArray2,aiArray3,chiV,chiVSqr,cM,nUB2,nUB2a,nUB2b,nUB3,nUB3a,nUB3b,nUB,u1UB,u2UB,u3UB,vSqr);

	if(poldegree(bigR1)!=0,
		print("bad R1 in do_generate_case(): R1=",bigR1);
		return([]);
	);
	if(poldegree(bigS1)!=1,
		print("bad S1 in do_generate_case(): S1=",bigS1);
		return([]);
	);
	if(poldegree(bigT1)!=1,
		print("bad T1 in eg2_do_generate_case(): T1=",bigT1);
		return([]);
	);
	
	if(dbg!=0,
		print("\nstarting degenerate case");
		printf("(C1) and (C2) bounds: n<max(R1,R2)=%7d\n",max(bigR1,bigR2));
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
		printf("chiV=%9.6f*log(x)\n",polcoef(chiV,1));
		printf("initial cM=%9.6f*log(x)\n",polcoef(cM,1));
	);
	if(polcoef(chiV,1)>polcoef(cM,1),
		cM=chiV;
	);
	if(dbg!=0,
		printf("cM=%9.6f*log(x)\n",polcoef(cM,1));
	);

	u1UB=(bigS1+1)*(bigT1+1);
	u1UB=subst(u1UB,logX,logXLB)/logXLB^2*logX^2;
	if(subst(bigS1,logX,logXLB)<subst(bigT1,logX,logXLB),
		u1UB=u1UB/(cM-bigT1);
	);
	if(subst(bigT1,logX,logXLB)<subst(bigS1,logX,logXLB),
		u1UB=u1UB/(cM-bigS1);
	);
	
	\\ we assume that S_1>R_1 here for u3UB calcs
	if(subst(bigS1,logX,logXLB)<bigR1,
		return([]);
	);
	u3UB=(bigR1+1)*(bigS1+1);
	u3UB=subst(u3UB,logX,logXLB)/logXLB*logX;
	u3UB=polcoef(u3UB/(cM-bigS1),0); \\ need this to make u3UB real, not a poly
	u3UB=floor(u3UB);

	\\ we assume that T_1>R_1 here for u2UB calcs
	if(subst(bigT1,logX,logXLB)<bigR1,
		return([]);
	);
	u2UB=(bigR1+1)*(bigT1+1);
	u2UB=polcoef(u2UB/(cM-bigT1),0); \\ need this to make u2UB real, not a poly
	u2UB=floor(u2UB);

	if(dbg!=0,
		printf("u1UB=%10.6f*logX\n",polcoef(u1UB,1)); \\ bR
		printf("u2UB=%10.6f: must have logX>%10.6f(check value)\n",u2UB,u3UB/polcoef(u1UB,1)); \\ bS
		printf("u3UB=%10.6f\n",u3UB); \\ bT
	);
	
	nUB1=0;
	if(type(u1UB)!="t_POL",
		if(dbg!=0,
			print("\nalpha1_do_degenerate_case(): about to check eliminating b1");
		);
		aiArray1=get_A1A2_from_b1(d,absLogA1,hgtA1,absLogA2,hgtA2,absLogA3,hgtA3,logXLB,u1UB,u2UB,u3UB,dbg);
		nUB1a=get_degen_nUB([aiArray1[1],aiArray1[2]],d,nLB,lamUB0,lamUB1,u1UB,logXLB,isComplex,dbg);
		\\ nUB1b is the value by eliminating b1 when u1=0
		nUB1b=get_degen_nUB([aiArray1[3],aiArray1[4]],d,nLB,lamUB0,lamUB1,u1UB,logXLB,isComplex,dbg);
		nUB1=max(nUB1a,nUB1b);
		if(dbg!=0,
			printf("alpha1_do_degenerate_case(): nUB1=%9.6e, nUB1a=%9.6e, nUB1b=%9.6e\n",nUB1,nUB1a,nUB1b);
		);
	);

	nUB2=0;
	if(type(u2UB)!="t_POL",
		if(dbg!=0,
			print("\nalpha1_do_degenerate_case(): about to check eliminating b2");
		);
		aiArray2=get_A1A2_from_b2(d,absLogA1,hgtA1,absLogA2,hgtA2,absLogA3,hgtA3,logXLB,u1UB,u2UB,u3UB,dbg);
		nUB2a=get_degen_nUB([aiArray2[1],aiArray2[2]],d,nLB,lamUB0,lamUB1,u2UB,logXLB,isComplex,dbg);
		\\ nUB2b is the value by eliminating b2 when u2=0
		nUB2b=get_degen_nUB([aiArray2[3],aiArray2[4]],d,nLB,lamUB0,lamUB1,u2UB,logXLB,isComplex,dbg);
		nUB2=max(nUB2a,nUB2b);
		if(dbg!=0,
			printf("alpha1_do_degenerate_case(): nUB2=%9.6e, nUB2a=%9.6e, nUB2b=%9.6e\n",nUB2,nUB2a,nUB2b);
		);
	);
	
	nUB3=0;
	if(type(u3UB)!="t_POL",
		if(dbg!=0,
			print("\nalpha1_do_degenerate_case(): about to check eliminating b3");
			print("u3UB=",u3UB,", type(u3UB)=",type(u3UB));
		);
		aiArray3=get_A1A2_from_b3(d,absLogA1,hgtA1,absLogA2,hgtA2,absLogA3,hgtA3,logXLB,u1UB,u2UB,u3UB,dbg);
		nUB3a=get_degen_nUB([aiArray3[1],aiArray3[2]],d,nLB,lamUB0,lamUB1,u3UB,logXLB,isComplex,dbg);
		\\ nUB3b is the value by eliminating b3 when u3=0
		nUB3b=get_degen_nUB([aiArray3[3],aiArray3[4]],d,nLB,lamUB0,lamUB1,u3UB,logXLB,isComplex,dbg);
		\\ value from trying to eliminate b3
		nUB3=max(nUB3a,nUB3b);
		if(dbg!=0,
			printf("alpha1_do_degenerate_case(): nUB3=%9.6e, nUB3a=%9.6e, nUB3b=%9.6e\n",nUB3,nUB3a,nUB3b);
		);
	);

	\\ check our guess that using the b_i for which u_iUB is smallest is best
	if(u2UB<u3UB && nUB3<nUB2,
		printf("BAD: in alpha3_do_degenerate_case(), u2UB=%10.6f<u3UB=%10.6f, but nUB2=%10.6f>nUB3=%10.6f\n",u2UB,u3UB,nUB2,nUB3);
		1/0;
	);
	if(u3UB<u2UB && nUB2<nUB3,
		printf("BAD: in alpha3_do_degenerate_case(), u3UB=%10.6f<u2UB=%10.6f, but nUB3=%10.6f>nUB2=%10.6f\n",u3UB,u2UB,nUB3,nUB2);
		1/0;
	);
	if(dbg!=0,
		printf("nUB1=%9.6e\n",nUB1);
		printf("nUB2=%9.6e\n",nUB2);
		printf("nUB3=%9.6e\n",nUB3);
		print("alpha1_do_degenerate_case(): END\n");
	);
	return([nUB1,nUB2,nUB3]);
}
