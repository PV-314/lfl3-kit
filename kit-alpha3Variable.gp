\\ \r lfl3\kit-alpha3Variable.gp

read("lfl3\\utils-step2.gp");
read("lfl3\\utils-step3.gp");
read("lfl3\\utils-step4.gp");

\\ assume that \alpha_3 is the variable \alpha_i (i.e., \alpha_1 and \alpha_2 are fixed numbers)

\\ uses:
\\ used for mignotte-eg1 and the Fibonacci perfect powers linear form

\\ returns step3Result=[K, L, m, rho(3logs), chi, R1, R2, S1, T1, T2, nonDegenNUB]
\\ default code where we assume that d1 and d2 in defined in our Notation section are both 1
\\ i.e., no info about gcd(b1,b3)=d1 or gcd(b2,b3)=d2
\\ 11 Dec 2021
alpha3_do_step3(step3Result,d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,chi,bigL,m,rho,nUB,logXLB,nLB,lamUB1,lamUB0,dbg)={
	local(d1,d2);
	
	d1=1;
	d2=1;
	return(alpha3_do_step3_with_d1d2(step3Result,d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,d1,d2,chi,bigL,m,rho,nUB,logXLB,nLB,lamUB1,lamUB0,dbg));
}

\\ returns step3Result=[K, L, m, rho(3logs), chi, R1, R2, S1, T1, T2, nonDegenNUB]
\\ here we use knowledge of d1 and d2. Lower bounds for either of them can be used here
\\ 29 June 2022
alpha3_do_step3_with_d1d2(step3Result,d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,d1,d2,chi,bigL,m,rho,nUB,logXLB,nLB,lamUB1,lamUB0,dbg=0)={
	my(newStep3Result,val);
	
	val=alpha3_do_step3_calcs(d,type(al1),a1,absLogA1,hgtA1,type(al2),a2,absLogA2,hgtA2,type(al3),a3,absLogA3,hgtA3,b1,b2,b3,d1,d2,chi,bigL,m,rho,nUB,logXLB,nLB,lamUB1,lamUB0,dbg);
	newStep3Result=step3_update_min(val,bigL,m,rho,chi,step3Result,dbg);
	return(newStep3Result);
}

\\ returns val=[bigK,bigR1,bigR2,bigS1,bigT1,bigT2,newNonDegenNUB] or else it is empty
\\ here we use knowledge of d1 and d2. Lower bounds for either of them can be used here
\\ 29 June 2022
alpha3_do_step3_calcs(d,al1Type,a1,absLogA1,hgtA1,al2Type,a2,absLogA2,hgtA2,al3Type,a3,absLogA3,hgtA3,b1,b2,b3,d1,d2,chi,bigL,m,rho,nUB,logXLB,nLB,lamUB1,lamUB0,dbg=0)={
	my(a,aP,areMultIndep,bigK,bigR,bigR1,bigR2,bigR3,bigS,bigS1,bigS2,bigS3,bigT,bigT1,bigT2,bigT3,c1,c2,c3,chiV,chiVSqr,cM,eqn42,g,isComplex,isZeroEstOK,lambdaFactor,logBUB,logLambdaLB,rt,vSqr);

	if(al3Type!="t_POL",
		error("ERROR in alpha3_do_step3_calcs(): al3Type=",al3Type," must be t_POL");
	);
	if(al1Type!="t_INT" && al1Type!="t_FRAC" && al1Type!="t_REAL" && al1Type!="t_COMPLEX",
		error("ERROR in alpha3_do_step3_calcs(): al1Type=",al1Type," must be t_INT, t_FRAC, t_COMPLEX or t_REAL");
	);
	if(al2Type!="t_INT" && al2Type!="t_FRAC" && al2Type!="t_REAL" && al2Type!="t_COMPLEX",
		error("ERROR in alpha3_do_step3_calcs(): al2Type=",al2Type," must be t_INT, t_FRAC, t_COMPLEX or t_REAL");
	);
	if(type(b3)!="t_POL",
		error("ERROR in alpha3_do_step3_calcs(): type(b3)=",type(b3)," must be t_POL for the kit to succeed here.");
	);
	if(nLB<1,
		error("ERROR in alpha3_do_step3_calcs(): nLB=",nLB," must be a positive integer");
	);

	\\ al3 not used again, al1 and al2 only used for c2
	a=min(a1,a2);  \\ the assumption here is that a3 will be c*logX and logX large, so a1 or a2, the constant alphas, must be smaller
	aP=max(a1,a2); \\ aPrime=second largest
	if(dbg>0,
		\\print("a1=",a1,", a3=",a3,", a=",a);
		printf("alpha3_do_step3_calcs(): a=%10.6f\n",a);
		printf("alpha3_do_step3_calcs(): aPrime=%10.6f\n",aP);
	);
	
	areMultIndep=are_multiplicatively_independent(al1Type,al2Type);
	c1=get_c1(a,bigL,chi,m);
	c2=get_c2(areMultIndep,a,bigL,m);
	c3=(3*m*m)^(1/3)*bigL;
	if(dbg>0,
		printf("alpha3_do_step3_calcs(): c1=%10.6f\n",c1);
		printf("alpha3_do_step3_calcs(): c2=%10.6f\n",c2);
		printf("alpha3_do_step3_calcs(): c3=%10.6f\n",c3);
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
	if(dbg>0,
		printf("alpha3_do_step3_calcs(): L=%s, m=%s, a1=%s, a2=%s, a3=%s, bigK=%s\n",bigL,m,a1,a2,a3,bigK);
	);
	bigK=(polcoef(bigK,1)+polcoef(bigK,0)/logXLB)*logX;
	if(bigK==0,
		printf("ERROR in alpha3_do_step3_calcs(): K=%s cannot be 0, L=%4d, m=%9.6f, a1=%s, a2=%s, a3=%s\n",bigK,bigL,m,a1,a2,a3);
		error();
	);
	
	\\ use the expression of \Lambda' at the end of our Theorem 2.1
	\\ if a*x*exp(a*x)>A, then log(x)>log(A)-log(a)-log(a*x).
	\\ to be safe, we have a=L*max(R/2, S/2, T/2), x=\Lambda and assume that ax<1
	lambdaFactor=bigL*max(max(subst(bigR,logX,logXLB)/2,subst(bigS,logX,logXLB)/2),subst(bigT,logX,logXLB)/2);
	logLambdaLB=-bigK*bigL*log(rho)-log(lambdaFactor)/logXLB*logX-1;
	logLambdaLB=subst(logLambdaLB,logX,logXLB)/logXLB*logX;
	
	if(dbg>0,
		printf("alpha3_do_step3_calcs(): nUB=%8.6e\n",nUB);
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
	if(dbg>0,
		isZeroEstOK=check_zero_estimate(bigK,bigL,bigR1,bigR2,bigR3,bigS1,bigS2,bigS3,bigT1,bigT2,bigT3,chi,logXLB,areMultIndep,dbg);
		\\if(isZeroEstOK==0,
		\\	return([]);
		\\);
	);

	logBUB=get_logB_UB(bigK,bigR,bigS,bigT,b1,b2,b3,d1,d2,logXLB,nUB,dbg);
	if(dbg>0,
		printf("log(b)<%9.6f\n",logBUB);
	);

	eqn42=get_eqn42(a1,a2,a3,bigK,bigL,bigR,bigS,bigT,d,rho,logBUB,nUB,logXLB,dbg);
	rt=polrootsreal(eqn42)[1];
	if(dbg>0,
		printf("eqn42(x)=0 at x=%9.6f\n",rt);
		printf("subst(eqn42,logX,logXLB)=%9.6f\n",subst(eqn42,logX,logXLB));
	);

	if(rt<logXLB && subst(eqn42,logX,logXLB)>0,
		\\ we only handle lower bounds for lambda that are linear in logX at the moment (May 2022)
		if(poldegree(logLambdaLB,logX)==1 && polcoef(logLambdaLB,0)==0,
			newNonDegenNUB=logLambdaLB/lamUB1;
			newNonDegenNUB=polcoef(newNonDegenNUB,0); \\ need this to make newNonDegenNUB of Pari type real, not a poly
			if(dbg>0,
				printf("newNonDegenNUB=%9.6e\n",newNonDegenNUB);
				print("about to return bigK=",bigK,", bigR1=",bigR1,", bigR2=",bigR2);
				print("bigS1=",bigS1,", bigT1=",bigT1,", bigT2=",bigT2);
			);
			return([bigK,bigR1,bigR2,bigS1,bigT1,bigT2,newNonDegenNUB]);
		);
		if(poldegree(logLambdaLB,logX)!=1 || polcoef(logLambdaLB,0)!=0,
			print("BAD in alpha3_do_step3_calcs(): logLambdaLB must be of form c*logX, but=",logLambdaLB);
			return([]);
		);
	);
	if(dbg>0,
		print("BAD in alpha3_do_step3_calcs(): eqn42 does not hold.");
		print("   eqn42=",eqn42);
		print("   rt=",rt," must be <logXLB=",logXLB);
		print("   and subst(eqn42,logX,logXLB)=",subst(eqn42,logX,logXLB)," must be positive");
		print("   Stopping here.");
	);
	return([]);
}

alpha3_do_step4(step3Result,minNUB,d,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,chi,rho2LB,rho2UB,muLB,muUB,logXLB,nLB,lamUB1,lamUB0,dbg=0)={
	my(newMinNUB,step4Result);
	
	step4Result=alpha3_do_step4_calcs(step3Result,d,type(al1),absLogA1,hgtA1,type(al2),absLogA2,hgtA2,type(al3),absLogA3,hgtA3,chi,rho2LB,rho2UB,muLB,muUB,logXLB,nLB,lamUB1,lamUB0,dbg);
	newMinNUB=step4_update_minNUB(step3Result,minNUB,step4Result,dbg);
	return(newMinNUB);
}

\\ assume the T_i's are constant and R_i's and S_i's are both linear in logZ
\\ also needs to have the alpha_i's
\\ to get a linear form in two logs from the original linear form in three logs
\\ it returns [nUB1,rhoB1,muB1,nUB2,rhoB2,muB2,nUB3,rhoB3,muB3]
\\ one triple of (nUB, rho, mu) for each of the possibilities of eliminating a b_i term
\\ 21 Nov 2021
alpha3_do_step4_calcs(step3Result,d,al1Type,absLogA1,hgtA1,al2Type,absLogA2,hgtA2,al3Type,absLogA3,hgtA3,chi,rho2LB,rho2UB,muLB,muUB,logXLB,nLB,lamUB1,lamUB0,dbg=0)={
	my(a,b1aResult,b1bResult,b2aResult,b2bResult,b3aResult,b3bResult,bigR1,bigS1,bigT1,bigT2,chiV,chiVSqr,cM,isComplex,logAiArray1,logAiArray2,logAiArray3,muB1,muB2,muB3,nUB1,nUB1a,nUB1b,nUB2,nUB2a,nUB2b,nUB3,nUB3a,nUB3b,rhoB1,rhoB2,rhoB3,u1UB,u2UB,u3UB,vSqr);

	if(al3Type!="t_POL",
		print("ERROR in alpha3_do_step4_calcs(): al3Type=",al3Type," must be t_POL");
		error();
	);
	if(al1Type!="t_INT" && al1Type!="t_FRAC" && al1Type!="t_REAL" && al1Type!="t_COMPLEX",
		print("ERROR in alpha3_do_step4_calcs(): al1Type=",al1Type," must be t_INT, t_FRAC, t_COMPLEX or t_REAL");
		error();
	);
	if(al2Type!="t_INT" && al2Type!="t_FRAC" && al2Type!="t_REAL" && al2Type!="t_COMPLEX",
		print("ERROR in alpha3_do_step4_calcs(): al2Type=",al2Type," must be t_INT, t_FRAC, t_COMPLEX or t_REAL");
		error();
	);

	bigR1=step3Result[6];
	bigS1=step3Result[8];
	bigT1=step3Result[9];
	bigT2=step3Result[10];
	if(poldegree(bigR1)!=1,
		print("bad R1 in alpha3_do_step4_calcs(): R1=",bigR1);
		return([]);
	);
	if(poldegree(bigS1)!=1,
		print("bad S1 in alpha3_do_step4_calcs(): S1=",bigS1);
		return([]);
	);
	if(poldegree(bigT1)!=0,
		print("bad T1 in alpha3_do_step4_calcs(): T1=",bigT1);
		return([]);
	);
	if(poldegree(bigT2)!=0,
		print("bad T2 in alpha3_do_step4_calcs(): T2=",bigT2);
		return([]);
	);
	
	if(dbg>0,
		print("\nalpha3_do_step4_calcs(): START");
		print("bigR1=",bigR1,", bigS1=",bigS1,", bigT1=",bigT1,", bigT2=",bigT2);
		printf("(C1) and (C2) bounds: n<max(T1,T2)=%7d\n",max(bigT1,bigT2));
	);

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
	if(dbg>0,
		printf("alpha3_do_step4_calcs(): chi*cV=%9.6f*log(x)\n",polcoef(chiV,1));
		printf("alpha3_do_step4_calcs(): initial cM=%9.6f*log(x)\n",polcoef(cM,1));
	);
	if(polcoef(chiV,1)>polcoef(cM,1),
		cM=chiV;
	);
	if(dbg>0,
		printf("alpha3_do_step4_calcs(): cM=%9.6f*log(x)\n",polcoef(cM,1));
	);

	\\ u1UB
	\\ we assume that S_1>T_1 there
	if(subst(bigS1,logX,logXLB)<bigT1,
		print("ERROR in alpha3_do_step4_calcs(): need S1=",bigS1,">T1=",bigT1,". Increase logXLB.");
		return([]);
	);
	u1UB=(bigS1+1)*(bigT1+1);
	u1UB=subst(u1UB,logX,logXLB)/logXLB*logX;
	u1UB=polcoef(u1UB/(cM-bigS1),0); \\ need this to make bR real, not a poly

	\\ doing u2UB calculation: we assume that R_1>T_1 there
	if(subst(bigR1,logX,logXLB)<bigT1,
		print("ERROR in alpha3_do_step4_calcs(): need R1=",bigR1,">R1=",bigR1,". Increase logXLB.");
		return([]);
	);
	u2UB=(bigR1+1)*(bigT1+1);
	u2UB=subst(u2UB,logX,logXLB)/logXLB*logX;
	u2UB=polcoef(u2UB/(cM-bigR1),0); \\ need this to make bT real, not a poly (otherwise, floor() and other functions do not work)

	\\ u3UB calc:
	u3UB=(bigR1+1)*(bigS1+1);
	u3UB=subst(u3UB,logX,logXLB)/logXLB^2*logX^2;
	if(subst(bigR1,logX,logXLB)<subst(bigS1,logX,logXLB),
		u3UB=u3UB/(cM-bigS1);
	);
	if(subst(bigS1,logX,logXLB)<subst(bigR1,logX,logXLB),
		u3UB=u3UB/(cM-bigR1);
	);

	if(dbg!=0,
		printf("alpha3_do_step4_calcs(): before rounding down, u1UB=%10.6f\n",u1UB); \\ was B_R, bR
		printf("alpha3_do_step4_calcs(): before rounding down, u2UB=%10.6f\n",u2UB); \\ was B_S, bS
	);

	u1UB=floor(u1UB);
	u2UB=floor(u2UB);
	if(dbg!=0,
		printf("alpha3_do_step4_calcs(): u1UB=%10.6f\n",u1UB); \\ was B_R, bR
		printf("alpha3_do_step4_calcs(): u2UB=%10.6f\n",u2UB); \\ was B_S, bS
		printf("alpha3_do_step4_calcs(): u3UB=%10.6f*logX\n",polcoef(u3UB,1)); \\ was B_S, bS
	);

	isComplex=is_complex(al1Type,al2Type,al3Type);
	nUB1=0;
	if(type(u1UB)!="t_POL",
		if(dbg>0,
			print("\nalpha3_do_step4_calcs(): about to check eliminating b1");
		);
		\\ it returns the array [hgtNewA1,logNewA1,hgtNewA2,logNewA2,hgtNewA1_0,logNewA1_0,hgtNewA2_0,logNewA2_0]
		logAiArray1=get_logA1A2_from_b1(d,absLogA1,hgtA1,absLogA2,hgtA2,absLogA3,hgtA3,isComplex,logXLB,u1UB,u2UB,u3UB,dbg);
		\\ get_degen_nUB() returns [rho2Log,max(nUBMultDep,nUBMultIndep)]
		if(dbg>0,
			print("\nabout to call get_degen_nUB() for eliminating b1 when u1!=0");
		);
		b1aResult=get_degen_nUB(logAiArray1[1],logAiArray1[2],logAiArray1[3],logAiArray1[4],d,nLB,lamUB0,lamUB1,u1UB,rho2LB,rho2UB,muLB,muUB,logXLB,isComplex,dbg);
		\\ b1bResult is the value by eliminating b1 when u1=0
		if(dbg>0,
			print("\nabout to call get_degen_nUB() for eliminating b1 when u1=0");
		);
		b1bResult=get_degen_nUB(logAiArray1[5],logAiArray1[6],logAiArray1[7],logAiArray1[8],d,nLB,lamUB0,lamUB1,u1UB,rho2LB,rho2UB,muLB,muUB,logXLB,isComplex,dbg);
		nUB1=b1aResult[1];
		rhoB1=b1aResult[2];
		if(length(b1aResult)>2,
			muB1=b1aResult[3];
		);
		if(b1bResult[1]>nUB1,
			nUB1=b1bResult[1];
			rhoB1=b1bResult[2];
			if(length(b1bResult)>2,
				muB1=b1bResult[3];
			);
		);
		if(dbg>0,
			printf("alpha3_do_step4_calcs(): nUB1=%9.6e, rhoB1=%9.6f, muB1=%9.6f, nUB1a=%9.6e, nUB1b=%9.6e\n",nUB1,rhoB1,muB1,b1aResult[1],b1bResult[1]);
		);
	);

	nUB2=0;
	if(type(u2UB)!="t_POL",
		if(dbg>0,
			print("\nalpha3_do_step4_calcs(): about to check eliminating b2");
		);
		logAiArray2=get_logA1A2_from_b2(d,absLogA1,hgtA1,absLogA2,hgtA2,absLogA3,hgtA3,isComplex,logXLB,u1UB,u2UB,u3UB,dbg);
		\\ get_degen_nUB() returns [rho2Log,max(nUBMultDep,nUBMultIndep)]
		if(dbg>0,
			print("\nabout to call get_degen_nUB() for eliminating b2 when u2!=0");
		);
		b2aResult=get_degen_nUB(logAiArray2[1],logAiArray2[2],logAiArray2[3],logAiArray2[4],d,nLB,lamUB0,lamUB1,u2UB,rho2LB,rho2UB,muLB,muUB,logXLB,isComplex,dbg);
		\\ b2bResult is the value by eliminating b2 when u2=0
		if(dbg>0,
			print("\nabout to call get_degen_nUB() for eliminating b2 when u2=0");
		);
		b2bResult=get_degen_nUB(logAiArray2[5],logAiArray2[6],logAiArray2[7],logAiArray2[8],d,nLB,lamUB0,lamUB1,u2UB,rho2LB,rho2UB,muLB,muUB,logXLB,isComplex,dbg);
		rhoB2=b2aResult[2];
		nUB2=b2aResult[1];
		if(length(b2aResult)>2,
			muB2=b2aResult[3];
		);
		if(b2bResult[1]>nUB2,
			rhoB2=b2bResult[2];
			nUB2=b2bResult[1];
			if(length(b2bResult)>2,
				muB2=b2bResult[3];
			);
		);
		if(dbg>0,
			printf("alpha3_do_step4_calcs(): nUB2=%9.6e, rhoB2=%9.6f, muB2=%9.6f, nUB2a=%9.6e, nUB2b=%9.6e\n",nUB2,rhoB2,muB2,b2aResult[1],b2bResult[1]);
		);
	);
	
	nUB3=0;
	if(type(u3UB)!="t_POL",
		if(dbg>0,
			print("\nalpha3_do_step4_calcs(): about to check eliminating b3");
			print("u3UB=",u3UB,", type(u3UB)=",type(u3UB));
		);
		logAiArray3=get_logA1A2_from_b3(d,absLogA1,hgtA1,absLogA2,hgtA2,absLogA3,hgtA3,isComplex,logXLB,u1UB,u2UB,u3UB,dbg);
		\\ get_degen_nUB() returns [rho2Log,max(nUBMultDep,nUBMultIndep)]
		if(dbg>0,
			print("\nabout to call get_degen_nUB() for eliminating b3 when u3!=0");
		);
		b3aResult=get_degen_nUB(logAiArray3[1],logAiArray3[2],logAiArray3[3],logAiArray3[4],d,nLB,lamUB0,lamUB1,u3UB,rho2LB,rho2UB,muLB,muUB,logXLB,isComplex,dbg);
		\\ b3bResult is the value by eliminating b3 when u3=0
		if(dbg>0,
			print("\nabout to call get_degen_nUB() for eliminating b3 when u3=0");
		);
		b3bResult=get_degen_nUB(logAiArray3[5],logAiArray3[6],logAiArray3[7],logAiArray3[8],d,nLB,lamUB0,lamUB1,u3UB,rho2LB,rho2UB,muLB,muUB,logXLB,isComplex,dbg);
		nUB3=b3aResult[1];
		rhoB3=b3aResult[2];
		if(length(b3aResult)>2,
			muB3=b3aResult[3];
		);
		if(b3bResult[1]>nUB3,
			nUB3=b3bResult[1];
			rhoB3=b3bResult[2];
			if(length(b3bResult)>2,
				muB3=b3bResult[3];
			);
		);
		if(dbg>0,
			printf("alpha3_do_step4_calcs(): nUB3=%9.6e, rhoB3=%9.6f, muB3=%9.6f, nUB3a=%9.6e, nUB3b=%9.6e\n",nUB3,rhoB3,muB3,b3aResult[1],b3bResult[1]);
		);
	);

	\\ check our guess that using the b_i for which u_iUB is smallest is best
	if(dbg>0 && u2UB<u1UB && 1.001*nUB1<nUB2,
		printf("ODD in alpha3_do_step4_calcs(): u2UB=%10.6f<u1UB=%10.6f, but nUB2=%10.6f>nUB1=%10.6f\n",u2UB,u1UB,nUB2,nUB1);
	);
	if(dbg>0 && u1UB<u2UB && 1.001*nUB2<nUB1,
		printf("ODD in alpha3_do_step4_calcs(): u1UB=%10.6f<u2UB=%10.6f, but nUB1=%10.6f>nUB2=%10.6f\n",u1UB,u2UB,nUB1,nUB2);
	);
	if(dbg>0,
		printf("nUB1=%9.6e\n",nUB1);
		printf("nUB2=%9.6e\n",nUB2);
		printf("nUB3=%9.6e\n",nUB3);
		print("alpha3_do_step4_calcs(): END\n");
	);
	return([nUB1,rhoB1,muB1,nUB2,rhoB2,muB2,nUB3,rhoB3,muB3]);
}