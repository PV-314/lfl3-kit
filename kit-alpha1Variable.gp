\\ \r lfl3\kit-alpha1Variable.gp

read("lfl3\\utils-step2.gp");
read("lfl3\\utils-step3.gp");
read("lfl3\\utils-step4.gp");

\\ assume that \alpha_1 is the variable \alpha_i (i.e., \alpha_2 and \alpha_3 are fixed numbers)

\\ uses:
\\ be used for mignotte-eg2

\\ returns step3Result=[K, L, m, rho(3logs), chi, R1, R2, S1, T1, T2, nonDegenNUB]
\\ default code where we assume that d1 and d2 in defined in our Notation section are both 1
\\ i.e., no info about gcd(b1,b3)=d1 or gcd(b2,b3)=d2
\\ 4 Jan 2022
alpha1_do_step3(step3Result,d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,chi,bigL,m,rho,nUB,logXLB,nLB,lamUB1,lamUB0,dbg)={
	local(d1,d2);
	
	d1=1;
	d2=1;
	return(alpha1_do_step3_with_d1d2(step3Result,d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,d1,d2,chi,bigL,m,rho,nUB,logXLB,nLB,lamUB1,lamUB0,dbg));
}

\\ returns step3Result=[K, L, m, rho(3logs), chi, R1, R2, S1, T1, T2, nonDegenNUB]
\\ here we use knowledge of d1 and d2. Lower bounds for either of them can be used here
\\ 29 June 2022
alpha1_do_step3_with_d1d2(step3Result,d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,d1,d2,chi,bigL,m,rho,nUB,logXLB,nLB,lamUB1,lamUB0,dbg=0)={
	my(newStep3Result,val);
	
	val=alpha1_do_step3_calcs(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,d1,d2,chi,bigL,m,rho,nUB,logXLB,nLB,lamUB1,lamUB0,dbg);
	newStep3Result=step3_update_min(val,bigL,m,rho,chi,step3Result,dbg);
	return(newStep3Result);
}

\\ returns val=[bigK,bigR1,bigR2,bigS1,bigT1,bigT2,newNonDegenNUB] or else it is empty
\\ here we use knowledge of d1 and d2. Lower bounds for either of them can be used here
\\ 29 June 2022
alpha1_do_step3_calcs(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,d1,d2,chi,bigL,m,rho,nUB,logXLB,nLB,lamUB1,lamUB0,dbg=0)={
	my(a,aP,areMultIndep,bigK,bigR,bigR1,bigR2,bigR3,bigS,bigS1,bigS2,bigS3,bigT,bigT1,bigT2,bigT3,c1,c2,c3,chiV,chiVSqr,cM,eqn42,g,isComplex,isZeroEstOK,logBUB,logLambdaLB,rt,vSqr);

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
	if(type(b1)!="t_POL",
		print("ERROR: type(b1)=",type(b1)," must be t_POL for the kit to succeed here.");
		return([]);
	);
	if(nLB<1,
		print("ERROR: nLB=",nLB," must be a positive integer");
		return([]);
	);
	
	\\ al1 not used again, al2 and al3 only used for c2
	a=min(a2,a3); \\ the assumption here is that a1 will be c*logX and logX large, so a2 or a3, the constant alphas, must be smaller
	aP=max(a2,a3); \\ aPrime=second largest
	if(dbg>0,
		\\print("a1=",a1,", a3=",a3,", a=",a);
		printf("alpha1_do_step3_calcs(): a=%10.6f\n",a);
		printf("alpha1_do_step3_calcs(): aPrime=%10.6f\n",aP);
	);
	
	areMultIndep=are_multiplicatively_independent(al2,al3);
	c1=get_c1(a,bigL,chi,m);
	c2=get_c2(areMultIndep,a,bigL,m);
	c3=(3*m*m)^(1/3)*bigL;
	if(dbg>0,
		printf("alpha1_do_step3_calcs(): c1=%10.6f\n",c1);
		printf("alpha1_do_step3_calcs(): c2=%10.6f\n",c2);
		printf("alpha1_do_step3_calcs(): c3=%10.6f\n",c3);
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
	if(dbg>0,
		printf("alpha1_do_step3_calcs(): L=%s, m=%s, a1=%s, a2=%s, a3=%s, bigK=%s\n",bigL,m,a1,a2,a3,bigK);
	);
	bigK=(polcoef(bigK,1)+polcoef(bigK,0)/logXLB)*logX;
	if(bigK==0,
		printf("alpha1_do_step3_calcs(): L=%s, m=%s, a1=%s, a2=%s, a3=%s, bigK=%s\n",bigL,m,a1,a2,a3,bigK);
		1/0;
	);
	logLambdaLB=-bigK*bigL*log(rho);
	if(dbg>0,
		printf("alpha1_do_step3_calcs(): nUB=%8.6e\n",nUB);
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
		printf("bigK=%10.6f*logX\n",polcoef(bigK,1));
		printf("bigL=%6d\n",bigL);

		printf("bigR=%10.6f\n",bigR);
		printf("bigS=%10.6f*logX\n",polcoef(bigS,1));
		printf("bigT=%10.6f*logX\n",polcoef(bigT,1));
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
			print("BAD in alpha1_do_step3_calcs(): logLambdaLB must be of form c*logX, but="+logLambdaLB);
			return([]);
		);
	);
	if(dbg>0,
		print("\nBAD in alpha1_do_step3_calcs(): eqn42 does not hold.");
		print("   eqn42=",eqn42);
		print("   rt=",rt," must be <logXLB=",logXLB);
		print("   and subst(eqn42,logX,logXLB)=",subst(eqn42,logX,logXLB)," must be positive");
		print("   Stopping here.");
	);
	return([]);
}

alpha1_do_step4(step3Result,minNUB,d,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,chi,rho2LB,rho2UB,muLB,muUB,logXLB,nLB,lamUB1,lamUB0,dbg=0)={
	my(newMinNUB,step4Result);
	
	step4Result=alpha1_do_step4_calcs(step3Result,d,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,chi,rho2LB,rho2UB,muLB,muUB,logXLB,nLB,lamUB1,lamUB0,dbg);
	newMinNUB=step4_update_minNUB(step3Result,minNUB,step4Result,dbg);
	return(newMinNUB);
}

\\ assume the T_i's are constant and R_i's and S_i's are both linear in logZ
\\ also needs to have the alpha_i's
\\ to get a linear form in two logs from the original linear form in three logs
\\ it returns [nUB1,rhoB1,muB1,nUB2,rhoB2,muB2,nUB3,rhoB3,muB3]
\\ one triple of (nUB, rho, mu) for each of the possibilities of eliminating a b_i term
\\ 21 Nov 2021
alpha1_do_step4_calcs(step3Result,d,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,chi,rho2LB,rho2UB,muLB,muUB,logXLB,nLB,lamUB1,lamUB0,dbg=0)={
	my(a,b1aResult,b1bResult,b2aResult,b2bResult,b3aResult,b3bResult,bigR1,bigR2,bigS1,bigT1,chiV,chiVSqr,cM,isComplex,logAiArray1,logAiArray2,logAiArray3,muB1,muB2,muB3,nUB1,nUB1a,nUB1b,nUB2,nUB2a,nUB2b,nUB3,nUB3a,nUB3b,rhoB1,rhoB2,rhoB3,u1UB,u2UB,u3UB,vSqr);

	if(type(al1)!="t_POL",
		print("ERROR in alpha1_do_step4_calcs(): type(al1)=",type(al1)," must be t_POL");
		return([]);
	);
	if(type(al2)!="t_INT" && type(al2)!="t_FRAC" && type(al2)!="t_REAL" && type(al2)!="t_COMPLEX",
		print("ERROR in alpha1_do_step4_calcs(): type(al2)=",type(al2)," must be t_INT, t_FRAC, t_COMPLEX or t_REAL");
		return([]);
	);
	if(type(al3)!="t_INT" && type(al3)!="t_FRAC" && type(al3)!="t_REAL" && type(al3)!="t_COMPLEX",
		print("ERROR in alpha1_do_step4_calcs(): type(al3)=",type(al3)," must be t_INT, t_FRAC, t_COMPLEX or t_REAL");
		return([]);
	);

	bigR1=step3Result[6];
	bigR2=step3Result[7];
	bigS1=step3Result[8];
	bigT1=step3Result[9];
	if(poldegree(bigR1)!=0,
		print("bad R1 in alpha1_do_step4_calcs(): R1=",bigR1);
		return([]);
	);
	if(poldegree(bigR2)!=0,
		print("bad R2 in alpha1_do_step4_calcs(): R2=",bigR2);
		return([]);
	);
	if(poldegree(bigS1)!=1,
		print("bad S1 in alpha1_do_step4_calcs(): S1=",bigS1);
		return([]);
	);
	if(poldegree(bigT1)!=1,
		print("bad T1 in alpha1_do_step4_calcs(): T1=",bigT1);
		return([]);
	);
	
	if(dbg>0,
		print("\nalpha1_do_step4_calcs(): START");
		print("bigR1=",bigR1,", bigR2=",bigR2,", bigS1=",bigS1,", bigT1=",bigT1);
		printf("(C1) and (C2) bounds: n<max(R1,R2)=%7d\n",max(bigR1,bigR2));
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
		printf("alpha1_do_step4_calcs(): chi*cV=%9.6f*log(x)\n",polcoef(chiV,1));
		printf("alpha1_do_step4_calcs(): initial cM=%9.6f*log(x)\n",polcoef(cM,1));
	);
	if(polcoef(chiV,1)>polcoef(cM,1),
		cM=chiV;
	);
	if(dbg>0,
		printf("alpha1_do_step4_calcs(): cM=%9.6f*log(x)\n",polcoef(cM,1));
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

	if(dbg>0,
		printf("alpha1_do_step4_calcs(): u1UB=%10.6f*logX\n",polcoef(u1UB,1)); \\ bR
		printf("alpha1_do_step4_calcs(): u2UB=%10.6f: must have logX>%10.6f(check value)\n",u2UB,u3UB/polcoef(u1UB,1)); \\ bS
		printf("alpha1_do_step4_calcs(): u3UB=%10.6f\n",u3UB); \\ bT
	);
	
	isComplex=is_complex(al1,al2,al3);
	nUB1=0;
	if(type(u1UB)!="t_POL",
		if(dbg>0,
			print("\nalpha1_do_step4_calcs(): about to check eliminating b1");
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
		if(b1bResult[2]>nUB1,
			nUB1=b1bResult[1];
			rhoB1=b1bResult[2];
			if(length(b1bResult)>2,
				muB1=b1bResult[3];
			);
		);
		if(dbg>0,
			printf("alpha1_do_step4_calcs(): nUB1=%9.6e, rhoB1=%9.6f, muB1=%9.6f, nUB1a=%9.6e, nUB1b=%9.6e\n",nUB1,rhoB1,muB1,b1aResult[2],b1bResult[2]);
		);
	);

	nUB2=0;
	if(type(u2UB)!="t_POL",
		if(dbg>0,
			print("\nalpha1_do_step4_calcs(): about to check eliminating b2");
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
			printf("alpha1_do_step4_calcs(): nUB2=%9.6e, rhoB2=%9.6f, muB2=%9.6f, nUB2a=%9.6e, nUB2b=%9.6e\n",nUB2,rhoB2,muB2,b2aResult[1],b2bResult[1]);
		);
	);
	
	nUB3=0;
	if(type(u3UB)!="t_POL",
		if(dbg>0,
			print("\nalpha1_do_step4_calcs(): about to check eliminating b3");
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
			printf("alpha1_do_step4_calcs(): nUB3=%9.6e, rhoB3=%9.6f, muB3=%9.6f, nUB3a=%9.6e, nUB3b=%9.6e\n",nUB3,rhoB3,muB3,b3aResult[1],b3bResult[1]);
		);
	);

	\\ check our guess that using the b_i for which u_iUB is smallest is best
	if(u2UB<u3UB && 1.001*nUB3<nUB2,
		printf("BAD in alpha1_do_step4_calcs(): u2UB=%10.6f<u3UB=%10.6f, but nUB2=%10.6f>nUB3=%10.6f\n",u2UB,u3UB,nUB2,nUB3);
		\\1/0;
	);
	if(u3UB<u2UB && 1.001*nUB2<nUB3,
		printf("BAD in alpha1_do_step4_calcs(): u3UB=%10.6f<u2UB=%10.6f, but nUB3=%10.6f>nUB2=%10.6f\n",u3UB,u2UB,nUB3,nUB2);
		\\1/0;
	);
	if(dbg>0,
		printf("nUB1=%9.6e\n",nUB1);
		printf("nUB2=%9.6e\n",nUB2);
		printf("nUB3=%9.6e\n",nUB3);
		print("alpha1_do_step4_calcs(): END\n");
	);
	return([nUB1,rhoB1,muB1,nUB2,rhoB2,muB2,nUB3,rhoB3,muB3]);
}