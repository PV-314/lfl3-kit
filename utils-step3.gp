\\ \r lfl3\utils-step3.gp

read("lfl3\\utils-general.gp");

\\ returns 1 if the zero estimate conditions are satisfied
\\ returns 0 otherwise
check_zero_estimate(bigK,bigL,bigR1,bigR2,bigR3,bigS1,bigS2,bigS3,bigT1,bigT2,bigT3,chi,logXLB,areMultIndep,dbg=0)={
	my(iLHS,iRHSa,iRHSb,iRHSc,iRHSd,iiLHS,iiRHS,iiiLHS,iiiRHS,ivLHS,ivRHS,rhsDeg,vLHS,vRHS);
	
	\\ check condition (i)
	iLHS=(bigR1+1)*(bigS1+1)*(bigT1+1);
	iRHSa=bigR1+bigS1+1;
	rhsDeg=poldegree(iRHSa,logX);
	if(rhsDeg>1,
		print("check_zero_estimate() FAIL condition (i): bigR1+bigS1+1=",(bigR1+bigS1+1)," is of degree ",rhsDeg,": too large");
		return(0);
	);
	iRHSa=subst(iRHSa,logX,logXLB)/logXLB*logX;

	iRHSb=bigR1+bigT1+1;
	rhsDeg=poldegree(iRHSb,logX);
	if(rhsDeg>1,
		print("check_zero_estimate() FAIL condition (i): bigR1+bigT1+1=",(bigR1+bigT1+1)," is of degree ",rhsDeg,": too large");
		return(0);
	);
	iRHSb=subst(iRHSb,logX,logXLB)/logXLB*logX;

	iRHSc=bigS1+bigT1+1;
	rhsDeg=poldegree(iRHSc,logX);
	if(rhsDeg>1,
		print("check_zero_estimate() FAIL condition (i): bigS1+bigT1+1=",(bigS1+bigT1+1)," is of degree ",rhsDeg,": too large");
		return(0);
	);
	iRHSc=subst(iRHSc,logX,logXLB)/logXLB*logX;
	
	iRHSd=(bigR1+1)*(bigS1+1)*(bigT1+1);
	rhsDeg=poldegree(iRHSd,logX);
	if(rhsDeg!=2,
		print("check_zero_estimate() FAIL condition (i): (bigR1+1)*(bigS1+1)*(bigT1+1)=",((bigR1+1)*(bigS1+1)*(bigT1+1))," is of wrong degree ",rhsDeg);
		return(0);
	);
	iRHSd=subst(iRHSd,logX,logXLB)/logXLB^2*logX^2;
	iRHSd=sqrt(polcoef(iRHSd,2))*chi*logX;
	printf("check_zero_estimate(): condition  (i), LHS=%9.6e*logX^2\n",polcoef(iLHS,2));
	printf("    bigR1+bigS1+1=%9.6e*logX^2, bigR1+bigT1+1=%9.6e*logX^2, bigS1+bigT1+1=%9.6e*logX^2, chi*cV=%9.6e*logX^2\n",polcoef(bigK*iRHSa,2),polcoef(bigK*iRHSb,2),polcoef(bigK*iRHSc,2),polcoef(bigK*iRHSd,2));
	
	\\ check condition (ii). Actually checking (C.i.2)
	if(poldegree(bigR1)==1 && polcoef(bigR1,1)>polcoef(bigS1,1) && polcoef(bigR1,1)>polcoef(bigT1,1),
		iiLHS=2*(bigS1+1)*(bigT1+1);
	);
	if(poldegree(bigS1)==1 && polcoef(bigS1,1)>polcoef(bigR1,1) && polcoef(bigS1,1)>polcoef(bigT1,1),
		iiLHS=2*(bigR1+1)*(bigT1+1);
	);
	if(poldegree(bigT1)==1 && polcoef(bigT1,1)>polcoef(bigR1,1) && polcoef(bigT1,1)>polcoef(bigS1,1),
		iiLHS=2*(bigR1+1)*(bigS1+1);
	);
	if(areMultIndep,
		iiLHS=(bigR1+1)*(bigS1+1)*(bigT1+1);
	);
	if(poldegree(iiLHS)==1,
		iiLHS=subst(iiLHS,logX,logXLB)/logXLB*logX;
	);
	if(poldegree(iiLHS)==2,
		iiLHS=subst(iiLHS,logX,logXLB)/logXLB^2*logX^2;
	);

	iiRHS=bigL;
	if(poldegree(iiLHS)==1,
		printf("check_zero_estimate(): condition  (ii), LHS=%9.6e*log(x), while RHS=%9.6f\n",polcoef(iiLHS,1),iiRHS);
	);
	if(poldegree(iiLHS)==2,
		printf("check_zero_estimate(): condition  (ii), LHS=%9.6e*log(x)^2, while RHS=%9.6f\n",polcoef(iiLHS,2),iiRHS);
	);
	if(poldegree(iiLHS)<1 || poldegree(iiLHS)>2,
		print("check_zero_estimate() FAIL condition (ii): LHS is not linear or quadratic=",iiLHS);
		return(0);
	);
	if(poldegree(iiRHS)!=0,
		print("check_zero_estimate() FAIL condition (ii): RHS is not a constant=",iiRHS);
		return(0);
	);
	if(subst(iiLHS,logX,logXLB)<iiRHS,
		printf("check_zero_estimate() FAIL condition (ii): LHS=%9.6e*log(x)^2, while RHS=%9.6f\n",polcoef(iiLHS,2),iiRHS);
		return(0);
	);

	\\ check condition (iii). Actually checking (C.ii.1)
	if(poldegree(bigR2)==1 && polcoef(bigR2,1)>polcoef(bigS2,1) && polcoef(bigR2,1)>polcoef(bigT2,1),
		iiiLHS=(bigS2+1)*(bigT2+1);
	);
	if(poldegree(bigS2)==1 && polcoef(bigS2,1)>polcoef(bigR2,1) && polcoef(bigS2,1)>polcoef(bigT2,1),
		iiiLHS=(bigR2+1)*(bigT2+1);
	);
	if(poldegree(bigT2)==1 && polcoef(bigT2,1)>polcoef(bigR2,1) && polcoef(bigT2,1)>polcoef(bigS2,1),
		iiiLHS=(bigR2+1)*(bigS2+1);
	);
	if(areMultIndep,
		iiiLHS=(bigR2+1)*(bigS2+1)*(bigT2+1)/2;
		if(poldegree(iiiLHS)==2,
			iiiLHS=subst(iiiLHS,logX,logXLB)/logXLB*logX;
		);
	);
	iiiRHS=bigK*bigL;
	if(poldegree(iiiLHS)!=1,
		print("check_zero_estimate() FAIL condition (iii): LHS is not a poly of degree 1=",iiiLHS);
		return(0);
	);
	iiiLHS=subst(iiiLHS,logX,logXLB)/logXLB*logX;
	if(poldegree(iiiRHS)!=1,
		print("check_zero_estimate() FAIL condition (iii): RHS is not a poly of degree 1=",iiiRHS);
		return(0);
	);
	iiiRHS=subst(iiiRHS,logX,logXLB)/logXLB*logX;
	iiiLHS=polcoef(iiiLHS,1);
	iiiRHS=polcoef(iiiRHS,1);
	printf("check_zero_estimate(): condition (iii), LHS=%9.6e*log(x), while RHS=%9.6e*log(x)\n",iiiLHS,iiiRHS);
	if(iiiLHS<iiiRHS,
		printf("check_zero_estimate() FAIL condition (iii): LHS=%9.6f*log(x), while RHS=%9.6f*log(x)\n",iiiLHS,iiiRHS);
		return(0);
	);

	\\ check condition (iv)
	ivLHS=(bigR2+1)*(bigS2+1)*(bigT2+1);
	ivRHS=bigK*bigK;
	if(poldegree(ivLHS)!=2,
		print("check_zero_estimate() FAIL condition (iv): LHS is not a poly of degree 2=",ivLHS);
		return(0);
	);
	if(poldegree(ivRHS)!=2,
		print("check_zero_estimate() FAIL condition (iv): RHS is not a poly of degree 2=",ivRHS);
		return(0);
	);
	ivLHS=polcoef(ivLHS,2);
	ivRHS=polcoef(ivRHS,2);
	printf("check_zero_estimate(): condition  (iv), LHS=%9.6e*log(x)^2, while RHS=%9.6e*log(x)^2\n",ivLHS,ivRHS);
	if(ivLHS<ivRHS,
		printf("check_zero_estimate() FAIL condition (iv): LHS=%9.6f*log(x)^2, while RHS=%9.6f*log(x)^2\n",ivLHS,ivRHS);
		return(0);
	);

	\\ check condition (v)
	vLHS=(bigR3+1)*(bigS3+1)*(bigT3+1);
	vRHS=3*bigK*bigK*bigL;
	if(poldegree(vLHS)!=2,
		print("check_zero_estimate() FAIL condition (v): LHS is not a poly of degree 2=",vLHS);
		return(0);
	);
	if(poldegree(vRHS)!=2,
		print("check_zero_estimate() FAIL condition (v): RHS is not a poly of degree 2=",vRHS);
		return(0);
	);
	vLHS=polcoef(vLHS,2);
	vRHS=polcoef(vRHS,2);
	printf("check_zero_estimate(): condition   (v), LHS=%9.6e*log(x)^2, while RHS=%9.6e*log(x)^2\n",vLHS,vRHS);
	if(vLHS<vRHS,
		printf("check_zero_estimate() FAIL condition (v): LHS=%9.6e*log(x)^2, while RHS=%9.6e*log(x)^2\n",vLHS,vRHS);
		return(0);
	);
	print();
	return(1);
}

\\ pull out common code for alpha1Variable and alpha3Variable cases
\\ 6 July 2022
get_eqn42(a1,a2,a3,bigK,bigL,bigR,bigS,bigT,d,rho,logBUB,nUB,logXLB,dbg=0)={
	my(eqn42,eqn42LHS,eqn42Rem,eqn42RHS,eqn42RHS1,eqn42RHS2,eqn42RHS3,eqn42X0,eqn42X1,gDenom,gNumer);
	
	\\ miw matrix lhs with K_0=K-1:
	eqn42LHS=(bigK*bigL/2+bigL/2-0.37*bigK-2)*log(rho);

	gDenom=12*bigR*bigS*bigT/logX/logX;
	gDenom=subst(gDenom,logX,logXLB)*logX*logX;
	gNumer=bigK*(bigK+1)*bigL/2;
	gNumer=polcoef(gNumer,2,logX)*logX*logX;
	g=1/4-gNumer/gDenom;
	if(dbg>0,
		printf("get_eqn42(): D=%4d, g=%8.6f\n",d,g);
		\\printf("get_eqn42(): bigK=%s\n",bigK);
		\\printf("get_eqn42(): bigL=%s\n",bigL);
		\\printf("get_eqn42(): bigR=%s\n",bigR);
		\\printf("get_eqn42(): bigS=%s\n",bigS);
		\\printf("get_eqn42(): bigT=%s\n",bigT);
	);
	
	\\ we break up eqn42RHS for display purposes below
	\\ first the (cD+1)*log(N) term: log(N)=log(K(K+1)*L/2)~log((k*log(X))^2*L)
	\\ and we write it as a constant times logX:
	\\ following two lineas give the following for Catalan:
	\\ 0.027631021115928548208215897456212370491*logX + 86.75096956874051628333651063594036501
	\\eqn42RHS1=log(polcoef(bigK,1)*polcoef(bigK+1,1)*bigL/2)+2*log(logXLB)/logXLB*logX;
	\\eqn42RHS1=(d+1)*eqn42RHS1;
	
	\\ constant slightly bigger here
	\\ 0.027631021115928548208215897456212370491*logX + 86.75096956874733508633243174550927413
	eqn42RHS1=subst(bigK*(bigK+1)*bigL/2,logX,logXLB)/logXLB/logXLB;
	eqn42RHS1=log(eqn42RHS1)+2*log(logXLB)/logXLB*logX;
	eqn42RHS1=(d+1)*eqn42RHS1;

	\\ next the gL(a1*R+a2*S+a3*T) term
	eqn42RHS2=g*bigL*(a1*bigR+a2*bigS+a3*bigT);
	\\ lastly the 2*D(K-1)*log(b)/3 term
	eqn42RHS3=2*d*(bigK-1)*logBUB/3;
	eqn42RHS=eqn42RHS1+eqn42RHS2+eqn42RHS3;
	eqn42RHS=subst(eqn42RHS,logN,log(nUB));
	eqn42=eqn42LHS-eqn42RHS;
	eqn42X1=pollead(polcoef(eqn42,1))*logX;
	if(dbg>0,
		printf("get_eqn42(): before simplifying: deg(eqn42)=%4d, eqn42=%s\n",poldegree(eqn42),eqn42);
		print("pollead(polcoef(eqn42,0))=",pollead(polcoef(eqn42,0)));
		print("subst(eqn42-eqn42X1,logX,logXLB)=",subst(eqn42-eqn42X1,logX,logXLB));
	);
	eqn42X0=pollead(polcoef(eqn42,0));
	eqn42Rem=max(0, pollead(subst(eqn42-eqn42X1-eqn42X0,logX,logXLB)));
	eqn42=eqn42X1+eqn42X0+eqn42Rem;
	if(dbg>0,
		printf("eqn42LHS=                        %s\n",eqn42LHS);
		printf("eqn42RHS (cD+1)*log(N) term=     %s\n",eqn42RHS1);
		printf("eqn42RHS gL(a1*R+a2*S+a3*T) term=%s\n",eqn42RHS2);
		printf("eqn42RHS 2D(K-1)*log(b)/3 term=  %s\n",eqn42RHS3);
		printf("eqn42RHS=                        %s\n",eqn42RHS);
		printf("eqn42LHS-eqn42RHS=               %s\n",eqn42);
		printf("logBUB=%9.6f\n",logBUB);
	);
	return(eqn42);
}

\\ to save the same duplicated code.
\\ it returns minResult
\\          return([bigK,bigR1,bigR2,bigS1,bigT1,bigT2,newNonDegenNUB]);
\\ recall that val=[bigK,bigR1,bigR2,bigS1,bigT1,bigT2,newNonDegenNUB] or else it is empty
\\ also minResult=[minBigK, minBigL, minM, minRho, minChi, minBigR1, minBigR2, minBigS1, minBigT1, minBigT2, minNonDegenNUB]
\\ 3 March 2022
step3_update_min(val,bigL,m,rho,chi,minResult,dbg=0)={
	my(bigK,localMinResult,nNonDegenUB);

	localMinResult=minResult;
	if(length(val)!=7,
		\\printf("ERROR in step3_update_min(): val must have 7 elements. It has %2d for bigL=%4d, m=%9.6f, rho=%9.6f, chi=%9.6f\n",length(val),bigL,m,rho,chi);
		return(minResult);
	);
	if(length(minResult)!=11,
		printf("ERROR in step3_update_min(): minResult must have 11 elements. It has %2d for bigL=%4d, m=%9.6f, rho=%9.6f, chi=%9.6f\n",length(minResult),bigL,m,rho,chi);
		return([]);
	);
	
	nNonDegenUB=val[7]; \\ newNonDegenNUB
	if(dbg>0,
		print("nNonDegenUB=",nNonDegenUB,", type=",type(nNonDegenUB));
		print("localMinResult[11]=",localMinResult[11]);
	);
	if(nNonDegenUB>0 && nNonDegenUB<localMinResult[11],
		bigK=val[1];
		if(abs(polcoef(bigK,0))>0.0001 || abs(polcoef(bigK,1))<0.0001,
			printf("ERROR: K=%s is in incorrect form\n",bigK);
			error();
		);
		localMinResult[1]=bigK;
		localMinResult[2]=bigL;
		localMinResult[3]=m;
		localMinResult[4]=rho;
		localMinResult[5]=chi;
		localMinResult[6]=val[2]; \\ R_1
		localMinResult[7]=val[3]; \\ R_2
		localMinResult[8]=val[4]; \\ S_1
		if(dbg>0,
			print("bigS1=",val[4]);
		);
		localMinResult[9]=val[5]; \\ T_1
		localMinResult[10]=val[6]; \\ T_2
		localMinResult[11]=nNonDegenUB;
		\\printf("L=%4d, m=%8.4f, rho=%8.4f, chi=%6.4f, K=%12.3f*logX, nonDegen log|Lambda|>%9.6e*logX, nonDegenNUB=%10.6e\n",bigL,m,rho,chi,polcoef(bigK,1),-polcoef(bigK,1)*bigL*log(rho),nNonDegenUB);
		return(localMinResult);
	);
	\\if(isComplex,
	\\	printf("step3_update_min(): %s*n+%9.6f>log |u_i*Lambda|>-%9.6f*logX*(log(b')+0.21)^2=-%9.6f*logX*(log(n)%9.6f)^2\n",lamUB1,lamUB0+log(lamMul),polcoef(lambdaLB,1,logX),polcoef(lambdaLB,1,logX),bPCnstUB);
	\\);
	return(minResult);
}

\\ we estimate (from above) each of three parts of expression for b
\\ in terms of logX and n
\\ where n is the quantity we are trying to bound, and often an upper bound for |b_1|, |b_2| and |b_3|
\\ we use the definition of b and Lemma 3.4(a) to bound log(b) from above here
\\ assuming that K=c*log(x)
\\ 21 Nov 2021
get_logB_UB(bigK,bigR,bigS,bigT,b1,b2,b3,d1,d2,logXLB,nUB,dbg=0)={
	my(logB1,logB2,logBUB);
	
	\\ here we break the expression for b (defined in the notation section)
	\\ into three parts:
	\\ 1) B1=b3p*eta0=(b3/d1)*eta0;
	\\ 2) B2=b3pp*zeta0=(b3/d2)*zeta0;
	\\ 3) the factorial part
	\\ also note that despite the name, it is not actually a log of anything yet
	\\ (but that does happen at the very end of handling logB1 here. Same with logB2)
	logB1=internal_get_logB1(bigR,bigT,b1,b3,d1,logXLB,nUB,dbg);
	logB2=internal_get_logB2(bigS,bigT,b2,b3,d2,logXLB,nUB,dbg);

	logBUB=logB1+logB2;
	logBUB=logBUB-2*(log(polcoef(bigK,1))+logLogX)+11/3;
	return(logBUB);
}

\\ for b3'*eta_0 term in expression for b (equation (3.5) on 3 March 2022)
\\ here b3'=b3/d1
\\ 15 Dec 2021
internal_get_logB1(bigR,bigT,b1,b3,d1,logXLB,nUB,dbg=0)={
	my(logB1,logB1a,logB1b,logB1c,logB1d);

	if(dbg>0,
		print("internal_get_logB1(): b1=",b1,", b3=",b3);
	);
	logB1=((bigR-1)*b3+(bigT-1)*b1)/2/d1;
	logB1a=polcoef(polcoef(logB1,1,n),1,logX);
	if(logB1a<0,
		printf("ERROR in internal_get_logB1(): logB1 has negative lead coefficient: %s\n",logB1);
	);
	logB1b=polcoef(polcoef(logB1,1,n),0,logX);
	logB1b=max(logB1b,0);
	logB1c=polcoef(polcoef(logB1,0,n),1,logX);
	logB1c=max(logB1c,0);
	logB1d=polcoef(polcoef(logB1,0,n),0,logX);
	logB1d=max(logB1d,0);
	if(dbg>0,
		printf("internal_get_logB1(): nUB=%9.6e\n",nUB);
		printf("initial b3'*eta_0 value=%s\n",logB1);
		printf("        n*logX coeff=%12.5f, type(n*logX coeff)=%s\n",logB1a,type(logB1a));
		printf("        n      coeff=%12.5f, type(n      coeff)=%s\n",logB1b,type(logB1b));
		printf("        logX   coeff=%12.5f, type(logX   coeff)=%s\n",logB1c,type(logB1c));
		printf("        cnst   coeff=%12.5f, type(cnst   coeff)=%s\n",logB1d,type(logB1d));
	);
	
	logB1=subst(logB1a*n*logX+logB1b*n+logB1c*logX+logB1d,n,nUB);
	
	logB1a=polcoef(logB1,1,logX);
	logB1b=polcoef(logB1,0,logX);
	logB1b=max(logB1b,0);
	if(dbg>0,
		printf("after substituting nUB for n: b3'*eta_0 value<%s\n",logB1);
		printf("        logX term  =%9.6e\n",logB1a);
		printf("        const term =%9.6e\n",logB1b);
	);
	logB1=log(logB1a+logB1b/logXLB)+logLogX;
	if(dbg>0,
		printf("actual log(b3'*eta_0) value<%s\n",logB1);
	);
	return(logB1);
}

\\ here we break the expression for b (defined in the notation section)
\\ into three parts:
\\ 1) B1=b3p*eta0=(b3/d1)*eta0 -- in fact, we use an upper bound for |b3|
\\ 2) B2=b3pp*zeta0=(b3/d2)*zeta0 -- in fact, we use an upper bound for |b3|
\\ 3) the factorial part
\\ also note that despite the name, it is not actually a log of anything yet
\\ (but that does happen at the very end of handling logB1 here. Same with logB2)
\\ 15 Dec 2021
get_logB2(bigM,bigS,bigT,b2UB,b3UB,d2,logXLB,nUB,dbg=0)={
	my(logB2);
	
	logB2=internal_get_logB2(bigS,bigT,b2UB,b3UB,d2,logXLB,nUB,dbg);
	logB2=logB2-(log(polcoef(bigM,1))+logLogX-3/2);
	return(logB2);
}

\\ for b3''*zeta_0 term in expression for b (equation (3.5) on 3 March 2022)
\\ 15 Dec 2021
internal_get_logB2(bigS,bigT,b2,b3,d2,logXLB,nUB,dbg=0)={
	my(logB2,logB2a,logB2b,logB2c,logB2d);
	
	logB2=((bigS-1)*b3+(bigT-1)*b2)/2/d2;
	logB2a=polcoef(polcoef(logB2,1,n),1,logX);
	if(logB2a<0,
		printf("ERROR: logB2 has negative lead coefficient: %s\n",logB2);
	);
	logB2b=polcoef(polcoef(logB2,1,n),0,logX);
	logB2b=max(logB2b,0);
	logB2c=polcoef(polcoef(logB2,0,n),1,logX);
	logB2c=max(logB2c,0);
	logB2d=polcoef(polcoef(logB2,0,n),0,logX);
	logB2d=max(logB2d,0);
	if(dbg>0,
		printf("\ninitial b3''*zeta_0 value=%s\n",logB2);
		printf("        n*logX coeff=%12.5f\n",logB2a);
		printf("        n      coeff=%12.5f\n",logB2b);
		printf("        logX   coeff=%12.5f\n",logB2c);
		printf("        cnst   coeff=%12.5f\n",logB2d);
	);
	
	logB2=subst(logB2a*n*logX+logB2b*n+logB2c*logX+logB2d,n,nUB);
	logB2a=polcoef(logB2,1,logX);
	logB2b=polcoef(logB2,0,logX);
	logB2b=max(logB2b,0);
	if(dbg>0,
		printf("after substituting nUB for n: b3''*zeta_0 value<%s\n",logB2);
		printf("        logX term  =%9.6e\n",logB2a);
		printf("        const term =%9.6e\n",logB2b);
	);
	logB2=log(logB2a+logB2b/logXLB)+logLogX;
	if(dbg>0,
		printf("actual log(b3''*zeta_0) value<%s\n",logB2);
	);
	return(logB2);
}