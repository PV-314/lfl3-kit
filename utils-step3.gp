\\ \r lfl3\step3-utils.gp

read("lfl3\\lfl-utils-general.gp");

\\ pull out common code for alpha1Variable and alpha3Variable cases
\\ 6 July 2022
get_eqn42(a1,a2,a3,bigK,bigL,bigR,bigS,bigT,d,rho,logBUB,nUB,logXLB,dbg=0)={
	my(eqn42,eqn42LHS,eqn42Rem,eqn42RHS,eqn42X0,eqn42X1,gDenom,gNumer);
	
	\\ with K_0=2(K-1):
	\\eqn42LHS=(bigK*bigL/2+bigL/4+bigL/8/bigK-2*bigK/3/bigL-1)*log(rho);
	\\ with K_0=K-1, rather than K_0=2(K-1)
	eqn42LHS=(bigK*bigL-bigK-bigK/3/bigL)*log(rho);

	gDenom=12*bigR*bigS*bigT/logX/logX;
	gDenom=subst(gDenom,logX,logXLB)*logX*logX;
	gNumer=bigK*bigK*bigL;
	gNumer=polcoef(gNumer,2,logX)*logX*logX;
	g=1/4-gNumer/gDenom;
	if(dbg!=0,
		printf("get_eqn42(): D=%4d, g=%8.6f\n",d,g);
		\\printf("get_eqn42(): bigK=%s\n",bigK);
		\\printf("get_eqn42(): bigL=%s\n",bigL);
		\\printf("get_eqn42(): bigR=%s\n",bigR);
		\\printf("get_eqn42(): bigS=%s\n",bigS);
		\\printf("get_eqn42(): bigT=%s\n",bigT);
	);
	
	\\ first the (cD+1)*log(N) term
	eqn42RHS=(d+1)*log(polcoef(bigK,1)*polcoef(bigK,1)*bigL)+(d+1)*log(logXLB)/logXLB*logX;
	\\ next the gL(a1*R+a2*S+a3*T) term
	eqn42RHS=eqn42RHS+g*bigL*(a1*bigR+a2*bigS+a3*bigT);
	\\ lastly the D(K-1)*log(b) term
	eqn42RHS=eqn42RHS+d*(bigK-1)*logBUB;
	eqn42RHS=subst(eqn42RHS,logN,log(nUB));
	eqn42=eqn42LHS-eqn42RHS;
	eqn42X1=pollead(polcoef(eqn42,1))*logX;
	if(dbg!=0,
		printf("get_eqn42(): before simplifying: deg(eqn42)=%4d, eqn42=%s\n",poldegree(eqn42),eqn42);
		print("pollead(polcoef(eqn42,0))=",pollead(polcoef(eqn42,0)));
		print("subst(eqn42-eqn42X1,logX,logXLB)=",subst(eqn42-eqn42X1,logX,logXLB));
	);
	eqn42X0=pollead(polcoef(eqn42,0));
	eqn42Rem=max(0, pollead(subst(eqn42-eqn42X1-eqn42X0,logX,logXLB)));
	eqn42=eqn42X1+eqn42X0+eqn42Rem;
	if(dbg!=0,
		printf("eqn42LHS=%s\n",eqn42LHS);
		printf("eqn42RHS D(K-1)*log(b) term=%s\n",d*(bigK-1)*logBUB);
		printf("eqn42RHS gL(a1*R+a2*S+a3*T) term=%s\n",g*bigL*(a1*bigR+a2*bigS+a3*bigT));
		printf("eqn42RHS (cD+1)*log(N) term=%s\n",(d+1)*log(polcoef(bigK,1)*polcoef(bigK,1)*bigL)+(d+1)*log(logXLB)/logXLB*logX);
		printf("eqn42RHS=%s\n",eqn42RHS);
		printf("eqn42LHS-eqn42RHS=%s\n",eqn42);
	);
	return(eqn42);
}

\\ to save the same duplicated code
\\ recall that val=[bigK,bigR1,bigR2,bigS1,bigT1,bigT2,newNonDegenNUB] or else it is empty
\\ also minResult=[minBigK, minBigL, minM, minRho, minChi, minBigR1, minBigR2, minBigS1, minBigT1, minBigT2, minNonDegenNUB]
\\ 3 March 2022
XXX_step3_update_min(val,bigL,m,rho,chi,minResult,dbg=0)={
	my(bigK,localMinResult,nNonDegenUB);
	
	localMinResult=minResult;
	if(length(val)==7,
		nNonDegenUB=val[7]; \\ newNonDegenNUB
		if(dbg!=0,
			print("nNonDegenUB=",nNonDegenUB,", type=",type(nNonDegenUB));
			print("localMinResult[11]=",localMinResult[11]);
		);
		if(nNonDegenUB>0 && nNonDegenUB<localMinResult[11],
			bigK=val[1];
			if(abs(polcoef(bigK,0))>0.0001 || abs(polcoef(bigK,1))<0.0001,
				printf("ERROR: K=%s is in incorrect form\n",bigK);
				1/0;
			);
			localMinResult[1]=bigK;
			localMinResult[2]=bigL;
			localMinResult[3]=m;
			localMinResult[4]=rho;
			localMinResult[5]=chi;
			localMinResult[6]=val[2]; \\ R_1
			localMinResult[7]=val[3]; \\ R_2
			localMinResult[8]=val[4]; \\ S_1
			localMinResult[9]=val[5]; \\ T_1
			localMinResult[10]=val[6]; \\ T_2
			localMinResult[11]=nNonDegenUB;
			\\printf("L=%4d, m=%8.4f, rho=%8.4f, chi=%6.4f, K=%12.3f*logX, nonDegen log|Lambda|>%9.6e*logX, nonDegenNUB=%10.6e\n",bigL,m,rho,chi,polcoef(bigK,1),-polcoef(bigK,1)*bigL*log(rho),nNonDegenUB);
		);
		if(isComplex,
			printf("in get_degen_nUB(): %s*n+%9.6f>log |u_i*Lambda|>-%9.6f*logX*(log(b')+0.21)^2=-%9.6f*logX*(log(n)%9.6f)^2\n",lamUB1,lamUB0+log(lamMul),polcoef(lambdaLB,1,logX),polcoef(lambdaLB,1,logX),bPCnstUB);
		);
	);
	return(localMinResult);
}

\\ we estimate (from above) each of three parts of expression for b'
\\ in terms of logX and n
\\ where n is the quantity we are trying to bound, and often an upper bound for |b_1|, |b_2| and |b_3|
\\ we use the definition of b and Lemma 3.4(a) to bound log(b) from above here
\\ assuming that K=c*log(x)
\\ 21 Nov 2021
get_logBPrime_UB(bigK,bigR,bigS,bigT,b1,b2,b3,d1,d2,logXLB,nUB,dbg=0)={
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
	logBUB=logBUB-2*(log(polcoef(bigK,1))+logLogX-3/2);
	return(logBUB);
}

\\ for b3'*eta_0 term in expression for b (equation (3.5) on 3 March 2022)
\\ here b3'=b3/d1
\\ 15 Dec 2021
internal_get_logB1(bigR,bigT,b1,b3,d1,logXLB,nUB,dbg=0)={
	my(logB1,logB1a,logB1b,logB1c,logB1d);

	if(dbg!=0,
		print("internal_get_logB1(): b1=",b1,", b3=",b3);
	);
	logB1=((bigR-1)*b3+(bigT-1)*b1)/2/d1;
	logB1a=polcoef(polcoef(logB1,1,n),1,logX);
	if(logB1a<0,
		printf("ERROR: in internal_get_logB1(), logB1 has negative lead coefficient: %s\n",logB1);
	);
	logB1b=polcoef(polcoef(logB1,1,n),0,logX);
	logB1b=max(logB1b,0);
	logB1c=polcoef(polcoef(logB1,0,n),1,logX);
	logB1c=max(logB1c,0);
	logB1d=polcoef(polcoef(logB1,0,n),0,logX);
	logB1d=max(logB1d,0);
	if(dbg!=0,
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
	if(dbg!=0,
		printf("after substituting for n: b3'*eta_0 value=%s\n",logB1);
		printf("        logX term  =%9.6e\n",logB1a);
		printf("        const term =%9.6e\n",logB1b);
	);
	logB1=log(logB1a+logB1b/logXLB)+logLogX;
	if(dbg!=0,
		printf("actual b3'*eta_0 value=%s\n",logB1);
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
	if(dbg!=0,
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
	if(dbg!=0,
		printf("after substituting for n: b3''*zeta_0 value=%s\n",logB2);
		printf("        logX term  =%9.6e\n",logB2a);
		printf("        const term =%9.6e\n",logB2b);
	);
	logB2=log(logB2a+logB2b/logXLB)+logLogX;
	if(dbg!=0,
		printf("actual b3''*zeta_0 value=%s\n",logB2);
	);
	return(logB2);
}
