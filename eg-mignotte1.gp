\\ \r lfl3\eg-mignotte1.gp

read("lfl3\\kit-alpha1Variable.gp");

\\ search for parameter choice for Example 1 in Mignotte's original kit:
\\ for the diophantine equation x^n-2^\alpha*5^\beta*y^n = \pm 1
\\ where n \geq 3 is prime, \alpha=1,2,3, 0 \leq \beta<n
\\ \Lamnda = n \log (x/y)-\alpha \log(2)-\beta \log(5)

\\ L=  81, m= 10.5000, rho(3logs)=  8.5000, chi= 0.8000, K=  179708.542*logX, nonDegen log|Lambda|>-3.115164 e7*logX,
\\ nonDegenNUB=3.115164 e7, rho(2logs)= 7.800000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=3.074700 e7, degenNUB3=6.337995 e7, nUB=3.115164 e7, transB-b2
\\ 4 Jan 2022
eg1_check_it1()={
	my(bigL,chi,m,mu,rho2Logs,rho3Logs);
	
	bigL=81;
	m=10.5;
	rho2Logs=7.8;
	rho3Logs=8.5;
	chi=0.8;
	mu=0.61;

	actMinNUB=eg1_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,,2);
	expMinNUB=31.15164*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg1_check_it1(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  53, m= 10.0000, rho(3logs)=  9.0000, chi= 1.1250, K=  124475.071*logX,
\\ nonDegen log|Lambda|>-1.449548 e7*logX, nonDegenNUB=1.449548 e7, rho(2logs)= 7.800000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=1.407796 e7, degenNUB3=2.895016 e7, nUB=1.449548 e7, eliminate-b2
\\ 4 Jan 2022 (updated with above on 27 Feb 2022)
eg1_check_it2()={
	my(bigL,chi,m,mu,nUB,rho2Logs,rho3Logs);

	bigL=53;
	m=10.0;
	rho2Logs=7.8;
	rho3Logs=9.0;
	chi=1.125;
	mu=0.61;
	nUB=32*10^6;

	actMinNUB=eg1_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=14.49548*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg1_check_it2(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  53, m= 11.2500, rho(3logs)=  8.2500, chi= 1.1500, K=  119255.328*logX, nonDegen log|Lambda|>-1.333767 e7*logX,
\\ nonDegenNUB=1.333767 e7, rho(2logs)= 7.800000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=1.290824 e7, degenNUB3=2.652192 e7, nUB=1.333767 e7, transB-b2
\\ 4 Jan 2022
eg1_check_it3()={
	my(bigL,chi,m,mu,nUB,rho2Logs,rho3Logs);

	bigL=53;
	m=11.25;
	rho2Logs=7.8;
	rho3Logs=8.25;
	chi=1.15;
	mu=0.61;
	nUB=15*10^6;

	actMinNUB=eg1_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=13.33767*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg1_check_it3(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ 4 Jan 2022
eg1_search_it1(dbg=0)={
my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigD,bigK,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,m,matveevChi,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,rho,val,w);
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	bigLLB=50;
	bigLUB=200;
	mLB=5;
	mUB=mLB+10;
	rho3LB=5;
	rho3UB=rho3LB+10;
	chiLB=0.5;
	chiUB=chiLB+0.5;
	rho2LB=6;
	rho2UB=rho2LB+4;
	eg1_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,,dbg);
}

\\ 4 Jan 2022 (27 Feb 2022: should use nUBInit=32*10^6 -- upper bound from iteration 1)
eg1_search_it2(nUBInit,dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	bigLLB=30;
	bigLUB=150;
	mLB=5;
	mUB=mLB+10;
	rho3LB=5;
	rho3UB=rho3LB+10;
	chiLB=0.75;
	chiUB=chiLB+0.5;
	rho2LB=7;
	rho2UB=rho2LB+4;
	eg1_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,nUBInit,dbg);
}

\\ 4 Jan 2022 (should use nUBInit=15*10^6 -- upper bound from iteration 2)
eg1_search_it3(nUBInit,dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	bigLLB=30;
	bigLUB=150;
	mLB=8;
	mUB=mLB+5;
	rho3LB=5;
	rho3UB=rho3LB+5;
	chiLB=0.8;
	chiUB=chiLB+0.5;
	rho2LB=7;
	rho2UB=rho2LB+4;
	eg1_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,nUBInit,dbg);
}

\\ only set muLB and muUB when using "check" functions above, not when using "search" functions
\\ (they are bounds for the mu in Theorem 2 of Laurent's 2008 paper)
\\ 3 July 2022
eg1_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,muLB=0,muUB=0,nUBInit=0,dbg=0)={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,areBoundsOK,b1,b2,b3,bigD,bigK,chiStep,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,matveevChi,minNUB,mStep,nDegenUB,nLB,nNonDegenUB,nUB,rho3Step,startTime,step3Result);

	startTime=getwalltime();
	\\ bigD=[Q(al_1,al_2,al_3):Q] -- used for Matveev's bounds
	bigD=1;
	\\ matveevChi=[R(al_1,al_2,al_3):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2,al_3):Q]/[R(al_1,al_2,al_3):R]
	d=bigD/matveevChi; \\ d is the "degree" value used in the kit

	al1=x; \\ x is just because it is not known.
	al2=2;
	al3=5;

	b1=n;
	b2=3; \\ alpha=1, 2 or 3
	b3=n-1; \\ 0 \leq beta<n
	logXLB=floor(nLB/2600-1.5);
	if(dbg>0,
		printf("calced logXLB=%5d\n",logXLB);
	);
	logXLB=380;
	nLB=10^6; \\ assumption at the start of proof of Theorem 6.3

	hgtA1=logX;
	absLogA1=5.0001; \\ from Step 3 after defns of a_i's
	hgtA2=log(2);
	absLogA2=abs(log(al2));
	hgtA3=log(5);
	absLogA3=abs(log(al3));
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX;
	lamUB0=log(2);

	if(nUBInit==0,
		\\print("bigD=",bigD,", matveevChi=",matveevChi,", al1=",al1,", absLogA1=",absLogA1,", hgtA1=",hgtA1,", al2=",al2);
		\\print("absLogA2=",absLogA2,", hgtA2=",hgtA2,", al3=",al3,", absLogA3=",absLogA3,", hgtA3=",hgtA3);
		\\print("b1=",b1,", b2=",b2,b3,", logXLB=",logXLB,", nLB=",nLB,", lamUB1=",lamUB1,", lamUB0=",lamUB0,", dbg=",dbg);
		nUBInit=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,dbg);
		printf("nUBInit=%9.6e\n",nUBInit);
		printf("eg1_search_general(): calculated nUBInit=%9.6e\n",nUBInit);
	);
	minNUB=nUBInit;
	printf("used nUBInit=%9.6e\n",nUBInit);

	areBoundsOK=check_bounds(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB);
	if(areBoundsOK==0,
		return();
	);
	chiStep=get_step(chiLB,chiUB);
	mStep=get_step(mLB,mUB);
	rho3Step=get_step(rho3LB,rho3UB);

	\\ need chi as the outer loop, as we need to get step 3 and step 4 values for each chi
	forstep(chi=chiLB,chiUB,chiStep,
		\\ step3Result=[minBigK, minBigL, minM, minRho, minChi, minBigR1, minBigR2, minBigS1, minBigT1, minBigT2, minNonDegenNUB]
		step3Result=[0,0,0,0,0,0,0,0,0,0,minNUB];
		for(bigL=bigLLB,bigLUB, \\ L=5 is the lower bound in Theorem 4.1
			\\if(bigL%10==0,print("L=",bigL));
			forstep(m=mLB,mUB,mStep,
				forstep(rho3=rho3LB,rho3UB,rho3Step,
					a1=(rho3-1)*absLogA1+2*logX;
					a2=(rho3+1)*log(al2);
					a3=(rho3+1)*log(al3);
					step3Result=alpha1_do_step3(step3Result,d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,chi,bigL,m,rho3,nUBInit,logXLB,nLB,lamUB1,lamUB0,dbg);
				);
			);
		);
		if(length(step3Result)>0 && step3Result[11]<minNUB,
			printf("for chi=%9.6f, minNonDegenNUB=%9.6e\n",chi,step3Result[11]);
			minNUB=alpha1_do_step4(step3Result,minNUB,d,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,chi,rho2LB,rho2UB,muLB,muUB,logXLB,nLB,lamUB1,lamUB0,dbg);
		);
	);
	print("time taken=",(getwalltime()-startTime));	
	return(minNUB);
}