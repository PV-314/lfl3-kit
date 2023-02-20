\\ \r lfl3\eg-mignotte2.gp

read("lfl3\\kit-alpha1Variable.gp");

\\ search for parameter choice for Example 2 in Mignotte's original kit:
\\ for the diophantine equation 2^\alpha*x^n-5^\beta*y^n = \pm 1
\\ where n \geq 3 is prime, 0 \leq \alpha, \beta<n
\\ \Lamnda = \alpha \log(p)-n \log (y/x)-\beta \log(q)
\\ where p=2 and q=5

\\ L= 112, m= 25.0000, rho(3logs)=  5.0000, chi= 0.9000, K=  596929.615*logX, nonDegen log|Lambda|>-1.076008 e8*logX,
\\ nonDegenNUB=1.076008 e8, rho(2logs)= 5.500000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=2.111392 e8, degenNUB3=1.018226 e8, nUB=1.076008 e8, transB-b3
\\ 25 Nov 2021
eg2_check_it1()={
	my(bigL,chi,m,mu,rho2Logs,rho3Logs);

	bigL=112;
	m=25;
	rho2Logs=5.5;
	rho3Logs=5;
	chi=0.9;
	mu=0.61;

	actMinNUB=eg2_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,,2);
	expMinNUB=107.6008*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg2_check_it1(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  74, m= 30.0000, rho(3logs)=  5.0000, chi= 1.2000, K=  473279.909*logX, nonDegen log|Lambda|>-5.636688 e7*logX,
\\ nonDegenNUB=5.636688 e7, rho(2logs)= 6.000000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=1.181038 e8, degenNUB3=5.691803 e7, nUB=5.691803 e7, transB-b3
\\ 9 July 2022
eg2_check_it2()={
	my(bigL,chi,m,mu,nUB,rho2Logs,rho3Logs);

	bigL=74;
	m=30;
	rho2Logs=6;
	rho3Logs=5;
	chi=1.2;
	mu=0.61;
	nUB=108*10^6;

	actMinNUB=eg2_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=56.91803*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg2_check_it2(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  72, m= 30.0000, rho(3logs)=  5.0000, chi= 1.2500, K=  460488.560*logX, nonDegen log|Lambda|>-5.336120 e7*logX,
\\ nonDegenNUB=5.336120 e7, rho(2logs)= 6.000000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=1.087008 e8, degenNUB3=5.256999 e7, nUB=5.336120 e7, transB-b3
\\ 9 July 2022
eg2_check_it3()={
	my(bigL,chi,m,mu,nUB,rho2Logs,rho3Logs);

	bigL=72;
	m=30;
	rho2Logs=6;
	rho3Logs=5;
	chi=1.25;
	mu=0.61;
	nUB=57*10^6;

	actMinNUB=eg2_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=53.36120*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg2_check_it3(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ 25 Nov 2021
eg2_search_it1(dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	bigLLB=50;
	bigLUB=150;
	mLB=10;
	mUB=mLB+20;
	rho3LB=2;
	rho3UB=rho3LB+10;
	chiLB=0.5;
	chiUB=chiLB+1;
	rho2LB=2;
	rho2UB=rho2LB+10;
	eg2_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,,dbg);
}

\\ 4 Jan 2022 (should use nUBInit=108*10^6 -- upper bound from iteration 1)
eg2_search_it2(nUBInit,dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);

	bigLLB=50;
	bigLUB=150;
	mLB=15;
	mUB=mLB+20;
	rho3LB=2;
	rho3UB=rho3LB+10;
	chiLB=0.5;
	chiUB=chiLB+2.0;
	rho2LB=6;
	rho2UB=rho2LB+5;
	eg2_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,nUBInit,dbg);
}

\\ 5 Jan 2022 (should use nUBInit=57*10^6 -- upper bound from iteration 2)
eg2_search_it3(nUBInit,dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);

	bigLLB=50;
	bigLUB=150;
	mLB=15;
	mUB=mLB+20;
	rho3LB=3;
	rho3UB=rho3LB+5;
	chiLB=0.7;
	chiUB=chiLB+1.0;
	rho2LB=6;
	rho2UB=rho2LB+5;
	eg2_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,nUBInit,dbg);
}

\\ only set muLB and muUB when using "check" functions above, not when using "search" functions
\\ (they are bounds for the mu in Theorem 2 of Laurent's 2008 paper)
\\ 9 July 2022
eg2_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,muLB=0,muUB=0,nUBInit=0,dbg=0)={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,areBoundsOK,b1,b2,b3,bigD,bigK,chiStep,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,matveevChi,minNUB,mStep,nDegenUB,nLB,nNonDegenUB,nUB,rho3Step,startTime,step3Result);

	startTime=getwalltime();
	\\ bigD=[Q(al_1,al_2,al_3):Q] -- used for Matveev's bounds
	bigD=1;
	\\ matveevChi=[R(al_1,al_2,al_3):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2,al_3):Q]/[R(al_1,al_2,al_3):R]
	d=bigD/matveevChi; \\ d is the "degree" value used in the kit

	al1=x; \\ x is just because it is not known. Here it is actually al_1=x/y (or y/x)
	al2=5;
	al3=2;
	b1=n; \\ 0 \leq alpha < n
	b2=n;
	b3=n-1; \\ 0 \leq beta < n
	logXLB=log(7);
	nLB=5*10^6; \\ assumption at the start of proof of Theorem 6.3
	
	hgtA1=logX;
	absLogA1=logX;
	hgtA2=log(al2);
	absLogA2=abs(log(al2));
	hgtA3=log(al3);
	absLogA3=abs(log(al3));
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX;
	lamUB0=0;
	
	if(nUBInit==0,
		\\print("bigD=",bigD,", matveevChi=",matveevChi,", al1=",al1,", absLogA1=",absLogA1,", hgtA1=",hgtA1,", al2=",al2);
		\\print("absLogA2=",absLogA2,", hgtA2=",hgtA2,", al3=",al3,", absLogA3=",absLogA3,", hgtA3=",hgtA3);
		\\print("b1=",b1,", b2=",b2,b3,", logXLB=",logXLB,", nLB=",nLB,", lamUB1=",lamUB1,", lamUB0=",lamUB0,", dbg=",dbg);
		nUBInit=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,dbg);
		printf("nUBInit=%9.6e\n",nUBInit);
		printf("eg2_search_general(): calculated nUBInit=%9.6e\n",nUBInit);
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
					a1=(rho3-1)*log(5.0001)+2*logX;
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