\\ \r lfl3\eg-bms-annals.gp

read("lfl3\\kit-alpha3Variable.gp");

\\ this is the linear form in the Annals paper of Yann, Maurice and Samir

\\ L= 167, m=  6.0000, rho(3logs)= 10.0000, chi= 0.7500, K=  221944.442*logX, nonDegen log|Lambda|>-8.534468 e7*logX,
\\ nonDegenNUB=4.267234 e7, rho(2logs)=10.000000, mu(2logs)= 0.610000, degenNUB1=3.883929 e7, degenNUB2=7.804339 e7, degenNUB3=    0.e-19, nUB=4.267234 e7, transB-b1
\\ 10 May 2022 (after improving constant to use from Laurent's 2008 paper)
bms1_check_it1()={
	my(bigL,chi,m,mu,nUB,rho2Logs,rho3Logs);
	
	bigL=167;
	m=6.0;
	rho2Logs=10.0;
	rho3Logs=10.0;
	chi=0.75;
	mu=0.61;

	actMinNUB=bms1_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,,2);
	expMinNUB=4.267234*10^7;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: bms1_check_it1(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L= 105, m=  7.2500, rho(3logs)=  9.7500, chi= 1.0300, K=  161616.787*logX, nonDegen log|Lambda|>-3.864469 e7*logX,
\\ nonDegenNUB=1.932234 e7, rho(2logs)=10.000000, mu(2logs)= 0.610000, degenNUB1=1.892864 e7, degenNUB2=3.780783 e7, degenNUB3=    0.e-19, nUB=1.932234 e7, transB-b1
\\ second iteration values
\\ 11 May 2022
bms1_check_it2()={
	my(bigL,chi,m,mu,nUB,rho2Logs,rho3Logs);
	
	bigL=105;
	m=7.25;
	rho2Logs=10.0;
	rho3Logs=9.75;
	chi=1.03;
	mu=0.61;
	nUB=43*10^6;

	actMinNUB=bms1_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=19.32234*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: bms1_check_it2(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L= 100, m=  7.0000, rho(3logs)= 10.0000, chi= 1.0500, K=  155051.007*logX, nonDegen log|Lambda|>-3.570181 e7*logX,
\\ nonDegenNUB=1.785091 e7, rho(2logs)=10.000000, mu(2logs)= 0.610000, degenNUB1=1.784192 e7, degenNUB2=3.562176 e7, degenNUB3=    0.e-19, nUB=1.785091 e7, eliminate-b1
\\ third iteration values
\\ 13 Dec 2021
bms1_check_it3()={
	my(bigL,chi,m,mu,nUB,rho2Logs,rho3Logs);
	
	bigL=100;
	m=7.0;
	rho2Logs=10;
	rho3Logs=10;
	chi=1.05;
	mu=0.61;
	nUB=19.4*10^6;

	actMinNUB=bms1_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=17.85091*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: bms1_check_it3(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ for iteration 1
\\ 11 Dec 2021 (changed to use "general" function on 3 July 2022)
bms1_search_it1(dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	bigLLB=100;
	bigLUB=250;
	mLB=4;
	mUB=mLB+5;
	rho3LB=7;
	rho3UB=rho3LB+5;
	chiLB=0.5;
	chiUB=chiLB+1;
	rho2LB=7;
	rho2UB=rho2LB+4;
	bms1_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,,dbg);
}

\\ for iteration 2
\\ 13 Dec 2021
bms1_search_it2(nUBInit,dbg=0)={ \\use nUBInit=8.37*10^6
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	bigLLB=50;
	bigLUB=200;
	mLB=4;
	mUB=mLB+5;
	rho3LB=7;
	rho3UB=rho3LB+5;
	chiLB=0.65;
	chiUB=chiLB+0.4;
	rho2LB=7;
	rho2UB=rho2LB+4;
	bms1_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,nUBInit,dbg);
}

\\ for iteration 3
\\ 13 Dec 2021
bms1_search_it3(nUBInit,dbg=0)={ \\ use nUBInit=3.87*10^6
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);

	bigLLB=60;
	bigLUB=120;
	mLB=5;
	mUB=mLB+4;
	rho3LB=7;
	rho3UB=rho3LB+4;
	chiLB=0.9;
	chiUB=chiLB+0.2;
	rho2LB=7;
	rho2UB=rho2LB+4;
	bms1_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,nUBInit,dbg);
}

\\ only set muLB and muUB when using "check" functions above, not when using "search" functions
\\ (they are bounds for the mu in Theorem 2 of Laurent's 2008 paper)
\\ 3 July 2022
bms1_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,muLB=0,muUB=0,nUBInit=0,dbg=0)={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,areBoundsOK,b1,b2,b3,bigD,bigK,chiStep,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,matveevChi,minNUB,mStep,nDegenUB,nLB,nNonDegenUB,nUB,rho3Step,startTime,step3Result,w);

	startTime=getwalltime();
	\\ bigD=[Q(al_1,al_2,al_3):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2,al_3):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2,al_3):Q]/[R(al_1,al_2,al_3):R]
	d=bigD/matveevChi; \\ d is the "degree" value used in the kit

	w=(1+sqrt(5))/2;
	logW=log(w);
	al1=w;
	al2=sqrt(5);
	al3=x; \\ x is just because it is not known. Here it is al_3=omega^k/y
	b1=n-1;
	b2=1;
	b3=n;
	logXLB=10^20;
	nLB=10*10^6;
	
	hgtA1=logW/2;
	absLogA1=abs(log(al1));
	hgtA2=log(sqrt(5));
	absLogA2=abs(log(al2));
	hgtA3=logW/2+logX+10^(-6); \\ from equation (6.5)
	absLogA3=logW+10^(-6); \\ from equation (6.4)
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-2*logX;
	lamUB0=1;

	if(nUBInit==0,
		\\print("bigD=",bigD,", matveevChi=",matveevChi,", al1=",al1,", absLogA1=",absLogA1,", hgtA1=",hgtA1,", al2=",al2);
		\\print("absLogA2=",absLogA2,", hgtA2=",hgtA2,", al3=",al3,", absLogA3=",absLogA3,", hgtA3=",hgtA3);
		\\print("b1=",b1,", b2=",b2,b3,", logXLB=",logXLB,", nLB=",nLB,", lamUB1=",lamUB1,", lamUB0=",lamUB0,", dbg=",dbg);
		nUBInit=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,dbg);
		printf("nUBInit=%9.6e\n",nUBInit);
		printf("bms1_search_general(): calculated nUBInit=%9.6e\n",nUBInit);
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
					a1=(rho3+1)*log(al1);
					a2=(rho3+3)*log(al2);
					a3=(rho3+1)*logW+4*logX+(rho3+3)*10^(-6); \\ here logX=log(y)
					step3Result=alpha3_do_step3(step3Result,d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,chi,bigL,m,rho3,nUBInit,logXLB,nLB,lamUB1,lamUB0,dbg);
				);
			);
		);
		if(length(step3Result)>0 && step3Result[11]<minNUB,
			printf("\nfor chi=%9.6f, minNonDegenNUB=%9.6e\n",chi,step3Result[11]);
			minNUB=alpha3_do_step4(step3Result,minNUB,d,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,chi,rho2LB,rho2UB,muLB,muUB,logXLB,nLB,lamUB1,lamUB0,dbg);
		);
	);
	print("time taken=",(getwalltime()-startTime));	
	return(minNUB);
}