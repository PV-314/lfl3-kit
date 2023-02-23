\\ \r lfl3\eg-bms-compositio.gp

read("lfl3\\kit-alpha1Variable.gp");

\\ this is the Case I linear form in the 2006 Compositio paper of Yann, Maurice and Samir
\\ but we have switched alpha_1 and alpha_2

\\ L= 108, m= 12.0000, rho(3logs)=  6.5000, chi= 0.0700, K=  229320.984*logX, nonDegen log|Lambda|>-4.635830 e7*logX,
\\ nonDegenNUB=9.271660 e7, rho(2logs)=200.000000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=6.946463 e7, degenNUB3=6.946463 e7, nUB=9.271660 e7, eliminate-b3
\\ 12 Feb 2023
bms2_check_it1()={
	my(bigL,chi,m,mu,rho2Logs,rho3Logs);
	
	bigL=108;
	m=12.0;
	rho2Logs=200;
	rho3Logs=6.5;
	chi=0.07;
	mu=0.61;

	actMinNUB=bms2_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,,2);
	expMinNUB=92.71660*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: bms2_check_it1(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  61, m= 13.0000, rho(3logs)=  6.5000, chi= 0.0850, K=  140317.547*logX, nonDegen log|Lambda|>-1.602145 e7*logX,
\\ nonDegenNUB=3.204290 e7, rho(2logs)=200.000000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=3.186990 e7, degenNUB3=3.186990 e7, nUB=3.204290 e7, eliminate-b3
\\ second iteration values
\\ 12 Feb 2023
bms2_check_it2()={
	my(bigL,chi,m,mu,nUB,rho2Logs,rho3Logs);
	
	bigL=61;
	m=13.0;
	rho2Logs=200;
	rho3Logs=6.5;
	chi=0.085;
	mu=0.61;
	nUB=93*10^6;

	actMinNUB=bms2_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=32.04290*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: bms2_check_it1(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  57, m= 16.5000, rho(3logs)=  6.0000, chi= 0.0875, K=  139152.681*logX, nonDegen log|Lambda|>-1.421170 e7*logX,
\\ nonDegenNUB=2.842341 e7, rho(2logs)=200.000000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=2.791274 e7, degenNUB3=2.791274 e7, nUB=2.842341 e7, eliminate-b3
\\ third iteration values
\\ 12 Feb 2023
bms2_check_it3()={
	my(bigL,chi,m,mu,nUB,rho2Logs,rho3Logs);
	
	bigL=57;
	m=16.5;
	rho2Logs=200;
	rho3Logs=6.0;
	chi=0.0875;
	mu=0.61;
	nUB=33*10^6;

	actMinNUB=bms2_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=28.42341*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: bms2_check_it3(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  55, m= 15.5000, rho(3logs)=  6.2500, chi= 0.0875, K=  138145.397*logX, nonDegen log|Lambda|>-1.392395 e7*logX,
\\ nonDegenNUB=2.784790 e7, rho(2logs)=250.000000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=2.749874 e7, degenNUB3=2.749874 e7, nUB=2.784790 e7, eliminate-b3
\\ fourth iteration values
\\ 12 Feb 2023
bms2_check_it4()={
	my(actMinNUB,bigL,chi,expMinNUB,m,mu,nUB,rho2Logs,rho3Logs);
	
	bigL=55;
	m=15.5;
	rho2Logs=250;
	rho3Logs=6.25;
	chi=0.0875;
	mu=0.61;
	nUB=28.5*10^6;

	actMinNUB=bms2_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=27.84790*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: bms2_check_it4(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ for iteration 1
\\ 11 May 2022
bms2_search_it1(dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	bigLLB=50;
	bigLUB=200;
	mLB=5;
	mUB=mLB+20;
	rho3LB=3;
	rho3UB=rho3LB+10;
	chiLB=0.02;
	chiUB=chiLB+0.1;
	rho2LB=100;
	rho2UB=rho2LB+100;

	bms2_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,,dbg);
}

\\ for iteration 2
\\ 11 May 2022
bms2_search_it2(nUBInit,dbg=0)={ \\use nUBInit=2.03*10^7
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);

	bigLLB=40;
	bigLUB=200;
	mLB=5;
	mUB=mLB+20;
	rho3LB=3;
	rho3UB=rho3LB+10;
	chiLB=0.02;
	chiUB=chiLB+0.1;
	rho2LB=100;
	rho2UB=rho2LB+100;

	bms2_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,nUBInit,dbg);
}

\\ to be used for iteration 3
\\ 11 May 2022
bms2_search_it3(nUBInit,dbg=0)={ \\use nUBInit=7.98*10^6
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);

	bigLLB=40;
	bigLUB=150;
	mLB=9;
	mUB=mLB+10;
	rho3LB=4;
	rho3UB=rho3LB+5;
	chiLB=0.05;
	chiUB=chiLB+0.05;
	rho2LB=100;
	rho2UB=rho2LB+100;

	bms2_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,nUBInit,dbg);
}

\\ to be used for iteration 4
\\ 11 May 2022
bms2_search_it4(nUBInit,dbg=0)={ \\use nUBInit=7.28*10^6
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);

	bigLLB=40;
	bigLUB=150;
	mLB=10;
	mUB=mLB+10;
	rho3LB=4;
	rho3UB=rho3LB+5;
	chiLB=0.05;
	chiUB=chiLB+0.05;
	rho2LB=150;
	rho2UB=rho2LB+100;

	bms2_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,nUBInit,dbg);
}

\\ only set muLB and muUB when using "check" functions above, not when using "search" functions
\\ (they are bounds for the mu in Theorem 2 of Laurent's 2008 paper)
\\ 3 July 2022
bms2_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,muLB=0,muUB=0,nUBInit=0,dbg=0)={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA2a,absLogA2b,absLogA3,al1,al2,al3,areBoundsOK,b1,b2,b3,bigD,chiStep,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,matveevChi,minNUB,mStep,nLB,rho3Step,startTime,step3Result);

	startTime=getwalltime();
	\\ bigD=[Q(al_1,al_2,al_3):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2,al_3):R] -- used for Matveev's bounds
	matveevChi=2;
	\\ remember that d=[Q(al_1,al_2,al_3):Q]/[R(al_1,al_2,al_3):R]
	d=bigD/matveevChi; \\ d is the "degree" value used in the kit

	al1=x; \\ this is just a placeholder so that al1=alpha_1 is considered to be a polynomial
	hgtA1=logX/2; \\ this must be correct though and consistent with logX(=log(x)) usages that follow
	absLogA1=Pi/2;

	al2=(1+sqrt(-7))/(1-sqrt(-7));
	hgtA2=log(2)/2;
	absLogA2a=abs(log(al2));
	absLogA2b=abs(log(-al2));
	absLogA2=min(absLogA2a,absLogA2b);
	print("absLogA2=",absLogA2);
	
	al3=-1;
	hgtA3=0;
	absLogA3=Pi;

	b1=n;
	b2=2;
	b3=n-1; \\ since b3=q with |q|<p
	nLB=20*10^6;
	logXLB=log((sqrt(nLB)-1)^2);
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX/2;
	lamUB0=log(2.2*sqrt(7));

	if(nUBInit==0,
		\\print("bigD=",bigD,", matveevChi=",matveevChi,", al1=",al1,", absLogA1=",absLogA1,", hgtA1=",hgtA1,", al2=",al2);
		\\print("absLogA2=",absLogA2,", hgtA2=",hgtA2,", al3=",al3,", absLogA3=",absLogA3,", hgtA3=",hgtA3);
		\\print("b1=",b1,", b2=",b2,b3,", logXLB=",logXLB,", nLB=",nLB,", lamUB1=",lamUB1,", lamUB0=",lamUB0,", dbg=",dbg);
		nUBInit=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,dbg);
		printf("nUBInit=%9.6e\n",nUBInit);
		printf("bms2_search_general(): calculated nUBInit=%9.6e\n",nUBInit);
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
		\\printf("chi=%9.6f\n",chi);
		for(bigL=bigLLB,bigLUB, \\ L=5 is the lower bound in Theorem 4.1
			\\if(bigL%10==0,print("L=",bigL));
			forstep(m=mLB,mUB,mStep,
				forstep(rho3=rho3LB,rho3UB,rho3Step,
					a1=rho3*absLogA1+2*d*hgtA1;
					a2=rho3*absLogA2+2*d*hgtA2;
					a3=rho3*absLogA3+2*d*hgtA3;
					step3Result=alpha1_do_step3(step3Result,d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,chi,bigL,m,rho3,nUBInit,logXLB,nLB,lamUB1,lamUB0,dbg);
				);
			);
		);
		if(length(step3Result)>0 && step3Result[11]<minNUB,
			printf("\nfor chi=%9.6f, minNonDegenNUB=%9.6e\n",chi,step3Result[11]);
			minNUB=alpha1_do_step4(step3Result,minNUB,d,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,chi,rho2LB,rho2UB,muLB,muUB,logXLB,nLB,lamUB1,lamUB0,dbg);
		);
	);
	print("time taken=",(getwalltime()-startTime));	
	return(minNUB);
}