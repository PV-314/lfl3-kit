\\ \r lfl3\eg-bms-compositio.gp

read("lfl3\\kit-alpha1Variable.gp");

\\ this is the Case I linear form in the 2006 Compositio paper of Yann, Maurice and Samir
\\ but we have switched alpha_1 and alpha_2

\\ L= 106, m= 21.0000, rho(3logs)=  5.5000, chi= 0.0800, K=  231548.100*logX, nonDegen log|Lambda|>-4.184151 e7*logX,
\\ nonDegenNUB=8.368301 e7, rho(2logs)=101.100000, mu(2logs)= 0.604000, degenNUB1=    0.e-19, degenNUB2=6.622400 e7, degenNUB3=6.622400 e7, nUB=8.368301 e7, transB-b3
\\ 12 Feb 2023
bms2_check_it1()={
	my(bigL,chi,m,mu,rho2Logs,rho3Logs);
	
	bigL=106;
	m=21.0;
	rho2Logs=180;
	rho3Logs=5.5;
	chi=0.08;
	mu=0.61;

	actMinNUB=bms2_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,,2);
	expMinNUB=83.68301*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: bms2_check_it1(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  59, m= 18.0000, rho(3logs)=  6.0000, chi= 0.1000, K=  135524.397*logX, nonDegen log|Lambda|>-1.432680 e7*logX,
\\ nonDegenNUB=2.865360 e7, rho(2logs)=101.100000, mu(2logs)= 0.604000, degenNUB1=    0.e-19, degenNUB2=2.640695 e7, degenNUB3=2.640695 e7, nUB=2.865360 e7, transB-b3
\\ second iteration values
\\ 12 Feb 2023
bms2_check_it2()={
	my(bigL,chi,m,mu,nUB,rho2Logs,rho3Logs);
	
	bigL=59;
	m=18.0;
	rho2Logs=180;
	rho3Logs=6.0;
	chi=0.1;
	mu=0.61;
	nUB=84*10^6;

	actMinNUB=bms2_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=28.65360*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: bms2_check_it1(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  59, m= 18.0000, rho(3logs)=  5.7500, chi= 0.1000, K=  122602.952*logX, nonDegen log|Lambda|>-1.265297 e7*logX,
\\ nonDegenNUB=2.530593 e7, rho(2logs)=150.000000, mu(2logs)= 0.604000, degenNUB1=    0.e-19, degenNUB2=2.346161 e7, degenNUB3=2.346161 e7, nUB=2.530593 e7, transB-b3
\\ third iteration values
\\ 12 Feb 2023
bms2_check_it3()={
	my(bigL,chi,m,mu,nUB,rho2Logs,rho3Logs);
	
	bigL=59;
	m=18.0;
	rho2Logs=180;
	rho3Logs=5.75;
	chi=0.0975;
	mu=0.61;
	nUB=29*10^6;

	actMinNUB=bms2_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=25.30593*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: bms2_check_it3(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  57, m= 19.0000, rho(3logs)=  5.7500, chi= 0.1000, K=  125027.304*logX, nonDegen log|Lambda|>-1.246577 e7*logX,
\\ nonDegenNUB=2.493154 e7, rho(2logs)=150.000000, mu(2logs)= 0.604000, degenNUB1=    0.e-19, degenNUB2=2.346161 e7, degenNUB3=2.346161 e7, nUB=2.493154 e7, transB-b3
\\ fourth iteration values
\\ 12 Feb 2023
bms2_check_it4()={
	my(actMinNUB,bigL,chi,expMinNUB,m,mu,nUB,rho2Logs,rho3Logs);
	
	bigL=57;
	m=19.0;
	rho2Logs=180;
	rho3Logs=5.75;
	chi=0.1;
	nUB=25.4*10^6;

	actMinNUB=bms2_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=24.93154*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: bms2_check_it4(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ for iteration 1
\\ 11 May 2022
bms2_search_it1(dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	bigLLB=30;
	bigLUB=200;
	mLB=10;
	mUB=mLB+20;
	rho3LB=3;
	rho3UB=rho3LB+10;
	chiLB=0.04;
	chiUB=chiLB+0.2;
	rho2LB=80;
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

	bigLLB=30;
	bigLUB=200;
	mLB=5;
	mUB=mLB+20;
	rho3LB=2;
	rho3UB=rho3LB+10;
	chiLB=0.05;
	chiUB=chiLB+0.2;
	rho2LB=80;
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

	bigLLB=30;
	bigLUB=120;
	mLB=11;
	mUB=mLB+10;
	rho3LB=4;
	rho3UB=rho3LB+5;
	chiLB=0.07;
	chiUB=chiLB+0.05;
	rho2LB=80;
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

	bigLLB=30;
	bigLUB=120;
	mLB=11;
	mUB=mLB+10;
	rho3LB=4;
	rho3UB=rho3LB+5;
	chiLB=0.07;
	chiUB=chiLB+0.05;
	rho2LB=80;
	rho2UB=rho2LB+100;

	bms2_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,nUBInit,dbg);
}

\\ only set muLB and muUB when using "check" functions above, not when using "search" functions
\\ (they are bounds for the mu in Theorem 2 of Laurent's 2008 paper)
\\ 3 July 2022
bms2_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,muLB=0,muUB=0,nUBInit=0,dbg=0)={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,areBoundsOK,b1,b2,b3,bigD,bigK,chiStep,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,matveevChi,minNUB,mStep,nDegenUB,nLB,nNonDegenUB,nUB,rho3Step,startTime,step3Result);

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
					a1=rho3*Pi/2+logX;
					a2=0.723*rho3;
					a3=rho3*Pi;
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