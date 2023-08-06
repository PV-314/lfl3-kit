\\ \r lfl3\eg-bdms.gp

read("lfl3\\kit-alpha3Variable.gp");

\\ this example comes from after Theorem 4 in Bennett, Dahmen, Mignotte and Siksek's paper
\\ Shifted powers in binary recurrence sequences
\\ in Math Proc Cambridge Phil Soc, Volume 158, March 2015, pp. 305--329.
\\ they write \Lambda =p*log |x| - log(sqrt(2)-k*log(1+sqrt(2)) (they use k*log(eps), where eps=1+sqrt(2))
\\ to put it in our form, we will write
\\ \Lambda =log(sqrt(2)+k*log(1+sqrt(2))-p*log |x|

\\ for chi= 0.750000, minNonDegenNUB=1.841827 e8
\\ L= 231, m= 12.5000, rho(3logs)= 5.7500, chi= 0.7500, K=  455824.190*logX,
\\ nonDegen log|Lambda|>-1.841827 e8*logX, nonDegenNUB=1.841827 e8,
\\ rho(2logs)= 7.0000, mu(2logs)= 0.6100,
\\ degenNUB1=3.726314 e8, degenNUB2=1.816269 e8, degenNUB3=    0.e-19,
\\ nUB=1.841827 e8, eliminate-b2\\ 5 August 2023
eg_bdms_check_it1()={
	my(bigL,chi,m,mu,rho2Logs,rho3Logs);
	
	bigL=231;
	m=12.5;
	rho2Logs=7.0;
	rho3Logs=5.75;
	chi=0.75;
	mu=0.61;

	actMinNUB=eg_bdms_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,,2);
	expMinNUB=184.1827*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg_bdms_check_it1(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ for chi= 1.020000, minNonDegenNUB=9.372721 e7
\\ L= 153, m= 14.5000, rho(3logs)= 5.7500, chi= 1.0200, K=  350215.053*logX,
\\ nonDegen log|Lambda|>-9.372721 e7*logX, nonDegenNUB=9.372721 e7,
\\ rho(2logs)= 7.0000, mu(2logs)= 0.6100,
\\ degenNUB1=1.927057 e8, degenNUB2=9.392528 e7, degenNUB3=    0.e-19,
\\ nUB=9.392528 e7, eliminate-b2
eg_bdms_check_it2()={
	my(bigL,chi,m,mu,nUB,rho2Logs,rho3Logs);

	bigL=153;
	m=14.5;
	rho2Logs=7.0;
	rho3Logs=5.75;
	chi=1.02;
	mu=0.61;
	nUB=185*10^6;

	actMinNUB=eg_bdms_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=93.92528*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg29_check_it2(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

eg_bdms_check_paper()={
	my(bigL,chi,m,mu,rho2Logs,rho3Logs);
	
	bigL=545;
	m=25.0;
	rho2Logs=283.0;
	rho3Logs=5;
	chi=2;
	mu=0.6;

	actMinNUB=eg_bdms_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,,2);
	expMinNUB=6.499120*10^9;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg_bdms_check_it1(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

eg_bdms_search_it1(dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	bigLLB=150;
	bigLUB=300;
	mLB=10;
	mUB=mLB+5;
	rho3LB=3;
	rho3UB=rho3LB+5;
	chiLB=0.65;
	chiUB=chiLB+0.2;
	rho2LB=5;
	rho2UB=rho2LB+4;
	eg_bdms_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,,dbg);
}

\\ use nUBInit=185*10^6
eg_bdms_search_it2(nUBInit,dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	bigLLB=100;
	bigLUB=200;
	mLB=10;
	mUB=mLB+5;
	rho3LB=3;
	rho3UB=rho3LB+5;
	chiLB=0.7;
	chiUB=chiLB+0.4;
	rho2LB=5;
	rho2UB=rho2LB+4;
	eg_bdms_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,nUBInit,dbg);
}

\\ only set muLB and muUB when using "check" functions above, not when using "search" functions
\\ (they are bounds for the mu in Theorem 2 of Laurent's 2008 paper)
eg_bdms_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,muLB=0,muUB=0,nUBInit=0,dbg=0)={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,areBoundsOK,b1,b2,b3,bigD,bigK,chiStep,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,matveevChi,minNUB,mStep,nDegenUB,nLB,nNonDegenUB,nUB,rho3Step,startTime,step3Result);

	startTime=getwalltime();
	\\ bigD=[Q(al_1,al_2,al_3):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2,al_3):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2,al_3):Q]/[R(al_1,al_2,al_3):R]
	d=bigD/matveevChi; \\ d is the "degree" value used in the kit

	\\ \Lambda =log(sqrt(2)) + k*log(1+sqrt(2)) - p*log |x|

	al1=1+sqrt(2);
	al2=sqrt(2);
	al3=x; \\ x is just because it is not known.

	b1=n*logX/log(al1);
	b2=1;
	b3=n;
	nLB=20*10^6;
	logXLB=log(7);
	
	if(dbg>0,
		printf("calced logXLB=%9.6f\n",logXLB);
	);

	hgtA1=log(al1)/2;
	absLogA1=abs(log(al1));
	hgtA2=2*log(al2)/2;
	absLogA2=abs(log(al2));
	hgtA3=logX;
	absLogA3=logX; \\ from Step 3 after defns of a_i's
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	\\ from equation (12) of the paper
	lamUB1=-logX;
	lamUB0=2;

	if(nUBInit==0,
		\\print("bigD=",bigD,", matveevChi=",matveevChi,", al1=",al1,", absLogA1=",absLogA1,", hgtA1=",hgtA1,", al2=",al2);
		\\print("absLogA2=",absLogA2,", hgtA2=",hgtA2,", al3=",al3,", absLogA3=",absLogA3,", hgtA3=",hgtA3);
		\\print("b1=",b1,", b2=",b2,b3,", logXLB=",logXLB,", nLB=",nLB,", lamUB1=",lamUB1,", lamUB0=",lamUB0,", dbg=",dbg);
		nUBInit=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,dbg);
		printf("nUBInit=%9.6e\n",nUBInit);
		printf("eg_bdms_search_general(): calculated nUBInit=%9.6e\n",nUBInit);
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
					a3=(rho3+3)*logX;
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