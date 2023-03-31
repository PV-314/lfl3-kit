\\ \r lfl3\eg-bgmp-eqn29.gp

\\read("lfl3\\kit-alpha1Variable.gp");

\\ this example comes from equation (2.9) in Bennett, Gyory, Mignotte and Pinter
\\ Comp. Math. (2006)
\\ for the diophantine equation 2^\alpha*x^n-5^\beta*y^n = \pm 1
\\ where n \geq 3 is prime, 0 \leq \alpha, \beta<n
\\ \Lamnda = \alpha \log(p)-n \log (y/x)-\beta \log(q)
\\ where p=2 and q=5
\\ \Lambda = n \log (x/y)-\alpha \log(2)-\beta \log(5)

\\ L=  85, m= 13.0000, rho(3logs)=  7.5000, chi= 0.7750, K=  235935.394*logX,
\\ nonDegen log|Lambda|>-4.040789 e7*logX, nonDegenNUB=4.040789 e7,
\\ rho(2logs)= 7.000000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=4.002593 e7, degenNUB3=8.168069 e7,
\\ nUB=4.040789 e7, eliminate-b2
eg29_check_it1()={
	my(bigL,chi,m,mu,rho2Logs,rho3Logs);
	
	bigL=85;
	m=13;
	rho2Logs=7.0;
	rho3Logs=7.5;
	chi=0.775;
	mu=0.61;

	actMinNUB=eg29_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,,2);
	expMinNUB=40.40789*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg29_check_it1(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  56, m= 12.0000, rho(3logs)=  8.0000, chi= 1.0750, K=  163891.508*logX,
\\ nonDegen log|Lambda|>-1.908496 e7*logX, nonDegenNUB=1.908496 e7,
\\ rho(2logs)= 7.000000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=1.863687 e7, degenNUB3=3.826044 e7,
\\ nUB=1.908496 e7, eliminate-b2
eg29_check_it2()={
	my(bigL,chi,m,mu,nUB,rho2Logs,rho3Logs);

	bigL=56;
	m=12.0;
	rho2Logs=7.0;
	rho3Logs=8.0;
	chi=1.075;
	mu=0.61;
	nUB=41*10^6;

	actMinNUB=eg29_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=19.08496*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg29_check_it2(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  59, m= 17.0000, rho(3logs)=  6.5000, chi= 1.1000, K=  160446.094*logX,
\\ nonDegen log|Lambda|>-1.771908 e7*logX, nonDegenNUB=1.771908 e7,
\\ rho(2logs)= 7.000000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=1.754319 e7, degenNUB3=3.597471 e7,
\\ nUB=1.771908 e7, eliminate-b2
eg29_check_it3()={
	my(bigL,chi,m,mu,nUB,rho2Logs,rho3Logs);

	bigL=59;
	m=17.0;
	rho2Logs=7.0;
	rho3Logs=6.5;
	chi=1.1;
	mu=0.61;
	nUB=20*10^6;

	actMinNUB=eg29_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=17.71908*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg29_check_it3(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

eg29_search_it1(dbg=0)={
my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigD,bigK,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,m,matveevChi,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,rho,val,w);
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	bigLLB=50;
	bigLUB=150;
	mLB=7;
	mUB=mLB+10;
	rho3LB=5;
	rho3UB=rho3LB+10;
	chiLB=0.5;
	chiUB=chiLB+1.0;
	rho2LB=6;
	rho2UB=rho2LB+4;
	eg29_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,,dbg);
}

eg29_search_it2(nUBInit,dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	bigLLB=30;
	bigLUB=150;
	mLB=8;
	mUB=mLB+10;
	rho3LB=5;
	rho3UB=rho3LB+10;
	chiLB=0.7;
	chiUB=chiLB+0.5;
	rho2LB=5;
	rho2UB=rho2LB+5;
	eg29_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,nUBInit,dbg);
}

eg29_search_it3(nUBInit,dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	bigLLB=30;
	bigLUB=100;
	mLB=8;
	mUB=mLB+10;
	rho3LB=5;
	rho3UB=rho3LB+10;
	chiLB=0.7;
	chiUB=chiLB+0.5;
	rho2LB=5;
	rho2UB=rho2LB+5;
	eg29_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,nUBInit,dbg);
}

\\ only set muLB and muUB when using "check" functions above, not when using "search" functions
\\ (they are bounds for the mu in Theorem 2 of Laurent's 2008 paper)
eg29_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,muLB=0,muUB=0,nUBInit=0,dbg=0)={
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
	nLB=10*10^6;
	logXLB=log(nLB-2*sqrt(nLB)+1); \\ from Lemma 3.4 of Comp. Math. paper
	
	if(dbg>0,
		printf("calced logXLB=%5d\n",logXLB);
	);

	hgtA1=logX;
	absLogA1=log(5); \\ from Step 3 after defns of a_i's
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
		printf("eg29_search_general(): calculated nUBInit=%9.6e\n",nUBInit);
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
			printf("\nfor chi=%9.6f, minNonDegenNUB=%9.6e\n",chi,step3Result[11]);
			minNUB=alpha1_do_step4(step3Result,minNUB,d,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,chi,rho2LB,rho2UB,muLB,muUB,logXLB,nLB,lamUB1,lamUB0,dbg);
		);
	);
	print("time taken=",(getwalltime()-startTime));	
	return(minNUB);
}