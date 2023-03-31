\\ \r lfl3\eg-mignotte2.gp

read("lfl3\\kit-alpha1Variable.gp");

\\ search for parameter choice for Example 2 in Mignotte's original kit:
\\ for the diophantine equation 2^\alpha*x^n-5^\beta*y^n = \pm 1
\\ where n \geq 3 is prime, 0 \leq \alpha, \beta<n
\\ \Lamnda = \alpha \log(p)-n \log (y/x)-\beta \log(q)
\\ where p=2 and q=5
\\ this example is a generalisation of equation (2.9) in Bennett, Gyory, Mignotte and Pinter
\\ Comp. Math. (2006)
\\ here alpha<n, but in the Comp Math paper, 2 \leq \alpha \leq 3
\\ see the proof of their Prop 3.5.
\\ but see eg-bgmp-eqn29.gp for code for equation (2.9) itself.

\\ L=  90, m= 16.0000, rho(3logs)=  7.0000, chi= 1.2500, K=  267222.598*logX,
\\ nonDegen log|Lambda|>-4.679920 e7*logX, nonDegenNUB=4.679920 e7,
\\ rho(2logs)= 5.500000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=9.296262 e7, degenNUB3=4.553920 e7,
\\ nUB=4.679920 e7, eliminate-b3
eg2_check_it1()={
	my(bigL,chi,m,mu,rho2Logs,rho3Logs);

	bigL=90;
	m=16;
	rho2Logs=5.5;
	rho3Logs=7.0;
	chi=1.25;
	mu=0.61;

	actMinNUB=eg2_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,,2);
	expMinNUB=46.79920*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg2_check_it1(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  60, m= 18.0000, rho(3logs)=  7.0000, chi= 1.7000, K=  200416.948*logX,
\\ nonDegen log|Lambda|>-2.339960 e7*logX, nonDegenNUB=2.339960 e7,
\\ rho(2logs)= 5.400000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=4.774301 e7, degenNUB3=2.336776 e7,
\\ nUB=2.339960 e7, eliminate-b3
eg2_check_it2()={
	my(bigL,chi,m,mu,nUB,rho2Logs,rho3Logs);

	bigL=60;
	m=18;
	rho2Logs=5.4;
	rho3Logs=7;
	chi=1.7;
	mu=0.61;
	nUB=47*10^6;

	actMinNUB=eg2_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=23.39960*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg2_check_it2(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  61, m= 17.5000, rho(3logs)=  6.8000, chi= 1.7000, K=  186869.249*logX,
\\ nonDegen log|Lambda|>-2.185105 e7*logX, nonDegenNUB=2.185105 e7,
\\ rho(2logs)= 5.400000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=4.427621 e7, degenNUB3=2.187998 e7,
\\ nUB=2.187998 e7, eliminate-b3
eg2_check_it3()={
	my(bigL,chi,m,mu,nUB,rho2Logs,rho3Logs);

	bigL=61;
	m=17.5;
	rho2Logs=5.4;
	rho3Logs=6.8;
	chi=1.7;
	mu=0.61;
	nUB=24*10^6;

	actMinNUB=eg2_search_general(bigL,bigL,m,m,rho3Logs,rho3Logs,chi,chi,rho2Logs,rho2Logs,mu,mu,nUB,2);
	expMinNUB=21.87998*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg2_check_it3(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

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

eg2_search_it2(nUBInit,dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);

	bigLLB=50;
	bigLUB=150;
	mLB=12;
	mUB=mLB+10;
	rho3LB=5;
	rho3UB=rho3LB+5;
	chiLB=0.8;
	chiUB=chiLB+1.0;
	rho2LB=4;
	rho2UB=rho2LB+4;
	eg2_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,nUBInit,dbg);
}

eg2_search_it3(nUBInit,dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rho2LB,rho2UB,rho3LB,rho3UB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);

	bigLLB=50;
	bigLUB=150;
	mLB=12;
	mUB=mLB+10;
	rho3LB=5;
	rho3UB=rho3LB+4;
	chiLB=1.0;
	chiUB=chiLB+1.0;
	rho2LB=4;
	rho2UB=rho2LB+4;
	eg2_search_general(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB,,,nUBInit,dbg);
}

\\ only set muLB and muUB when using "check" functions above, not when using "search" functions
\\ (they are bounds for the mu in Theorem 2 of Laurent's 2008 paper)
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
	nLB=10*10^6; \\ assumption at the start of proof of Theorem 6.3
	logXLB=log(nLB-2*sqrt(nLB)+1); \\ from Lemma 3.4 of Comp. Math. paper
	
	if(dbg>0,
		printf("calced logXLB=%5d\n",logXLB);
	);

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
			printf("\nfor chi=%9.6f, minNonDegenNUB=%9.6e\n",chi,step3Result[11]);
			minNUB=alpha1_do_step4(step3Result,minNUB,d,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,chi,rho2LB,rho2UB,muLB,muUB,logXLB,nLB,lamUB1,lamUB0,dbg);
		);
	);
	print("time taken=",(getwalltime()-startTime));	
	return(minNUB);
}