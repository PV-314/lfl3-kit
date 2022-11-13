\\ \r lfl3\param-search-mignotte-eg2.gp

read("lfl3\\lfl-utils-alpha1Variable.gp");

\\ search for parameter choice for Example 2 in Mignotte's original kit:
\\ for the diophantine equation 2^\alpha*x^n-5^\beta*y^n = \pm 1
\\ where n \geq 3 is prime, 0 \leq \alpha, \beta<n
\\ \Lamnda = \alpha \log(p)-n \log (y/x)-\beta \log(q)
\\ where p=2 and q=5

\\ Mignotte's values
\\ 25 Nov 2021
eg2_check1()={
	my(bigL,chi,m,nUB,rho);

	bigL=110;
	m=41.28955;
	rho=7.0;
	chi=1.0;
	nUB=5.7*10^11;

	actMinNUB=eg2_search_general(bigL,bigL,m,m,rho,rho,chi,chi,nUB,1);
	expMinNUB=4.832797*10^8;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg2_check1(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  86, m=  7.0000, rho=  5.0000, chi= 1.8000, K=  128339.867*logX, nonDegen log|Lambda|>-1.776373 e7*logX, nonDegenNUB=1.776373 e7, degenNUB1=    0.e-19, degenNUB2=3.862020 e7, degenNUB3=1.755363 e7, nUB=1.776373 e7, transB-b1
\\ 25 Nov 2021
eg2_check_it1()={
	my(bigL,chi,m,nUB,rho);

	bigL=86;
	m=7.0;
	rho=5.0;
	chi=1.8;
	\\nUB=5.7*10^11;

	actMinNUB=eg2_search_general(bigL,bigL,m,m,rho,rho,chi,chi,,1);
	expMinNUB=1.776373*10^7;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg2_check_it1(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  62, m=  9.0000, rho=  4.5000, chi= 2.6000, K=   92171.776*logX, nonDegen log|Lambda|>-8.595276 e6*logX, nonDegenNUB=8.595276 e6, degenNUB1=    0.e-19, degenNUB2=1.752490 e7, degenNUB3=7.928707 e6, nUB=8.595276 e6, transB-b1
\\ 9 July 2022
eg2_check_it2()={
	my(bigL,chi,m,nUB,rho);

	bigL=62;
	m=9.0;
	rho=4.5;
	chi=2.6;
	nUB=1.78*10^7;

	actMinNUB=eg2_search_general(bigL,bigL,m,m,rho,rho,chi,chi,nUB,1);
	expMinNUB=8.595276*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg2_check_it2(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  58, m=  9.5000, rho=  4.5000, chi= 2.6000, K=   91015.499*logX, nonDegen log|Lambda|>-7.939873 e6*logX, nonDegenNUB=7.939873 e6, degenNUB1=    0.e-19, degenNUB2=1.748558 e7, degenNUB3=7.910565 e6, nUB=7.939873 e6, transB-b1
\\ 9 July 2022
eg2_check_it3()={
	my(bigL,chi,m,nUB,rho);

	bigL=58;
	m=9.5;
	rho=4.5;
	chi=2.6;
	nUB=8.6*10^6;

	actMinNUB=eg2_search_general(bigL,bigL,m,m,rho,rho,chi,chi,nUB,1);
	expMinNUB=7.939873*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg2_check_it3(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ same as it3: L=  58, m=  9.5000, rho=  4.5000, chi= 2.6000, K=   91015.499*logX, nonDegen log|Lambda|>-7.939873 e6*logX, nonDegenNUB=7.939873 e6, degenNUB1=    0.e-19, degenNUB2=1.748558 e7, degenNUB3=7.910565 e6, nUB=7.939873 e6, transB-b1
\\ 9 July 2022
eg2_check_it4()={
	my(bigL,chi,m,nUB,rho);

	bigL=58;
	m=9.5;
	rho=4.5;
	chi=2.6;
	nUB=8.6*10^6;

	actMinNUB=eg2_search_general(bigL,bigL,m,m,rho,rho,chi,chi,nUB,1);
	expMinNUB=7.939873*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg2_check_it4(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  86, m=  7.0000, rho=  5.0000, chi= 1.8000, K=  128339.867*logX, nonDegen log|Lambda|>-1.776373 e7*logX, nonDegenNUB=1.776373 e7, degenNUB1=    0.e-19, degenNUB2=3.862020 e7, degenNUB3=1.755363 e7, nUB=1.776373 e7, transB-b1
\\ 25 Nov 2021
eg2_search_it1(dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rhoLB,rhoUB);

	bigLLB=50;
	bigLUB=350;
	mLB=2;
	mUB=22;
	rhoLB=2;
	rhoUB=22;
	chiLB=1.0;
	chiUB=3.0;
	eg2_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,,dbg);
}

\\ 4 Jan 2022 (should use nUBInit=1.78*10^7 -- upper bound from iteration 1)
eg2_search_it2(nUBInit,dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rhoLB,rhoUB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);

	bigLLB=50;
	bigLUB=350;
	mLB=2;
	mUB=22;
	rhoLB=2;
	rhoUB=12;
	chiLB=1.0;
	chiUB=5.0;
	eg2_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,nUBInit,dbg);
}

\\ 5 Jan 2022 (should use nUBInit=8.6*10^6 -- upper bound from iteration 2)
eg2_search_it3(nUBInit,dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rhoLB,rhoUB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);

	bigLLB=40;
	bigLUB=350;
	mLB=2;
	mUB=12;
	rhoLB=2;
	rhoUB=12;
	chiLB=1.0;
	chiUB=5.0;
	eg2_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,nUBInit,dbg);
}

\\ 5 Jan 2022 (should use nUBInit=7.94*10^6 -- upper bound from iteration 3)
eg2_search_it4(nUBInit,dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rhoLB,rhoUB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);

	bigLLB=40;
	bigLUB=350;
	mLB=2;
	mUB=12;
	rhoLB=2;
	rhoUB=12;
	chiLB=1.0;
	chiUB=5.0;
	eg2_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,nUBInit,dbg);
}

\\ 9 July June 2022
eg2_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,nUBInit=0,dbg=0)={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigD,bigK,chiStep,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,matveevChi,minNUB,mStep,nDegenUB,nLB,nNonDegenUB,nUB,rhoStep,val,w);

	\\ bigD=[Q(al_1,al_2,al_3):Q] -- used for Matveev's bounds
	bigD=1;
	\\ matveevChi=[R(al_1,al_2,al_3):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2,al_3):Q]/[R(al_1,al_2,al_3):R]
	d=bigD/matveevChi; \\ d is the "degree" value used in the kit

	al1=x; \\ actually x/y (or y/x)
	al2=5;
	al3=2;
	b1=n; \\ 0 \leq alpha < n
	b2=n;
	b3=n; \\ 0 \leq beta < n
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
		nUBInit=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,1);
		printf("nUBInit=%9.6e\n",nUBInit);
		printf("eg2_search_general(): calculated nUBInit=%9.6e\n",nUBInit);
		\\printf("                 but reverting to nUBInit=%9.6e\n\n",nUBInit);
	);
	minNUB=nUBInit;

	areBoundsOK=check_bounds(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB);
	if(areBoundsOK==0,
		return();
	);
	chiStep=get_step(chiLB,chiUB);
	mStep=get_step(mLB,mUB);
	rhoStep=get_step(rhoLB,rhoUB);

	for(bigL=bigLLB,bigLUB, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%10==0,print("L=",bigL));
	forstep(m=mLB,mUB,mStep,
		forstep(rho=rhoLB,rhoUB,rhoStep,
			a1=(rho-1)*log(5.0001)+2*logX;
			a2=(rho+1)*log(al2);
			a3=(rho+1)*log(al3);
			forstep(chi=chiLB,chiUB,chiStep,
				val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
				minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
			);
		);
	);
	);
	return(minNUB);
}
