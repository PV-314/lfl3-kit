\\ \r lfl3\param-search-mignotte-eg1.gp

read("lfl3\\lfl-utils-alpha1Variable.gp");

\\ search for parameter choice for Example 1 in Mignotte's original kit:
\\ for the diophantine equation x^n-2^\alpha*5^\beta*y^n = \pm 1
\\ where n \geq 3 is prime, \alpha=1,2,3, 0 \leq \beta<n
\\ \Lamnda = n \log (x/y)-\alpha \log(2)-\beta \log(5)

\\ L=  69, m=  3.5000, rho=  8.0000, chi= 3.0000, K=   45654.728*logX, nonDegen log|Lambda|>-6.550607 e6*logX, nonDegenNUB=6.550607 e6, degenNUB1=    0.e-19, degenNUB2=6.570613 e6, degenNUB3=1.296015 e7, nUB=6.570613 e6, transB-b1
\\ 4 Jan 2022
eg1_check_it1()={
	my(actMinNUB,bigL,chi,expMinNUB,m,rho);

	bigL=69;
	m=3.5;
	rho=8.0;
	chi=3.0;

	actMinNUB=eg1_search_general(bigL,bigL,m,m,rho,rho,chi,chi,,1);
	expMinNUB=6.570613*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg1_check_it1(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  46, m=  5.0000, rho=  7.0000, chi= 4.2500, K=   34139.041*logX, nonDegen log|Lambda|>-3.055849 e6*logX, nonDegenNUB=3.055849 e6, degenNUB1=    0.e-19, degenNUB2=3.066622 e6, degenNUB3=5.916155 e6, nUB=3.066622 e6, transB-b1
\\ 4 Jan 2022 (updated with above on 27 Feb 2022)
eg1_check_it2()={
	my(bigL,chi,m,nUB,rho);

	bigL=46;
	m=5.0;
	rho=7.0;
	chi=4.25;
	nUB=6.58*10^6;

	actMinNUB=eg1_search_general(bigL,bigL,m,m,rho,rho,chi,chi,nUB,1);
	expMinNUB=3.066622*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg1_check_it2(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  46, m=  4.0000, rho=  7.5000, chi= 4.2500, K=   30929.393*logX, nonDegen log|Lambda|>-2.866707 e6*logX, nonDegenNUB=2.866707 e6, degenNUB1=    0.e-19, degenNUB2=2.821734 e6, degenNUB3=5.534157 e6, nUB=2.866707 e6, transB-b1
\\ 4 Jan 2022
eg1_check_it3()={
	my(bigL,chi,m,nUB,rho);

	bigL=46;
	m=4.0;
	rho=7.5;
	chi=4.25;
	nUB=3.07*10^6;

	actMinNUB=eg1_search_general(bigL,bigL,m,m,rho,rho,chi,chi,nUB,1);
	expMinNUB=2.866707*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg1_check_it3(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  49, m=  4.0000, rho=  7.0000, chi= 4.0000, K=   29092.401*logX, nonDegen log|Lambda|>-2.773949 e6*logX, nonDegenNUB=2.773949 e6, degenNUB1=    0.e-19, degenNUB2=2.862710 e6, degenNUB3=5.614504 e6, nUB=2.862710 e6, transB-b1
\\ 4 Jan 2022
eg1_check_it4()={
	my(bigL,chi,m,nUB,rho);

	bigL=49;
	m=4.0;
	rho=7.0;
	chi=4.0;
	nUB=2.87*10^6;

	actMinNUB=eg1_search_general(bigL,bigL,m,m,rho,rho,chi,chi,nUB,1);
	expMinNUB=2.862710*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg1_check_it4(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ 4 Jan 2022
eg1_search_it1(dbg=0)={
my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigD,bigK,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,m,matveevChi,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,rho,val,w);

	bigLLB=60;
	bigLUB=350;
	mLB=1;
	mUB=6;
	rhoLB=2;
	rhoUB=12;
	chiLB=1;
	chiUB=6;
	eg1_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB);
}

\\ 4 Jan 2022 (27 Feb 2022: should use nUBInit=6.58*10^6 -- upper bound from iteration 1)
eg1_search_it2(nUBInit,dbg=0)={
	my(al1,al2,al3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,degenNUB,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,minNUB,nLB,nonDegenNUB,nUB,rho,val);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	bigLLB=40;
	bigLUB=100;
	mLB=1;
	mUB=11;
	rhoLB=4;
	rhoUB=14;
	chiLB=1;
	chiUB=6;
	eg1_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,nUBInit);
}

\\ 4 Jan 2022 (should use nUBInit=3.07*10^6 -- upper bound from iteration 2)
eg1_search_it3(nUBInit,dbg=0)={
	my(al1,al2,al3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,degenNUB,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,minNUB,nLB,nonDegenNUB,nUB,rho,val);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	bigLLB=30;
	bigLUB=100;
	mLB=1;
	mUB=11;
	rhoLB=1;
	rhoUB=11;
	chiLB=2.0;
	chiUB=7.0;
	eg1_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,nUBInit);
}

\\ 4 Jan 2022 (should use nUBInit=2.87*10^6 -- upper bound from iteration 3)
eg1_search_it4(nUBInit,dbg=0)={
	my(al1,al2,al3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,degenNUB,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,minNUB,nLB,nonDegenNUB,nUB,rho,val);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);

	bigLLB=30;
	bigLUB=100;
	mLB=1;
	mUB=11;
	rhoLB=1;
	rhoUB=11;
	chiLB=2.0;
	chiUB=7.0;
	eg1_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,nUBInit);
}

\\ 30 June 2022
eg1_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,nUBInit=0,dbg=0)={
my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigD,bigK,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,matveevChi,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,val,w);

	\\ bigD=[Q(al_1,al_2,al_3):Q] -- used for Matveev's bounds
	bigD=1;
	\\ matveevChi=[R(al_1,al_2,al_3):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2,al_3):Q]/[R(al_1,al_2,al_3):R]
	d=bigD/matveevChi; \\ d is the "degree" value used in the kit

	al1=x;
	al2=2;
	al3=5;

	b1=n;
	b2=3; \\ alpha=1, 2 or 3
	b3=n; \\ 0 \leq beta<n
	nLB=10^6; \\ assumption at the start of proof of Theorem 6.3
	logXLB=floor(nLB/2600-1.5);
	if(dbg!=0,
		printf("calced logXLB=%5d\n",logXLB);
	);
	logXLB=380;

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
		nUBInit=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,1);
		printf("nUBInit=%9.6e\n",nUBInit);
		printf("eg1_search_general(): calculated nUBInit=%9.6e\n",nUBInit);
		\\printf("                 but reverting to nUBInit=%9.6e\n\n",nUBInit);
	);
	minNUB=nUBInit;

	areBoundsOK=check_bounds(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB);
	if(areBoundsOK==0,
		return();
	);
	mStep=(mUB-mLB)/20.0;
	if(mStep==0,mStep=0.000001);
	rhoStep=(rhoUB-rhoLB)/20.0;
	if(rhoStep==0,rhoStep=0.000001);
	chiStep=(chiUB-chiLB)/20.0;
	if(chiStep==0,chiStep=0.000001);

	for(bigL=bigLLB,bigLUB,
	if(bigL%10==0,print("L=",bigL));
	forstep(m=mLB,mUB,mStep,
		forstep(rho=rhoLB,rhoUB,rhoStep,
			a1=(rho-1)*absLogA1+2*logX;
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
