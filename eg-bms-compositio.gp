\\ \r lfl3\param-search-bms-compositio-eg.gp

read("lfl3\\lfl-utils-alpha1Variable.gp");

\\ this is the Case I linear form in the 2006 Compositio paper of Yann, Maurice and Samir
\\ but we have switched alpha_1 and alpha_2

\\ L=  98, m=  8.0000, rho=  5.0000, chi= 3.0000, K=   64242.779*logX, nonDegen log|Lambda|>-1.013269 e7*logX, nonDegenNUB=2.026537 e7, degenNUB1=    0.e-19, degenNUB2=1.733578 e7, degenNUB3=3.748719 e7, nUB=2.026537 e7, transB-b1
\\ 3 July 2022
eg4_check_it1()={
	my(bigL,chi,m,rho);
	
	bigL=98;
	m=8.0;
	rho=5.0;
	chi=3.0;

	actMinNUB=eg4_search_general(bigL,bigL,m,m,rho,rho,chi,chi,,1);
	expMinNUB=2.026537*10^7;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg4_check_it1(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  61, m=  8.0000, rho=  5.0000, chi= 4.0000, K=   39987.852*logX, nonDegen log|Lambda|>-3.925836 e6*logX, nonDegenNUB=7.851672 e6, degenNUB1=    0.e-19, degenNUB2=7.970023 e6, degenNUB3=1.709112 e7, nUB=7.970023 e6, transB-b1
\\ second iteration values
\\ 11 May 2022
eg4_check_it2()={
	my(bigL,chi,m,nUB,rho);
	
	bigL=61;
	m=8.0;
	rho=5.0;
	chi=4.0;
	nUB=2.03*10^7;

	actMinNUB=eg4_search_general(bigL,bigL,m,m,rho,rho,chi,chi,nUB,1);
	expMinNUB=7.970023*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg4_check_it2(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  50, m= 11.0000, rho=  5.0000, chi= 4.5000, K=   45068.276*logX, nonDegen log|Lambda|>-3.626730 e6*logX, nonDegenNUB=7.253459 e6, degenNUB1=    0.e-19, degenNUB2=7.277166 e6, degenNUB3=1.519463 e7, nUB=7.277166 e6, transB-b1
\\ third iteration values
\\ 11 May 2022
eg4_check_it3()={
	my(bigL,chi,m,nUB,rho);
	
	bigL=50;
	m=11.0;
	rho=5.0;
	chi=4.5;
	nUB=7.98*10^6;

	actMinNUB=eg4_search_general(bigL,bigL,m,m,rho,rho,chi,chi,nUB,1);
	expMinNUB=7.277166*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg4_check_it3(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  50, m= 11.0000, rho=  5.0000, chi= 4.5000, K=   45068.276*logX, nonDegen log|Lambda|>-3.626730 e6*logX, nonDegenNUB=7.253459 e6, degenNUB1=    0.e-19, degenNUB2=7.277166 e6, degenNUB3=1.519463 e7, nUB=7.277166 e6, transB-b1
\\ fourth iteration values
\\ 11 May 2022
eg4_check_it4()={
	my(bigL,chi,m,nUB,rho);
	
	bigL=50;
	m=11.0;
	rho=5.0;
	chi=4.5;
	nUB=7.28*10^6;

	actMinNUB=eg4_search_general(bigL,bigL,m,m,rho,rho,chi,chi,nUB,1);
	expMinNUB=7.277166*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg4_check_it4(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ for iteration 1
\\ 11 May 2022
eg4_search_it1(dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rhoLB,rhoUB);

	bigLLB=30;
	bigLUB=350;
	mLB=1;
	mUB=21;
	rhoLB=2;
	rhoUB=22;
	chiLB=1;
	chiUB=11;

	eg4_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,,dbg);
}

\\ 11 May 2022
eg4_search_it2(nUBInit,dbg=0)={ \\use nUBInit=2.03*10^7
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rhoLB,rhoUB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);

	bigLLB=30;
	bigLUB=350;
	mLB=1;
	mUB=21;
	rhoLB=2;
	rhoUB=22;
	chiLB=1;
	chiUB=11;

	eg4_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,nUBInit,dbg);
}

\\ to be used for iteration 3
\\ 11 May 2022
eg4_search_it3(nUBInit,dbg=0)={ \\use nUBInit=7.98*10^6
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rhoLB,rhoUB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);

	bigLLB=30;
	bigLUB=350;
	mLB=1;
	mUB=21;
	rhoLB=2;
	rhoUB=22;
	chiLB=1;
	chiUB=11;

	eg4_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,nUBInit,dbg);
}

\\ to be used for iteration 4
\\ 11 May 2022
eg4_search_it4(nUBInit,dbg=0)={ \\use nUBInit=7.28*10^6
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rhoLB,rhoUB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);

	bigLLB=30;
	bigLUB=350;
	mLB=1;
	mUB=21;
	rhoLB=2;
	rhoUB=22;
	chiLB=1;
	chiUB=11;

	eg4_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,nUBInit,dbg);
}

\\ 3 July 2022
eg4_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,nUBInit=0,dbg=0)={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigD,bigK,chiStep,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,matveevChi,minNUB,mStep,nDegenUB,nLB,nNonDegenUB,nUB,rhoStep,val,w);

	\\ bigD=[Q(al_1,al_2,al_3):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2,al_3):R] -- used for Matveev's bounds
	matveevChi=2;
	\\ remember that d=[Q(al_1,al_2,al_3):Q]/[R(al_1,al_2,al_3):R]
	d=bigD/matveevChi; \\ d is the "degree" value used in the kit

	al1=x; \\ this is just a placeholder so that al1=alpha_1 is considered to be a polynomial
	hgtA1=logX/2; \\ this must be correct though and consistent with logX usages that follow
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
	b3=n; \\ since b3=q with |q|<p
	nLB=50*10^6;
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
		printf("eg4_search_general(): calculated nUBInit=%9.6e\n",nUBInit);
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
			a1=rho*Pi/2+logX;
			a2=0.723*rho;
			a3=rho*Pi;
			forstep(chi=chiLB,chiUB,chiStep,
				val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
				minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
			);
		);
	);
	);
	return(minNUB);
}
