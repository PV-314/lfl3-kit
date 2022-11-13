\\ \r lfl3\param-search-bms-annals-eg.gp

read("lfl3\\lfl-utils-alpha3Variable.gp");

\\ this is the linear form in the Annals paper of Yann, Maurice and Samir

\\ L= 147, m=  2.0000, rho=  9.0000, chi= 3.0000, K=   54647.348*logX, nonDegen log|Lambda|>-1.765066 e7*logX, nonDegenNUB=8.825328 e6, degenNUB1=8.994871 e6, degenNUB2=2.372594 e7, degenNUB3=    0.e-19, nUB=8.994871 e6, transB-b1
\\ 10 May 2022 (after improving constant to use from Laurent's 2008 paper)
eg3_check_it1()={
	my(bigL,chi,m,nUB,rho);
	
	bigL=147;
	m=2.0;
	rho=9.0;
	chi=3.0;

	actMinNUB=eg3_search_general(bigL,bigL,m,m,rho,rho,chi,chi,,1);
	expMinNUB=8.994871*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg3_check_it1(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  94, m=  2.5000, rho=  8.5000, chi= 4.2500, K=   39767.640*logX, nonDegen log|Lambda|>-7.999906 e6*logX, nonDegenNUB=3.999953 e6, degenNUB1=3.891046 e6, degenNUB2=1.034027 e7, degenNUB3=    0.e-19, nUB=3.999953 e6, transB-b1
\\ second iteration values
\\ 11 May 2022
eg3_check_it2()={
	my(bigL,chi,m,nUB,rho);
	
	bigL=94;
	m=2.5;
	rho=8.5;
	chi=4.25;
	nUB=9*10^6;

	actMinNUB=eg3_search_general(bigL,bigL,m,m,rho,rho,chi,chi,nUB,1);
	expMinNUB=3.999953*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg3_check_it2(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ L=  81, m=  2.5000, rho=  9.5000, chi= 4.5000, K=   41168.482*logX, nonDegen log|Lambda|>-7.507264 e6*logX, nonDegenNUB=3.753632 e6, degenNUB1=3.788629 e6, degenNUB2=1.006819 e7, degenNUB3=    0.e-19, nUB=3.788629 e6, transB-b1
\\ third iteration values
\\ 13 Dec 2021
eg3_check_it3()={
	my(bigL,chi,m,nUB,rho);
	
	bigL=81;
	m=2.5;
	rho=9.5;
	chi=4.5;
	nUB=4*10^6;

	actMinNUB=eg3_search_general(bigL,bigL,m,m,rho,rho,chi,chi,nUB,1);
	expMinNUB=3.788629*10^6;
	if(abs(actMinNUB/expMinNUB-1)>0.0001,
		printf("FAIL: eg3_check_it3(), actMinNUB=%9.6e, expMinNUB=%9.6e\n",actMinNUB,expMinNUB);
	);
}

\\ for iteration 1
\\ 11 Dec 2021 (changed to use "general" function on 3 July 2022)
eg3_search_it1(nUBInit=0,dbg=0)={
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rhoLB,rhoUB);

	bigLLB=60;
	bigLUB=350;
	mLB=1;
	mUB=11;
	rhoLB=2;
	rhoUB=22;
	chiLB=1;
	chiUB=11;

	eg3_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,nUBInit,dbg);
}

\\ for iteration 2
\\ 13 Dec 2021
eg3_search_it2(nUBInit,dbg=0)={ \\use nUBInit=8.37*10^6
	my(bigLLB,bigLUB,chiLB,chiUB,mLB,mUB,rhoLB,rhoUB);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	bigLLB=40;
	bigLUB=500;
	mLB=1;
	mUB=6;
	rhoLB=2;
	rhoUB=12;
	chiLB=1;
	chiUB=6;
	eg3_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,nUBInit,dbg);
}

\\ for iteration 3
\\ 13 Dec 2021
eg3_search_it3(nUBInit,dbg=0)={ \\ use nUBInit=3.87*10^6
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,hgtA1,hgtA2,hgtA3,logXLB,m,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,rho,val,w);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);

	bigLLB=40;
	bigLUB=500;
	mLB=1;
	mUB=6;
	rhoLB=2;
	rhoUB=12;
	chiLB=1;
	chiUB=6;
	eg3_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,nUBInit,dbg);
}

\\ 3 July 2022
eg3_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,nUBInit=0,dbg=0)={ \\ should use nUBInit=8.73*10^11
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigD,bigK,chi,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,m,matveevChi,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,rho,val,w);

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
	b1=n;
	b2=1;
	b3=n;
	logXLB=10^20;
	nLB=2*10^6;
	
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
		nUBInit=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,1);
		printf("nUBInit=%9.6e\n",nUBInit);
		printf("eg3_search_general(): calculated nUBInit=%9.6e\n",nUBInit);
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

	for(bigL=bigLLB,bigLUB, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%10==0,print("L=",bigL));
	forstep(m=mLB,mUB,mStep,
		forstep(rho=rhoLB,rhoUB,rhoStep,
			a1=(rho+1)*log(al1);
			a2=(rho+3)*log(al2);
			a3=(rho+1)*logW+4*logX+(rho+3)*10^(-6); \\ here logX=log(y)
			forstep(chi=chiLB,chiUB,chiStep,
				val=alpha3_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
				minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
			);
		);
	);
	);
	return(minNUB);
}
