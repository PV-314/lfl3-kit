\\ \r lfl3\param-search-eg3.gp

read("lfl3\\lfl-utils-alpha3Variable.gp");

\\ this is the linear form in the Annals paper of Yann, Maurice and Samir

\\ 13 Dec 2021
eg3_check_it1_old1()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigD,bigL,chi,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,matveevChi,nLB,nUB,rho,w);
	
	bigL=128;
	m=16.0;
	rho=8.0;
	chi=1.1;
	logXLB=10^20;
	nLB=2*10^6;
	nUB=4.7*10^11;

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi;
	w=(1+sqrt(5))/2;
	al1=w;
	al2=sqrt(5);
	al3=x; \\ x is just because it is not known. Here it is al_3=omega^k/y
	b1=n;
	b2=1;
	b3=n;
	a1=(rho+1)*log(al1);
	a2=(rho+3)*log(al2);
	a3=(rho+1)*log(w)+4*logX+(rho+3)*10^(-6); \\ here logX=log(y)
	
	hgtA1=log(w)/2;
	absLogA1=abs(log(al1));
	hgtA2=log(sqrt(5));
	absLogA2=abs(log(al2));
	hgtA3=log(w)/2+logX+10^(-6); \\ from equation (6.5)
	absLogA3=log(w)+10^(-6); \\ from equation (6.4)
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-2*logX;
	lamUB0=1;

	\\nUB=get_matveev_ubnd(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,lamUB1,lamUB0,1);
	nUB=4.7*10^11;
	nUB=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,1);
	printf("eg3_check_it1_old1(): nUB=%9.6e\n",nUB);

	val=alpha3_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\ 3 Jan 2022
eg3_check_it1_old2()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigL,chi,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,nLB,nUB,rho);
	
	bigL=144;
	m=2.0;
	rho=9.0;
	chi=3.2;
	logXLB=10^20;
	nLB=2*10^6;
	nUB=1.9*10^12;

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi;
	w=(1+sqrt(5))/2;
	al1=w;
	al2=sqrt(5);
	al3=x; \\ x is just because it is not known. Here it is al_3=omega^k/y
	b1=n;
	b2=1;
	b3=n;
	a1=(rho+1)*log(al1);
	a2=(rho+3)*log(al2);
	a3=(rho+1)*log(w)+4*logX+(rho+3)*10^(-6); \\ here logX=log(y)
	
	hgtA1=log(w)/2;
	absLogA1=abs(log(al1));
	hgtA2=log(sqrt(5));
	absLogA2=abs(log(al2));
	hgtA3=log(w)/2+logX+10^(-6); \\ from equation (6.5)
	absLogA3=log(w)+10^(-6); \\ from equation (6.4)
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-2*logX;
	lamUB0=1;

	val=alpha3_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\ 3 Jan 2022 (after degen bug fix)
eg3_check_it1_old3()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigL,chi,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,nLB,nUB,rho);
	
	bigL=155;
	m=3.0;
	rho=8.0;
	chi=6.2;
	logXLB=10^20;
	nLB=2*10^6;
	nUB=1.9*10^12;

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi;
	w=(1+sqrt(5))/2;
	al1=w;
	al2=sqrt(5);
	al3=x; \\ x is just because it is not known. Here it is al_3=omega^k/y
	b1=n;
	b2=1;
	b3=n;
	a1=(rho+1)*log(al1);
	a2=(rho+3)*log(al2);
	a3=(rho+1)*log(w)+4*logX+(rho+3)*10^(-6); \\ here logX=log(y)
	
	hgtA1=log(w)/2;
	absLogA1=abs(log(al1));
	hgtA2=log(sqrt(5));
	absLogA2=abs(log(al2));
	hgtA3=log(w)/2+logX+10^(-6); \\ from equation (6.5)
	absLogA3=log(w)+10^(-6); \\ from equation (6.4)
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-2*logX;
	lamUB0=1;

	val=alpha3_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\ L= 164, m=  2.0000, rho=  8.0000, chi=3.5000, K=   50297.865*logX, nonDegen log|Lambda|>-1.715300 e7*logX, nonDegenNUB=8.576501 e6, degenNUB1=8.028405 e6, degenNUB3=2.234751 e7, nUB=8.576501 e6, transB-b1
\\ 3 Mar 2022 (after checking both degen possibilities)
eg3_check_it1_old4()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigD,bigL,chi,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,matveevChi,nLB,nUB,rho);
	
	bigL=164;
	m=2.0;
	rho=8.0;
	chi=3.5;
	logXLB=10^20;
	nLB=2*10^6;
	nUB=8.73*10^11;

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi;
	w=(1+sqrt(5))/2;
	al1=w;
	al2=sqrt(5);
	al3=x; \\ x is just because it is not known. Here it is al_3=omega^k/y
	b1=n;
	b2=1;
	b3=n;
	a1=(rho+1)*log(al1);
	a2=(rho+3)*log(al2);
	a3=(rho+1)*log(w)+4*logX+(rho+3)*10^(-6); \\ here logX=log(y)
	
	hgtA1=log(w)/2;
	absLogA1=abs(log(al1));
	hgtA2=log(sqrt(5));
	absLogA2=abs(log(al2));
	hgtA3=log(w)/2+logX+10^(-6); \\ from equation (6.5)
	absLogA3=log(w)+10^(-6); \\ from equation (6.4)
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-2*logX;
	lamUB0=1;

	nUB=1.9*10^12;
	nUB=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,1);
	printf("eg3_check_it1(): nUB=%9.6e\n\n",nUB);
	nUB=8.73*10^11;
	val=alpha3_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\L= 149, m=  2.0000, rho=  9.0000, chi=3.8000, K=   55390.849*logX, nonDegen log|Lambda|>-1.813421 e7*logX, nonDegenNUB=9.067107 e6, degenNUB1=8.941538 e6, degenNUB2=2.285469 e7, degenNUB3=   0.e-19, nUB=9.067107 e6, transB-b1
\\ 9 May 2022 (after fixing nUB)
eg3_check_it1_old5()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigD,bigL,chi,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,matveevChi,nLB,nUB,rho);
	
	bigL=149;
	m=2.0;
	rho=9.0;
	chi=3.8;
	logXLB=10^20;
	nLB=2*10^6;
	nUB=1.9*10^12;

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi;
	w=(1+sqrt(5))/2;
	al1=w;
	al2=sqrt(5);
	al3=x; \\ x is just because it is not known. Here it is al_3=omega^k/y
	b1=n;
	b2=1;
	b3=n;
	a1=(rho+1)*log(al1);
	a2=(rho+3)*log(al2);
	a3=(rho+1)*log(w)+4*logX+(rho+3)*10^(-6); \\ here logX=log(y)
	
	hgtA1=log(w)/2;
	absLogA1=abs(log(al1));
	hgtA2=log(sqrt(5));
	absLogA2=abs(log(al2));
	hgtA3=log(w)/2+logX+10^(-6); \\ from equation (6.5)
	absLogA3=log(w)+10^(-6); \\ from equation (6.4)
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-2*logX;
	lamUB0=1;
	val=alpha3_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\L= 143, m=  2.0000, rho=  9.0000, chi=3.1000, K=   53160.345*logX, nonDegen log|Lambda|>-1.670315 e7*logX, nonDegenNUB=8.351573 e6, degenNUB1=8.238661 e6, degenNUB2=2.162236 e7, degenNUB3=   0.e-19, nUB=8.351573 e6, transB-b1
\\ 10 May 2022 (after improving constant to use from Laurent's 2008 paper)
eg3_check_it1()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigD,bigL,chi,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,matveevChi,nLB,nUB,rho);
	
	bigL=143;
	m=2.0;
	rho=9.0;
	chi=3.1;
	logXLB=10^20;
	nLB=2*10^6;
	nUB=1.9*10^12;

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi;
	w=(1+sqrt(5))/2;
	al1=w;
	al2=sqrt(5);
	al3=x; \\ x is just because it is not known. Here it is alpha_3=omega^k/y
	b1=n;
	b2=1;
	b3=n;
	a1=(rho+1)*log(al1);
	a2=(rho+3)*log(al2);
	a3=(rho+1)*log(w)+4*logX+(rho+3)*10^(-6); \\ here logX=log(y)
	
	hgtA1=log(w)/2;
	absLogA1=abs(log(al1));
	hgtA2=log(sqrt(5));
	absLogA2=abs(log(al2));
	hgtA3=log(w)/2+logX+10^(-6); \\ from equation (6.5)
	absLogA3=log(w)+10^(-6); \\ from equation (6.4)
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-2*logX;
	lamUB0=1;

	nUB=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,1);
	printf("eg3_check_it1(): calculated nUB=%9.6e\n",nUB);
	nUB=1.9*10^12;
	printf("                 but reverting to nUB=%9.6e\n\n",nUB);

	val=alpha3_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\ second iteration values
\\ 13 Dec 2021
eg3_check_it2_old1()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigL,chi,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,nLB,nUB,rho);
	
	bigL=90;
	m=17.0;
	rho=8.0;
	chi=1.5;
	logXLB=10^20;
	nLB=2*10^6;
	nUB=99.77*10^6;
	
	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi;
	w=(1+sqrt(5))/2;
	al1=w;
	al2=sqrt(5);
	al3=x; \\ x is just because it is not known. Here it is al_3=omega^k/y
	b1=n;
	b2=1;
	b3=n;
	a1=(rho+1)*log(al1);
	a2=(rho+3)*log(al2);
	a3=(rho+1)*log(w)+4*logX+(rho+3)*10^(-6); \\ here logX=log(y)
	
	hgtA1=log(w)/2;
	absLogA1=abs(log(al1));
	hgtA2=log(sqrt(5));
	absLogA2=abs(log(al2));
	hgtA3=log(w)/2+logX+10^(-6); \\ from equation (6.5)
	absLogA3=log(w)+10^(-6); \\ from equation (6.4)
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-2*logX;
	lamUB0=1;

	val=alpha3_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\L=  90, m=  2.3000, rho=  9.0000, chi=4.4000, K=   38476.194*logX, nonDegen log|Lambda|>-7.608676 e6*logX, nonDegenNUB=3.804338 e6, degenNUB1=3.814004 e6, degenNUB2=1.013560 e7, degenNUB3=   0.e-19, nUB=3.814004 e6, transB-b1
\\ second iteration values
\\ 11 May 2022
eg3_check_it2()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigL,chi,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,nLB,nUB,rho);
	
	bigL=90;
	m=2.3;
	rho=9.0;
	chi=4.4;
	logXLB=10^20;
	nLB=2*10^6;
	nUB=8.36*10^6;
	
	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi;
	w=(1+sqrt(5))/2;
	al1=w;
	al2=sqrt(5);
	al3=x; \\ x is just because it is not known. Here it is al_3=omega^k/y
	b1=n;
	b2=1;
	b3=n;
	a1=(rho+1)*log(al1);
	a2=(rho+3)*log(al2);
	a3=(rho+1)*log(w)+4*logX+(rho+3)*10^(-6); \\ here logX=log(y)
	
	hgtA1=log(w)/2;
	absLogA1=abs(log(al1));
	hgtA2=log(sqrt(5));
	absLogA2=abs(log(al2));
	hgtA3=log(w)/2+logX+10^(-6); \\ from equation (6.5)
	absLogA3=log(w)+10^(-6); \\ from equation (6.4)
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-2*logX;
	lamUB0=1;

	val=alpha3_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\ L=  93, m=  2.8000, rho=  7.7000, chi=4.4000, K=   37547.802*logX, nonDegen log|Lambda|>-7.127830 e6*logX, nonDegenNUB=3.563915 e6, degenNUB1=3.567114 e6, degenNUB2=9.159273 e6, degenNUB3=   0.e-19, nUB=3.567114 e6, transB-b1
\\ third iteration values
\\ 13 Dec 2021
eg3_check_it3()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigL,chi,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,nLB,nUB,rho);
	
	bigL=93;
	m=2.8;
	rho=7.7;
	chi=4.4;
	logXLB=10^20;
	nLB=2*10^6;
	nUB=3.82*10^6;
	
	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi;
	w=(1+sqrt(5))/2;
	al1=w;
	al2=sqrt(5);
	al3=x; \\ x is just because it is not known. Here it is al_3=omega^k/y
	b1=n;
	b2=1;
	b3=n;
	a1=(rho+1)*log(al1);
	a2=(rho+3)*log(al2);
	a3=(rho+1)*log(w)+4*logX+(rho+3)*10^(-6); \\ here logX=log(y)
	
	hgtA1=log(w)/2;
	absLogA1=abs(log(al1));
	hgtA2=log(sqrt(5));
	absLogA2=abs(log(al2));
	hgtA3=log(w)/2+logX+10^(-6); \\ from equation (6.5)
	absLogA3=log(w)+10^(-6); \\ from equation (6.4)
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-2*logX;
	lamUB0=1;

	val=alpha3_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\ for iteration 1
\\ 11 Dec 2021
eg3_search_it1(nUBInit,dbg=0)={ \\ should use nUBInit=8.73*10^11
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigD,bigK,chi,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,m,matveevChi,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,rho,val,w);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	minNUB=nUBInit;

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
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
	nUBInit=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,1);
	printf("nUBInit=%9.6e\n",nUBInit);
	minNUB=nUBInit;
	
	for(bigL=50,350, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%10==0,print("L=",bigL));
	for(mv=1,25, \\ m \geq 1
		m=mv/1.0;
		for(pv=2,40, \\ rho \geq 2
		rho=pv/1.0;
		a1=(rho+1)*log(al1);
		a2=(rho+3)*log(al2);
		a3=(rho+1)*logW+4*logX+(rho+3)*10^(-6); \\ here logX=log(y)
		for(cv=1,100,
			chi=cv/10.0;
			val=alpha3_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
			minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
		);
		);
	);
	);
}

\\ difference with eg3_search_it1() is that we step by 0.1 for m
\\ 13 Dec 2021
eg3_search_it2(nUBInit,dbg=0)={ \\use 8.576501 e6 (3 Mar 2022)
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,rho,val,w);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	minNUB=nUBInit;

	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=2;
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

	for(bigL=50,350, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%10==0,print("L=",bigL));
	for(mv=10,50, \\ m \geq 1
		m=mv/10.0;
		for(pv=2,15, \\ rho \geq 2
		rho=pv/1.0;
		a1=(rho+1)*log(al1);
		a2=(rho+3)*log(al2);
		a3=(rho+1)*logW+4*logX+(rho+3)*10^(-6); \\ here logX=log(y)
		for(cv=1,100,
			chi=cv/10.0;
			val=alpha3_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
			minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
		);
		);
	);
	);
}

\\ to be used for iteration 3
\\ difference with eg3_search_it2() is that we step by 0.1 for rho
\\ 13 Dec 2021
eg3_search_it3(nUBInit,dbg=0)={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,hgtA1,hgtA2,hgtA3,logXLB,m,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,rho,val,w);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	minNUB=nUBInit;

	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=2;
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

	for(bigL=70,150, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%1==0,print("L=",bigL));
	for(mv=10,60, \\ m \geq 1
		m=mv/10.0;
		for(pv=50,120,
		rho=pv/10.0;
		a1=(rho+1)*log(al1);
		a2=(rho+3)*log(al2);
		a3=(rho+1)*logW+4*logX+(rho+3)*10^(-6); \\ here logX=log(y)
		\\print("in eg3_search_it3(): a1=",a1,", a3=",a3);
		for(cv=10,90,
			chi=cv/10.0;
			val=alpha3_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
			minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
		);
		);
	);
	);
}

\\ to be used for iteration 4
\\ difference with eg3_search_it3() is that we step by 0.01 for m, rho, chi
\\ 13 Dec 2021
eg3_search_it4(nUBInit,dbg=0)={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,hgtA1,hgtA2,hgtA3,logXLB,m,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,rho,val,w);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	minNUB=nUBInit;

	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=2;
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

	for(bigL=80,150, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%1==0,print("L=",bigL));
	for(mv=1400,1600, \\ m \geq 1
		m=mv/100.0;
		\\for(pv=20,300,
		for(pv=700,900,
		rho=pv/100.0;
		a1=(rho+1)*log(al1);
		a2=(rho+3)*log(al2);
		a3=(rho+1)*logW+4*logX+(rho+3)*10^(-6); \\ here logX=log(y)
		for(cv=140,155,
			chi=cv/100.0;
			val=alpha3_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
			minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
		);
		);
	);
	);
}