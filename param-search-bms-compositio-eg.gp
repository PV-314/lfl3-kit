\\ \r lfl3\param-search-bms-compositio-eg.gp

read("lfl3\\lfl-utils-alpha1Variable.gp");

\\ this is the Case I linear form in the 2006 Compositio paper of Yann, Maurice and Samir
\\ but we have switched alpha_1 and alpha_2

\\ L= 104, m=  7.0000, rho=  5.0000, chi=3.3000, K=   63724.743*logX, nonDegen log|Lambda|>-1.066635 e7*logX, nonDegenNUB=2.133269 e7, degenNUB1=   0.e-19, degenNUB2=2.106323 e7, degenNUB3=4.158053 e7, nUB=2.133269 e7, transB-b1
\\ this was before improving the get_logA() function
\\ 11 May 2022
eg4_check_it1_old1()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigD,bigL,chi,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,matveevChi,nLB,nUB,rho);
	
	bigL=104;
	m=7.0;
	rho=5.0;
	chi=3.3;
	logXLB=log(1.99*10^6);
	nLB=2*10^6;
	nUB=2.9*10^13;

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=2;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi; \\ d is the "degree" value used in the kit

	al1=x; \\ this is just a placeholder so that al1=alpha_1 is considered to be a polynomial
	hgtA1=logX/2; \\ this must be correct though and consistent with logX usages that follow
	absLogA1=Pi/2;
	a1=rho*Pi/2+logX;

	al2=(1+sqrt(-7))/(1-sqrt(-7));
	hgtA2=log(2)/2;
	absLogA2a=abs(log(al2));
	absLogA2b=abs(log(-al2));
	absLogA2=min(absLogA2a,absLogA2b);
	a2=0.73*rho;
	
	al3=-1;
	hgtA3=0;
	absLogA3=Pi;
	a3=rho*Pi;

	b1=2;
	b2=n;
	b3=n; \\ since b3=q with |q|<p
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX/2;
	lamUB0=log(2.2*sqrt(7));

	nUB=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,1);
	printf("eg4_check_it1_old1(): calculated nUB=%9.6e\n",nUB);
	nUB=2.9*10^13;
	printf("                 but reverting to nUB=%9.6e\n\n",nUB);

	val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\ L= 101, m=  7.0000, rho=  5.0000, chi=2.8000, K=   61886.529*logX, nonDegen log|Lambda|>-1.005986 e7*logX, nonDegenNUB=2.011971 e7, degenNUB1=   0.e-19, degenNUB2=2.039189 e7, degenNUB3=4.125844 e7, nUB=2.039189 e7, transB-b1
\\ this was after improving the get_logA() function
\\ 17 May 2022
eg4_check_it1()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigD,bigL,chi,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,matveevChi,nLB,nUB,rho);
	
	bigL=101;
	m=7.0;
	rho=5.0;
	chi=2.8;
	logXLB=log(1.99*10^6);
	nLB=2*10^6;
	nUB=2.9*10^13;

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=2;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi; \\ d is the "degree" value used in the kit

	al1=x; \\ this is just a placeholder so that al1=alpha_1 is considered to be a polynomial
	hgtA1=logX/2; \\ this must be correct though and consistent with logX usages that follow
	absLogA1=Pi/2;
	a1=rho*Pi/2+logX;

	al2=(1+sqrt(-7))/(1-sqrt(-7));
	hgtA2=log(2)/2;
	absLogA2a=abs(log(al2));
	absLogA2b=abs(log(-al2));
	absLogA2=min(absLogA2a,absLogA2b);
	a2=0.73*rho;
	
	al3=-1;
	hgtA3=0;
	absLogA3=Pi;
	a3=rho*Pi;

	b1=2;
	b2=n;
	b3=n; \\ since b3=q with |q|<p
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX/2;
	lamUB0=log(2.2*sqrt(7));

	nUB=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,1);
	printf("eg4_check_it1(): calculated nUB=%9.6e\n",nUB);
	nUB=2.9*10^13;
	printf("                 but reverting to nUB=%9.6e\n\n",nUB);

	val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\ second iteration values
\\ 11 May 2022
eg4_check_it2()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigL,chi,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,nLB,nUB,rho);
	
	bigL=60;
	m=8.1;
	rho=5.0;
	chi=4.3;
	logXLB=log(1.99*10^6);
	nLB=2*10^6;
	nUB=20.4*10^6;

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=2;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi; \\ d is the "degree" value used in the kit

	al1=x; \\ this is just a placeholder so that al1=alpha_1 is considered to be a polynomial
	hgtA1=logX/2; \\ this must be correct though and consistent with logX usages that follow
	absLogA1=Pi/2;
	a1=rho*Pi/2+logX;

	al2=(1+sqrt(-7))/(1-sqrt(-7));
	hgtA2=log(2)/2;
	absLogA2a=abs(log(al2));
	absLogA2b=abs(log(-al2));
	absLogA2=min(absLogA2a,absLogA2b);
	a2=0.723*rho;
	
	al3=-1;
	hgtA3=0;
	absLogA3=Pi;
	a3=rho*Pi;

	b1=n;
	b2=2;
	b3=n; \\ since b3=q with |q|<p
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX/2;
	lamUB0=log(2.2*sqrt(7));

	val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\ L=  58, m= 10.3000, rho=  4.6000, chi=4.5000, K=   43016.769*logX, nonDegen log|Lambda|>-3.807469 e6*logX, nonDegenNUB=7.614937 e6, degenNUB1=   0.e-19, degenNUB2=7.421633 e6, degenNUB3=1.498226 e7, nUB=7.614937 e6, transB-b1
\\ third iteration values
\\ 11 May 2022
eg4_check_it3()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigL,chi,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,nLB,nUB,rho);
	
	bigL=58;
	m=10.3;
	rho=4.6;
	chi=4.5;
	logXLB=log(1.99*10^6);
	nLB=2*10^6;
	nUB=8.22*10^6;

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=2;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi; \\ d is the "degree" value used in the kit

	al1=x; \\ this is just a placeholder so that al1=alpha_1 is considered to be a polynomial
	hgtA1=logX/2; \\ this must be correct though and consistent with logX usages that follow
	absLogA1=Pi/2;
	a1=rho*Pi/2+logX;

	al2=(1+sqrt(-7))/(1-sqrt(-7));
	hgtA2=log(2)/2;
	absLogA2a=abs(log(al2));
	absLogA2b=abs(log(-al2));
	absLogA2=min(absLogA2a,absLogA2b);
	a2=0.723*rho;
	
	al3=-1;
	hgtA3=0;
	absLogA3=Pi;
	a3=rho*Pi;

	b1=n;
	b2=2;
	b3=n; \\ since b3=q with |q|<p
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX/2;
	lamUB0=log(2.2*sqrt(7));

	val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\ for iteration 1
\\ 11 May 2022
eg4_search_it1(dbg=0)={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigD,bigK,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,m,matveevChi,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,rho,val,w);

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=2;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
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
	logXLB=log(1.99*10^6);
	nLB=2*10^6;
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX/2;
	lamUB0=log(2.2*sqrt(7));
	nUBInit=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,1);
	printf("nUBInit=%9.6e\n",nUBInit);
	printf("eg4_search_it1(): calculated nUBInit=%9.6e\n",nUBInit);
	nUBInit=2.9*10^13;
	printf("                 but reverting to nUBInit=%9.6e\n\n",nUBInit);
	minNUB=nUBInit;
	
	for(bigL=40,350, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%10==0,print("L=",bigL));
	for(mv=1,25, \\ m \geq 1
		m=mv/1.0;
		for(pv=2,20, \\ rho \geq 2
		rho=pv/1.0;
		a1=rho*Pi/2+logX;
		a2=0.723*rho;
		a3=rho*Pi;
		for(cv=1,100,
			chi=cv/10.0;
			val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
			minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
		);
		);
	);
	);
}

\\ difference with eg4_search_it1() is that we step by 0.1 for m
\\ 11 May 2022
eg4_search_it2(nUBInit,dbg=0)={ \\use 8.576501 e6 (3 Mar 2022)
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,rho,val,w);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	minNUB=nUBInit;

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=2;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi; \\ d is the "degree" value used in the kit

	al1=x;
	hgtA1=logX/2;
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
	logXLB=log(1.99*10^6);
	nLB=2*10^6;
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX/2;
	lamUB0=log(2.2*sqrt(7));

	for(bigL=30,350, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%10==0,print("L=",bigL));
	for(mv=10,100, \\ m \geq 1
		m=mv/10.0;
		for(pv=2,15, \\ rho \geq 2
		rho=pv/1.0;
		a1=rho*Pi/2+logX;
		a2=0.723*rho;
		a3=rho*Pi;
		for(cv=1,100,
			chi=cv/10.0;
			val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
			minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
		);
		);
	);
	);
}

\\ to be used for iteration 3
\\ difference with eg4_search_it2() is that we step by 0.1 for rho
\\ 11 May 2022
eg4_search_it3(nUBInit,dbg=0)={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,hgtA1,hgtA2,hgtA3,logXLB,m,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,rho,val,w);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	minNUB=nUBInit;

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=2;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi; \\ d is the "degree" value used in the kit

	al1=x;
	hgtA1=logX/2;
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
	logXLB=log(1.99*10^6);
	nLB=2*10^6;
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX/2;
	lamUB0=log(2.2*sqrt(7));

	for(bigL=40,150, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%1==0,print("L=",bigL));
	for(mv=50,130, \\ m \geq 1
		m=mv/10.0;
		for(pv=30,80,
		rho=pv/10.0;
		a1=rho*Pi/2+logX;
		a2=0.723*rho;
		a3=rho*Pi;
		\\print("in eg4_search_it3(): a1=",a1,", a3=",a3);
		for(cv=30,80,
			chi=cv/10.0;
			val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
			minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
		);
		);
	);
	);
}

\\ to be used for iteration 4
\\ difference with eg4_search_it3() is that we step by 0.01 for m, rho, chi
\\ 11 May 2022
eg4_search_it4(nUBInit,dbg=0)={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,hgtA1,hgtA2,hgtA3,logXLB,m,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,rho,val,w);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	minNUB=nUBInit;

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=2;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=2;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi; \\ d is the "degree" value used in the kit

	al1=x;
	hgtA1=logX/2;
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
	logXLB=log(1.99*10^6);
	nLB=2*10^6;
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX/2;
	lamUB0=log(2.2*sqrt(7));

	for(bigL=80,150, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%1==0,print("L=",bigL));
	for(mv=1400,1600, \\ m \geq 1
		m=mv/100.0;
		\\for(pv=20,300,
		for(pv=700,900,
		rho=pv/100.0;
		a1=rho*Pi/2+logX;
		a2=0.723*rho;
		a3=rho*Pi;
		for(cv=140,155,
			chi=cv/100.0;
			val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
			minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
		);
		);
	);
	);
}
