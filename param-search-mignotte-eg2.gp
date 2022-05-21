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
	my(al1,al2,al3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,degenNUB,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,minNUB,nLB,nonDegenNUB,nUB,rho,val);

	d=1;
	al1=x; \\ actually x/y (or y/x)
	al2=5;
	al3=2;
	b1=n; \\ 0 \leq alpha < n
	b2=n;
	b3=n; \\ 0 \leq beta < n
	logXLB=log(7);
	nLB=10*10^6; \\ assumption at the start of proof of Theorem 6.3
	
	hgtA1=logX;
	absLogA1=logX;
	hgtA2=log(al2);
	absLogA2=abs(log(al2));
	hgtA3=log(al3);
	absLogA3=abs(log(al3));
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX;
	lamUB0=0;
	
	bigL=110;
	m=41.28955;
	rho=7.0;
	chi=1.0;
	nUB=5.7*10^11;

	a1=(rho-1)*log(5.0001)+2*logX;
	a2=(rho+1)*log(al2);
	a3=(rho+1)*log(al3);

	val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
	if(length(val)==0,
		print("invalid param choice. Step (3) fails");
	);
}

\\ L=  93, m=  6.0000, rho=  5.0000, chi=1.5000, K=  118959.545*logX, nonDegen log|Lambda|>-1.780559 e7*logX, nonDegenNUB=1.780559 e7, degenNUB1=4.868621 e7, degenNUB3=1.747417 e7, nUB=1.780559 e7, transB-b3
\\ 27 Jan 2022
eg2_check_it1_old1()={
	my(al1,al2,al3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,degenNUB,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,minNUB,nLB,nonDegenNUB,nUB,rho,val);

	d=1;
	al1=x; \\ actually x/y (or y/x)
	al2=5;
	al3=2;
	b1=n;
	b2=n; \\ 0 \leq beta < n
	b3=n; \\ 0 \leq alpha < n
	logXLB=log(7);
	nLB=10*10^6; \\ assumption at the start of proof of Theorem 6.3
	
	hgtA1=logX;
	absLogA1=logX;
	hgtA2=log(al2);
	absLogA2=abs(log(al2));
	hgtA3=log(al3);
	absLogA3=abs(log(al3));
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX;
	lamUB0=0;
	
	bigL=93;
	m=6.0;
	rho=5.0;
	chi=1.5;
	nUB=5.7*10^11;

	a1=(rho-1)*log(5.0001)+2*logX;
	a2=(rho+1)*log(al2);
	a3=(rho+1)*log(al3);

	val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
	if(length(val)==0,
		print("invalid param choice. Step (3) fails");
	);
}

\\L=  82, m=  7.0000, rho=  5.0000, chi=1.6000, K=  122370.571*logX, nonDegen log|Lambda|>-1.614972 e7*logX, nonDegenNUB=1.614972 e7, degenNUB1=4.513984 e7, degenNUB3=1.614451 e7, nUB=1.614972 e7, transB-b3
\\ 28 Jan 2022 (after correcting expression for a going into de Weger-Petho lemma)
eg2_check_it1_old2()={
	my(al1,al2,al3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,degenNUB,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,minNUB,nLB,nonDegenNUB,nUB,rho,val);

	d=1;
	al1=x; \\ actually x/y (or y/x)
	al2=5;
	al3=2;
	b1=n;
	b2=n; \\ 0 \leq beta < n
	b3=n; \\ 0 \leq alpha < n
	logXLB=log(7);
	nLB=10*10^6; \\ assumption at the start of proof of Theorem 6.3
	
	hgtA1=logX;
	absLogA1=logX;
	hgtA2=log(al2);
	absLogA2=abs(log(al2));
	hgtA3=log(al3);
	absLogA3=abs(log(al3));
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX;
	lamUB0=0;
	
	bigL=82;
	m=7.0;
	rho=5.0;
	chi=1.6;
	nUB=5.7*10^11;

	a1=(rho-1)*log(5.0001)+2*logX;
	a2=(rho+1)*log(al2);
	a3=(rho+1)*log(al3);

	val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
	if(length(val)==0,
		print("invalid param choice. Step (3) fails");
	);
}

\\L=  90, m=  6.0000, rho=  5.0000, chi=1.8000, K=  115122.140*logX, nonDegen log|Lambda|>-1.667537 e7*logX, nonDegenNUB=1.667537 e7, degenNUB1=3.636226 e7, degenNUB3=1.621235 e7, nUB=1.667537 e7, transB-b3
\\ 28 Jan 2022 (after correcting expression for a going into de Weger-Petho lemma)
eg2_check_it1_old3()={
	my(al1,al2,al3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,degenNUB,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,minNUB,nLB,nonDegenNUB,nUB,rho,val);

	d=1;
	al1=x; \\ actually x/y (or y/x)
	al2=5;
	al3=2;
	b1=n;
	b2=n; \\ 0 \leq beta < n
	b3=n; \\ 0 \leq alpha < n
	logXLB=log(7);
	nLB=10*10^6; \\ assumption at the start of proof of Theorem 6.3
	
	hgtA1=logX;
	absLogA1=logX;
	hgtA2=log(al2);
	absLogA2=abs(log(al2));
	hgtA3=log(al3);
	absLogA3=abs(log(al3));
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX;
	lamUB0=0;
	
	bigL=90;
	m=6.0;
	rho=5.0;
	chi=1.8;
	nUB=5.7*10^11;

	a1=(rho-1)*log(5.0001)+2*logX;
	a2=(rho+1)*log(al2);
	a3=(rho+1)*log(al3);

	val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
	if(length(val)==0,
		print("invalid param choice. Step (3) fails");
	);
}

\\ L= 155, m= 33.0000, rho=  5.0000, chi=0.1000, K= 1090462.493*logZ, minBnd=-2.720299 e8*logZ, nNonDegenUB=2.720299 e8, nDegenUB=2.720299 e8, nUB=2.720299 e8
\\ new 26 Nov 2021: L= 162, m= 39.0000, rho=  5.0000, chi=1.2000, K= 1346929.039*logZ, minBnd=-3.511834 e8*logZ, nNonDegenUB=3.511834 e8, nDegenUB=3.505548 e8, nUB=3.511834 e8
\\ 25 Nov 2021
eg2_check_it1()={
	my(al1,al2,al3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,degenNUB,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,minNUB,nLB,nonDegenNUB,nUB,rho,val);

	d=1;
	al1=x; \\ actually x/y (or y/x)
	al2=5;
	al3=2;
	b1=n; \\ 0 \leq alpha < n
	b2=n;
	b3=n; \\ 0 \leq beta < n
	logXLB=log(7);
	nLB=10*10^6; \\ assumption at the start of proof of Theorem 6.3
	
	hgtA1=logX;
	absLogA1=logX;
	hgtA2=log(al2);
	absLogA2=abs(log(al2));
	hgtA3=log(al3);
	absLogA3=abs(log(al3));
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX;
	lamUB0=0;
	
	bigL=162;
	m=39.0;
	rho=5.0;
	chi=1.2;
	nUB=5.7*10^11;

	a1=(rho-1)*log(5.0001)+2*logX;
	a2=(rho+1)*log(al2);
	a3=(rho+1)*log(al3);

	val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
	if(length(val)==0,
		print("invalid param choice. Step (3) fails");
	);
}

\\ L=  91, m=  6.0000, rho=  5.0000, chi=1.8000, K=  116401.275*logX, nonDegen log|Lambda|>-1.704800 e7*logX, nonDegenNUB=1.704800 e7, degenNUB1=   0.e-19, degenNUB2=3.719930 e7, degenNUB3=1.674490 e7, nUB=1.704800 e7, transB-b1
\\ 25 Nov 2021
eg2_search_it1(dbg=0)={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigD,bigK,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,m,matveevChi,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,rho,val,w);

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=1;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
	d=bigD/matveevChi; \\ d is the "degree" value used in the kit

	al1=x; \\ actually x/y (or y/x)
	al2=5;
	al3=2;
	b1=n;
	b2=n; \\ 0 \leq beta < n
	b3=n; \\ 0 \leq alpha < n
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

	nUBInit=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,1);
	printf("nUBInit=%9.6e\n",nUBInit);
	printf("eg2_search_it1(): calculated nUBInit=%9.6e\n",nUBInit);
	nUBInit=5.4*10^11;
	printf("                 but reverting to nUBInit=%9.6e\n\n",nUBInit);
	minNUB=nUBInit;

	for(bigL=70,200, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%10==0,print("L=",bigL));
	for(mv=2,30, \\ m \geq 1
		m=mv/1.0;
		for(pv=2,10,
			rho=pv/1.0;
			a1=(rho-1)*log(5.0001)+2*logX;
			a2=(rho+1)*log(al2);
			a3=(rho+1)*log(al3);
			for(cv=1,30,
				chi=cv/10.0;
				val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
				minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
			);
		);
	);
	);
}

\\ difference with eg2_search_it1() is that we step by 0.1 for m
\\ 4 Jan 2022
eg2_search_it2(nUBInit,dbg=0)={
	my(al1,al2,al3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,degenNUB,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,minNUB,nLB,nonDegenNUB,nUB,rho,val);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	minNUB=nUBInit;

	d=1;
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

	for(bigL=30,100, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%1==0,print("L=",bigL));
	for(mv=10,50, \\ m \geq 1
		m=mv/10.0;
		for(pv=2,30, \\ rho \geq 2
			rho=pv/1.0;
			a1=(rho-1)*log(5.0001)+2*logX;
			a2=(rho+1)*log(al2);
			a3=(rho+1)*log(al3);
			for(cv=1,50,
				chi=cv/10.0;
				val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
					minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
			);
		);
	);
	);
}

\\ to be used for iteration 3
\\ difference with eg2_search_it2() is that we step by 0.1 for rho
\\ 5 Jan 2022
eg2_search_it3(nUBInit,dbg=0)={
	my(al1,al2,al3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,degenNUB,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,minNUB,nLB,nonDegenNUB,nUB,rho,val);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	minNUB=nUBInit;

	d=1;
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

	for(bigL=30,100, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%1==0,print("L=",bigL));
	for(mv=10,50, \\ m \geq 1
		m=mv/10.0;
		for(pv=20,150, \\ rho \geq 2
			rho=pv/10.0;
			a1=(rho-1)*log(5.0001)+2*logX;
			a2=(rho+1)*log(al2);
			a3=(rho+1)*log(al3);
			for(cv=1,50,
				chi=cv/10.0;
				val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
					minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
			);
		);
	);
	);
}

\\ to be used for iteration 4
\\ difference with eg2_search_it3() is that we step by 0.01 for m, rho, chi
\\ 5 Jan 2022
eg2_search_it4(nUBInit,dbg=0)={
	my(al1,al2,al3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,degenNUB,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,minNUB,nLB,nonDegenNUB,nUB,rho,val);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	minNUB=nUBInit;

	d=1;
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

	for(bigL=30,100, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%1==0,print("L=",bigL));
	for(mv=10,50, \\ m \geq 1
		m=mv/10.0;
		for(pv=20,150, \\ rho \geq 2
			rho=pv/10.0;
			a1=(rho-1)*log(5.0001)+2*logX;
			a2=(rho+1)*log(al2);
			a3=(rho+1)*log(al3);
			for(cv=1,50,
				chi=cv/10.0;
				val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
				minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
			);
		);
	);
	);
}
