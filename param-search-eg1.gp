\\ \r lfl3\param-search-eg1.gp

read("lfl3\\lfl-utils-alpha1Variable.gp");

\\ search for parameter choice for Example 1 in Mignotte's original kit:
\\ for the diophantine equation x^n-2^\alpha*5^\beta*y^n = \pm 1
\\ where n \geq 3 is prime, \alpha=1,2,3, 0 \leq \beta<n
\\ \Lamnda = n \log (x/y)-\alpha \log(2)-\beta \log(5)

\\ L=  70, m=  3.0000, rho=  8.0000, chi=2.4100, K=   39699.763*logX, nonDegen log|Lambda|>-5.778734 e6*logX, nonDegenNUB=5.778734 e6, degenNUB1=5.921966 e6, degenNUB3=1.661151 e7, nUB=5.921966 e6, transB-b1
\\ L=  70, m=  3.0000, rho=  8.0000, chi=2.5000, K=   39699.763*logX, nonDegen log|Lambda|>-5.778734 e6*logX, nonDegenNUB=5.778734 e6, degenNUB1=5.826089 e6, degenNUB3=1.601925 e7, nUB=5.826089 e6, transB-b1
\\ 4 Jan 2022
eg1_check_it1()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigL,chi,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,nUB,rho);

	d=1;
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
	
	bigL=70;
	m=3.0;
	rho=8.0;
	chi=2.41;
	nUB=5.6*10^11;

	a1=(rho+1)*log(al1);
	a3=(rho+1)*log(al3);
	a2=(rho-1)*absLogA2+2*logX;

	val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\ L=  51, m=  5.0000, rho=  6.0000, chi=3.3000, K=   28795.345*logX, nonDegen log|Lambda|>-2.631311 e6*logX, nonDegenNUB=2.631311 e6, degenNUB1=2.493316 e6, degenNUB3=7.193410 e6, nUB=2.631311 e6, transB-b1
\\ 4 Jan 2022 (updated with above on 27 Feb 2022)
eg1_check_it2()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigL,chi,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,nUB,rho);

	d=1;
	al1=2;
	al2=x;
	al3=5;
	b1=3; \\ alpha=1, 2 or 3
	b2=n;
	b3=n; \\ 0 \leq beta<n
	nLB=10^6; \\ assumption at the start of proof of Theorem 6.3
	logXLB=floor(nLB/2600-1.5);
	if(dbg!=0,
		printf("calced logXLB=%5d\n",logXLB);
	);
	logXLB=380;

	hgtA1=log(2);
	absLogA1=abs(log(al1));
	hgtA2=logX;
	absLogA2=5.0001; \\ from Step 3 after defns of a_i's
	hgtA3=log(5);
	absLogA3=abs(log(al3));
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX;
	lamUB0=log(2);
	
	bigL=51;
	m=5.0;
	rho=6.0;
	chi=3.3;
	nUB=5.83*10^6;
	a1=(rho+1)*log(al1);
	a3=(rho+1)*log(al3);
	a2=(rho-1)*absLogA2+2*logX;
	val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\ 4 Jan 2022
eg1_check_it3()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigL,chi,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,nUB,rho);

	d=1;
	al1=2;
	al2=x;
	al3=5;
	b1=3; \\ alpha=1, 2 or 3
	b2=n;
	b3=n; \\ 0 \leq beta<n
	nLB=10^6; \\ assumption at the start of proof of Theorem 6.3
	logXLB=floor(nLB/2600-1.5);
	if(dbg!=0,
		printf("calced logXLB=%5d\n",logXLB);
	);
	logXLB=380;

	hgtA1=log(2);
	absLogA1=abs(log(al1));
	hgtA2=logX;
	absLogA2=5.0001; \\ from Step 3 after defns of a_i's
	hgtA3=log(5);
	absLogA3=abs(log(al3));
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX;
	lamUB0=log(2);
	
	bigL=70;
	m=3.0;
	rho=8.0;
	chi=2.41;
	nUB=5.6*10^11;
	a1=(rho+1)*log(al1);
	a3=(rho+1)*log(al3);
	a2=(rho-1)*absLogA2+2*logX;
	val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\ 4 Jan 2022
eg1_check_it4()={
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,al1,al2,al3,b1,b2,b3,bigL,chi,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,nUB,rho);

	d=1;
	al1=2;
	al2=x;
	al3=5;
	b1=3; \\ alpha=1, 2 or 3
	b2=n;
	b3=n; \\ 0 \leq beta<n
	nLB=10^6; \\ assumption at the start of proof of Theorem 6.3
	logXLB=floor(nLB/2600-1.5);
	if(dbg!=0,
		printf("calced logXLB=%5d\n",logXLB);
	);
	logXLB=380;

	hgtA1=log(2);
	absLogA1=abs(log(al1));
	hgtA2=logX;
	absLogA2=5.0001; \\ from Step 3 after defns of a_i's
	hgtA3=log(5);
	absLogA3=abs(log(al3));
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX;
	lamUB0=log(2);
	
	bigL=70;
	m=3.0;
	rho=8.0;
	chi=2.41;
	nUB=5.6*10^11;
	a1=(rho+1)*log(al1);
	a3=(rho+1)*log(al3);
	a2=(rho-1)*absLogA2+2*logX;
	val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUB,lamUB1,lamUB0,1);
}

\\ nUBInit should be 5.66*10^11 here
\\ 4 Jan 2022
eg1_search_it1(dbg=0)={
my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigD,bigK,d,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,m,matveevChi,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,rho,val,w);

	\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	bigD=1;
	\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	matveevChi=1;
	\\ remember that d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R]
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

	nUBInit=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,1);
	printf("nUBInit=%9.6e\n",nUBInit);
	printf("eg2_search_it1(): calculated nUBInit=%9.6e\n",nUBInit);
	nUBInit=5.66*10^11;
	printf("                 but reverting to nUBInit=%9.6e\n\n",nUBInit);
	minNUB=nUBInit;

	for(bigL=40,350, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%10==0,print("L=",bigL));
	for(mv=1,20, \\ m \geq 1
		m=mv/1.0;
		for(pv=2,20, \\ rho \geq 2
			rho=pv/1.0;
			a1=(rho-1)*absLogA1+2*logX;
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

\\ 4 Jan 2022 (27 Feb 2022: should use nUBInit=5.83*10^6)
eg1_search_it2(nUBInit,dbg=0)={
	my(al1,al2,al3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,degenNUB,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,minNUB,nLB,nonDegenNUB,nUB,rho,val);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	minNUB=nUBInit;
	d=1;
	al1=2;
	al2=x;
	al3=5;
	b1=3; \\ alpha=1, 2 or 3
	b2=n;
	b3=n; \\ 0 \leq beta<n
	nLB=10^6; \\ assumption at the start of proof of Theorem 6.3
	logXLB=floor(nLB/2600-1.5);
	if(dbg!=0,
		printf("calced logXLB=%5d\n",logXLB);
	);
	logXLB=380;

	hgtA1=log(2);
	absLogA1=abs(log(al1));
	hgtA2=logX;
	absLogA2=5.0001; \\ from Step 3 after defns of a_i's
	hgtA3=log(5);
	absLogA3=abs(log(al3));
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX;
	lamUB0=log(2);

	for(bigL=20,100, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%1==0,print("L=",bigL));
	for(mv=1,10, \\ m \geq 1
		m=mv/1.0;
		for(pv=20,150, \\ rho \geq 2
			rho=pv/10.0;
			a1=(rho-1)*absLogA1+2*logX;
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

\\ 4 Jan 2022
eg1_search_it3(nUBInit,dbg=0)={
	my(al1,al2,al3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,degenNUB,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,minNUB,nLB,nonDegenNUB,nUB,rho,val);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	minNUB=nUBInit;
	d=1;
	al1=2;
	al2=x;
	al3=5;
	b1=3; \\ alpha=1, 2 or 3
	b2=n;
	b3=n; \\ 0 \leq beta<n
	nLB=10^6; \\ assumption at the start of proof of Theorem 6.3
	logXLB=floor(nLB/2600-1.5);
	if(dbg!=0,
		printf("calced logXLB=%5d\n",logXLB);
	);
	logXLB=380;

	hgtA1=log(2);
	absLogA1=abs(log(al1));
	hgtA2=logX;
	absLogA2=5.0001; \\ from Step 3 after defns of a_i's
	hgtA3=log(5);
	absLogA3=abs(log(al3));
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX;
	lamUB0=log(2);

	for(bigL=45,150, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%1==0,print("L=",bigL));
	for(mv=140,190, \\ m \geq 1
		m=mv/10.0;
		\\for(pv=20,300,
		for(pv=60,90, \\ rho \geq 2
			rho=pv/10.0;
			a1=(rho-1)*absLogA1+2*logX;
			a2=(rho+1)*log(al2);
			a3=(rho+1)*log(al3);
			for(cv=10,20,
				chi=cv/10.0;
				val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
				minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
			);
		);
	);
	);
}

\\ 4 Jan 2022
eg1_search_it4(nUBInit,dbg=0)={
	my(al1,al2,al3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigK,chi,d,degenNUB,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logXLB,m,minNUB,nLB,nonDegenNUB,nUB,rho,val);

	if(nUBInit<0.00001,
		printf("ERROR: nUBInit=%9.6f must be a positive real number\n",nUBInit);
		return();
	);
	minNUB=nUBInit;
	d=1;
	al1=2;
	al2=x;
	al3=5;
	b1=3; \\ alpha=1, 2 or 3
	b2=n;
	b3=n; \\ 0 \leq beta<n
	nLB=10^6; \\ assumption at the start of proof of Theorem 6.3
	logXLB=floor(nLB/2600-1.5);
	if(dbg!=0,
		printf("calced logXLB=%5d\n",logXLB);
	);
	logXLB=380;

	hgtA1=log(2);
	absLogA1=abs(log(al1));
	hgtA2=logX;
	absLogA2=5.0001; \\ from Step 3 after defns of a_i's
	hgtA3=log(5);
	absLogA3=abs(log(al3));
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX;
	lamUB0=log(2);

	for(bigL=36,150, \\ L=5 is the lower bound in Theorem 4.1
	if(bigL%1==0,print("L=",bigL));
	for(mv=280,500, \\ m \geq 1
		m=mv/100.0;
		for(pv=650,850,
			rho=pv/100.0;
			a1=(rho-1)*absLogA1+2*logX;
			a2=(rho+1)*log(al2);
			a3=(rho+1)*log(al3);
			for(cv=300,400,
				chi=cv/100.0;
				val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
				minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
			);
		);
	);
	);
}