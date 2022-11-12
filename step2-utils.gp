\\ \r lfl3\step2-utils.gp

\\ read("lfl3\\lfl-utils-general.gp");

\\ 17 Mar 2022
get_matveev_ubnd(bigD,chi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,dbg=0)={
	my(a,actB1,actB2,actB3,b1gA1,b1gA2,b1gA3,bigB,bigBCnst,c1,lamLB,nUB,rts,t);

	if(type(al1)!="t_POL",
		if(dbg!=0,
			print("get_matveev_ubnd(): swapping al1 and al2 to try to make al1 a polynomial");
		);
		t=al1;
		al1=al2;
		al2=t;

		t=absLogA1;
		absLogA1=absLogA2;
		absLogA2=t;

		t=hgtA1;
		hgtA1=hgtA2;
		hgtA2=t;

		t=b1;
		b1=b2;
		b2=t;
	);
	if(type(al1)!="t_POL",
		if(dbg!=0,
			print("get_matveev_ubnd(): swapping al1 and al3 to try to make al1 a polynomial");
		);
		t=al1;
		al1=al3;
		al3=t;

		t=absLogA1;
		absLogA1=absLogA3;
		absLogA3=t;

		t=hgtA1;
		hgtA1=hgtA3;
		hgtA3=t;

		t=b1;
		b1=b3;
		b3=t;
	);
	if(type(al1)!="t_POL",
		print("ERROR in get_matveev_ubnd(): type(al1)=",type(al1)," must be t_POL");
		return([]);
	);
	if(type(al2)!="t_INT" && type(al2)!="t_FRAC" && type(al2)!="t_REAL" && type(al2)!="t_COMPLEX",
		print("ERROR in get_matveev_ubnd(): type(al2)=",type(al2)," must be t_INT, t_FRAC, t_COMPLEX or t_REAL");
		return([]);
	);
	if(type(al3)!="t_INT" && type(al3)!="t_FRAC" && type(al3)!="t_REAL" && type(al3)!="t_COMPLEX",
		print("ERROR in get_matveev_ubnd(): type(al3)=",type(al3)," must be t_INT, t_FRAC, t_COMPLEX or t_REAL");
		return([]);
	);

	if(dbg!=0,
		print("get_matveev_ubnd(): bigD*hgtA1-absLogA1=",bigD*hgtA1-absLogA1);
	);
	if(bigD*hgtA1-absLogA1!=0,
		rts=polrootsreal(bigD*hgtA1-absLogA1);
		if(length(rts)>0 && rts[1]<logXLB,
			bigA1=bigD*hgtA1;
		);
		if(length(rts)>0 && rts[1]>logXLB,
			bigA1=absLogA1;
		);
	);
	if(bigD*hgtA1-absLogA1==0,
		bigA1=bigD*hgtA1; \\ they are both equal, so just pick one
	);

	bigA2=max(bigD*hgtA2,absLogA2);

	bigA3=max(bigD*hgtA3,absLogA3);
	if(dbg!=0,
		print("bigA1=",bigA1);
		print("bigA2=",bigA2);
		print("bigA3=",bigA3);
	);
	
	\\ need to re-order the alpha_i's, so bigA1 is largest
	actB1=b1;
	actB2=b2;
	actB3=b3;
	rts=polrootsreal(bigA1-bigA2);
	if(length(rts)>0 && rts[1]>logXLB,
		t=bigA1;
		bigA1=bigA2;
		bigA2=t;
		t=actB1;
		actB1=actB2;
		actB2=t;
	);
	rts=polrootsreal(bigA1-bigA3);
	if(length(rts)>0 && rts[1]>logXLB,
		t=bigA1;
		bigA1=bigA3;
		bigA3=t;
		t=actB1;
		actB1=actB3;
		actB3=t;
	);
	if(dbg!=0,
		print("get_matveev_ubnd(): actB1=",actB1,", actB2=",actB2,", actB3=",actB3);
		\\print("b1-b2*bigA2/bigA1=",actB1-actB2*bigA2/bigA1);
	);
	
	rts=polrootsreal(actB1-subst(actB2*bigA2/bigA1,logX,logXLB));
	\\print("b1 rts=",rts);
	if(length(rts)>0 && rts[1]<nLB,
		bigB=actB1;
		\\print("bigB=b1=",bigB);
	);
	if(length(rts)>0 && rts[1]>nLB,
		bigB=actB2*bigA2/bigA1;
		\\print("bigB=b2");
	);
	\\print("bigB=",bigB);

	rts=polrootsreal(bigB-subst(actB3*bigA3/bigA1,logX,logXLB));
	if(length(rts)>0 && rts[1]>nLB,
		bigB=actB3*bigA3/bigA1;
		\\print("bigB=b3");
	);
	\\print("bigB=",bigB);

	\\ make bigB of the form constant*n
	\\
	bigBCnst=subst(bigB,n,nLB)/nLB;
	bigB=bigBCnst*n;
	if(dbg!=0,
		print("get_matveev_ubnd(): bigB=",bigB);
	);
	
	c1=5*16^5/6/chi*exp(3.0)*(7+2*chi)*(3*exp(1)/2)^chi*(26.25+log(bigD*bigD*log(exp(1)*bigD)));
	lamLB=-c1*bigD*bigD*bigA1*bigA2*bigA3*(log(1.5*exp(1.0)*bigD*log(exp(1)*bigD))+log(bigBCnst)+logN);
	if(dbg!=0,
		print("get_matveev_ubnd(): lamLB=",lamLB);
	);

	\\ for solve n=a+b*(log n)^h
	a=(polcoef(lamLB,0,logN)-lamUB0)/lamUB1;
	a=subst(a,logX,logXLB);
	b=polcoef(lamLB,1,logN)/lamUB1;
	if(dbg!=0,
		print("get_matveev_ubnd(): a=",a,", b=",b,", polcoef(lamLB,1,logN)=",polcoef(lamLB,1,logN),", lamUB1=",lamUB1);
	);
	if(poldegree(numerator(b),logX)>poldegree(denominator(b),logX),
		print("BAD in get_matveev_ubnd(): numerator of b has larger degree than denominator. b=",b);
	);
	b=subst(b,logX,logXLB);

	h=1;
	nUB=get_solnUB(a,b,h,dbg);
	if(dbg!=0,
		printf("get_matveev_ubnd(): nUB=%12.6e\n\n",nUB);
	);
	return(nUB);
}

\\ test using the Compositio example
\\ 9 May 2022
eg4_test_matveev()={
	my(absLogA1,absLogA2,absLogA2a,absLogA2b,absLogA3,al1,al2,al3,b1,b2,b3,bigD,d,hgtA1,hgtA2,hgtA3,logXLB,nLB,lamUB1,lamUB0,matveevChi,nUB);
	
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
	logXLB=log(1.99*10^6);
	nLB=2*10^6;
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1=-logX/2;
	lamUB0=log(2.2*sqrt(7));
	nUB=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,1);
}