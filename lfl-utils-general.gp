\\ \r lfl3\lfl-utils-general.gp

\\ 17 Mar 2022
get_matveev_ubnd(bigD,chi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,dbg=0)={
	my(a,actB1,actB2,actB3,b1gA1,b1gA2,b1gA3,bigB,bigBCnst,c1,lamLB,nUB,rts,t);

	if(type(al1)!="t_POL",
		if(dbg!=0,
			print("in get_matveev_ubnd(): swapping al1 and al2 to try to make al1 a polynomial");
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
			print("in get_matveev_ubnd(): swapping al1 and al3 to try to make al1 a polynomial");
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
		print("in get_matveev_ubnd(): bigD*hgtA1-absLogA1=",bigD*hgtA1-absLogA1);
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
		print("in get_matveev_ubnd(): actB1=",actB1,", actB2=",actB2,", actB3=",actB3);
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
		print("in get_matveev_ubnd(): bigB=",bigB);
	);
	
	c1=5*16^5/6/chi*exp(3.0)*(7+2*chi)*(3*exp(1)/2)^chi*(26.25+log(bigD*bigD*log(exp(1)*bigD)));
	lamLB=-c1*bigD*bigD*bigA1*bigA2*bigA3*(log(1.5*exp(1.0)*bigD*log(exp(1)*bigD))+log(bigBCnst)+logN);
	if(dbg!=0,
		print("in get_matveev_ubnd(): lamLB=",lamLB);
	);

	\\ for solve n=a+b*(log n)^h
	a=(polcoef(lamLB,0,logN)-lamUB0)/lamUB1;
	a=subst(a,logX,logXLB);
	b=polcoef(lamLB,1,logN)/lamUB1;
	if(dbg!=0,
		print("in get_matveev_ubnd(): a=",a,", b=",b,", polcoef(lamLB,1,logN)=",polcoef(lamLB,1,logN),", lamUB1=",lamUB1);
	);
	if(poldegree(numerator(b),logX)>poldegree(denominator(b),logX),
		print("BAD: numerator of b has larger degree than denominator. b=",b);
	);
	b=subst(b,logX,logXLB);

	h=1;
	nUB=get_solnUB(a,b,h,dbg);
	if(dbg!=0,
		printf("in get_matveev_ubnd(): nUB=%12.6e\n\n",nUB);
	);
	return(nUB);
}

\\ pull out common code for alpha1Variable and alpha3Variable cases
\\ 6 July 2022
get_eqn42(a1,a2,a3,bigK,bigL,bigR,bigS,bigT,d,rho,logBUB,nUB,logXLB,dbg=0)={
	my(eqn42,eqn42LHS,eqn42Rem,eqn42RHS,eqn42X0,eqn42X1,gDenom,gNumer);
	
	\\ with K_0=2(K-1):
	eqn42LHS=(bigK*bigL/2+bigL/4+bigL/8/bigK-2*bigK/3/bigL-1)*log(rho);
	\\ with K_0=K-1, rather than K_0=2(K-1)
	\\eqn42LHS=(bigK*bigL-bigK-bigK/3/bigL)*log(rho);

	gDenom=12*bigR*bigS*bigT/logX/logX;
	gDenom=subst(gDenom,logX,logXLB)*logX*logX;
	gNumer=bigK*bigK*bigL;
	gNumer=polcoef(gNumer,2,logX)*logX*logX;
	g=1/4-gNumer/gDenom;
	if(dbg!=0,
		printf("g=%8.6f\n",g);
	);
	
	\\ first the (cD+1)*log(N) term
	eqn42RHS=(d+1)*log(polcoef(bigK,1)*polcoef(bigK,1)*bigL)+(d+1)*log(logXLB)/logXLB*logX;
	\\ next the gL(a1*R+a2*S+a3*T) term
	eqn42RHS=eqn42RHS+g*bigL*(a1*bigR+a2*bigS+a3*bigT);
	\\ lastly the D(K-1)*log(b) term
	eqn42RHS=eqn42RHS+d*(bigK-1)*logBUB;
	eqn42RHS=subst(eqn42RHS,logN,log(nUB));
	eqn42=eqn42LHS-eqn42RHS;
	eqn42X1=pollead(polcoef(eqn42,1))*logX;
	if(dbg!=0,
		printf("before simplifying: deg=%4d, eqn42=%s\n",poldegree(eqn42),eqn42);
		print("pollead(polcoef(eqn42,0))=",pollead(polcoef(eqn42,0)));
		print("subst(eqn42-eqn42X1,logX,logXLB)=",subst(eqn42-eqn42X1,logX,logXLB));
	);
	eqn42X0=pollead(polcoef(eqn42,0));
	eqn42Rem=max(0, pollead(subst(eqn42-eqn42X1-eqn42X0,logX,logXLB)));
	eqn42=eqn42X1+eqn42X0+eqn42Rem;
	if(dbg!=0,
		printf("eqn42LHS=%s\n",eqn42LHS);
		printf("eqn42RHS D(K-1)*log(b) term=%s\n",d*(bigK-1)*logBUB);
		printf("eqn42RHS gL(a1*R+a2*S+a3*T) term=%s\n",g*bigL*(a1*bigR+a2*bigS+a3*bigT));
		printf("eqn42RHS (cD+1)*log(N) term=%s\n",(d+1)*log(polcoef(bigK,1)*polcoef(bigK,1)*bigL)+(d+1)*log(logXLB)/logXLB*logX);
		printf("eqn42RHS=%s\n",eqn42RHS);
		printf("eqn42LHS-eqn42RHS=%s\n",eqn42);
	);
	return(eqn42);
}

\\ to save the same duplicated code
\\ 3 March 2022
update_minNUB(val,bigL,m,rho,chi,minNUB,dbg=0)={
	my(bigK,localMinNUB,minB,minNDegenUB,nDegenUB,nNonDegenUB,nUB);
	
	localMinNUB=minNUB;
	if(length(val)==3,
		bigK=val[1];
		nNonDegenUB=val[2];
		nDegenUB=val[3];
		if(dbg!=0,
			print("nNonDegenUB=",nNonDegenUB,", type=",type(nNonDegenUB));
			print("nDegenUB=",nDegenUB,", type=",type(nDegenUB));
		);
		\\ initialise minNDegenUB
		minNDegenUB=0;
		if(nDegenUB[1]>0,minNDegenUB=nDegenUB[1]);
		if(nDegenUB[2]>0 && minNDegenUB==0,minNDegenUB=nDegenUB[2]);
		if(nDegenUB[3]>0 && minNDegenUB==0,minNDegenUB=nDegenUB[3]);
		if(nDegenUB[1]>0,minNDegenUB=min(minNDegenUB,nDegenUB[1]));
		if(nDegenUB[2]>0,minNDegenUB=min(minNDegenUB,nDegenUB[2]));
		if(nDegenUB[3]>0,minNDegenUB=min(minNDegenUB,nDegenUB[3]));
		nUB=max(minNDegenUB,nNonDegenUB);
		if(minNDegenUB>0 && nNonDegenUB>0 && nUB<minNUB,
			localMinNUB=nUB;
			minB="b1";
			if(nDegenUB[2]<nDegenUB[1],
				minB="b3";
			);
			printf("L=%4d, m=%8.4f, rho=%8.4f, chi=%7.4f, K=%12.3f*logX, nonDegen log|Lambda|>%9.6e*logX, nonDegenNUB=%10.6e, degenNUB1=%10.6e, degenNUB2=%10.6e, degenNUB3=%10.6e, nUB=%10.6e, transB-%s\n",bigL,m,rho,chi,bigK,-bigK*bigL*log(rho),nNonDegenUB,nDegenUB[1],nDegenUB[2],nDegenUB[3],nUB,minB);
		);
	);
	return(localMinNUB);
}

\\ on 25 Nov 2021
round_down(num,places)={
	if(num<0 || num==0,
		print("num=",num," must be positive");
		return();
	);
	e=floor(log(num)/log(10))-(places-1);
	\\print("num=",num,", e=",e);
	v=floor(num/10.0^e)*10.0^e;
	\\print("v=",v);
	return(v);
}

\\ on 25 Nov 2021
round_up(num,places)={
	if(num<0 || num==0,
		print("num=",num," must be positive");
		return();
	);
	e=floor(log(num)/log(10))-(places-1);
	v=ceil(num/10.0^e)*10.0^e;
	return(v);
}

\\ for solve x=a+b*(log x)^h
\\ uses proof of Lemma 2.2 from the Petho-de Weger paper
\\
\\ the rescalingFactor *was* used for solving
\\ (x/rescalingFactor)=a+b*(log (x/rescalingFactor))^h
\\ 21 Nov 2021
\\get_solnUB(a,b,h,rescalingFactor,dbg=0)={
get_solnUB(a,b,h,dbg=0)={
	my(c,logC,xHUB);
	
	c=h*b^(1/h);
	if(dbg!=0,
		printf("in get_solnUB(): a=%9.6f, b=%9.6f, c=%9.6f\n",a,b,c);
	);
	logC=log(c);
	xHUB=c*logC+logC/(logC-1)*(a^(1/h)+c*log(logC));
	if(dbg!=0,
		printf("in get_solnUB(): a=%9.6f, b=%9.6f, c=%9.6f, xHUB=%9.6f\n",a,b,c,xHUB);
	);
	\\return(xHUB^h*rescalingFactor);
	return(xHUB^h);
}

\\ from Table in Laurent's 2008 paper
\\ 11 May 2022
get_lflComplexCnst(bPLB,d)={
	my(lflCnst);
	
	if(d*bPLB<10,
		\\printf("in get_lflCnst(): bPLB=%9.6f is too small.\n",bPLB);
		lflCnst=1000000000.0;
	);
	if(d*bPLB>10,
		lflCnst=32.3;
	);
	if(d*bPLB>12,
		lflCnst=29.9;
	);
	if(d*bPLB>14,
		lflCnst=28.2;
	);
	if(d*bPLB>16,
		lflCnst=26.9;
	);
	if(d*bPLB>18,
		lflCnst=26.0;
	);
	if(d*bPLB>20,
		lflCnst=25.2;
	);
	if(d*bPLB>22,
		lflCnst=24.5;
	);
	if(d*bPLB>24,
		lflCnst=24.0;
	);
	if(d*bPLB>26,
		lflCnst=23.5;
	);
	if(d*bPLB>28,
		lflCnst=23.1;
	);
	if(d*bPLB>30,
		lflCnst=22.8;
	);
	return(lflCnst);
}

\\ from Table in Laurent's 2008 paper
\\ 14 Dec 2021 (pulled out of above and more cases considered
get_lflRealCnst(bPLB,d)={
	my(lflCnst);
	
	if(d*bPLB<10,
		\\printf("in get_lflCnst(): bPLB=%9.6f is too small.\n",bPLB);
		lflCnst=1000000000.0;
	);
	if(d*bPLB>10,
		lflCnst=25.2;
	);
	if(d*bPLB>12,
		lflCnst=23.4;
	);
	if(d*bPLB>14,
		lflCnst=22.1;
	);
	if(d*bPLB>16,
		lflCnst=21.1;
	);
	if(d*bPLB>18,
		lflCnst=20.3;
	);
	if(d*bPLB>20,
		lflCnst=19.7;
	);
	if(d*bPLB>22,
		lflCnst=19.2;
	);
	if(d*bPLB>24,
		lflCnst=18.8;
	);
	if(d*bPLB>26,
		lflCnst=18.4;
	);
	if(d*bPLB>28,
		lflCnst=18.1;
	);
	if(d*bPLB>30,
		lflCnst=17.9;
	);
	return(lflCnst);
}

\\ here we eliminate the b_1*log(alpha1) term (hence function name)
\\ b1*log(a1)+b2*log(a2)-b3*log(a3) and u1*b1+u2*b2+u3*b3=0, so
\\ u1*b1*log(a1)+u1*b2*log(a2)-u1*b3*log(a3)
\\ =-(u2*b2+u3*b3)*log(a1)+u1*b2*log(a2)-u1*b3*log(a3)
\\ =b2*log(a1^(-u2)*a2^u1)-b3*log(a1^u3*a3^u1)
\\
\\ if u1=0, then u2*b2+u3*b3=0, so we eliminate either b2 or b3 depending on which one is constant
\\ if we eliminate b2, then
\\ u2*b1*log(a1)+u2*b2*log(a2)-u2*b3*log(a3)=u2*b1*log(a1)-u3*b3*log(a2)-u2*b3*log(a3)
\\ =b1*log(a1^u2)-b3*log(a2^u3*a3^u2)
\\
\\ if we eliminate b3, then
\\ u3*b1*log(a1)+u3*b2*log(a2)-u3*b3*log(a3)=u3*b1*log(a1)+u3*b2*log(a2)+u2*b2*log(a3)
\\ =b1*log(a1^u3)-b2*log(a2^(-u3)*a3^(-u2))
\\
\\ we will always return [log(A_1), log(A_2)]
\\ such that b_i*log(alpha_1)-b_j*log(alpha_2) with i<j
\\ 3 Jan 2022
get_A1A2_from_b1(d,absLogA1,hgtA1,absLogA2,hgtA2,absLogA3,hgtA3,logXLB,u1UB,u2UB,u3UB,dbg=0)={
	my(hgtNewA1,hgtNewA2,logA1,logA2,logNewA1,logNewA2);
	
	\\ "New" as it is for A1 and A2 in linear form in two logs
	hgtNewA1=u1UB*hgtA2+u2UB*hgtA1;
	if(poldegree(hgtNewA1)==1,
		hgtNewA1=subst(hgtNewA1,logX,logXLB)/logXLB*logX;
	);
	logNewA1=u1UB*absLogA2+u2UB*absLogA1;
	if(poldegree(logNewA1)==1,
		logNewA1=subst(logNewA1,logX,logXLB)/logXLB*logX;
	);
	logA1=get_logA(hgtNewA1,logNewA1,logXLB,d,dbg);
	if(dbg!=0,
		\\printf("\nget_A1A2_from_b1(): hgtA1=%s\n",hgtA1);
		\\printf("get_A1A2_from_b1(): hgtA2=%s\n",hgtA2);
		\\print("d=",d,", absLogA1=",absLogA1,", hgtA1=",hgtA1);
		\\print("absLogA2=",absLogA2,", hgtA2=",hgtA2);
		\\print("absLogA3=",absLogA3,", hgtA3=",hgtA3);
		\\print("u1UB=",u1UB,", u2UB=",u2UB,", u3UB=",u3UB);
		printf("get_A1A2_from_b1(): hgtNewA1=%s\n",hgtNewA1);
		printf("get_A1A2_from_b1(): logNewA1=%s\n",logNewA1);
		printf("get_A1A2_from_b1(): logA1=%s\n",logA1);
	);

	hgtNewA2=u3UB*hgtA1+u1UB*hgtA3;
	if(poldegree(hgtNewA2)==1,
		hgtNewA2=subst(hgtNewA2,logX,logXLB)/logXLB*logX;
	);
	logNewA2=u3UB*absLogA1+u1UB*absLogA3;
	if(poldegree(logNewA2)==1,
		logNewA2=subst(logNewA2,logX,logXLB)/logXLB*logX;
	);
	logA2=get_logA(hgtNewA2,logNewA2,logXLB,d,dbg);
	if(dbg!=0,
		printf("get_A1A2_from_b1(): hgtNewA2=%s\n",hgtNewA2);
		printf("get_A1A2_from_b1(): logNewA2=%s\n",logNewA2);
		printf("get_A1A2_from_b1(): logA2=%s\n",logA2);
	);

	\\ get values if u_1=0
	if(type(u2UB)!="t_POL",
		\\ linear form is: b1*log(a1^u2)-b3*log(a2^u3*a3^u2)
		if(dbg!=0,
			print("get_A1A2_from_b1(): considering u_1=0 and eliminating b_2");
		);
		hgtNewA1_0=u2UB*hgtA1;
		logNewA1_0=u2UB*absLogA1;
		logA1_0=get_logA(hgtNewA1_0,logNewA1_0,logXLB,d,dbg);
		hgtNewA2_0=u3UB*hgtA2+u2UB*hgtA3;
		logNewA2_0=u3UB*absLogA2+u2UB*absLogA3;
		logA2_0=get_logA(hgtNewA2_0,logNewA2_0,logXLB,d,dbg);
	);
	if(type(u2UB)=="t_POL" && type(u3UB)!="t_POL",
		\\ linear form is: b1*log(a1^u3)-b2*log(a2^(-u3)*a3^(-u2))

		if(dbg!=0,
			print("get_a1a2_from_b1(): considering u_1=0 and eliminating b_3");
		);
		hgtNewA1_0=u3UB*hgtA1;
		logNewA1_0=u3UB*absLogA1;
		hgtNewA2_0=u3UB*hgtA2+u2UB*hgtA3;
		logNewA2_0=u3UB*absLogA2+u2UB*absLogA3;
	);
	if(type(u2UB)=="t_POL" && type(u3UB)=="t_POL",
		print("ERROR in get_A1A2_from_b1(): u2=",u2UB," and u3=",u3UB," are both polynomials");
		return([]);
	);
	if(dbg!=0,
		printf("get_A1A2_from_b1(): hgtNewA1_0=%s\n",hgtNewA1_0);
		printf("get_A1A2_from_b1(): logNewA1_0=%s\n",logNewA1_0);
		printf("get_A1A2_from_b1(): logA1_0=%s\n",logA1_0);
		printf("get_A1A2_from_b1(): hgtNewA2_0=%s\n",hgtNewA2_0);
		printf("get_A1A2_from_b1(): logNewA2_0=%s\n",logNewA2_0);
		printf("get_A1A2_from_b1(): logA2_0=%s\n",logA2_0);
	);
	return([logA1,logA2,logA1_0,logA2_0]);
}

\\ here we eliminate the b_2*log(alpha2) term (hence function name)
\\ b1*log(a1)+b2*log(a2)-b3*log(a3) and u1*b1+u2*b2+u3*b3=0, so
\\ u2*b1*log(a1)+u2*b2*log(a2)-u2*b3*log(a3)=u2*b1*log(a1)-(u1*b1+u3*b3)*log(a2)-u2*b3*log(a3)
\\ =b1*log(a1^u2*a2^(-u1)) - b3*log(a2^u3*a3^u2)
\\
\\ if u2=0, then u1*b1+u3*b3=0, so we eliminate either b1 or b3 depending on which one is constant
\\ if we eliminate b3, then
\\ u3*b1*log(a1)+u3*b2*log(a2)-u3*b3*log(a3)=u3*b1*log(a1)+u3*b2*log(a2)+u1*b1*log(a3)
\\ =b1*log(a1^u3*a3^u1)-b2*log(a2^(-u3))
\\
\\ if we eliminate b1, then
\\ u1*b1*log(a1)+u1*b2*log(a2)-u1*b3*log(a3)=-u3*b3*log(a1)+u1*b2*log(a2)-u1*b3*log(a3)
\\ =b2*log(a2^u1)-b3*log(a1^u3*a3^u1)
\\
\\ we will always return [log(A_1), log(A_2)]
\\ such that b_i*log(alpha_1)-b_j*log(alpha_2) with i<j
\\ 3 Jan 2022
get_A1A2_from_b2(d,absLogA1,hgtA1,absLogA2,hgtA2,absLogA3,hgtA3,logXLB,u1UB,u2UB,u3UB,dbg=0)={
	my(hgtNewA1,hgtNewA2,logA1,logA2,logNewA1,logNewA2);
	
	\\ "New" as it is for A1 and A2 in linear form in two logs
	hgtNewA1=u1UB*hgtA2+u2UB*hgtA1;
	logNewA1=u1UB*absLogA2+u2UB*absLogA1;
	logA1=get_logA(hgtNewA1,logNewA1,logXLB,d,dbg);
	if(dbg!=0,
		printf("get_A1A2_from_b2(): hgtNewA1=%s\n",hgtNewA1);
		printf("get_A1A2_from_b2(): logNewA1=%s\n",logNewA1);
		printf("get_A1A2_from_b2(): logA1=%s\n",logA1);
	);

	hgtNewA2=u3UB*hgtA2+u2UB*hgtA3;
	logNewA2=u3UB*absLogA2+u2UB*absLogA3;
	logA2=get_logA(hgtNewA2,logNewA2,logXLB,d,dbg);
	if(dbg!=0,
		printf("get_A1A2_from_b2(): hgtNewA2=%s\n",hgtNewA2);
		printf("get_A1A2_from_b2(): logNewA2=%s\n",logNewA2);
		printf("get_A1A2_from_b2(): logA2=%s\n",logA2);
	);
	
	\\ get values if u_2=0
	if(type(u3UB)!="t_POL",
		\\ linear form is: b1*log(a1^u3*a3^u1)-b2*log(a2^(-u3))
		if(dbg!=0,
			print("get_A1A2_from_b2(): considering u_2=0 and eliminating b_3");
		);
		hgtNewA1_0=u3UB*hgtA1+u1UB*hgtA3;
		logNewA1_0=u3UB*absLogA1+u1UB*absLogA3;
		logA1_0=get_logA(hgtNewA1_0,logNewA1_0,logXLB,d,dbg);
		hgtNewA2_0=u3UB*hgtA2;
		logNewA2_0=u3UB*absLogA2;
		logA2_0=get_logA(hgtNewA2_0,logNewA2_0,logXLB,d,dbg);
	);
	if(type(u3UB)=="t_POL" && type(u1UB)!="t_POL",
		\\ linear form is: b2*log(a2^u1)-b3*log(a1^u3*a3^u1)
		if(dbg!=0,
			print("get_A1A2_from_b2(): considering u_2=0 and eliminating b_1");
		);
		hgtNewA1_0=u1UB*hgtA2;
		logNewA1_0=u1UB*absLogA2;
		logA1_0=get_logA(hgtNewA1_0,logNewA1_0,logXLB,d,dbg);
		hgtNewA2_0=u3UB*hgtA1+u1UB*hgtA3;
		logNewA2_0=u3UB*absLogA1+u1UB*absLogA3;
		logA2_0=get_logA(hgtNewA2_0,logNewA2_0,logXLB,d,dbg);
	);
	if(type(u1UB)=="t_POL" && type(u3UB)=="t_POL",
		print("ERROR in get_A1A2_from_b2(): u1=",u1UB," and u3=",u3UB," are both polynomials");
		return([]);
	);

	if(dbg!=0,
		printf("get_A1A2_from_b2(): hgtNewA1_0=%s\n",hgtNewA1_0);
		printf("get_A1A2_from_b2(): logNewA1_0=%s\n",logNewA1_0);
		printf("get_A1A2_from_b2(): logA1_0=%s\n",logA1_0);
		printf("get_A1A2_from_b2(): hgtNewA2_0=%s\n",hgtNewA2_0);
		printf("get_A1A2_from_b2(): logNewA2_0=%s\n",logNewA2_0);
		printf("get_A1A2_from_b2(): logA2_0=%s\n",logA2_0);
	);
	return([logA1,logA2,logA1_0,logA2_0]);
}

\\ here we eliminate the b_3*log(alpha3) term (hence function name)
\\ b1*log(a1)+b2*log(a2)-b3*log(a3) and u1*b1+u2*b2+uu3*b3=0, so
\\ u3*b1*log(a1)+u3*b2*log(a2)+u3*b3*log(a3)=u3*b1*log(a1)+u3*b2*log(a2)+(u1*b1+u2*b2)*log(a3)
\\ =b1*log(a1^u3*a3^u1)-b2*log(a2^(-u3)*a3^(-u2))
\\
\\ if u3=0, then u1*b1+u2*b2=0, so we eliminate either b1 or b2 depending on which one is constant
\\ if we eliminate b2, then
\\ u2*b1*log(a1)+u2*b2*log(a2)-u2*b3*log(a3)=u2*b1*log(a1)-(u1*b1)*log(a2)-u2*b3*log(a3)
\\ =b1*log(a1^u2*a2^(-u1))-b3*log(a3^u2)
\\
\\ if we eliminate b1, then
\\ u1*b1*log(a1)+u1*b2*log(a2)-u1*b3*log(a3)=-u2*b2*log(a1)+u1*b2*log(a2)-u1*b3*log(a3)
\\ =b2*log(a1^(-u2)*a2^u1)-b3*log(a3^(-u1))
\\
\\ we will always return [log(A_1), log(A_2)]
\\ such that b_i*log(alpha_1)-b_j*log(alpha_2) with i<j
\\
\\ 3 Jan 2022
get_A1A2_from_b3(d,absLogA1,hgtA1,absLogA2,hgtA2,absLogA3,hgtA3,logXLB,u1UB,u2UB,u3UB,dbg=0)={
	my(hgtNewA1,hgtNewA1_0,hgtNewA2,hgtNewA2_0,logA1,logA1_0,logA2,logA2_0,logNewA1,logNewA1_0,logNewA2,logNewA2_0);
	
	\\ "New" as it is for A1 and A2 in linear form in two logs
	\\ linear form = b1*log(a1^u3*a3^u1)-b2*log(a2^(-u3)*a3^(-u2))
	hgtNewA1=u3UB*hgtA1+u1UB*hgtA3;
	logNewA1=u3UB*absLogA1+u1UB*absLogA3;
	logA1=get_logA(hgtNewA1,logNewA1,logXLB,d,dbg);
	if(dbg!=0,
		printf("get_A1A2_from_b3(): hgtNewA1=%s\n",hgtNewA1);
		printf("get_A1A2_from_b3(): logNewA1=%s\n",logNewA1);
		printf("get_A1A2_from_b3(): logA1=%s\n",logA1);
	);

	hgtNewA2=u3UB*hgtA2+u2UB*hgtA3;
	logNewA2=u3UB*absLogA2+u2UB*absLogA3;
	logA2=get_logA(hgtNewA2,logNewA2,logXLB,d,dbg);
	if(dbg!=0,
		printf("get_A1A2_from_b3(): hgtNewA2=%s\n",hgtNewA2);
		printf("get_A1A2_from_b3(): logNewA2=%s\n",logNewA2);
		printf("get_A1A2_from_b3(): logA2=%s\n",logA2);
	);

	\\ get values if u_3=0
	if(type(u1UB)!="t_POL",
		\\ linear form is: b2*log(a1^(-u2)*a2^u1)-b3*log(a3^(-u1))
		if(dbg!=0,
			print("get_A1A2_from_b3(): considering u_3=0 and eliminating b_1");
		);
		hgtNewA1_0=u2UB*hgtA1+u1UB*hgtA2;
		logNewA1_0=u2UB*absLogA1+u1UB*absLogA2;
		logA1_0=get_logA(hgtNewA1_0,logNewA1_0,logXLB,d,dbg);
		hgtNewA2_0=u1UB*hgtA3;
		logNewA2_0=u1UB*absLogA1;
		logA2_0=get_logA(hgtNewA2_0,logNewA2_0,logXLB,d,dbg);
	);
	if(type(u1UB)=="t_POL" && type(u2UB)!="t_POL",
		\\ linear form is: b1*log(a1^u2*a2^(-u1))-b3*log(a3^u2)
		if(dbg!=0,
			print("get_A1A2_from_b3(): considering u_3=0 and eliminating b_2");
		);
		hgtNewA1_0=u2UB*hgtA1+u1UB*hgtA2;
		logNewA1_0=u2UB*absLogA1+u1UB*absLogA2;
		logA1_0=get_logA(hgtNewA1_0,logNewA1_0,logXLB,d,dbg);
		hgtNewA2_0=u2UB*hgtA3;
		logNewA2_0=u2UB*absLogA3;
		logA2_0=get_logA(hgtNewA2_0,logNewA2_0,logXLB,d,dbg);
	);
	if(type(u1UB)=="t_POL" && type(u2UB)=="t_POL",
		print("ERROR in get_A1A2_from_b3(): u1=",u1UB," and u2=",u2UB," are both polynomials");
		return([]);
	);

	if(dbg!=0,
		printf("get_A1A2_from_b3(): hgtNewA1_0=%s\n",hgtNewA1_0);
		printf("get_A1A2_from_b3(): logNewA1_0=%s\n",logNewA1_0);
		printf("get_A1A2_from_b3(): logA1_0=%s\n",logA1_0);
		printf("get_A1A2_from_b3(): hgtNewA2_0=%s\n",hgtNewA2_0);
		printf("get_A1A2_from_b3(): logNewA2_0=%s\n",logNewA2_0);
		printf("get_A1A2_from_b3(): logA2_0=%s\n",logA2_0);
	);
	return([logA1,logA2,logA1_0,logA2_0]);
}

\\ returns max(hgt(alpha_i), log |alpha_i|/D, 1/D)
\\ make sure that the logNew argument is log |alpha_i|, not log |alpha_i|/d etc
\\ (we divide logNew by d in this function)
\\ 28 Jan 2022
get_logA(hgtNew,logNew,logXLB,d,dbg=0)={
	my(actHgtNew,actLowNew,aDeg,cf,logA);
	
	logA=0;
	if(type(logNew)=="t_POL" || type(hgtNew)=="t_POL",
		\\ make arguments into monomials
		actHgtNew=subst(hgtNew,logX,logXLB);
		actLogNew=subst(logNew,logX,logXLB)/d; \\ divide by d here
		if(max(actHgtNew,actLogNew)<1.00001/d,
			logA=1/d;
		);
		if(max(actHgtNew,actLogNew)>1/d,
			aDeg=max(poldegree(hgtNew,logX),poldegree(logNew,logX));
			actHgtNew=actHgtNew*(logX/logXLB)^aDeg;
			actLogNew=actLogNew*(logX/logXLB)^aDeg;
			if(dbg!=0,
				print("in get_logA(): hgtNew=",actHgtNew,", logNew=",actLogNew,", aDeg=",aDeg);
			);
			logA=max(polcoef(actHgtNew,aDeg,logX),polcoef(actLowNew,aDeg,logX))*logX^aDeg;
		);
	);
	if(type(logNew)!="t_POL" && type(hgtNew)!="t_POL",
		logA=max(logNew/d,hgtNew); \\ and divide by d here too
		logA=max(1/d,logA);
	);
	return(logA);
}

\\ some checking code that is not presently used (but can be)
\\ 15 Dec 2021
passes_zero_estimate(chi,bigK,bigL,bigM,bigR1,bigR2,bigR3,bigS1,bigS2,bigS3,bigT1,bigT2,bigT3,logXLB,dbg=0)={
	my(calVUB,check1RHS,check1RHS_t1,check1RHS_t2,check2RHS,check3RHS,check5RHS);
	
	calVUB=(bigR1+1)*(bigS1+1)*(bigT1+1);
	calVUB=subst(calVUB,logX,logXLB)/logXLB/logXLB;
	calVUB=sqrt(calVUB);
	
	check1RHS_t1=bigR1+bigS1+1;
	check1RHS_t1=subst(check1RHS_t1,logX,logXLB)/logXLB;
	check1RHS_t2=bigS1+bigT1+1;
	check1RHS_t2=subst(check1RHS_t2,logX,logXLB)/logXLB;
	check1RHS=max(check1RHS_t1,check1RHS_t2);

	check1RHS_t1=bigR1+bigT1+1;
	check1RHS_t1=subst(check1RHS_t1,logX,logXLB)/logXLB;
	check1RHS=max(check1RHS,check1RHS_t1);
	check1RHS=max(check1RHS,chi*calVUB);
	check1RHS_t1=subst(bigK,logX,logXLB)/logXLB;
	check1RHS_t2=subst(bigM,logX,logXLB)/logXLB;
	check1RHS=max(check1RHS_t1,check1RHS_t2)*check1RHS;
	if(dbg!=0,
		printf("(bigR1+1)*(bigS1+1)*(bigT1+1)=%s\n",(bigR1+1)*(bigS1+1)*(bigT1+1));
		printf("check1RHS=%s*logX^2\n",check1RHS);
	);
	if(polcoef((bigR1+1)*(bigS1+1)*(bigT1+1),2)<check1RHS,
		return(0);
	);

	check2RHS=bigL;
	check2RHS=subst(check2RHS,logX,logXLB)/logXLB;
	if(polcoef((bigR1+1)*(bigS1+1)*(bigT1+1),2)<check2RHS,
		if(dbg!=0,
			printf("(bigR1+1)*(bigS1+1)*(bigT1+1)=%s\n",(bigR1+1)*(bigS1+1)*(bigT1+1));
			printf("check2RHS=%s*logX\n",check2RHS);
		);
		return(0);
	);
	
	check3RHS=2*bigK*bigM;
	check3RHS=subst(check3RHS,logX,logXLB)/logXLB/logXLB;
	if(polcoef((bigR2+1)*(bigS2+1)*(bigT2+1),2)<check3RHS,
		if(dbg!=0,
			printf("(bigR2+1)*(bigS2+1)*(bigT2+1)=%s\n",(bigR2+1)*(bigS2+1)*(bigT2+1));
			printf("check3RHS=%s*logX^2\n",check3RHS);
		);
		return(0);
	);

	\\ TODO: check 4:
	
	check5RHS=6*bigK*bigL*bigM;
	check5RHS=subst(check5RHS,logX,logXLB)/logXLB/logXLB/logXLB;
	if(polcoef((bigR3+1)*(bigS3+1)*(bigT3+1),2)<check5RHS,
		if(dbg!=0,
			printf("(bigR3+1)*(bigS3+1)*(bigT3+1)=%s\n",(bigR3+1)*(bigS3+1)*(bigT3+1));
			printf("check5RHS=%s*logX^3\n",check5RHS);
		);
		return(0);
	);
	return(1);
}

\\ lamMul is the factor we have multiplied lambda by to get a linear form in 2 logs
\\ it will be an upper bound for one of the u_i's in the linear relationship between the b_i's
\\ we use Corollary 1 or 2 of Laurent's 2008 paper here
\\ aiArray is an array of [logA1,logA2] for using in Laurent's linear forms results,
\\ upper bounds for logA1 and logA2 in fact. Not an array of the alpha_i values
\\
\\ 3 Jan 2022
get_degen_nUB(aiArray,d,nLB,lamUB0,lamUB1,lamMul,logXLB,isComplex,dbg=0)={
	my(a,absBPLB,b,bDenom,bPCnst,bPCnstUB,h,lambdaLB,lflCnst,logA1,logA2,nUB);

	if(dbg!=0,
		print("\nget_degen_nUB(): START");
	);
	
	if(length(aiArray)!=2,
		print("ERROR: in get_degen_nUB(), length(aiArray)=",length(aiArray),", but must have length 2");
		1/0;
	);
	logA1=aiArray[1];
	logA2=aiArray[2];
	if(type(logA1)=="t_POL",
		bPCnstLB=1/logA2/d;
	);
	if(type(logA2)=="t_POL",
		bPCnstLB=1/logA1/d;
	);
	if(type(logA1)=="t_POL" && type(logA2)=="t_POL",
		print("ERROR in get_degen_nUB(): logA1 and logA2 are both polynomials");
		return([]);
	);
	if(type(logA1)!="t_POL" && type(logA2)!="t_POL",
		bPCnstLB=1/max(logA1,logA2)/d;
	);
	\\if(dbg!=0,
	\\	printf("in get_degen_nUB(): bPCnstLB=%9.6f\n",bPCnstLB);
	\\);
	bPCnstUB=1/logA2/d+1/logA1/d; \\ b'<n*bPCnstUB
	bPCnstUB=subst(bPCnstUB,logX,logXLB);
	if(!isComplex,
	 	bPCnstUB=log(bPCnstUB)+0.38; \\ now log(b')+0.38=bPCnstUB+log(n)
		absBPLB=log(bPCnstLB*nLB)+0.38; \\ an absolute lower bound for b' in Corollary 2 (hence name, bP=b prime)
		\\ using Corollary 2 of Laurent 2008
		lflCnst=get_lflRealCnst(absBPLB,d);
	);
	if(isComplex,
	 	bPCnstUB=log(bPCnstUB)+0.21; \\ now log(b')+0.21=bPCnstUB+log(n)
		absBPLB=log(bPCnstLB*nLB)+0.21; \\ an absolute lower bound for b' in Corollary 2 (hence name, bP=b prime)
		\\ using Corollary 1 of Laurent 2008
		lflCnst=get_lflComplexCnst(absBPLB,d);
	);

	\\ log |lamMul*Lambda|>lambdaLB*(log(b')+0.38,m/D,1/D)^2
	lambdaLB=lflCnst*d*d*d*d*logA1*logA2;
	if(poldegree(lambdaLB)==1,
		lambdaLB=subst(lambdaLB,logX,logXLB)/logXLB*logX;
	);
	if(polcoef(lambdaLB,0)<0,
		print("BAD: constant term of lambdaLB is negative, lambdaLB=",lambdaLB);
	);
	
	if(dbg!=0,
		printf("in get_degen_nUB(): d=%6d\n",d);
		printf("in get_degen_nUB(): logA1=%s\n",logA1);
		printf("in get_degen_nUB(): logA2=%s\n",logA2);
		printf("in get_degen_nUB(): log(b')>log(n)%9.6f\n",log(bPCnstLB)); \\ no sign before %9.6f, as log(bPCnst) is negative
		if(!isComplex,
			printf("in get_degen_nUB(): log(n)%9.6f>log(b')+0.38>log(n)%9.6f>%9.6f\n",bPCnstUB,log(bPCnstLB)+0.38,absBPLB);
		);
		if(isComplex,
			printf("in get_degen_nUB(): log(n)%9.6f>log(b')+0.21>log(n)%9.6f>%9.6f\n",bPCnstUB,log(bPCnstLB)+0.21,absBPLB);
		);
		printf("in get_degen_nUB(): lflCnst=%9.6f\n",lflCnst);
		printf("in get_degen_nUB(): lamMul=%9.6f\n",lamMul);
		printf("in get_degen_nUB(): lambdaLB=%9.6f*logX\n",polcoef(lambdaLB,1,logX));
		printf("in get_degen_nUB(): recall |u_i*Lambda|>lambdaLB*(log(b')+0.38,m/D,1/D)^2\n");
		if(!isComplex,
			printf("in get_degen_nUB(): %s*n+%9.6f>log |u_i*Lambda|>-%9.6f*logX*(log(b')+0.38)^2=-%9.6f*logX*(log(n)%9.6f)^2\n",lamUB1,lamUB0+log(lamMul),polcoef(lambdaLB,1,logX),polcoef(lambdaLB,1,logX),bPCnstUB);
		);
		if(isComplex,
			printf("in get_degen_nUB(): %s*n+%9.6f>log |u_i*Lambda|>-%9.6f*logX*(log(b')+0.21)^2=-%9.6f*logX*(log(n)%9.6f)^2\n",lamUB1,lamUB0+log(lamMul),polcoef(lambdaLB,1,logX),polcoef(lambdaLB,1,logX),bPCnstUB);
		);
	);
	
	\\ log(b')+0.38=log(n/bDenom): used for rescaling in a momment
	bDenom=exp(-bPCnstUB);
	if(type(lamMul)=="t_POL",
		\\a=(-log(polcoef(lamMul,1,logX))-lamUB0)/lamUB1/bDenom-log(logXLB)/subst(lamUB1,logX,logXLB)/bDenom;
		print("BAD lamMul=",lamMul);
		1/0;
	);
	a=-(lamUB0+log(lamMul))/bDenom/lamUB1;
	if(poldegree(numerator(a),logX)>poldegree(denominator(a),logX),
		print("BAD: numerator of a has larger degree than denominator. a=",a);
	);
	a=subst(a,logX,logXLB);
	b=-lambdaLB/lamUB1/bDenom;
	if(poldegree(numerator(b),logX)>poldegree(denominator(b),logX),
		print("BAD: numerator of b has larger degree than denominator. b=",b);
	);
	b=subst(b,logX,logXLB);

	h=2;
	if(dbg!=0,
		printf("for our de Weger-Petho lemma: a=%9.6e, b=%9.6f, bDenom=%9.6f\n",a,b,bDenom);
	);
	nUB=get_solnUB(a,b,h,dbg); \\in fact this is an upper bound for n/bDenom
	nUB=nUB*bDenom;
	if(dbg!=0,
		\\print("logA1=",logA1,", logA2=",logA2);
		\\printf("b'=%8.6f*n\n",bPCnst);
		printf("degen nUB=%9.6f\n",nUB);
		print("get_degen_nUB(): END\n");
	);
	return(nUB);
}

\\ we estimate (from above) each of three parts of expression for b'
\\ in terms of logX and n
\\ where we assume that n is an upper bound for |b_1|, |b_2| and |b_3|
\\ we use the definition of b and Lemma 3.4(a) to bound log(b) from above here
\\ assuming that K=c*log(x)
\\ 21 Nov 2021
get_logBPrime_UB(bigK,bigR,bigS,bigT,b1,b2,b3,d1,d2,logXLB,nUB,dbg=0)={
	my(logB1,logB2,logBUB);
	
	\\ here we break the expression for b (defined in the notation section)
	\\ into three parts: 1) B1=b3'*eta0=(b3/d1)*eta0; 2) B2=b3''*zeta0=(b3/d2)*zeta0;
	\\ 3) the factorial part
	\\ also note that despite the name, it is not actually a log of anything yet
	\\ (but that does happen at the very end of handling logB1 here. Same with logB2)
	logB1=internal_get_logB1(bigR,bigT,b1,b3,d1,logXLB,nUB,dbg);
	logB2=internal_get_logB2(bigS,bigT,b2,b3,d2,logXLB,nUB,dbg);

	logBUB=logB1+logB2;
	logBUB=logBUB-2*(log(polcoef(bigK,1))+logLogX-3/2);
	return(logBUB);
}

\\ here we break the expression for b (defined in the notation section)
\\ into three parts:
\\ 1) B1=b3p*eta0=(b3/d1)*eta0;
\\ 2) B2=b3pp*zeta0=(b3/d2)*zeta0;
\\ 3) the factorial part
\\ also note that despite the name, it is not actually a log of anything yet
\\ (but that does happen at the very end of handling logB1 here. Same with logB2)
\\ 15 Dec 2021
get_logB1(bigK,bigR,bigT,b1,b3,d1,logXLB,nUB,dbg=0)={
	my(logB1);
	
	logB1=internal_get_logB1(bigR,bigT,b1,b3,d1,logXLB,nUB,dbg);
	logB1=logB1-(log(polcoef(bigK,1))+logLogX-3/2);
	return(logB1);
}

\\ for b3'*eta_0 term in expression for b (equation (3.5) on 3 March 2022)
\\ here b3'=b3/d1
\\ 15 Dec 2021
internal_get_logB1(bigR,bigT,b1,b3,d1,logXLB,nUB,dbg=0)={
	my(logB1,logB1a,logB1b,logB1c,logB1d);

	logB1=((bigR-1)*b3+(bigT-1)*b1)/2/d1;
	logB1a=polcoef(polcoef(logB1,1,n),1,logX);
	if(logB1a<0,
		printf("ERROR: in internal_get_logB1(), logB1 has negative lead coefficient: %s\n",logB1);
	);
	logB1b=polcoef(polcoef(logB1,1,n),0,logX);
	logB1b=max(logB1b,0);
	logB1c=polcoef(polcoef(logB1,0,n),1,logX);
	logB1c=max(logB1c,0);
	logB1d=polcoef(polcoef(logB1,0,n),0,logX);
	logB1d=max(logB1d,0);
	if(dbg!=0,
		printf("in internal_get_logB1(): nUB=%9.6e\n",nUB);
		printf("initial b3'*eta_0 value=%s\n",logB1);
		printf("        n*logX coeff=%12.5f, type(n*logX coeff)=%s\n",logB1a,type(logB1a));
		printf("        n      coeff=%12.5f, type(n      coeff)=%s\n",logB1b,type(logB1b));
		printf("        logX   coeff=%12.5f, type(logX   coeff)=%s\n",logB1c,type(logB1c));
		printf("        cnst   coeff=%12.5f, type(cnst   coeff)=%s\n",logB1d,type(logB1d));
	);
	
	logB1=subst(logB1a*n*logX+logB1b*n+logB1c*logX+logB1d,n,nUB);
	logB1a=polcoef(logB1,1,logX);
	logB1b=polcoef(logB1,0,logX);
	logB1b=max(logB1b,0);
	if(dbg!=0,
		printf("after substituting for n: b3'*eta_0 value=%s\n",logB1);
		printf("        logX term  =%9.6e\n",logB1a);
		printf("        const term =%9.6e\n",logB1b);
	);
	logB1=log(logB1a+logB1b/logXLB)+logLogX;
	if(dbg!=0,
		printf("actual b3'*eta_0 value=%s\n",logB1);
	);
	return(logB1);
}

\\ here we break the expression for b (defined in the notation section)
\\ into three parts:
\\ 1) B1=b3p*eta0=(b3/d1)*eta0;
\\ 2) B2=b3pp*zeta0=(b3/d2)*zeta0;
\\ 3) the factorial part
\\ also note that despite the name, it is not actually a log of anything yet
\\ (but that does happen at the very end of handling logB1 here. Same with logB2)
\\ 15 Dec 2021
get_logB2(bigM,bigS,bigT,b2,b3,d2,logXLB,nUB,dbg=0)={
	my(logB2);
	
	logB2=internal_get_logB2(bigS,bigT,b2,b3,d2,logXLB,nUB,dbg);
	logB2=logB2-(log(polcoef(bigM,1))+logLogX-3/2);
	return(logB2);
}

\\ for b3''*zeta_0 term in expression for b (equation (3.5) on 3 March 2022)
\\ 15 Dec 2021
internal_get_logB2(bigS,bigT,b2,b3,d2,logXLB,nUB,dbg=0)={
	my(logB2,logB2a,logB2b,logB2c,logB2d);
	
	logB2=((bigS-1)*b3+(bigT-1)*b2)/2/d2;
	logB2a=polcoef(polcoef(logB2,1,n),1,logX);
	if(logB2a<0,
		printf("ERROR: logB2 has negative lead coefficient: %s\n",logB2);
	);
	logB2b=polcoef(polcoef(logB2,1,n),0,logX);
	logB2b=max(logB2b,0);
	logB2c=polcoef(polcoef(logB2,0,n),1,logX);
	logB2c=max(logB2c,0);
	logB2d=polcoef(polcoef(logB2,0,n),0,logX);
	logB2d=max(logB2d,0);
	if(dbg!=0,
		printf("\ninitial b3''*zeta_0 value=%s\n",logB2);
		printf("        n*logX coeff=%12.5f\n",logB2a);
		printf("        n      coeff=%12.5f\n",logB2b);
		printf("        logX   coeff=%12.5f\n",logB2c);
		printf("        cnst   coeff=%12.5f\n",logB2d);
	);
	
	logB2=subst(logB2a*n*logX+logB2b*n+logB2c*logX+logB2d,n,nUB);
	logB2a=polcoef(logB2,1,logX);
	logB2b=polcoef(logB2,0,logX);
	logB2b=max(logB2b,0);
	if(dbg!=0,
		printf("after substituting for n: b3''*zeta_0 value=%s\n",logB2);
		printf("        logX term  =%9.6e\n",logB2a);
		printf("        const term =%9.6e\n",logB2b);
	);
	logB2=log(logB2a+logB2b/logXLB)+logLogX;
	if(dbg!=0,
		printf("actual b3''*zeta_0 value=%s\n",logB2);
	);
	return(logB2);
}

\\ 0 means bad, 1 means good
\\ 2 July 2022
check_bounds(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB)={
	if(bigLLB<5,
		printf("ERROR: bigLLB=%9.6f must be at least 5\n",bigLLB);
		return(0);
	);
	if(mUB-mLB<0,
		printf("ERROR: mLB=%9.6f must be at most mUB=9.6f\n",mLB,mUB);
		return(0);
	);
	if(mLB<1,
		printf("ERROR: mLB=%9.6f must be at least 1\n",mLB);
		return(0);
	);
	if(rhoUB<rhoLB<0,
		printf("ERROR: rhoLB=%9.6f must be at most rhoUB=9.6f\n",rhoLB,rhoUB);
		return(0);
	);
	if(rhoLB<1,
		printf("ERROR: rhoLB=%9.6f must be at least 2\n",rhoLB);
		return(0);
	);
	if(chiUB-chiLB<0,
		printf("ERROR: chiLB=%9.6f must be at most chiUB=9.6f\n",chiLB,chiUB);
		return(0);
	);
	if(chiLB<0.001,
		printf("ERROR: chiLB=%9.6f must be positve\n",chiLB);
		return(0);
	);
	return(1);
}
