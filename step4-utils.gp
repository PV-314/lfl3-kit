\\ \r lfl3\step4-utils.gp

\\ read("lfl3\\lfl-utils-general.gp");

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