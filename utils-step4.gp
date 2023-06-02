\\ \r lfl3\utils-step4.gp

read("lfl3\\utils-general.gp");

\\ here we eliminate the b_1*log(alpha1) term (hence function name)
\\ b1*log(a1)+b2*log(a2)-b3*log(a3) and u1*b1+u2*b2+u3*b3=0, so
\\ u1*Lambda
\\ =u1*b1*log(a1)+u1*b2*log(a2)-u1*b3*log(a3)
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
\\ it returns the array [hgtNewA1,logNewA1,hgtNewA2,logNewA2,hgtNewA1_0,logNewA1_0,hgtNewA2_0,logNewA2_0]
\\ A_1 and A_2 are for when u_1 is not zero, while A_1_0 and A_2_0 are for when u_1 is zero
\\ we order A_1 and A_2 (similarly A_1_0 and A_2_0)
\\ such that b_i*log(alpha_1)-b_j*log(alpha_2) with i<j
\\
\\ 3 Jan 2022
get_logA1A2_from_b1(d,absLogA1,hgtA1,absLogA2,hgtA2,absLogA3,hgtA3,isComplex,logXLB,u1UB,u2UB,u3UB,dbg=0)={
	my(hgtNewA1,hgtNewA1_0,hgtNewA2,hgtNewA2_0,logNewA1,logNewA1_0,logNewA2,logNewA2_0);
	
	\\ "New" as it is for A1 and A2 in linear form in two logs
	hgtNewA1=u1UB*hgtA2+u2UB*hgtA1;
	if(poldegree(hgtNewA1)==1,
		hgtNewA1=subst(hgtNewA1,logX,logXLB)/logXLB*logX;
	);
	logNewA1=u1UB*absLogA2+u2UB*absLogA1;
	logNewA1=update_logAi(logNewA1,isComplex,logXLB);
	if(dbg>0,
		\\printf("\nget_logA1A2_from_b1(): hgtA1=%s\n",hgtA1);
		\\printf("get_logA1A2_from_b1(): hgtA2=%s\n",hgtA2);
		\\print("d=",d,", absLogA1=",absLogA1,", hgtA1=",hgtA1);
		\\print("absLogA2=",absLogA2,", hgtA2=",hgtA2);
		\\print("absLogA3=",absLogA3,", hgtA3=",hgtA3);
		\\print("u1UB=",u1UB,", u2UB=",u2UB,", u3UB=",u3UB);
		printf("get_logA1A2_from_b1(): hgtNewA1=%s\n",hgtNewA1);
		printf("get_logA1A2_from_b1(): logNewA1=%s\n",logNewA1);
	);

	hgtNewA2=u3UB*hgtA1+u1UB*hgtA3;
	if(poldegree(hgtNewA2)==1,
		hgtNewA2=subst(hgtNewA2,logX,logXLB)/logXLB*logX;
	);
	logNewA2=u3UB*absLogA1+u1UB*absLogA3;
	logNewA2=update_logAi(logNewA2,isComplex,logXLB);

	if(dbg>0,
		printf("get_logA1A2_from_b1(): hgtNewA2=%s\n",hgtNewA2);
		printf("get_logA1A2_from_b1(): logNewA2=%s\n",logNewA2);
	);

	\\ get values if u_1=0
	if(type(u2UB)!="t_POL",
		\\ linear form is: b1*log(a1^u2)-b3*log(a2^u3*a3^u2)
		if(dbg>0,
			print("get_logA1A2_from_b1(): considering u_1=0 and eliminating b_2");
		);
		hgtNewA1_0=u2UB*hgtA1;
		logNewA1_0=u2UB*absLogA1;
		logNewA1_0=update_logAi(logNewA1_0,isComplex,logXLB);
		hgtNewA2_0=u3UB*hgtA2+u2UB*hgtA3;
		logNewA2_0=u3UB*absLogA2+u2UB*absLogA3;
		logNewA2_0=update_logAi(logNewA2_0,isComplex,logXLB);
	);
	if(type(u2UB)=="t_POL" && type(u3UB)!="t_POL",
		\\ linear form is: b1*log(a1^u3)-b2*log(a2^(-u3)*a3^(-u2))

		if(dbg>0,
			print("get_logA1A2_from_b1(): considering u_1=0 and eliminating b_3");
		);
		hgtNewA1_0=u3UB*hgtA1;
		logNewA1_0=u3UB*absLogA1;
		logNewA1_0=update_logAi(logNewA1_0,isComplex,logXLB);
		hgtNewA2_0=u3UB*hgtA2+u2UB*hgtA3;
		logNewA2_0=u3UB*absLogA2+u2UB*absLogA3;
		logNewA2_0=update_logAi(logNewA2_0,isComplex,logXLB);
	);
	if(type(u2UB)=="t_POL" && type(u3UB)=="t_POL",
		print("ERROR in get_logA1A2_from_b1(): u2=",u2UB," and u3=",u3UB," are both polynomials");
		return([]);
	);
	if(dbg>0,
		printf("get_logA1A2_from_b1(): hgtNewA1_0=%s\n",hgtNewA1_0);
		printf("get_logA1A2_from_b1(): logNewA1_0=%s\n",logNewA1_0);
		printf("get_logA1A2_from_b1(): hgtNewA2_0=%s\n",hgtNewA2_0);
		printf("get_logA1A2_from_b1(): logNewA2_0=%s\n",logNewA2_0);
	);
	return([hgtNewA1,logNewA1,hgtNewA2,logNewA2,hgtNewA1_0,logNewA1_0,hgtNewA2_0,logNewA2_0]);
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
\\ it returns the array [hgtNewA1,logNewA1,hgtNewA2,logNewA2,hgtNewA1_0,logNewA1_0,hgtNewA2_0,logNewA2_0]
\\ A_1 and A_2 are for when u_2 is not zero, while A_1_0 and A_2_0 are for when u_2 is zero
\\ we order A_1 and A_2 (similarly A_1_0 and A_2_0)
\\ such that b_i*log(alpha_1)-b_j*log(alpha_2) with i<j
\\
\\ 3 Jan 2022
get_logA1A2_from_b2(d,absLogA1,hgtA1,absLogA2,hgtA2,absLogA3,hgtA3,isComplex,logXLB,u1UB,u2UB,u3UB,dbg=0)={
	my(hgtNewA1,hgtNewA1_0,hgtNewA2,hgtNewA2_0,logNewA1,logNewA1_0,logNewA2,logNewA2_0);
	
	\\ "New" as it is for A1 and A2 in linear form in two logs
	hgtNewA1=u1UB*hgtA2+u2UB*hgtA1;
	logNewA1=u1UB*absLogA2+u2UB*absLogA1;
	logNewA1=update_logAi(logNewA1,isComplex,logXLB);
	if(dbg>0,
		printf("get_logA1A2_from_b2(): hgtNewA1=%s\n",hgtNewA1);
		printf("get_logA1A2_from_b2(): logNewA1=%s\n",logNewA1);
	);

	hgtNewA2=u3UB*hgtA2+u2UB*hgtA3;
	logNewA2=u3UB*absLogA2+u2UB*absLogA3;
	logNewA2=update_logAi(logNewA2,isComplex,logXLB);
	if(dbg>0,
		printf("get_logA1A2_from_b2(): hgtNewA2=%s\n",hgtNewA2);
		printf("get_logA1A2_from_b2(): logNewA2=%s\n",logNewA2);
	);
	
	\\ get values if u_2=0
	if(type(u3UB)!="t_POL",
		\\ linear form is: b1*log(a1^u3*a3^u1)-b2*log(a2^(-u3))
		if(dbg>0,
			print("get_logA1A2_from_b2(): considering u_2=0 and eliminating b_3");
			print("u3UB=",u3UB,", absLogA2=",absLogA2,", hgtA2=",hgtA2);
		);
		hgtNewA1_0=u3UB*hgtA1+u1UB*hgtA3;
		logNewA1_0=u3UB*absLogA1+u1UB*absLogA3;
		logNewA1_0=update_logAi(logNewA1_0,isComplex,logXLB);
		hgtNewA2_0=u3UB*hgtA2;
		logNewA2_0=u3UB*absLogA2;
		logNewA2_0=update_logAi(logNewA2_0,isComplex,logXLB);
	);
	if(type(u3UB)=="t_POL" && type(u1UB)!="t_POL",
		\\ linear form is: b2*log(a2^u1)-b3*log(a1^u3*a3^u1)
		if(dbg>0,
			print("get_logA1A2_from_b2(): considering u_2=0 and eliminating b_1");
		);
		hgtNewA1_0=u1UB*hgtA2;
		logNewA1_0=u1UB*absLogA2;
		logNewA1_0=update_logAi(logNewA1_0,isComplex,logXLB);
		hgtNewA2_0=u3UB*hgtA1+u1UB*hgtA3;
		logNewA2_0=u3UB*absLogA1+u1UB*absLogA3;
		logNewA2_0=update_logAi(logNewA2_0,isComplex,logXLB);
	);
	if(type(u1UB)=="t_POL" && type(u3UB)=="t_POL",
		print("ERROR in get_logA1A2_from_b2(): u1=",u1UB," and u3=",u3UB," are both polynomials");
		return([]);
	);

	if(dbg>0,
		printf("get_logA1A2_from_b2(): hgtNewA1_0=%s\n",hgtNewA1_0);
		printf("get_logA1A2_from_b2(): logNewA1_0=%s\n",logNewA1_0);
		printf("get_logA1A2_from_b2(): hgtNewA2_0=%s\n",hgtNewA2_0);
		printf("get_logA1A2_from_b2(): logNewA2_0=%s\n",logNewA2_0);
	);
	return([hgtNewA1,logNewA1,hgtNewA2,logNewA2,hgtNewA1_0,logNewA1_0,hgtNewA2_0,logNewA2_0]);
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
\\ it returns the array [hgtNewA1,logNewA1,hgtNewA2,logNewA2,hgtNewA1_0,logNewA1_0,hgtNewA2_0,logNewA2_0]
\\ A_1 and A_2 are for when u_3 is not zero, while A_1_0 and A_2_0 are for when u_3 is zero
\\ we order A_1 and A_2 (similarly A_1_0 and A_2_0)
\\ such that b_i*log(alpha_1)-b_j*log(alpha_2) with i<j
\\
\\ 3 Jan 2022
get_logA1A2_from_b3(d,absLogA1,hgtA1,absLogA2,hgtA2,absLogA3,hgtA3,isComplex,logXLB,u1UB,u2UB,u3UB,dbg=0)={
	my(hgtNewA1,hgtNewA1_0,hgtNewA2,hgtNewA2_0,logNewA1,logNewA1_0,logNewA2,logNewA2_0);
	
	\\ "New" as it is for A1 and A2 in linear form in two logs
	\\ linear form = b1*log(a1^u3*a3^u1)-b2*log(a2^(-u3)*a3^(-u2))
	hgtNewA1=u3UB*hgtA1+u1UB*hgtA3;
	logNewA1=u3UB*absLogA1+u1UB*absLogA3;
	logNewA1=update_logAi(logNewA1,isComplex,logXLB);
	if(dbg>0,
		printf("get_logA1A2_from_b3(): hgtNewA1=%s\n",hgtNewA1);
		printf("get_logA1A2_from_b3(): logNewA1=%s\n",logNewA1);
	);

	hgtNewA2=u3UB*hgtA2+u2UB*hgtA3;
	logNewA2=u3UB*absLogA2+u2UB*absLogA3;
	logNewA2=update_logAi(logNewA2,isComplex,logXLB);
	if(dbg>0,
		printf("get_logA1A2_from_b3(): hgtNewA2=%s\n",hgtNewA2);
		printf("get_logA1A2_from_b3(): logNewA2=%s\n",logNewA2);
	);

	\\ get values if u_3=0
	if(type(u1UB)!="t_POL",
		\\ linear form is: b2*log(a1^(-u2)*a2^u1)-b3*log(a3^(-u1))
		if(dbg>0,
			print("get_logA1A2_from_b3(): considering u_3=0 and eliminating b_1");
		);
		hgtNewA1_0=u2UB*hgtA1+u1UB*hgtA2;
		logNewA1_0=u2UB*absLogA1+u1UB*absLogA2;
		logNewA1_0=update_logAi(logNewA1_0,isComplex,logXLB);
		hgtNewA2_0=u1UB*hgtA3;
		logNewA2_0=u1UB*absLogA1;
		logNewA2_0=update_logAi(logNewA2_0,isComplex,logXLB);
	);
	if(type(u1UB)=="t_POL" && type(u2UB)!="t_POL",
		\\ linear form is: b1*log(a1^u2*a2^(-u1))-b3*log(a3^u2)
		if(dbg>0,
			print("get_logA1A2_from_b3(): considering u_3=0 and eliminating b_2");
		);
		hgtNewA1_0=u2UB*hgtA1+u1UB*hgtA2;
		logNewA1_0=u2UB*absLogA1+u1UB*absLogA2;
		logNewA1_0=update_logAi(logNewA1_0,isComplex,logXLB);
		hgtNewA2_0=u2UB*hgtA3;
		logNewA2_0=u2UB*absLogA3;
		logNewA2_0=update_logAi(logNewA2_0,isComplex,logXLB);
	);
	if(type(u1UB)=="t_POL" && type(u2UB)=="t_POL",
		print("ERROR in get_logA1A2_from_b3(): u1=",u1UB," and u2=",u2UB," are both polynomials");
		return([]);
	);

	if(dbg>0,
		printf("get_logA1A2_from_b3(): hgtNewA1_0=%s\n",hgtNewA1_0);
		printf("get_logA1A2_from_b3(): logNewA1_0=%s\n",logNewA1_0);
		printf("get_logA1A2_from_b3(): hgtNewA2_0=%s\n",hgtNewA2_0);
		printf("get_logA1A2_from_b3(): logNewA2_0=%s\n",logNewA2_0);
	);
	return([hgtNewA1,logNewA1,hgtNewA2,logNewA2,hgtNewA1_0,logNewA1_0,hgtNewA2_0,logNewA2_0]);
}

\\ lamMul is the factor we have multiplied lambda by to get a linear form in 2 logs
\\ it will be an upper bound for one of the u_i's in the linear relationship between the b_i's
\\ we use Corollary 1 or 2 of Laurent's 2008 paper here
\\ logAsArray is an array of [logA1,logA2] in the notation on page 328 of Laurent's 2008 paper,
\\ upper bounds for logA1 and logA2 in fact.
\\ using LMN, returns [rho2Log,max(nUBMultDep,nUBMultIndep)]
\\ using Laurent 2008, returns [max(nUBMultDep,nUBMultIndep),rho2Log,muLog]
\\ 3 Jan 2022
get_degen_nUB(hgtAlpha1,absLogAlpha1,hgtAlpha2,absLogAlpha2,d,nLB,lamUB0,lamUB1,lamMul,rho2LB,rho2UB,muLB,muUB,logXLB,isComplex,dbg=0)={
	my(muMultIndep,nUBMultDep,nUBMultIndep,resultMultIndep,rhoMultIndep);

	if(dbg>0,
		print("get_degen_nUB(): START");
	);
	
	if((type(hgtAlpha1)=="t_POL" || type(absLogAlpha1)=="t_POL") && (type(hgtAlpha2)=="t_POL" || type(absLogAlpha2)=="t_POL"),
		print("ERROR in get_degen_nUB(): logA1 and logA2 are both polynomials");
		return([]);
	);

	\\ returns [nUBMin,rhoMin,muMin]
	resultMultIndep=calc_degen_mult_indep_nUB(absLogAlpha1,hgtAlpha1,absLogAlpha2,hgtAlpha2,isComplex,d,lamMul,lamUB1,lamUB0,nLB,rho2LB,rho2UB,muLB,muUB,logXLB,dbg);
	nUBMultIndep=resultMultIndep[1];
	rhoMultIndep=resultMultIndep[2];
	muMultIndep=resultMultIndep[3];

	nUBMultDep=0;
	if(isComplex && (abs(subst(hgtAlpha1,logX,logXLB))<0.00001 || abs(subst(hgtAlpha2,logX,logXLB))<0.00001),
		nUBMultDep=calc_degen_mult_dep_nUB(absLogAlpha1,hgtAlpha1,absLogAlpha2,hgtAlpha2,d,lamMul,lamUB1,lamUB0,nLB,logXLB,dbg);
	);
	if(dbg>0,
		printf("get_degen_nUB(): nUBMultDep  =%9.6e\n",nUBMultDep);
		printf("get_degen_nUB(): nUBMultIndep=%9.6e\n",nUBMultIndep);
	);
	if(dbg>0,
		print("get_degen_nUB(): END");
	);

	return([max(nUBMultDep,nUBMultIndep),rhoMultIndep,muMultIndep]);
}

\\ uses equation (28) of Laurent 2008 (based on Theorem 2)
\\ returns [rhoMin,nUBMin,muMin]
\\ 4 Feb 2023
calc_degen_mult_indep_nUB(absLogAlpha1,hgtAlpha1,absLogAlpha2,hgtAlpha2,isComplex,bigD,lamMul,lamUB1,lamUB0,nLB,rho2LB,rho2UB,muLB,muUB,logXLB,dbg=0)={
	my(a1,a1PlusA2UB,a2,bCnst,bigLambdaLB,hLB,hTerm1LB,hTerm1UB,hTerm2UB,hTerm3sUB,hTerm4a,hTerm4b,hUB,lambda,logLambda,mu,muMin,muStep,mySigma,nUB,nUBMin,rho2Step,rhoMin,sqrtA1A2UB);

	nUBMin=10^100;
	muMin=0;
	rhoMin=0;
	if(muLB<0.001 || muUB<1/3, \\ mu \geq 1/3 by Theorem 2 of Laurent 2008
		muLB=0.54;
		muUB=0.74;
	);
	muStep=get_step(muLB,muUB);
	rho2Step=get_step(rho2LB,rho2UB);
	forstep(mu=muLB,muUB,muStep,
	if(dbg>1,
		printf("calc_degen_mult_indep_nUB(): mu=%9.6f\n",mu);
	);
	forstep(rho=rho2LB,rho2UB,rho2Step,
		\\ used this to scrape values for test
		\\if(rho==3.1,
		\\	dbg=1;
		\\	print("absLogAlpha1=",absLogAlpha1,"\nhgtAlpha1=",hgtAlpha1,"\nabsLogAlpha2=",absLogAlpha2,"\n",hgtAlpha2,"\n",isComplex,"\n",d,"\n",lamMul,"\n",lamUB1,"\n",lamUB0,"\nnLB=",nLB,"\nlogXLB=",logXLB,"\n");
		\\);
		mySigma=(1+2*mu-mu*mu)/2;
		lambda=mySigma*log(rho);
		logLambda=log(lambda);
		
		aiLB=1;
		a1=rho*absLogAlpha1+2*bigD*hgtAlpha1;
		a2=rho*absLogAlpha2+2*bigD*hgtAlpha2;
		if(!isComplex,
			a1=(rho-1)*absLogAlpha1+2*bigD*hgtAlpha1;
			a2=(rho-1)*absLogAlpha2+2*bigD*hgtAlpha2;
		);
		if(type(a1)=="t_REAL",
			a1=max(aiLB,a1);
		);
		if(type(a2)=="t_REAL",
			a2=max(aiLB,a2);
		);
		if(subst(a1*a2,logX,logXLB)<lambda*lambda,
			printf("BAD: a1*a2=%9.6f<lambda^2=%9.6f\n",a1*a2,lambda*lambda);
			print("rho=",rho,", bigD=",bigD,", absLogAlpha1=",absLogAlpha1,", hgtAlpha1=",hgtAlpha1);
			print("rho=",rho,", bigD=",bigD,", absLogAlpha2=",absLogAlpha2,", hgtAlpha2=",hgtAlpha2);
			error();
		);
		a1LB=subst(a1,logX,logXLB);
		a2LB=subst(a2,logX,logXLB);

		\\ bCnstLB*n<b1/a2+b2/a1
		bCnst=1/a1+1/a2;
		bCnstUB=subst(bCnst,logX,logXLB);
		bCnstLB=subst(bCnst,logX,10^1000);
		
		if(dbg>1,
			printf("\ncalc_degen_mult_indep_nUB(): mu=%9.6f, rho=%9.6f, sigma=%9.6f, lambda=%9.6f\n",mu,rho,mySigma,lambda);
			printf("a1=%s, a2=%s\n",a1,a2);
			printf("calc_degen_mult_indep_nUB(): %9.6f*n<(b1/a2+b2/a1)<%9.6f*n\n",bCnstLB,bCnstUB);
		);

		\\ we will have h=d*(log(n)+c), for some absolute constant c
		\\ 0.693147=log(2) as in Laurent's equaion (3)
		hLBInit=max(bigD*0.693147/2,lambda);
		hUB=bigD*(log(bCnstUB)+logN+logLambda+1.75)+0.06;
		hUB=max(0,hLBInit-subst(hUB,logN,log(nLB)))+hUB;
		
		hLB=bigD*(log(bCnstLB)+logN+logLambda+1.75)+0.06;
		hLB=max(0,hLBInit-subst(hLB,logN,log(nLB)))+hLB;

		if(dbg>1,
			printf("calc_degen_mult_indep_nUB(): hLB=%s\n",hLB);
			printf("calc_degen_mult_indep_nUB(): hUB=%s\n",hUB);
		);

		\\ for calculating C'' (=bigCPP here) use lower bounds for a1, a2 and h now (as Laurent shows after his eqn (28))
		bigH=hLB/lambda+1/mySigma;
		bigH=subst(bigH,logN,log(nLB));
		\\print("bigH=",bigH);
		myOmega=2*(1+sqrt(1+1/4/bigH/bigH));
		\\print("myOmega=",myOmega);
		myTheta=sqrt(1+1/4/bigH/bigH)+1/2/bigH;
		\\print("myTheta=",myTheta);

		bigC=myOmega/6+1/2*sqrt(myOmega*myOmega/9);
		bigC=mu/lambda/lambda/lambda/mySigma*bigC*bigC;
		bigC=subst(bigC,logX,logXLB);
		bigC=subst(bigC,n,nLB);
		\\print("bigC=",bigC);

		bigCP=sqrt(bigC*mySigma*myOmega*myTheta/lambda/lambda/lambda/mu);
		bigCP=subst(bigCP,logX,logXLB);
		bigCP=subst(bigCP,logN,log(nLB));
		\\print("bigCP=",bigCP);

		bigCP2=(1+lambda/hLB/mySigma)*(1+lambda/hLB/mySigma);
		bigCP2=subst(bigCP2,logX,logXLB);
		bigCP2=subst(bigCP2,logN,log(nLB));
		\\print("bigCP2=",bigCP2);
		tmp1=hLB*myTheta;
		tmp1=subst(tmp1,logX,logXLB);
		tmp1=subst(tmp1,logN,log(nLB));
		\\print("tmp1=",tmp1);
		tmp2=bigCP*(hLB+lambda/mySigma)*(hLB+lambda/mySigma)*a1LB*a2LB;
		tmp2=subst(tmp2,logX,logXLB);
		tmp2=subst(tmp2,logN,log(nLB));
		\\print("tmp2=",tmp2);
		bigCP2=bigCP2*(bigC+sqrt(tmp1)/(hLB+lambda/mySigma)/a1LB/a2LB+log(tmp2)/(hLB+lambda/mySigma)/(hLB+lambda/mySigma)/a1LB/a2LB);
		bigCP2=subst(bigCP2,logX,logXLB);
		bigCP2=subst(bigCP2,logN,log(nLB));
		\\print("bigCP2=",bigCP2);

		bigLambdaLB=-bigCP2*hUB*hUB*a1*a2;
		nUB=calc_degen_nUB(bigLambdaLB,lamMul,lamUB1,lamUB0,logXLB,dbg);

		if(nUB<nUBMin,
			nUBMin=nUB;
			rhoMin=rho;
			muMin=mu;
		);
	);
	);
	if(dbg>0,
		printf("calc_degen_mult_indep_nUB(): returning mu=%9.6f, rho=%9.6f, nUBMin=%9.6e\n",muMin,rhoMin,nUBMin);
	);
	return([nUBMin,rhoMin,muMin]);
}

\\ bigLambdaCnst is from a lower bound from LFL for |u_i*Lambda| of the form bigLambdaCnst*H^2,
\\ where bigLambdaCnst depends on logX only, not on n
\\ uses Theoreme 3 of LMN
\\ 30 July 2022
calc_degen_mult_dep_nUB(absLogAlpha1,hgtAlpha1,absLogAlpha2,hgtAlpha2,d,lamMul,lamUB1,lamUB0,nLB,logXLB,dbg=0)={
	my(a,bCnst,bigH,bigLambdaLBCnst,cf1,hLB,hUB);

	if(abs(subst(hgtAlpha1,logX,logXLB))<0.0001,
		a=10.98*absLogAlpha2+d*hgtAlpha2;
	);
	if(abs(subst(hgtAlpha2,logX,logXLB))<0.0001,
		a=10.98*absLogAlpha1+d*hgtAlpha1;
	);
	a=a+max(0,20-subst(a,logX,logXLB));

	\\ bCnst*n>b1/(2*a)+b2/68.9
	bCnst=1/68.9+1/(2*a);
	bCnst=subst(bCnst,logX,logXLB);
	
	if(dbg>0,
		printf("calc_degen_mult_dep_nUB(): a=%s\n",a);
		printf("calc_degen_mult_dep_nUB(): (b1/(2*a)+b2/68.9)<%9.6f*n\n",bCnst);
	);
	
	\\ we will have h=d*(log(n)+c), for some absolute constant c
	bigH=d*(log(bCnstUB)+logN)+2.35*d+5.03;
	bigHInit=max(17,sqrt(d)/10);
	bigH=bigH+max(0,bigHInit-subst(bigH,logN,log(nLB)));

	bigLambdaCnst=-8.87*a;
	\\ now we normalise bigH so it is monic:
	if(poldegree(bigH,logN)==1,
		cf1=polcoef(bigH,1,logN);
		bigLambdaCnst=bigLambdaCnst*cf1*cf1;
		bigH=bigH/cf1;
	);
	if(poldegree(bigH,logN)!=1,
		print("ERROR in calc_degen_mult_dep_nUB(): bigH=",bigH," must be a linear polynomial in logN");
		error();
	);
	if(dbg>0,
		if(polcoef(bigH,0,logN)<0,
			printf("calc_degen_mult_dep_nUB(): after normalisation, bigH=%9.6f*logN%9.6f\n",polcoef(bigH,1,logN),polcoef(bigH,0,logN));
		);
		if(polcoef(bigH,0,logN)>-0.00001,
			printf("calc_degen_mult_dep_nUB(): after normalisation, bigH=%9.6f*logN+%9.6f\n",polcoef(bigH,1,logN),polcoef(bigH,0,logN));
		);
	);

	\\ make bigLambdaCnst=c*logX
	if(poldegree(bigLambdaCnst,logX)==1,
		\\print("logXLB=",logXLB);
		bigLambdaCnst=subst(bigLambdaCnst,logX,logXLB)/logXLB*logX;
	);
	if(polcoef(bigLambdaCnst,0,logX)>0,
		print("BAD: constant term of bigLambdaCnst is positive, bigLambdaCnst=",bigLambdaCnst);
	);

	bigLambdaLB=bigLambdaCnst*bigH*bigH;
	nUB=calc_degen_nUB(bigLambdaLB,lamMul,lamUB1,lamUB0,logXLB,dbg);
	return(nUB);
}

\\ common code for both mult dep and mult indep cases
\\ 18 Feb 2023
calc_degen_nUB(bigLambdaLB,lamMul,lamUB1,lamUB0,logXLB,dbg=0)={
	my(a,b,c1,c2,cfX1,nUB);
	
	\\ change it from lower bound for |u_i*Lambda| to a lower bound for |Lambda|
	bigLambdaLB=bigLambdaLB-log(lamMul);
	if(dbg>0,
		print("calc_degen_nUB(): after subtracting log(lamMul), bigLambdaLB=",bigLambdaLB);
		printf("lamMul=%9.6f\n",lamMul);
	);

	\\ now write log |Lambda|>f(logN)*logX
	bigLambdaLB=polcoef(bigLambdaLB,1,logX)*logX+polcoef(bigLambdaLB,0,logX)/logXLB*logX;
	if(dbg>0,
		print("calc_degen_nUB(): bigLambdaLB=",bigLambdaLB);
	);

	\\ now write log |Lambda|>c1*(logN+c2)^2*logX
	if(poldegree(bigLambdaLB,logN)==2 && poldegree(bigLambdaLB,logX)==1,
		cfX1=polcoef(bigLambdaLB,1,logX);
		\\print("cfX1=",cfX1);
		c1=polcoef(cfX1,2,logN);
		\\print("c1=",c1);
		c2=polcoef(cfX1/c1,1,logN)/2;
		\\print("c2=",c2);
		if(polcoef(cfX1/c1,0,logN)>0,
			if(c2>0,
				c2=max(c2,sqrt(polcoef(cfX1/c1,0,logN)));
			);
			if(c2<0,
				c2=min(c2,-sqrt(polcoef(cfX1/c1,0,logN)));
			);
		);
		if(dbg>0,
			printf("calc_degen_nUB(): log |Lambda|>%9.6e*(logN+%9.6e)^2*logX\n",c1,c2);
		);
		b=c1*logX*exp(c2)/lamUB1;
		b=polcoef(b,0);
		a=-lamUB0*exp(c2)/lamUB1;
		if(poldegree(a)>0,
			printf("ERROR in calc_degen_nUB(): a=%s should be of non-positive degree\n",a);
		);
		if(poldegree(a)==0,
			a=polcoef(a,0);
		);
		if(poldegree(a)<0,
			a=subst(a,logX,logXLB);
		);
		if(dbg>0,
			print("calc_degen_nUB(): lamUB0=",lamUB0,", lamUB1=",lamUB1);
			printf("calc_degen_nUB(): a=%9.6f, b=%9.6f\n",a,b);
		);
		nUB=get_solnUB(a,b,2,dbg);
		nUB=nUB/exp(c2);
		if(dbg>0,
			printf("calc_degen_nUB(): nUB=%9.6e\n",nUB);
		);
	);
	if(poldegree(bigLambdaLB,logN)!=2 || poldegree(bigLambdaLB,logX)!=1,
		print("ERROR in calc_degen_nUB(): bigLambda must be of degree 2 in logN and degree 1 in logX, bigLambdaLB=",bigLambdaLB);
		error();
	);
	if(dbg>0,
		printf("calc_degen_nUB(): returning nUB=%9.6e\n",nUB);
	);

	return(nUB);
}

\\ \\ step4Result=[rhoB1,nUB1,rhoB2,nUB2,rhoB3,nUB3]
\\ 28 Jan 2023
step4_update_minNUB(step3Result,minNUB,step4Result,dbg=0)={
	my(bigK,bigL,chi,localMinNUB,m,minB,minNDegenUB,minRho2Logs,nDegenUB,nNonDegenUB,nUB,rho3Logs);

	if(dbg>0,
		print("step3Result=",step3Result);
		print("step4Result=",step4Result);
	);
	bigK=step3Result[1];
	bigL=step3Result[2];
	m=step3Result[3];
	rho3Logs=step3Result[4];
	chi=step3Result[5];
	nNonDegenUB=step3Result[11];
	
	localMinNUB=minNUB;
	if(dbg>0,
		printf("step4_update_minNUB(): nNonDegenUB=%9.6e, type=%s\n",nNonDegenUB,type(nNonDegenUB));
		printf("step4_update_minNUB(): nDegenUB=[%9.6e, %9.6f, %9.6f, %9.6e, %9.6f, %9.6f, %9.6e, %9.6f, %9.6f]\n",step4Result[1],step4Result[2],step4Result[3],step4Result[4],step4Result[5],step4Result[6],step4Result[7],,step4Result[8],step4Result[9]);
	);
	minB="b1";
	\\ initialise minNDegenUB
	minNDegenUB=0;
	if(step4Result[1]>0,
		minNDegenUB=step4Result[1];
		minRho2Logs=step4Result[2];
		minMu=step4Result[3];
	);
	if(step4Result[4]>0 && minNDegenUB==0,
		minB="b2";
		minNDegenUB=step4Result[4];
		minRho2Logs=step4Result[5];
		minMu=step4Result[6];
	);
	if(step4Result[7]>0 && minNDegenUB==0,
		minB="b3";
		minNDegenUB=step4Result[7];
		minRho2Logs=step4Result[8];
		minMu=step4Result[9];
	);
	if(step4Result[1]>0 && step4Result[1]<minNDegenUB,
		minB="b1";
		minNDegenUB=step4Result[1];
		minRho2Logs=step4Result[2];
		minMu=step4Result[3];
	);
	if(step4Result[4]>0 && step4Result[4]<minNDegenUB,
		minB="b2";
		minNDegenUB=step4Result[4];
		minRho2Logs=step4Result[5];
		minMu=step4Result[6];
	);
	if(step4Result[7]>0 && step4Result[7]<minNDegenUB,
		minB="b3";
		minNDegenUB=step4Result[7];
		minRho2Logs=step4Result[8];
		minMu=step4Result[9];
	);
	nUB=max(minNDegenUB,nNonDegenUB);
	if(minNDegenUB>0 && nNonDegenUB>0 && nUB<minNUB,
		localMinNUB=nUB;
		printf("L=%4d, m=%8.4f, rho(3logs)=%7.4f, chi=%7.4f, K=%12.3f*logX, nonDegen log|Lambda|>%9.6e*logX, nonDegenNUB=%10.6e, rho(2logs)=%7.4f, mu(2logs)=%7.4f, degenNUB1=%10.6e, degenNUB2=%10.6e, degenNUB3=%10.6e, nUB=%10.6e, eliminate-%s\n",bigL,m,rho3Logs,chi,polcoef(bigK,1),-polcoef(bigK,1)*bigL*log(rho3Logs),nNonDegenUB,minRho2Logs,minMu,step4Result[1],step4Result[4],step4Result[7],nUB,minB);
	);
	return(localMinNUB);
}

\\ 9 Feb 2023
update_logAi(logAi,isComplex,logXLB)={
	my(newLogAi);
	
	newLogAi=logAi;
	if(poldegree(newLogAi)==1,
		newLogAi=subst(newLogAi,logX,logXLB)/logXLB*logX;
	);
	if(isComplex!=0,
		if(type(newLogAi)=="t_POL",
			newLogAi=Pi;
		);
		newLogAi=min(newLogAi,Pi);
	);
	return(newLogAi);
}