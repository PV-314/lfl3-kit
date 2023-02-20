\\ \r lfl3\utils-general.gp

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
	if(dbg>0,
		\\print("a=",a,", type(a)=",type(a),", b=",b,", type(b)=",type(b),", c=",c);
		printf("get_solnUB(): a=%9.6f, b=%9.6f, c=%9.6f\n",a,b,c);
	);
	logC=log(c);
	xHUB=c*logC+logC/(logC-1)*(a^(1/h)+c*log(logC));
	if(dbg>0,
		printf("get_solnUB(): a=%9.6f, b=%9.6f, c=%9.6f, xHUB=%9.6f\n",a,b,c,xHUB);
	);
	\\return(xHUB^h*rescalingFactor);
	return(xHUB^h);
}

\\ used for checking values of the lower and upper bounds on the parameters
\\ used the "examples" code
\\ 0 means bad, 1 means good
\\ 2 July 2022
check_bounds(bigLLB,bigLUB,mLB,mUB,rho3LB,rho3UB,chiLB,chiUB,rho2LB,rho2UB)={
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
	if(rho3UB<rho3LB,
		printf("ERROR: rho3LB=%9.6f must be at most rho3UB=9.6f\n",rho3LB,rho3UB);
		return(0);
	);
	if(rho3LB<2,
		printf("ERROR: rho3LB=%9.6f must be at least 2\n",rho3LB);
		return(0);
	);
	if(chiUB-chiLB<0,
		printf("ERROR: chiLB=%9.6f must be at most chiUB=9.6f\n",chiLB,chiUB);
		return(0);
	);
	if(chiLB<0.0001,
		printf("ERROR: chiLB=%9.6f must be positve\n",chiLB);
		return(0);
	);
	if(rho2UB<rho2LB,
		printf("ERROR: rho2LB=%9.6f must be at most rho2UB=9.6f\n",rho2LB,rho2UB);
		return(0);
	);
	if(rho2LB<1.000001,
		printf("ERROR: rho2LB=%9.6f must be larger than 1.000001\n",rho2LB);
		return(0);
	);
	return(1);
}

\\ 12 Nov 2022
get_step(vLB,vUB)={
	my(vStep);
	
	if(vLB>vUB,
		printf("ERROR in get_step(), vLB=%9.6f must be at most vUB=%9.6f\n");
		1/0;
	);
	vStep=(vUB-vLB)/20.0;
	if(vStep==0,vStep=0.000001);
	return(vStep);
}

\\ the following functions are a bit overkill, but want to reduce duplicated code
\\ if two are real, then all three are real and so they are multiplicatively independent (our assumption at the start of Section 3)
\\ 15 Jan 2023
are_multiplicatively_independent(al1,al2)={
	my(areMultIndep);
	
	areMultIndep=0;
	if(type(al1)=="t_REAL" && type(al2)=="t_REAL",
		areMultIndep=1;
	);
	return(areMultIndep);
}

\\ this is correct from our assumption at the start of Section 3
\\ 15 Jan 2023
is_complex(al1,al2,al3)={
	my(isComplex);
	isComplex = type(al1)=="t_COMPLEX" || type(al2)=="t_COMPLEX" || type(al3)=="t_COMPLEX";
	return(isComplex);
}

\\ 15 Jan 2023
get_c1(a,bigL,chi,m)={
	my(c1,two12);
	
	two13=1.2599210498948731647672106072782283506; \\ 2^(1/3) to save recalculation every time
	c1=two13;
	c1=max(c1,(chi*m*bigL)^(2/3));
	c1=max(c1,sqrt(2*m*bigL/a));
	return(c1);
}

\\ 15 Jan 2023
get_c2(areMultIndep,a,bigL,m)={
	my(c2);
	
	if(areMultIndep,
		c2=(m*bigL)^(2/3);
	);
	if(!areMultIndep,
		c2=max((m*bigL)^(2/3),sqrt(m/a)*bigL);
	);
	return(c2);
}