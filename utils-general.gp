\\ \r lfl3\lfl-utils-general.gp

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

\\ used for checking values of the lower and upper bounds on the parameters
\\ used the "examples" code
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
	if(rhoLB<2,
		printf("ERROR: rhoLB=%9.6f must be at least 2\n",rhoLB);
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
	return(1);
}

\\ 12 Nov 2022
get_step(vLB,vUB)={
	my(vStep);
	
	if(vLB>vUB,
		print("ERROR in get_step(), vLB=%9.6f must be at most vUB=%9.6f\n");
		1/0;
	);
	vStep=(vUB-vLB)/20.0;
	if(vStep==0,vStep=0.000001);
	return(vStep);
}
