INSTALLATION
============
(1) have Pari installed
(2) in the directory where Pari is installed, create a directory called lfl3
(3) put all the files in https://github.com/PV-314/lfl3-kit in a directory called lfl3

CONTENTS
========
(1) lfl-utils-alpha1Variable.gp
(2) lfl-utils-alpha3Variable.gp
(3) param-search-bms-annals-eg.gp
(4) param-search-bms-compositio-eg.gp
(5) param-search-mignotte-eg1.gp
(6) param-search-mignotte-eg2.gp
(7) lfl-utils-general.gp
(8) step2-utils.gp
(9) step3-utils.gp

Files (3)-(6) above are for running the kit against four separate examples.
param-search-bms-annals-eg.gp is for the linear form in the Annals paper of Bugeaud, Mignotte and Siksek.
param-search-bms-compositio-eg.gp is for the linear form with D=7 in Section 15 of the Compositio paper of Bugeaud, Mignotte and Siksek.
param-search-mignotte-eg1.gp is for Example 1 in Mignotte's original kit paper
param-search-mignotte-eg2.gp is for Example 2 in Mignotte's original kit paper

Each param-search-* file contains a collection of
egX_search_itY() functions where "it" stands for iteration.
	So egX_search_it1() is used for the first iteration of the kit and so on.
egX_check_itY() functions which typically contain the best example found for iteration Y using the egX_search_itY() function.

The code assumes that \Lambda is of the form
b_1 \log(\alpha_1) + b_2 \log(\alpha_2) - b_3 \log(\alpha_3)
where b_1, b_2 and b_3 are positive rational numbers

There is no point in using the kit unless at least one of the alpha_i's is not constant
(use Matveev + LLL instead!).
So we can assume that either alpha_1 is variable or alpha_3 is variable.

lfl-utils-alpha1Variable.gp and lfl-utils-alpha3Variable.gp, respectively contain
the code for each of these two cases.

Finally, lfl-utils-general.gp, step2-utils.gp and step4-utils.gp
contain common code used by both
lfl-utils-alpha1Variable.gp and lfl-utils-alpha3Variable.gp.

HOW TO USE
==========
We recommend that you create a copy of one of the
param-search-egX.gp files for the particular linear form that is of interest to you.

Create/copy a egX_search_general() function. This will contain all the parameters
for the linear form you are considering, as well as the code to call the functions
that do the actual calculations.

The variables that need to be assigned values are:
	bigD,matveevChi,d

	al1,hgtA1,absLogA1,
	al2,hgtA2,absLogA2,
	al3,hgtA3,absLogA3,

	b1,b2,b3
	
	logXLB,nLB
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1,lamUB0
	nUBInit

We assume that the quantity you are trying to bound from above with the kit is called n
We also assume that x is the non-constant alpha_i - this is important and used throughout.
Once again, just follow the examples.

Tip: in the egX_search_general() functions. Do not forget that a1, a2 and a3 need to be set
inside the for loop for rho, not before that. This is because a1, a2 and a3 will depend on rho.

Start by writing and running egX_search_it1().
All it requires are lower and upper bounds for bigL, chi, m and rho,
along with calling the function:
egX_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,,dbg);

It will produce output like the following, giving values of the parameters used and some other values.
The key value is nUB, where UB=upper bound, so this is giving you the upper bound on n (see above)
for those values of the parameters.
nUB will be decreasing, so the last value displayed will be the smallest value found for that iteration.
L=  70, m= 11.0000, rho=  5.0000, chi=2.0000, K=  164155.644*logX, nonDegen log|Lambda|>-1.849388 e7*logX, nonDegenNUB=1.849388 e7, degenNUB1=   0.e-19, degenNUB2=4.047859 e7, degenNUB3=1.816552 e7, nUB=1.849388 e7, transB-b1
L=  71, m=  6.0000, rho=  6.0000, chi=1.9000, K=  142874.437*logX, nonDegen log|Lambda|>-1.817576 e7*logX, nonDegenNUB=1.817576 e7, degenNUB1=   0.e-19, degenNUB2=3.845822 e7, degenNUB3=1.747862 e7, nUB=1.817576 e7
L=  76, m=  5.0000, rho=  6.0000, chi=1.8000, K=  127446.681*logX, nonDegen log|Lambda|>-1.735489 e7*logX, nonDegenNUB=1.735489 e7, degenNUB1=   0.e-19, degenNUB2=3.844486 e7, degenNUB3=1.747243 e7, nUB=1.747243 e7
L=  85, m=  7.0000, rho=  5.0000, chi=1.9000, K=  126847.543*logX, nonDegen log|Lambda|>-1.735303 e7*logX, nonDegenNUB=1.735303 e7, degenNUB1=   0.e-19, degenNUB2=3.543732 e7, degenNUB3=1.616577 e7, nUB=1.735303 e7
L=  91, m=  6.0000, rho=  5.0000, chi=1.8000, K=  116401.275*logX, nonDegen log|Lambda|>-1.704800 e7*logX, nonDegenNUB=1.704800 e7, degenNUB1=   0.e-19, degenNUB2=3.719930 e7, degenNUB3=1.674490 e7, nUB=1.704800 e7

With that best value for the iteration, use the values of the parameters to make
a egX_check_it1() function that checks this particular choice of the parameters.
Again, the easiest way is to copy an existing egX_check_it1() function from one of the
example files provided in github and just change the variables in that function, as required.

Here set the lower and upper bounds for bigL, chi, m and rho to both be equal to the
values of bigL, chi, m and rho that you want to check.

This egX_check_it1() function outputs lots more information and can be used for testing,
as well as outputting details that you may want to include in a paper for your description
of how you used the kit (this is helpful info for readers and for referees -- see the
exposition of the examples in Section 6 of the kit paper).

CONTACT
=======
If you have any questions, problems, need some help, find a bug,..., please contact PV
you can find PV's details in the arxiv preprint for the kit.
