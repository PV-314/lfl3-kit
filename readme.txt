INSTALLATION
============
(1) have Pari installed
(2) in the directory where Pari is installed, create a directory called lfl3
(3) put all the files in https://github.com/PV-314/lfl3-kit in the lfl3 directory

CONTENTS
========
 (1) eg-mignotte1.gp
 (2) eg-mignotte2.gp
 (3) eg-bms-annals.gp
 (4) eg-bms-compositio.gp
 (5) kit-alpha1Variable.gp
 (6) kit-alpha3Variable.gp
 (7) utils-general.gp
 (8) utils-step2.gp
 (9) utils-step3.gp
(10) utils-step4.gp

The first four files are for running the kit against four separate examples.
eg-mignotte1.gp is for Example 1 in Mignotte's original kit paper
eg-mignotte2.gp is for Example 2 in Mignotte's original kit paper
eg-bms-annals.gp is for the linear form in the Annals paper of Bugeaud, Mignotte and Siksek.
eg-bms-compositio.gp is for the linear form with D=7 in Section 15 of the Compositio paper of Bugeaud, Mignotte and Siksek.

Each of these files contains a collection of
egX_search_itY() functions where "it" stands for iteration.
	So egX_search_it1() is used for the first iteration of the kit and so on.
egX_check_itY() functions which typically contain the best example found for iteration Y using the egX_search_itY() function.

The code assumes that \Lambda is of the form
b_1 \log(\alpha_1) + b_2 \log(\alpha_2) - b_3 \log(\alpha_3)
where b_1, b_2 and b_3 are positive rational numbers

There is no point in using the kit unless at least one of the alpha_i's is not constant
(use Matveev + LLL instead!).
So we can assume that either alpha_1 is variable or alpha_3 is variable.

kit-alpha1Variable.gp and kit-alpha3Variable.gp, respectively contain
the functions for users for each of these two cases.

Finally, utils-general.gp, utils-step2.gp, utils-step3.gp and utils-step4.gp contains common code used by both
kit-alpha1Variable.gp and kit-alpha3Variable.gp.

HOW TO USE
==========
We recommend that you create a copy of one of the
eg-*.gp files for the particular linear form that is of interest to you.

The functions with names like egX_search_it1() are for setting the bounds on the search.
The functions with names like egX_search_general() are for setting the values for your actual problem
and running the search.

The variables that need to be assigned values in egX_search_general() are:
	bigD       = [Q(al_1,al_2,al_3):Q] -- used for Matveev's bounds
	matveevChi = [R(al_1,al_2,al_3):R] -- used for Matveev's bounds

	al1,hgtA1,absLogA1,a1 = quantities associated with alpha_1 (al1 is an expression for alpha_1 itself)
	al2,hgtA2,absLogA2,a2 = quantities associated with alpha_2 (al2 is an expression for alpha_2 itself)
	al3,hgtA3,absLogA3,a3 = quantities associated with alpha_3 (al3 is an expression for alpha_3 itself)

	b1,b2,b3 = upper bound for absolute values of b_1, b_2 and b_3
	
	logXLB = lower bound (hence "LB") for the log of the variable that one of the alpha_i's depends on
	nLB    = lower bound for variable to denote the largest of the b_i's. This must be less than your final upper bound from the kit, but also as large as possible to get best results
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1,lamUB0

VARIABLE NAMING CONVENTION:
We assume that the quantity you are trying to bound from above with the kit is called n.
We also assume that x is the non-constant alpha_i.
These two naming conventions are important and used throughout the code.
Once again, just follow the examples.

Start by writing and running egX_search_it1().
It will produce output like the following, giving values of the parameters used and some other values.
The key value is nUB, where UB=upper bound, so this is giving you the upper bound on n (see above)
for those values of the parameters.
nUB will be decreasing, so the last value displayed will be the smallest value found for that iteration.
for chi= 2.500000, minNonDegenNUB=3.977630 e9
L= 301, m= 28.0000, rho(3logs)=  4.2500, chi= 2.5000, K= 4566501.848*logX, nonDegen log|Lambda|>-1.988815 e9*logX, nonDegenNUB=3.977630 e9, rho(2logs)= 7.000000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=5.687211 e9, degenNUB3=7.893746 e9, nUB=5.687211 e9, eliminate-b2
for chi= 2.750000, minNonDegenNUB=4.110874 e9
L= 306, m= 28.0000, rho(3logs)=  4.2500, chi= 2.7500, K= 4642357.360*logX, nonDegen log|Lambda|>-2.055437 e9*logX, nonDegenNUB=4.110874 e9, rho(2logs)= 7.000000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=4.998104 e9, degenNUB3=6.956398 e9, nUB=4.998104 e9, eliminate-b2
for chi= 3.000000, minNonDegenNUB=4.229909 e9
L= 305, m= 29.0000, rho(3logs)=  4.2500, chi= 3.0000, K= 4792442.910*logX, nonDegen log|Lambda|>-2.114954 e9*logX, nonDegenNUB=4.229909 e9, rho(2logs)= 7.000000, mu(2logs)= 0.610000, degenNUB1=    0.e-19, degenNUB2=4.501157 e9, degenNUB3=6.268175 e9, nUB=4.501157 e9, eliminate-b2

With that best value for the iteration, use the values of the parameters to make
a egX_check_it1() function that checks this particular choice of the parameters.
Again, the easiest way is to copy an existing egX_check_it1() function from one of the
example files provided in github and just change the variables in that function, as required.

This egX_check_it1() function outputs lots more information and can be used for testing,
as well as outputting details that you may want to include in a paper for your description
of how you used the kit (this is helpful info for readers and for referees -- see the
exposition of the examples in Section 6 of the kit paper).

The purpose of egX_search_it2(), egX_search_it3(),... functions is to tune the search range for further iterations.
They take as an argument nUBInit = initial upper bound for n (the variable denoting the largest of the b_i's).
egX_search_it1() does not need nUBInit, as nUBInit gets calculated from Matveev's bound in this case.

choosing ranges:
Since m, rho3 (the value of rho in the non-degenerate case), chi and rho2 (the value of rho to use in the degenerate case)
take 20 equidistributed values in the ranges provided, a smaller range will provide finer step sizes between the values.
The ranges in the example files provided give a good indication of appropriate size of the ranges.
The endpoints of the ranges should be chosen so that the values of the parameters that get displayed (see above)
are roughly in the middle of the range (or at least away from the endpoints). So some experimentation and
adjustment of these endpoints will likely be needed when initially applying this code to your problem.

Tips:
(1) in the egX_search_general() functions, do not forget that a1, a2 and a3 need to be set
inside the for loop for rho3, not before that. This is because a1, a2 and a3 will depend on rho3.

(2) in the egX_search_general() functions, make sure you use the correct alphaX_do_step3() and alphaX_do_step4() functions
for your situation.
If alpha_1 is the algebraic number whose value is not constant, then use alpha1_do_step3() and alpha1_do_step4()
If alpha_3 is the algebraic number whose value is not constant, then use alpha3_do_step3() and alpha3_do_step4()
But the code will give you an error, if you use the incorrect one.

CONTACT
=======
If you have any questions, problems, need a hand, find a bug,..., please contact PV
you can find PV's details in the arxiv preprint for the kit or in the published
version in Math. Comp.