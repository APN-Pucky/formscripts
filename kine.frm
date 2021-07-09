dimension 4;
AutoDeclare Vector p;
Vectors p1,p2,p3,H,q,A,B;
CFunctions C,C0,C1,C2,C00,C11,C12,C22,Dim;
Symbols m1,...,m3,U,MQs,Q;
Indices mu,rho,nu,o,n,m,tm,tn;

Local M = Q*H(m)*
	(	- d_(o,n)*(pA(m)-pB(m)) 
		- d_(n,m)*(pB(o)-pC(o)) 
		- d_(o,m)*(pC(n)-pA(n)))*
	(p1(o)-q(o))*(q(n)-p1(n)-p2(n))
;


*Print M;Bracket H;.end
id pA = -q-p1-p2;
id pB = q+p1;
id pC = p2;



id p2.p2 = 0 ;
id H?.q*p?.q*Q = H(tm)*p(tn)*C(tm,tn);
id q.q*Q = C(rho,rho);
id p?.q*Q = p(m)*C(m);
id Q = C0;



repeat;
	id C(mu?,mu?)=
		 Dim(C00)
		+p1(mu)*p1(mu)*C11
		+(p2(mu)+p1(mu))*(p2(mu)+p1(mu))*C22
		+((p2(mu)+p1(mu))*p1(mu)+p1(mu)*(p2(mu)+p1(mu)))*C12;
	id C(mu?,nu?)=
		 d_(mu,nu)*C00
		+p1(mu)*p1(nu)*C11
		+(p2(mu)+p1(mu))*(p2(nu)+p1(nu))*C22
		+((p2(mu)+p1(mu))*p1(nu)+p1(mu)*(p2(nu)+p1(nu)))*C12;
	id C(mu?)=p1(mu)*C1+(p2(mu)+p1(mu))*C2;
*	id p3 = p1+p2;
	id p1.p2 = 1/2*(p3.p3-p1.p1-p2.p2);
	id p2.p2 = 0;
	id p1.p1 = U;
	id p3.p3 = MQs;
endrepeat;

Print M;Bracket H;.end