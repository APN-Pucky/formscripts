dimension 4;
AutoDeclare Vector p,l;
Vectors p1,H,A,B,pu,pb;
CFunctions C,C0,C1,C2,C00,C11,C12,C22,Dim;
Symbols m1,...,m3,U,MQs,Q;
Indices mu,rho,nu,o,n,m,tm,tn;

Local M = Q*H(m)*
	(	- d_(m,n)*(pA(o)-pB(o)) 
		- d_(n,o)*(pB(m)-pC(m)) 
		- d_(o,m)*(pC(n)-pA(n)))*
	(pu(o)-l(o))*(l(n)-p1(n))
;


*Print M;Bracket H;.end
id p1 = pu+pb;
id pA = pb;
id pB = -(l+pu+pb);
id pC = l+pu;

id pb.pb = 0 ;


id H?.l*p?.l*Q = H(tm)*p(tn)*C(tm,tn);
id l.l*Q = C(rho,rho);
id p?.l*Q = p(m)*C(m);
id Q = C0;


*Print M;Bracket H;.end

repeat;
	id C(mu?,mu?)=
		 Dim(C00)
		+l1(mu)*l1(mu)*C11
		+(l2(mu)+l1(mu))*(l2(mu)+l1(mu))*C22
		+((l2(mu)+l1(mu))*l1(mu)+l1(mu)*(l2(mu)+l1(mu)))*C12;
	id C(mu?,nu?)=
		 d_(mu,nu)*C00
		+l1(mu)*l1(nu)*C11
		+(l2(mu)+l1(mu))*(l2(nu)+l1(nu))*C22
		+((l2(mu)+l1(mu))*l1(nu)+l1(mu)*(l2(nu)+l1(nu)))*C12;
	id C(mu?)=l1(mu)*C1+(l2(mu)+l1(mu))*C2;
	id l1 = pu;
	id l2 = pb;
	id pu.pb = 1/2*(p1.p1-pu.pu-pb.pb);
	id pb.pb = 0;
	id pu.pu = U;
	id p1.p1 = MQs;
endrepeat;

Print M;Bracket H;.end