Vectors p1,p2,p3,H;
CFunctions C,C0,C1,C2,C00,C11,C12,C22;
Symbols m1,...,m3,U,MQ;
Indices mu,rho,nu;

Local M = H(mu)*(
	 p1(mu)*(C(rho,p1,p2,p3,m1,m2,m3)*(4*p1(rho)+2*p2(rho))-4*C(rho,rho,p1,p2,p3,m1,m2,m3)-C0(p1,p2,p3,m1,m2,m3)*p1.p2)
	+p2(mu)*(C(rho,p1,p2,p3,m1,m2,m3)*p1(rho)-2*C(rho,rho,p1,p2,p3,m1,m2,m3)+C0(p1,p2,p3,m1,m2,m3)*p1.p1)
	+C(mu,rho,p1,p2,p3,m1,m2,m3)*(4*p1(rho)+2*p2(rho))-C(mu,p1,p2,p3,m1,m2,m3)*(4*p1.p1+3*p1.p2)
	)
;
repeat;
	id C(mu?,p1?,p2?,p3?,m1?,m2?,m3?)=p1(mu)*C1(p1,p2,p3,m1,m2,m3)+p2(mu)*C2(p1,p2,p3,m1,m2,m3);
	id C(mu?,nu?,p1?,p2?,p3?,m1?,m2?,m3?)=
		 d_(mu,nu)*C00(p1,p2,p3,m1,m2,m3)
		+p1(mu)*p1(nu)*C11(p1,p2,p3,m1,m2,m3)
		+p2(mu)*p2(nu)*C22(p1,p2,p3,m1,m2,m3)
		+(p2(mu)*p1(nu)+p1(mu)*p2(nu))*C12(p1,p2,p3,m1,m2,m3);
	id p3 = p1+p2;
	id p1.p2 = 1/2*(p3.p3-p1.p1-p2.p2);
	id p2.p2 = 0;
	id p1.p1 = U;
	id p3.p3 = MQ*MQ;
endrepeat;

* simplify
repeat;
	id C0(p1,p2,p3,m1,m2,m3) = C0;
	id C1(p1,p2,p3,m1,m2,m3) = C1;
	id C2(p1,p2,p3,m1,m2,m3) = C2;
	id C00(p1,p2,p3,m1,m2,m3) = C00;
	id C11(p1,p2,p3,m1,m2,m3) = C11;
	id C12(p1,p2,p3,m1,m2,m3) = C12;
	id C22(p1,p2,p3,m1,m2,m3) = C22;
endrepeat;

Bracket H;

Print M;