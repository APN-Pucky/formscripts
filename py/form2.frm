
Symbols s,t,u,dim,m1,...,m3,U,MG,MGs,MU,MX,MXs,MUs,Q,gs,L#C,R#C,LG#C,RG#C,MG,Tr,Nc,Cf,Ca,dr;
dimension dim;
AutoDeclare Tensor tens;
AutoDeclare Vector p,ll,k,l;
AutoDeclare Index mu,nu,rho,kappa,gamma,sigma;
AutoDeclare Symbol sym;
Vectors p1,H,pu,pb,q,l1,l2;
Tensors f(antisymmetric),polsum(symmetric),Gamma;
Function PL,PR,df,da,VF;
CFunctions C,C0,C1,C2,C00,C11,C12,C22,T,A,B,D,Denom,DIMD;
Indices a,o,n,m,tm,tn,beta,b,betap,alphap,a,alpha,ind,delta,k,j,l,c,d;

    Local MMMM = 
	2*
	polsum(mu,nu)/96*
	i_*(R *PL(1)+L*PR(1))*
	i_*g_(1,nu1)*(pa(nu1)+pb(nu1))/q.q*
	VF(mu,a,gamma,beta)*
	g_(1,pa)*
	i_*gs*g_(1,nu)*T(a,beta,gamma)*
	(-i_)*g_(1,nu2)*(pa(nu2)+pb(nu2))/q.q*
	(-i_)*(L*PL(1)+R*PR(1))*
	(g_(1,p2)+MX*gi_(1))
;

id VF(mu?,a?,beta?,alpha?) = 
+C(tm)*(1/32*g_(1,pb,tm,mu)*T(a,beta,alpha)*i_*pi_^-2*gs^3*Ca-1/32*g_(1,pa,tm,mu)*T(a,beta,alpha)*i_*pi_^-2*gs^3*Ca+1/16*g_(1,tm)*T(a,beta,alpha)*pb(mu)*i_*pi_^-2*gs^3*Ca-1/32*g_(1,tm)*T(a,beta,alpha)*pb(mu)*i_*pi_^-2*dim*gs^3*Ca*dr+1/8*g_(1,tm)*T(a,beta,alpha)*pa(mu)*i_*pi_^-2*gs^3*Ca-1/16*g_(1,tm)*T(a,beta,alpha)*pa(mu)*i_*pi_^-2*dim*gs^3*Ca*dr-1/16*g_(1,mu,tm,pb)*T(a,beta,alpha)*i_*pi_^-2*gs^3*Ca-1/32*g_(1,mu,tm,pa)*T(a,beta,alpha)*i_*pi_^-2*gs^3*Ca)+C(tm,mu)*(1/8*g_(1,tm)*T(a,beta,alpha)*i_*pi_^-2*gs^3*Ca-1/16*g_(1,tm)*T(a,beta,alpha)*i_*pi_^-2*dim*gs^3*Ca*dr)+C(rho,rho)*(-1/16*g_(1,mu)*T(a,beta,alpha)*i_*pi_^-2*gs^3*Ca*dr);;
*id polsum(mu?,nu?) = -d_(mu,nu) + (pa(mu)*pb(nu)+pa(nu)*pb(mu))/pa.pb;
*id polsum(mu?,nu?) = -d_(mu,nu) + ax*((q(mu)*pb(nu)+q(nu)*pb(mu))/q.pb - q.q*pb(mu)*pb(nu)/q.pb/q.pb);
id polsum(mu?,nu?) = -d_(mu,nu) + ((q(mu)*pb(nu)+q(nu)*pb(mu))/q.pb - q.q*pb(mu)*pb(nu)/q.pb/q.pb);

repeat;
id PL(ind?) = g7_(ind)/2;
id PR(ind?) = g6_(ind)/2;
endrepeat;

**********************************************************
*                  COLOUR STRUCTURE SIMPLIFY             *
**********************************************************
    
repeat;
* remove df(k,j)
   id df(k?,l?)*df(l?,j?)=df(k,j);
   id T(a?,k?,l?)*df(k?,j?)=T(a,j,l);id T(a?,k?,l?)*df(l?,j?)=T(a,k,j);
* remove da(a,b)
   id da(a?,b?)*da(b?,c?)=da(a,c);
   id T(a?,k?,l?)*da(a?,b?)=T(b,k,l);
   id f(a?,b?,c?)*da(a?,d?)=f(d,b,c);
* length-three objects simplify:
   id T(b?,k?,j?)*T(a?,j?,c?)*T(b?,c?,l?)=(Cf-Nc*Tr)*T(a,k,l);
   id T(b?,j?,l?)*T(c?,l?,k?)*f(a?,b?,c?)=i_*Nc*Tr*T(a,j,k);
* length-two objects that give out df(k,j)
   id T(a?,k?,l?)*T(a?,l?,j?)=Cf*df(k,j);
* length-two objects that give out da(a,b)
   id T(a?,k?,l?)*T(b?,l?,k?)=Tr*da(a,b);
   id f(a?,b?,c?)*f(d?,b?,c?)=Nc*da(a,d); 
* simplify traces
   id T(b?,k?,k?)=0;
   id da(a?,a?) = Nc*Cf/Tr;
   id df(a?,a?) = Nc;
* simplify combination of factors
   id Nc^-2=2-Nc^2+Cf^2*Tr^-2;
   id Nc^2=1+Nc*Cf/Tr;
   id Tr=1/2;
   id Tr^-1=2;
endrepeat;
Trace4 1;
repeat;
	id C(mu?,nu?)=
		 d_(mu,nu)*C00
		+l1(mu)*l1(nu)*C11
		+(l2(mu)+l1(mu))*(l2(nu)+l1(nu))*C22
		+((l2(mu)+l1(mu))*l1(nu)+l1(mu)*(l2(nu)+l1(nu)))*C12;
	id C(mu?)=l1(mu)*C1+(l2(mu)+l1(mu))*C2;
	id l1 = pa;
	id l2 = pb;
endrepeat;
id C00*dr = DIMD(C00);
id C11*dr = (C11);
id C12*dr = (C12);
id C22*dr = (C22);

repeat;
id q = pa+pb;
id q.pb^-1 = pa.pb^-1;
id q.q^-1 = 1/2*pa.pb^-1;
id k1 = pa-p2;
id k2 = pa-p2;
id pa.pa =0;
id pb.pb = 0;
id p1.p1 = MUs;
id p2.p2 = MXs;

id p1 = pa+pb-p2;
endrepeat;

* Bracket Ca,Cf,gs,L,R;print MMMM;.end

repeat;
id pa.pb = s/2;
id pa.pb^-1 = 2/s;
id pb.p1 = (MUs-u)/2;
id pa.p2 = (MXs-u)/2;
id p2.pb = (MXs-t)/2;
id u = -t-s +MXs+MUs;
endrepeat;
Factorize;print MMMM;.sort;