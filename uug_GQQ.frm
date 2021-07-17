dimension Dim;
AutoDeclare Tensor tens;
AutoDeclare Vector p,ll;
Vectors p1,H,A,B,pu,pb;
Tensors f(antisymmetric);
Function PL,PR,df,da;
CFunctions C,C0,C1,C2,C00,C11,C12,C22,T;
Symbols m1,...,m3,U,MQs,Q,gs,L#C,R#C,MG,Tr,Nc,Cf,Ca;
Indices mu,rho,nu,o,n,m,tm,tn,beta,b,betap,alphap,a,alpha,ind,delta,k,j,l,c,d;

Local M = Q*H(mu)/16/pi_/pi_*
	i_*gs*(L*PL(1)+R*PR(1))* T(b,beta,betap)*
	i_*(-g_(1,ll)+MG)*
	i_*gs*(R*PL(1)+L*PR(1))* T(b,alphap,alpha)*
	i_*
	(-i_)*gs*(2*ll(mu)+2*pa(mu)+pb(mu))*T(a,betap,alphap)

;

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

*Print M;Bracket H;.end
*id p1 = pu+pb;
*id pA = pb;
*id pB = -(ll+pu+pb);
*id pC = l+pu;

id pb.pb = 0;
id pa.pa = 0;


id g_(l?,7_,ll)*p?.ll*Q = g_(l,7_,tm)*p(tn)*C(tm,tn);
id g_(l?,6_,ll)*p?.ll*Q = g_(l,6_,tm)*p(tn)*C(tm,tn);
id g_(l?,pu?,ll)*p?.ll*Q = g_(l,pu,tm)*p(tn)*C(tm,tn);
id H?.ll*p?.ll*Q = H(tm)*p(tn)*C(tm,tn);
id ll.ll*Q = C(rho,rho);

id g_(l?,7_,ll)*Q = g_(l,7_,tm)*C(tm);
id g_(l?,6_,ll)*Q = g_(l,6_,tm)*C(tm);
id g_(l?,pu?,ll)*Q = g_(l,pu,tm)*C(tm);
id p?.ll*Q = p(m)*C(m);
id Q = C0;
id L*R = 0;
id Nc = 2*(1/2/Ca+Cf);

Print M;Bracket H;.end

repeat;
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

*Format C;
Factorize ;
Print M;.end