Dimension 4;
AutoDeclare Vector p,k,l;
AutoDeclare Index mu,nu;
Vector q;
Tensor f(antisymmetric),polsum(symmetric);
Functions V,Ei,Eo,Ui,Uo,ui,uo,Leg,
	  df,PL,PR,VF,da,
;
CFunctions D,eps,epsbar,T,C,C0,C1,C2,C00,C11,C12,C22,g(symmetric),DIM,Denom;
Symbol MUs,MXs,MGs,MG,L#C,R#C,LG#C,RG#C,gs,MX,Ca,Cf,u,s,Nc,Tr,t,P1,P2,MG,dr,ax,dim;
Index a,gamma,beta,alpha,l,tm,tn,o,b,c,d,k,j;

Local MMMM = 
	2*
    polsum(mu,nu)/96*
    VF(mu,a,b)*
    (g_(1,tn)*(pa(tn)-p1(tn)) + MG)*Denom((pa(c)-p1(c))*(pa(c)-p1(c))-MGs)*
   	-i_*gs*(RG*PL(1)+LG*PR(1))*T(b,gamma,beta)*
*   i_*(R *PL(1)+L*PR(1))*
	g_(1,pa)*
	(-i_)*(R*PL(1)+L*PR(1))*
	(g_(1,p2)+MX)*
*	(-i_*gs)*(p1(mu)+pa(mu))*T(a,gamma,beta)*
*	(-i_*gs)*(p1(mu)*P1+pb(mu)*P2)*T(a,gamma,beta)*
	i_*gs*(pa(nu)+p2(nu))*T(a,beta,gamma)
;
	
id VF(mu,a,b) = 
	   + g_(1,7_,q)*q(mu) * ( 1/64*da(a,b)*C2*i_*pi_^-2*gs^2*L*RG + 1/32*da(a,
         b)*C22*i_*pi_^-2*gs^2*L*RG )

       + g_(1,7_,q)*pa(mu) * ( 1/64*da(a,b)*C2*i_*pi_^-2*gs^2*L*RG + 1/32*da(a
         ,b)*C12*i_*pi_^-2*gs^2*L*RG )

       + g_(1,6_,q)*q(mu) * ( 1/64*da(a,b)*C2*i_*pi_^-2*gs^2*R*LG + 1/32*da(a,
         b)*C22*i_*pi_^-2*gs^2*R*LG )

       + g_(1,6_,q)*pa(mu) * ( 1/64*da(a,b)*C2*i_*pi_^-2*gs^2*R*LG + 1/32*da(a
         ,b)*C12*i_*pi_^-2*gs^2*R*LG )

       + 1/32*g_(1,7_,mu)*da(a,b)*C00*i_*pi_^-2*gs^2*L*RG + 1/32*g_(1,6_,mu)*
         da(a,b)*C00*i_*pi_^-2*gs^2*R*LG;

*id polsum(mu?,nu?) = -d_(mu,nu);
*id polsum(mu?,nu?) = -d_(mu,nu) + (pa(mu)*pb(nu)+pa(nu)*pb(mu))/pa.pb;
*id polsum(mu?,nu?) = -d_(mu,nu) + ax*((q(mu)*pb(nu)+q(nu)*pb(mu))/q.pb - q.q*pb(mu)*pb(nu)/q.pb/q.pb);
id polsum(mu?,nu?) = -d_(mu,nu) + ((q(mu)*pb(nu)+q(nu)*pb(mu))/q.pb - q.q*pb(mu)*pb(nu)/q.pb/q.pb);
id PL(l?) = g7_(l)/2;
id PR(l?) = g6_(l)/2;

id T(a?,gamma?,beta?)*T(a?,beta?,alpha?)=Cf*df(alpha,gamma);
id df(alpha?,alpha?) = Ca;

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

* Bracket Ca,Cf,gs,L,R;print MMMM;.end

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

.sort
*Format C;
Factorize;
Bracket Ca,Cf,gs,L,R,pa.pb,dr;print MMMM;
#write<fun.f> "m *=%e",MMMM;
.end