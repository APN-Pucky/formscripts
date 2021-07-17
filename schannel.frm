Dimension d;
AutoDeclare Vector p,k,l;
AutoDeclare Index mu,nu;
Vector q;
Tensor f(antisymmetric),polsum;
Functions V,Ei,Eo,Ui,Uo,ui,uo,g,Leg,
	  df,PL,PR
;
CFunctions D,eps,epsbar,T;
Symbol b,c,MUs,MXs,L#C,R#C,gs,MX,Ca,Cf,u,s,t,P1,P2;
Index a,gamma,beta,alpha,l;

Local MMMM = 
	polsum(mu,nu)/96*
	i_*(R *PL(1)+L*PR(1))*
	i_*g_(1,nu1)*(pa(nu1)+pb(nu1))/q.q*
	-i_*gs*g_(1,mu)*T(a,gamma,beta)*
	g_(1,pa)*
	i_*gs*g_(1,nu)*T(a,beta,gamma)*
	(-i_)*g_(1,nu2)*(pa(nu2)+pb(nu2))/q.q*
	(-i_)*(L*PL(1)+R*PR(1))*
	(g_(1,p2)+MX)
;

*id polsum(mu?,nu?) = -d_(mu,nu);
*id polsum(mu?,nu?) = -d_(mu,nu) + (pa(mu)*pb(nu)+pa(nu)*pb(mu))/pa.pb;
id polsum(mu?,nu?) = -d_(mu,nu) + (q(mu)*pb(nu)+q(nu)*pb(mu))/q.pb - q.q*pb(mu)*pb(nu)/q.pb/q.pb;
id PL(l?) = g7_(l)/2;
id PR(l?) = g6_(l)/2;

id T(a?,gamma?,beta?)*T(a?,beta?,alpha?)=Cf*df(alpha,gamma);
id df(alpha?,alpha?) = Ca;

Trace4,1;

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
Format C;
Factorize;
Bracket Ca,Cf,gs,L,R,pa.pb;print MMMM;
#write<fun.f> "m *=%e",MMMM;
.end