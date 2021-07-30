
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

Local M = Q/16/pi_/pi_*i_*
	i_*gs*(LG*PL(1)+RG*PR(1))* T(b,beta,alphap)*
	i_*(-g_(1,o)*(ll(o)+pa(o)+pb(o))+MG*gi_(1))*
	(-gs)*f(b,a,c)*g_(1,mu)*
	i_*(-g_(1,rho)*(ll(rho)+pa(rho))+MG*gi_(1))*
	i_*gs*(RG*PL(1)+LG*PR(1))* T(c,alphap,alpha)*
	i_
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

id pa.pa = 0;
id pb.pb = 0;

**********************************************************
*                  Insert Loop tensors 			 *
**********************************************************

#procedure seti(name,N);
#redefine `name'si ",`name'1?,...,`name'`N'?"
#redefine `name'sf ",`name'1,...,`name'`N'"
#if ( `N' == 0)
#redefine `name'si ""
#redefine `name'sf ""
#endif
#endprocedure;


#define nusi ""
#define nusf ""
#define musi ""
#define musf ""
#define rhosi ""
#define rhosf ""

*TODO 3x l nom cases

* 2x l
#do M = {4,3,2,1,0}
#do N = {4,3,2,1,0}
#do R = {4,3,2,1,0}
#call seti(nu,`N')
#call seti(mu,`M')
#call seti(rho,`R')
id Q*g_(l? `nusi' , ll, mu?, ll)= -g_(l `nusf' , mu)*C(o,o)*dr + 2*g_(l `nusf' , o)*C(mu,o);
id Q*g_(l? `nusi' , ll `musi', ll `rhosi')= g_(l `nusf' , tm `musf', tn `rhosf')*C(tm,tn);
id Q*g_(l? `nusi' , ll `musi')*ll(tn?) = g_(l `nusf' , tm `musf')*C(tm,tn);
id Q*g_(l? `nusi' , ll `musi')*ll.p? =  g_(l `nusf' , tm `musf')*C(tm,tn)*p(tn);
#enddo
#enddo
#enddo
id ll(o?)*p?.ll*Q = p(tn)*C(o,tn);
id ll(o?)*ll(mu?)*Q = C(o,mu);
id ll.ll*Q = C(rho,rho);

* 1x l
#do M = {4,3,2,1,0}
#do N = {4,3,2,1,0}
#call seti(nu,`N')
#call seti(mu,`M')
id Q*g_(l? `nusi' , ll `musi') =  g_(l `nusf' , tm `musf')*C(tm);
*#message `nusi',`musi'
#enddo
#enddo
id p?.ll*Q = p(m)*C(m);

* 0x l
id Q = C0;

id LG*RG = 0;
id Nc = Ca;
Bracket+ A,B,C,D;.sort;print M;.sort;