import form
import re
import feyn
count = 0

init = '''
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
'''

pol = '''
repeat;
id PL(ind?) = g7_(ind)/2;
id PR(ind?) = g6_(ind)/2;
endrepeat;
'''

color = '''
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
'''
trace ="Trace4 1;"

polsum_feyn='''
id polsum(mu?,nu?) = -d_(mu,nu);
'''

polsum_phys='''
*id polsum(mu?,nu?) = -d_(mu,nu) + (pa(mu)*pb(nu)+pa(nu)*pb(mu))/pa.pb;
*id polsum(mu?,nu?) = -d_(mu,nu) + ax*((q(mu)*pb(nu)+q(nu)*pb(mu))/q.pb - q.q*pb(mu)*pb(nu)/q.pb/q.pb);
id polsum(mu?,nu?) = -d_(mu,nu) + ((q(mu)*pb(nu)+q(nu)*pb(mu))/q.pb - q.q*pb(mu)*pb(nu)/q.pb/q.pb);
'''

onshell='''
id pa.pa = 0;
id pb.pb = 0;
'''

insert_loops = '''
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

*****
* D gammas
*****

#do M = {0,1,2,3,4}
#do N = {0,1,2,3,4}
#do R = {4,3,2,1,0}
#call seti(nu,`N')
#call seti(mu,`M')
#call seti(rho,`R')
id g_(l? `nusi' , ll, mu?, ll `rhosi')= -g_(l `nusf' , mu  `rhosf')*ll.ll*dr + 2*g_(l `nusf' , o  `rhosf')*(ll(mu)*ll(o));
id g_(l? `nusi' , mu?, ll, mu?  `musi') =  (2-dim*dr)*g_(l `nusf' , ll `musf');
id g_(l? `nusi' , ll ,ll `rhosi')= g_(l `nusf' `rhosf')*ll.ll*dr;
*#message `nusi',`musi'
#enddo
#enddo
#enddo

*TODO 3x l nom cases

* 2x l
#do M = {4,3,2,1,0}
#do N = {4,3,2,1,0}
#do R = {4,3,2,1,0}
#call seti(nu,`N')
#call seti(mu,`M')
#call seti(rho,`R')
id Q*g_(l? `nusi' , ll ,ll `rhosi')= g_(l `nusf' `rhosf')*C(o,o)*dr;
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
#do M = {0,1,2,3,4}
#do N = {0,1,2,3,4}
#call seti(nu,`N')
#call seti(mu,`M')
id Q*g_(l? `nusi' , ll `musi') =  g_(l `nusf' , tm `musf')*C(tm);
*#message `nusi',`musi'
#enddo
#enddo
id p?.ll*Q = p(m)*C(m);

* 0x l
id Q = C0;
'''

reduce_loops='''
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
'''

kinematics='''
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
'''

def loop(s):
    return run(s + "Bracket+ A,B,C,D;.sort;")
def full(s):
    return run(s + "Factorize;")

def run(s):
   global count
   count =count+1
   with open("form" + str(count)+ ".frm","w") as frm:
        with form.open(keep_log=1000) as f:
            l = s.split("Local")[1].split("=")[0].strip()
            txt = s + 'print '  +l + ";.sort;"
            f.write(txt)
            frm.write(txt)
            r = f.read("" + l)
            r = re.sub(r'\+factor_\^?[0-9]*',r'',r).strip("*") + ";"
            print(r+ "\n")
            return r
# TODO debug print to file and run it manually

def toC(s):
    with form.open(keep_log=1000) as f:
        l = "TMP"
        txt = feyn.init + "L TMP = " + s + ";Factorize;Format C;.sort;"
        f.write(txt)
        r = f.read("" + l)
        r = re.sub(r'factor_',r'',r)
        r = re.sub(r';_\+=\+pow\(,\d*\)',r'',r)
        r = r.strip("+").strip("*")
        print(r)
        return r
