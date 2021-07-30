from feyn import *
'''
Local M = Q/16/pi_/pi_*i_*
	-i_*gs*g_(1,rho)*T(b,beta,alphap)*
	-i_*g_(1,ll)*
	-i_*gs*g_(1,kappa)*T(c,alphap,alpha)*
	-i_*d_(kappa,sigma)*
	-gs*f(c,b,a)*(
		- d_(mu,nu)*(pA(sigma)-pB(sigma)) 
		- d_(nu,sigma)*(pB(mu)-pC(mu)) 
		- d_(sigma,mu)*(pC(nu)-pA(nu)))*
	-i_*d_(nu,rho)
;
'''
def ggq():
    return loop( init+ '''
Local M = Q/16/pi_/pi_*i_*
	-i_*gs*g_(1,rho)*
	-i_*g_(1,ll)*
	-i_*gs*g_(1,kappa)*
	-i_*d_(kappa,sigma)*
	-gs*(
		- d_(mu,nu)*(pA(sigma)-pB(sigma)) 
		- d_(nu,sigma)*(pB(mu)-pC(mu)) 
		- d_(sigma,mu)*(pC(nu)-pA(nu)))*
	-i_*d_(nu,rho)*
	-i_/2*Ca*T(a,beta,alpha)
;
id pA = pb;
id pB = -(ll+pa+pb);
id pC = ll+pa;
'''
+ pol + color + onshell + insert_loops + 
'''
*id LG*RG = 0;
*id Nc = 2*(1/2/Ca+Cf);
'''
)
def qqg():
    return loop( init+ '''
Local M = Q/16/pi_/pi_*i_*
	-i_*gs*g_(1,nu)*T(b,beta,betap)*
	i_*(g_(1,o)*(ll(o)+pa(o)+pb(o)))*
	-i_*gs*g_(1,mu)*T(a,betap,alphap)*
	i_*(g_(1,rho)*(ll(rho)+pa(rho)))*
	-i_*gs*g_(1,sigma)*T(b,alphap,alpha)*
	-i_*d_(sigma,nu)
;
'''
+ pol + color + onshell + insert_loops + 
'''
*id LG*RG = 0;
*id Nc = 2*(1/2/Ca+Cf);
'''
)

def GQQ():
    return loop( init+ '''
Local M = Q/16/pi_/pi_*i_*
	i_*gs*(LG*PL(1)+RG*PR(1))* T(b,beta,betap)*
	i_*(-g_(1,ll)+MG*gi_(1))*
	i_*gs*(RG*PL(1)+LG*PR(1))* T(b,alphap,alpha)*
	i_*
	(-i_)*gs*(2*ll(mu)+2*pa(mu)+pb(mu))*T(a,betap,alphap)*
	i_
;
'''
+ pol + color + onshell + insert_loops + 
'''
id LG*RG = 0;
id Nc = 2*(1/2/Ca+Cf);
'''
)

def GGQ():
    return loop( init+ '''
Local M = Q/16/pi_/pi_*i_*
	i_*gs*(LG*PL(1)+RG*PR(1))* T(b,beta,alphap)*
	i_*(-g_(1,o)*(ll(o)+pa(o)+pb(o))+MG*gi_(1))*
	(-gs)*f(b,a,c)*g_(1,mu)*
	i_*(-g_(1,rho)*(ll(rho)+pa(rho))+MG*gi_(1))*
	i_*gs*(RG*PL(1)+LG*PR(1))* T(c,alphap,alpha)*
	i_
;
'''
+ pol + color + onshell + insert_loops + 
'''
id LG*RG = 0;
id Nc = Ca;
'''
)

def schannel(VF):
    return full(init+
	'''
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
'''
+
'''
id VF(mu?,a?,beta?,alpha?) = 
'''
+
VF() + ";"
+ polsum_phys
+pol
+color
+trace
+ reduce_loops
+ kinematics
)

