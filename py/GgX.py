from feyn import *
def QQq():
    return loop( init+ '''
    Local M = Q/16/pi_/pi_*i_*
	i_*(Lp*PL(1)+Rp*PR(1))*
	i_*(-g_(1,ll))*
	-i_*gs*(LG*PL(1)+RG*PR(1))* T(b,alphap,alpha)*
	i_*
	-i_*gs*(2*ll(mu)+2*pa(mu)-2*p1(mu)+pb(mu))*T(a,alpha,alphap)*
	i_
;
'''
+ pol + color + onshell + insert_loops + 
'''
id LG*RG = 0;
id Nc = Ca;
'''
)

def qqQ():
    return loop( init+ '''
    Local M = Q/16/pi_/pi_*i_*
	i_*(Lp*PL(1)+Rp*PR(1))*
	i_*(-g_(1,rho)*(ll(rho)+pt(rho)+pb(rho)))*
	-i_*gs*g_(1,mu)*T(a,alpha,alphap)*
	i_*(-g_(1,nu)*(ll(nu)+pt(nu)))*
	i_*gs*(LG*PL(1)+RG*PR(1))* T(b,alphap,alpha)*
	i_
;
'''
+ pol + color + onshell + insert_loops + 
'''
id pt = pa-p1;
id LG*RG = 0;
id Nc = Ca;
'''
)

def uchannel(VF):
    return full(init+
    '''
    Local MMMM = 
	2*
    polsum(mu,nu)/96*
    VF(mu,a,b)*
    (g_(1,tn)*(pa(tn)-p1(tn)) + MG)*Denom(t-MGs)*
    -i_*gs*(LGp*PL(1)+RGp*PR(1))*T(b,gamma,beta)*
    g_(1,pa)*
    (-i_)*(L*PL(1)+R*PR(1))*
    (g_(1,p2)+MX)*
    i_*gs*(pa(nu)-p2(nu)+p1(nu))*T(a,beta,gamma)*Denom(u-MUs)
;
'''
+
'''
id VF(mu?,a?,b) = 
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