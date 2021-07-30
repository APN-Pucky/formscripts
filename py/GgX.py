from feyn import *
def QQq():
    return loop( init+ '''
    Local M = Q/16/pi_/pi_*i_*
	i_*(L*PL(1)+R*PR(1))*
	i_*(-g_(1,ll))*
	-i_*gs*(LG*PL(1)+RG*PR(1))* T(b,alphap,alpha)*
	i_*
	-i_*gs*(2*ll(mu)+2*pa(mu)+pb(mu))*T(a,alpha,alphap)*
	i_
;
'''
+ pol + color + onshell + insert_loops + 
'''
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