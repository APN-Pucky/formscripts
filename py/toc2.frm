
Symbols s,t,u,dim,m1,...,m3,U,MG,MGs,MU,MX,MXs,MUs,Q,gs,L#C,R#C,LG#C,RG#C,Lp#C,Rp#C,LGp#C,RGp#C,MG,Tr,Nc,Cf,Ca,dr;
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
L TMP = C1*s^-1*MXs*MUs^3*R*RG*Lp*LGp + C1*s^-1*MXs*MUs^3*L*LG*Rp*RGp - C1*
         s^-1*t*MUs^3*R*RG*Lp*LGp - C1*s^-1*t*MUs^3*L*LG*Rp*RGp - 3*C1*s^-1*t*
         MXs*MUs^2*R*RG*Lp*LGp - 3*C1*s^-1*t*MXs*MUs^2*L*LG*Rp*RGp + 3*C1*s^-1
         *t^2*MUs^2*R*RG*Lp*LGp + 3*C1*s^-1*t^2*MUs^2*L*LG*Rp*RGp + 3*C1*s^-1*
         t^2*MXs*MUs*R*RG*Lp*LGp + 3*C1*s^-1*t^2*MXs*MUs*L*LG*Rp*RGp - 3*C1*
         s^-1*t^3*MUs*R*RG*Lp*LGp - 3*C1*s^-1*t^3*MUs*L*LG*Rp*RGp - C1*s^-1*
         t^3*MXs*R*RG*Lp*LGp - C1*s^-1*t^3*MXs*L*LG*Rp*RGp + C1*s^-1*t^4*R*RG*
         Lp*LGp + C1*s^-1*t^4*L*LG*Rp*RGp - C1*MXs*MUs^2*R*RG*Lp*LGp - C1*MXs*
         MUs^2*L*LG*Rp*RGp + 2*C1*t*MUs^2*R*RG*Lp*LGp + 2*C1*t*MUs^2*L*LG*Rp*
         RGp + 2*C1*t*MXs*MUs*R*RG*Lp*LGp + 2*C1*t*MXs*MUs*L*LG*Rp*RGp - 4*C1*
         t^2*MUs*R*RG*Lp*LGp - 4*C1*t^2*MUs*L*LG*Rp*RGp - C1*t^2*MXs*R*RG*Lp*
         LGp - C1*t^2*MXs*L*LG*Rp*RGp + 2*C1*t^3*R*RG*Lp*LGp + 2*C1*t^3*L*LG*
         Rp*RGp - C1*s*t*MUs*R*RG*Lp*LGp - C1*s*t*MUs*L*LG*Rp*RGp + C1*s*t^2*R
         *RG*Lp*LGp + C1*s*t^2*L*LG*Rp*RGp + C2*s^-1*MXs*MUs^3*R*RG*Lp*LGp + 
         C2*s^-1*MXs*MUs^3*L*LG*Rp*RGp - C2*s^-1*MXs^2*MUs^2*R*RG*Lp*LGp - C2*
         s^-1*MXs^2*MUs^2*L*LG*Rp*RGp - C2*s^-1*t*MUs^3*R*RG*Lp*LGp - C2*s^-1*
         t*MUs^3*L*LG*Rp*RGp - C2*s^-1*t*MXs*MUs^2*R*RG*Lp*LGp - C2*s^-1*t*MXs
         *MUs^2*L*LG*Rp*RGp + 2*C2*s^-1*t*MXs^2*MUs*R*RG*Lp*LGp + 2*C2*s^-1*t*
         MXs^2*MUs*L*LG*Rp*RGp + 2*C2*s^-1*t^2*MUs^2*R*RG*Lp*LGp + 2*C2*s^-1*
         t^2*MUs^2*L*LG*Rp*RGp - C2*s^-1*t^2*MXs*MUs*R*RG*Lp*LGp - C2*s^-1*t^2
         *MXs*MUs*L*LG*Rp*RGp - C2*s^-1*t^2*MXs^2*R*RG*Lp*LGp - C2*s^-1*t^2*
         MXs^2*L*LG*Rp*RGp - C2*s^-1*t^3*MUs*R*RG*Lp*LGp - C2*s^-1*t^3*MUs*L*
         LG*Rp*RGp + C2*s^-1*t^3*MXs*R*RG*Lp*LGp + C2*s^-1*t^3*MXs*L*LG*Rp*RGp
          - C2*MXs*MUs^2*R*RG*Lp*LGp - C2*MXs*MUs^2*L*LG*Rp*RGp + C2*MG*MX*MXs
         *MUs*R*LG*Rp*LGp + C2*MG*MX*MXs*MUs*L*RG*Lp*RGp + 2*C2*t*MUs^2*R*RG*
         Lp*LGp + 2*C2*t*MUs^2*L*LG*Rp*RGp - C2*t*MG*MX*MUs*R*LG*Rp*LGp - C2*t
         *MG*MX*MUs*L*RG*Lp*RGp - C2*t*MG*MX*MXs*R*LG*Rp*LGp - C2*t*MG*MX*MXs*
         L*RG*Lp*RGp - 2*C2*t^2*MUs*R*RG*Lp*LGp - 2*C2*t^2*MUs*L*LG*Rp*RGp + 
         C2*t^2*MXs*R*RG*Lp*LGp + C2*t^2*MXs*L*LG*Rp*RGp + C2*t^2*MG*MX*R*LG*
         Rp*LGp + C2*t^2*MG*MX*L*RG*Lp*RGp - C2*s*t*MUs*R*RG*Lp*LGp - C2*s*t*
         MUs*L*LG*Rp*RGp + C2*s*t*MG*MX*R*LG*Rp*LGp + C2*s*t*MG*MX*L*RG*Lp*RGp
          - 2*C00*s^-1*MXs*MUs^2*R*RG*Lp*LGp - 2*C00*s^-1*MXs*MUs^2*L*LG*Rp*
         RGp + 2*C00*s^-1*t*MUs^2*R*RG*Lp*LGp + 2*C00*s^-1*t*MUs^2*L*LG*Rp*RGp
          + 4*C00*s^-1*t*MXs*MUs*R*RG*Lp*LGp + 4*C00*s^-1*t*MXs*MUs*L*LG*Rp*
         RGp - 4*C00*s^-1*t^2*MUs*R*RG*Lp*LGp - 4*C00*s^-1*t^2*MUs*L*LG*Rp*RGp
          - 2*C00*s^-1*t^2*MXs*R*RG*Lp*LGp - 2*C00*s^-1*t^2*MXs*L*LG*Rp*RGp + 
         2*C00*s^-1*t^3*R*RG*Lp*LGp + 2*C00*s^-1*t^3*L*LG*Rp*RGp + C00*MXs*MUs
         *R*RG*Lp*LGp + C00*MXs*MUs*L*LG*Rp*RGp - 3*C00*t*MUs*R*RG*Lp*LGp - 3*
         C00*t*MUs*L*LG*Rp*RGp - C00*t*MXs*R*RG*Lp*LGp - C00*t*MXs*L*LG*Rp*RGp
          + 3*C00*t^2*R*RG*Lp*LGp + 3*C00*t^2*L*LG*Rp*RGp + C00*s*t*R*RG*Lp*
         LGp + C00*s*t*L*LG*Rp*RGp;

Bracket+ C00,C1,C2,;
Format C;Print TMP;.sort;
