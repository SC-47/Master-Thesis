function F = myfun_4_c(x,pb,pe,am,c_,k0,k1,k2,wf_L,wf_H,alpha,gamma)
%known parameter : pb pe am c_ k0 k1 k2 N_L N_H
%assumption paramter : rho 
%solve parameter: bd_L be_L b_H alpha gamma
F = [ (1+gamma+alpha)*pb*x(1)+((1+gamma)*k2+alpha)*wf_L*k1*(x(1))^(k2)-alpha*(wf_L*(1-k0)+am-c_);
      (1+gamma+alpha)*(pb+pe)*x(2)+((1+gamma)*k2+alpha)*wf_L*k1*(x(2))^(k2)-alpha*(wf_L*(1-k0)+am-c_);
      (1+gamma+alpha)*(pb+pe)*x(3)+((1+gamma)*k2+alpha)*wf_H*k1*(x(3))^(k2)-alpha*(wf_H*(1-k0)+am-c_);
      ];


