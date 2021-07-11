function F = myfun_5_c(x,pb,pe,am,c_,k0,k1,k2,wf_L,wf_H,phi,eta)
%known parameter : pb pe am c_ k0 k1 k2 N_L N_H
%assumption paramter : rho 
%solve parameter: bd_L be_L b_H 
F = [ (1+phi*(x(1)^(1-eta)))*pb*x(1)+(k2+phi*(x(1)^(1-eta)))*wf_L*k1*((x(1))^(k2))-(phi*(x(1)^(1-eta)))*(wf_L*(1-k0)+am-c_);
      (1+phi*(x(2)^(1-eta)))*(pb+pe)*x(2)+(k2+phi*(x(2)^(1-eta)))*wf_L*k1*((x(2))^(k2))-(phi*(x(2)^(1-eta)))*(wf_L*(1-k0)+am-c_);
      (1+phi*(x(3)^(1-eta)))*(pb+pe)*x(3)+(k2+phi*(x(3)^(1-eta)))*wf_H*k1*((x(3))^(k2))-(phi*(x(3)^(1-eta)))*(wf_H*(1-k0)+am-c_);
     ];



