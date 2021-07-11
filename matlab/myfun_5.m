function F = myfun_5(x,pb,pe,am,c_,k0,k1,k2,wf_L,wf_H,b_L,b_H,R)
%known parameter : pb pe am c_ k0 k1 k2 N_L N_H
%assumption paramter : rho 
%solve parameter: bd_L be_L b_H phi eta
F = [ (1+x(4)*(x(1)^(1-x(5))))*pb*x(1)+(k2+x(4)*(x(1)^(1-x(5))))*wf_L*k1*((x(1))^(k2))-(x(4)*(x(1)^(1-x(5))))*(wf_L*(1-k0)+am-c_);
      (1+x(4)*(x(2)^(1-x(5))))*(pb+pe)*x(2)+(k2+x(4)*(x(2)^(1-x(5))))*wf_L*k1*((x(2))^(k2))-(x(4)*(x(2)^(1-x(5))))*(wf_L*(1-k0)+am-c_);
      (1+x(4)*(x(3)^(1-x(5))))*(pb+pe)*x(3)+(k2+x(4)*(x(3)^(1-x(5))))*wf_H*k1*((x(3))^(k2))-(x(4)*(x(3)^(1-x(5))))*(wf_H*(1-k0)+am-c_);
      b_L-(R*x(1)+(1-R)*x(2));
      b_H-x(3);
     ];


