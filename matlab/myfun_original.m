function F = myfun_original (x,pb,pe,am,c_,k0,k1,k2,N_L,N_H,rho,b_total,tf_L,tf_H,wf_L,wf_H) % 
%parameter: pb pe am c_ k0 k1 k2 N_L N_H
%assume: rho
%data
%Target: b_total tf_L tf_H wf_L wf_H
%solve: bd_L be_L b_H R alpha gamma A pi 
F = [ (1+x(6)+x(5))*pb*x(1)+((1+x(6))*k2+x(5))*wf_L*k1*(x(1))^(k2)-x(5)*(wf_L*(1-k0)+am-c_);
      (1+x(6)+x(5))*(pb+pe)*x(2)+((1+x(6))*k2+x(5))*wf_L*k1*(x(2))^(k2)-x(5)*(wf_L*(1-k0)+am-c_);
      (1+x(6)+x(5))*(pb+pe)*x(3)+((1+x(6))*k2+x(5))*wf_H*k1*(x(3))^(k2)-x(5)*(wf_H*(1-k0)+am-c_);
      b_total-(N_L*(x(4)*x(1)+(1-x(4))*x(2))+N_H*x(3));
      tf_L-(x(4)*(1-k0-((x(6)*pb*x(1)+(x(6)*k2+x(5))*wf_L*k1*(x(1))^(k2))/(x(5)*wf_L)))+(1-x(4))*(1-k0-((x(6)*(pb+pe)*x(2)+(x(6)*k2+x(5))*wf_L*k1*(x(2))^(k2))/(x(5)*wf_L))));
      tf_H-(1-k0-((x(6)*(pb+pe)*x(3)+(x(6)*k2+x(5))*wf_H*k1*(x(3))^(k2))/(x(5)*wf_H)));
      wf_L-x(7)*x(8)*(N_L*(x(4)*(1-k0-((x(6)*pb*x(1)+(x(6)*k2+x(5))*wf_L*k1*(x(1))^(k2))/(x(5)*wf_L)))+(1-x(4))*(1-k0-((x(6)*(pb+pe)*x(2)+(x(6)*k2+x(5))*wf_L*k1*(x(2))^(k2))/(x(5)*wf_L)))))^(rho-1)*(x(8)*(N_L*(x(4)*(1-k0-((x(6)*pb*x(1)+(x(6)*k2+x(5))*wf_L*k1*(x(1))^(k2))/(x(5)*wf_L)))+(1-x(4))*(1-k0-((x(6)*(pb+pe)*x(2)+(x(6)*k2+x(5))*wf_L*k1*(x(2))^(k2))/(x(5)*wf_L)))))^rho+(1-x(8))*(N_H*(1-k0-((x(6)*(pb+pe)*x(3)+(x(6)*k2+x(5))*wf_H*k1*(x(3))^(k2))/(x(5)*wf_H))))^rho)^(1/rho-1);
      wf_H-x(7)*(1-x(8))*(N_H*(1-k0-((x(6)*(pb+pe)*x(3)+(x(6)*k2+x(5))*wf_H*k1*(x(3))^(k2))/(x(5)*wf_H))))^(rho-1)*(x(8)*(N_L*(x(4)*(1-k0-((x(6)*pb*x(1)+(x(6)*k2+x(5))*wf_L*k1*(x(1))^(k2))/(x(5)*wf_L)))+(1-x(4))*(1-k0-((x(6)*(pb+pe)*x(2)+(x(6)*k2+x(5))*wf_L*k1*(x(2))^(k2))/(x(5)*wf_L)))))^rho+(1-x(8))*(N_H*(1-k0-((x(6)*(pb+pe)*x(3)+(x(6)*k2+x(5))*wf_H*k1*(x(3))^(k2))/(x(5)*wf_H))))^rho)^(1/rho-1);
     ];

