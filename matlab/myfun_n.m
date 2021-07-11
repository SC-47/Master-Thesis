function F = myfun_n(x,pb,pe,am,c_,k0,k1,k2,wf_L,wf_H,tf_L,tf_H,b_L,b_H,N_L,N_H,rho) 
%known parameter : pb pe am c_ k0 k1 k2 N_L N_H 
%assumption paramter : rho 
%target bd_L be_L R gamma alpha eta beta A pi    
F = [ (1+x(4)+x(5))*pb*x(1)+((1+x(4))*k2+x(5))*wf_L*k1*(x(1))^(k2)-x(5)*(wf_L*(1-k0)+am-c_);
      (1+x(4)+x(5))*(pb+pe)*x(2)+((1+x(4))*k2+x(5))*wf_L*k1*(x(2))^(k2)-x(5)*(wf_L*(1-k0)+am-c_);
      (1+x(4)+x(5))*(pb+pe)*b_H+((1+x(4))*k2+x(5))*wf_H*k1*(b_H)^(k2)-x(5)*(wf_H*(1-k0)+am-c_);
      b_L-(x(3)*x(1)+(1-x(3))*x(2));
      (1+x(4))*log((pb*x(1)+wf_L*k1*k2*(x(1))^(k2))/((pb+pe)*x(2)+wf_L*k1*k2*(x(2))^(k2)))+x(5)*log(x(1)/x(2))+x(7)*log(wf_L/wf_H)-x(6);
      tf_L-(x(3)*(1-k0-((x(4)*pb*x(1)+(x(4)*k2+x(5))*wf_L*k1*x(1)^k2)/(x(5)*wf_L)))+(1-x(3))*(1-k0-((x(4)*(pb+pe)*x(2)+(x(4)*k2+x(5))*wf_L*k1*x(2)^k2)/(x(5)*wf_L))));
      tf_H-(1-k0-((x(4)*(pb+pe)*b_H+(x(4)*k2+x(5))*wf_H*k1*b_H^k2)/(x(5)*wf_H)));
      wf_L-x(8)*x(9)*(N_L*(x(3)*(1-k0-((x(4)*pb*x(1)+(x(4)*k2+x(5))*wf_L*k1*x(1)^k2)/(x(5)*wf_L)))+(1-x(3))*(1-k0-((x(4)*(pb+pe)*x(2)+(x(4)*k2+x(5))*wf_L*k1*x(2)^k2)/(x(5)*wf_L)))))^(rho-1)*(x(9)*(N_L*(x(3)*(1-k0-((x(4)*pb*x(1)+(x(4)*k2+x(5))*wf_L*k1*x(1)^k2)/(x(5)*wf_L)))+(1-x(3))*(1-k0-((x(4)*(pb+pe)*x(2)+(x(4)*k2+x(5))*wf_L*k1*x(2)^k2)/(x(5)*wf_L)))))^rho+(1-x(9))*(N_H*(1-k0-((x(4)*(pb+pe)*b_H+(x(4)*k2+x(5))*wf_H*k1*b_H^k2)/(x(5)*wf_H))))^rho)^(1/rho-1);
      wf_H-x(8)*(1-x(9))*(N_H*(1-k0-((x(4)*(pb+pe)*b_H+(x(4)*k2+x(5))*wf_H*k1*b_H^k2)/(x(5)*wf_H))))^(rho-1)*(x(9)*(N_L*(x(3)*(1-k0-((x(4)*pb*x(1)+(x(4)*k2+x(5))*wf_L*k1*x(1)^k2)/(x(5)*wf_L)))+(1-x(3))*(1-k0-((x(4)*(pb+pe)*x(2)+(x(4)*k2+x(5))*wf_L*k1*x(2)^k2)/(x(5)*wf_L)))))^rho+(1-x(9))*(N_H*(1-k0-((x(4)*(pb+pe)*b_H+(x(4)*k2+x(5))*wf_H*k1*b_H^k2)/(x(5)*wf_H))))^rho)^(1/rho-1);
    ];
