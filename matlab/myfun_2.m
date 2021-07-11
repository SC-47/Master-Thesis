function F = myfun_2(x,pb,pe,am,c_,k0,k1,k2,wf_L,wf_H,tf_L,tf_H,b_L,b_H,N_L,N_H,rho)
%known parameter : pb pe am c_ k0 k1 k2 N_L N_H pi R
%assumption paramter : rho 
%solve: bd_L be_L gamma R phi eta A pi
F = [ (1+x(3)+(x(5)*x(1)^(1-x(6))))*pb*x(1)+((1+x(3))*k2+(x(5)*x(1)^(1-x(6))))*wf_L*k1*(x(1))^(k2)-(x(5)*x(1)^(1-x(6)))*(wf_L*(1-k0)+am-c_);
      (1+x(3)+(x(5)*x(2)^(1-x(6))))*(pb+pe)*x(2)+((1+x(3))*k2+(x(5)*x(2)^(1-x(6))))*wf_L*k1*(x(2))^(k2)-(x(5)*x(2)^(1-x(6)))*(wf_L*(1-k0)+am-c_);
      (1+x(3)+(x(5)*b_H^(1-x(6))))*(pb+pe)*b_H+((1+x(3))*k2+(x(5)*b_H^(1-x(6))))*wf_H*k1*(b_H)^(k2)-(x(5)*b_H^(1-x(6)))*(wf_H*(1-k0)+am-c_);     
      tf_L-(x(4)*(1-k0-((x(3)*pb*x(1)+(x(3)*k2+(x(5)*x(1)^(1-x(6))))*wf_L*k1*(x(1))^(k2))/(wf_L*(x(5)*x(1)^(1-x(6))))))+(1-x(4))*(1-k0-((x(3)*(pb+pe)*x(2)+(x(3)*k2+(x(5)*x(2)^(1-x(6))))*wf_L*k1*(x(2))^(k2))/(wf_L*(x(5)*x(2)^(1-x(6)))))));
      tf_H-(1-k0-((x(3)*(pb+pe)*b_H+(x(3)*k2+(x(5)*b_H^(1-x(6))))*wf_H*k1*(b_H)^(k2))/(wf_H*(x(5)*b_H^(1-x(6))))));
      b_L-(x(4)*x(1)+(1-x(4))*x(2));
      wf_L-x(7)*(x(8))*(N_L*(x(4)*(1-k0-((x(3)*pb*x(1)+(x(3)*k2+(x(5)*x(1)^(1-x(6))))*wf_L*k1*(x(1))^(k2))/(wf_L*(x(5)*x(1)^(1-x(6))))))+(1-x(4))*(1-k0-((x(3)*(pb+pe)*x(2)+(x(3)*k2+(x(5)*x(2)^(1-x(6))))*wf_L*k1*(x(2))^(k2))/(wf_L*(x(5)*x(2)^(1-x(6))))))))^(rho-1)*(x(8)*(N_L*(x(4)*(1-k0-((x(3)*pb*x(1)+(x(3)*k2+(x(5)*x(1)^(1-x(6))))*wf_L*k1*(x(1))^(k2))/(wf_L*(x(5)*x(1)^(1-x(6))))))+(1-x(4))*(1-k0-((x(3)*(pb+pe)*x(2)+(x(3)*k2+(x(5)*x(2)^(1-x(6))))*wf_L*k1*(x(2))^(k2))/(wf_L*(x(5)*x(2)^(1-x(6))))))))^(rho)+(1-x(8))*(N_H*(1-k0-((x(3)*(pb+pe)*x(3)+(x(3)*k2+(x(5)*x(3)^(1-x(6))))*wf_H*k1*(x(3))^(k2))/(wf_H*(x(5)*x(3)^(1-x(6)))))))^rho)^(1/rho-1);
      wf_H-x(7)*(1-x(8))*(N_H*(1-k0-((x(3)*(pb+pe)*x(3)+(x(3)*k2+(x(5)*x(3)^(1-x(6))))*wf_H*k1*(x(3))^(k2))/(wf_H*(x(5)*x(3)^(1-x(6)))))))^(rho-1)*(x(8)*(N_L*(x(4)*(1-k0-((x(3)*pb*x(1)+(x(3)*k2+(x(5)*x(1)^(1-x(6))))*wf_L*k1*(x(1))^(k2))/(wf_L*(x(5)*x(1)^(1-x(6))))))+(1-x(4))*(1-k0-((x(3)*(pb+pe)*x(2)+(x(3)*k2+(x(5)*x(2)^(1-x(6))))*wf_L*k1*(x(2))^(k2))/(wf_L*(x(5)*x(2)^(1-x(6))))))))^(rho)+(1-x(8))*(N_H*(1-k0-((x(3)*(pb+pe)*x(3)+(x(3)*k2+(x(5)*x(3)^(1-x(6))))*wf_H*k1*(x(3))^(k2))/(wf_H*(x(5)*x(3)^(1-x(6)))))))^rho)^(1/rho-1);
    ];












