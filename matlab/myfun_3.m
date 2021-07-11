%function F = myfun_3(x,pb,pe,am,c_,k0,k1,k2,wf_L,wf_H,R,tf_L,tf_H,b_L) 
%function F = myfun_3(x,pb,pe,am,k0,k1,k2,wf_L,wf_H,R,tf_L,tf_H,b_H) 
function F = myfun_3(x,pb,pe,am,c_,k0,k1,k2,wf_L,wf_H,R,tf_L,tf_H,b_total,N_L,N_H,rho) 
%known parameter : pb pe am c_ k0 k1 k2 N_L N_H pi R
%assumption paramter : rho 
%solve: bd_L be_L b_H gamma phi eta A pi
F = [ (1+x(4)+(x(5)*x(1)^(1-x(6))))*pb*x(1)+((1+x(4))*k2+(x(5)*x(1)^(1-x(6))))*wf_L*k1*(x(1))^(k2)-(x(5)*x(1)^(1-x(6)))*(wf_L*(1-k0)+am-c_);
      (1+x(4)+(x(5)*x(2)^(1-x(6))))*(pb+pe)*x(2)+((1+x(4))*k2+(x(5)*x(2)^(1-x(6))))*wf_L*k1*(x(2))^(k2)-(x(5)*x(2)^(1-x(6)))*(wf_L*(1-k0)+am-c_);
      (1+x(4)+(x(5)*x(3)^(1-x(6))))*(pb+pe)*x(3)+((1+x(4))*k2+(x(5)*x(3)^(1-x(6))))*wf_H*k1*(x(3))^(k2)-(x(5)*x(3)^(1-x(6)))*(wf_H*(1-k0)+am-c_);     
      tf_L-(R*(1-k0-((x(4)*pb*x(1)+(x(4)*k2+(x(5)*x(1)^(1-x(6))))*wf_L*k1*(x(1))^(k2))/(wf_L*(x(5)*x(1)^(1-x(6))))))+(1-R)*(1-k0-((x(4)*(pb+pe)*x(2)+(x(4)*k2+(x(5)*x(2)^(1-x(6))))*wf_L*k1*(x(2))^(k2))/(wf_L*(x(5)*x(2)^(1-x(6)))))));
      tf_H-(1-k0-((x(4)*(pb+pe)*x(3)+(x(4)*k2+(x(5)*x(3)^(1-x(6))))*wf_H*k1*(x(3))^(k2))/(wf_H*(x(5)*x(3)^(1-x(6))))));
      %b_L-(R*x(1)+(1-R)*x(2));
      %b_H-x(3);
      b_total-(N_L*(R*x(1)+(1-R)*x(2))+N_H*x(3));
      wf_L-x(7)*x(8)*(((N_L*(R*(1-k0-((x(4)*pb*x(1)+(x(4)*k2+(x(5)*x(1)^(1-x(6))))*wf_L*k1*(x(1))^(k2))/(wf_L*(x(5)*x(1)^(1-x(6))))))+(1-R)*(1-k0-((x(4)*(pb+pe)*x(2)+(x(4)*k2+(x(5)*x(2)^(1-x(6))))*wf_L*k1*(x(2))^(k2))/(wf_L*(x(5)*x(2)^(1-x(6)))))))))^(rho-1))*(x(8)*((N_L*(R*(1-k0-((x(4)*pb*x(1)+(x(4)*k2+(x(5)*x(1)^(1-x(6))))*wf_L*k1*(x(1))^(k2))/(wf_L*(x(5)*x(1)^(1-x(6))))))+(1-R)*(1-k0-((x(4)*(pb+pe)*x(2)+(x(4)*k2+(x(5)*x(2)^(1-x(6))))*wf_L*k1*(x(2))^(k2))/(wf_L*(x(5)*x(2)^(1-x(6)))))))))^(rho)+(1-x(8))*((N_H*(1-k0-((x(4)*(pb+pe)*x(3)+(x(4)*k2+(x(5)*x(3)^(1-x(6))))*wf_H*k1*(x(3))^(k2))/(wf_H*(x(5)*x(3)^(1-x(6))))))))^(rho))^((1-rho)/rho) ;
      wf_H-x(7)*(1-x(8))*(((N_H*(1-k0-((x(4)*(pb+pe)*x(3)+(x(4)*k2+(x(5)*x(3)^(1-x(6))))*wf_H*k1*(x(3))^(k2))/(wf_H*(x(5)*x(3)^(1-x(6))))))))^(rho-1))*(x(8)*((N_L*(R*(1-k0-((x(4)*pb*x(1)+(x(4)*k2+(x(5)*x(1)^(1-x(6))))*wf_L*k1*(x(1))^(k2))/(wf_L*(x(5)*x(1)^(1-x(6))))))+(1-R)*(1-k0-((x(4)*(pb+pe)*x(2)+(x(4)*k2+(x(5)*x(2)^(1-x(6))))*wf_L*k1*(x(2))^(k2))/(wf_L*(x(5)*x(2)^(1-x(6)))))))))^(rho)+(1-x(8))*((N_H*(1-k0-((x(4)*(pb+pe)*x(3)+(x(4)*k2+(x(5)*x(3)^(1-x(6))))*wf_H*k1*(x(3))^(k2))/(wf_H*(x(5)*x(3)^(1-x(6))))))))^(rho))^((1-rho)/rho) ;
      ];