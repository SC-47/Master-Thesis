function F = myfun0(x,pb,pe,am,c_,k0,k1,k2,N_L,N_H,b_L,b_H,wf_L,wf_H,R) % x = [alpha beta gamma be_L pi rho]
%known parameter : pb pe am c_ k0 k1 k2 N_L N_H
%target : b_L b_H wf_L wf_H R              
%calibration parameter : alpha beta gamma pi rho
F = [ (1+x(3)+x(1))*pb*((b_L*(1-R)*x(4))/R)+((1+x(3))*k2+x(1))*wf_L*k1*(((b_L*(1-R)*x(4))/R))^(k2)-x(1)*(wf_L*(1-k0)+am-c_);
      (1+x(3)+x(1))*(pb+pe)*x(4)+((1+x(3))*k2+x(1))*wf_L*k1*(x(4))^(k2)-x(1)*(wf_L*(1-k0)+am-c_);
      (1+x(3)+x(1))*(pb+pe)*b_H+((1+x(3))*k2+x(1))*wf_H*k1*(b_H)^(k2)-x(1)*(wf_H*(1-k0)+am-c_);
      (1+x(3))*log((pb*((b_L*(1-R)*x(4))/R)+wf_L*k1*k2*(((b_L*(1-R)*x(4))/R))^(k2))/((pb+pe)*x(4)+wf_L*k1*k2*(x(4))^(k2)))+x(1)*log(((b_L*(1-R)*x(4))/R)/x(4))+x(2)*log(wf_L/wf_H)
      (wf_L/wf_H)-(x(5)/(1-x(5)))*((N_L*(R*(1-k0-((x(3)*pb*((b_L*(1-R)*x(4))/R)+(x(3)*k2+x(1))*wf_L*k1*(((b_L*(1-R)*x(4))/R))^(k2))/(x(1)*wf_L)))+(1-R)*(1-k0-((x(3)*(pb+pe)*x(4)+(x(3)*k2+x(1))*wf_L*k1*(x(4))^(k2))/(x(1)*wf_L)))))/(N_H*(1-k0-((x(3)*(pb+pe)*b_H+(x(3)*k2+x(1))*wf_H*k1*(b_H)^(k2))/(x(1)*wf_H)))))^(x(6)-1)
      wf_H-(1-x(5))*(N_H*(1-k0-((x(3)*(pb+pe)*b_H+(x(3)*k2+x(1))*wf_H*k1*(b_H)^(k2))/(x(1)*wf_H))))^(x(6)-1)*(x(5)*(N_L*(R*(1-k0-((x(3)*pb*((b_L*(1-R)*x(4))/R)+(x(3)*k2+x(1))*wf_L*k1*(((b_L*(1-R)*x(4))/R))^(k2))/(x(1)*wf_L)))+(1-R)*(1-k0-((x(3)*(pb+pe)*x(4)+(x(3)*k2+x(1))*wf_L*k1*(x(4))^(k2))/(x(1)*wf_L)))))+(1-x(5))*(N_H*(1-k0-((x(3)*(pb+pe)*b_H+(x(3)*k2+x(1))*wf_H*k1*(b_H)^(k2))/(x(1)*wf_H)))))
     ];
