function F = myfun0(x,pb,pe,am,c_,k0,k1,k2,N_L,N_H,rho,pi,alpha,beta,gamma,A) % 
%known parameter : pb pe am c_ k0 k1 k2 N_L N_H
%assumption paramter : rho
%caliabration parameter : pi alpha beta gamma A
%endogeneous variable: bd_L be_L b_H wf_L wf_H R
F = [ (1+gamma+alpha)*pb*x(1)+((1+gamma)*k2+alpha)*x(4)*k1*(x(1))^(k2)-alpha*(x(4)*(1-k0)+am-c_);
      (1+gamma+alpha)*(pb+pe)*x(2)+((1+gamma)*k2+alpha)*x(4)*k1*(x(2))^(k2)-alpha*(x(4)*(1-k0)+am-c_);
      (1+gamma+alpha)*(pb+pe)*x(3)+((1+gamma)*k2+alpha)*x(5)*k1*(x(3))^(k2)-alpha*(x(5)*(1-k0)+am-c_);
      (1+gamma)*log((pb*x(1)+x(4)*k1*k2*(x(1))^(k2))/((pb+pe)*x(2)+x(4)*k1*k2*(x(2))^(k2)))+alpha*log(x(1)/x(2))+beta*log(x(4)/x(5));
      (x(4)/x(5))-(pi/(1-pi))*((N_L*(x(6)*(1-k0-((gamma*pb*x(1)+(gamma*k2+alpha)*x(4)*k1*(x(1))^(k2))/(alpha*x(4))))+(1-x(6))*(1-k0-((gamma*(pb+pe)*x(2)+(gamma*k2+alpha)*x(4)*k1*(x(2))^(k2))/(alpha*x(4))))))/(N_H*(1-k0-((gamma*(pb+pe)*x(3)+(gamma*k2+alpha)*x(5)*k1*(x(3))^(k2))/(alpha*x(5))))))^(rho-1);
      x(5)-(A*(1-pi)*((N_H*(1-k0-((gamma*(pb+pe)*x(3)+(gamma*k2+alpha)*x(5)*k1*(x(3))^(k2))/(alpha*x(5)))))^(rho-1))*(pi*(N_L*(x(6)*(1-k0-((gamma*pb*x(1)+(gamma*k2+alpha)*x(4)*k1*(x(1))^(k2))/(alpha*x(4))))+(1-x(6))*(1-k0-((gamma*(pb+pe)*x(2)+(gamma*k2+alpha)*x(4)*k1*(x(2))^(k2))/(alpha*x(4))))))^(rho)+(1-pi)*(N_H*(1-k0-((gamma*(pb+pe)*x(3)+(gamma*k2+alpha)*x(5)*k1*(x(3))^(k2))/(alpha*x(5)))))^(rho))^((1-rho)/rho));
     ];
 
 


