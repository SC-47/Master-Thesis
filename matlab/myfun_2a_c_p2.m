function F = myfun_2a_c_p2 (x,pb,pe,am,c_,k0,k1,k2,N_L,N_H,rho,alpha_L,alpha_H,gamma,beta,A,pi,Q) 
%known parameter : pb pe am c_ k0 k1 k2 N_L N_H rho alpha_L alpha_H gamma beta A pi
%assumption paramter :rho
%target bd_L be_L b_H R wf_L wf_H             
F = [ (1+gamma+alpha_L)*(pb*(1-Q))*x(1)+((1+gamma)*k2+alpha_L)*x(5)*k1*(x(1))^(k2)-alpha_L*(x(5)*(1-k0)+(am-(N_L*((1-x(4))*x(1)+x(4)*x(2)))*pb*Q)-c_);
      (1+gamma+alpha_L)*(pb*(1-Q)+pe)*x(2)+((1+gamma)*k2+alpha_L)*x(5)*k1*(x(2))^(k2)-alpha_L*(x(5)*(1-k0)+(am-(N_L*((1-x(4))*x(1)+x(4)*x(2)))*pb*Q)-c_);
      (1+gamma+alpha_H)*(pb+pe)*x(3)+((1+gamma)*k2+alpha_H)*x(6)*k1*(x(3))^(k2)-alpha_H*(x(6)*(1-k0)+(am-(N_L*((1-x(4))*x(1)+x(4)*x(2)))*pb*Q)-c_);
      (1+gamma)*log((pb*(1-Q)*x(1)+x(5)*k1*k2*(x(1))^(k2))/((pb*(1-Q)+pe)*x(2)+x(5)*k1*k2*(x(2))^(k2)))+alpha_L*log(x(1)/x(2))+beta*log(x(5)/x(6))
      x(5)-(A*(pi*(N_L*((1-x(4))*(1-k0-((gamma*(pb*(1-Q))*x(1)+(gamma*k2+alpha_L)*x(5)*k1*x(1)^k2)/(alpha_L*x(5))))+x(4)*(1-k0-((gamma*((pb*(1-Q))+pe)*x(2)+(gamma*k2+alpha_L)*x(5)*k1*x(2)^k2)/(alpha_L*x(5))))))^rho+(1-pi)*(N_H*(1-k0-((gamma*(pb+pe)*x(3)+(gamma*k2+alpha_H)*x(6)*k1*x(3)^k2)/(alpha_H*x(6)))))^rho)^(1/rho-1)*pi*(N_L*((1-x(4))*(1-k0-((gamma*(pb*(1-Q))*x(1)+(gamma*k2+alpha_L)*x(5)*k1*x(1)^k2)/(alpha_L*x(5))))+x(4)*(1-k0-((gamma*((pb*(1-Q))+pe)*x(2)+(gamma*k2+alpha_L)*x(5)*k1*x(2)^k2)/(alpha_L*x(5))))))^(rho-1)); 
      x(6)-(A*(pi*(N_L*((1-x(4))*(1-k0-((gamma*(pb*(1-Q))*x(1)+(gamma*k2+alpha_L)*x(5)*k1*x(1)^k2)/(alpha_L*x(5))))+x(4)*(1-k0-((gamma*((pb*(1-Q))+pe)*x(2)+(gamma*k2+alpha_L)*x(5)*k1*x(2)^k2)/(alpha_L*x(5))))))^rho+(1-pi)*(N_H*(1-k0-((gamma*(pb+pe)*x(3)+(gamma*k2+alpha_H)*x(6)*k1*x(3)^k2)/(alpha_H*x(6)))))^rho)^(1/rho-1)*(1-pi)*(N_H*(1-k0-((gamma*(pb+pe)*x(3)+(gamma*k2+alpha_H)*x(6)*k1*x(3)^k2)/(alpha_H*x(6)))))^(rho-1))
     ];




