function F = myfun_n_c (x,pb,pe,am,c_,k0,k1,k2,N_L,N_H,rho,alpha,eta,gamma,beta,A,pi) % 
%known parameter : pb pe am c_ k0 k1 k2 N_L N_H
%assumption paramter : rho 
%caliabration parameter : alpha eta gamma beta A pi
%solve parameter: bd_L be_L b_H R wf_L wf_H     
F = [ (1+gamma+alpha)*pb*x(1)+((1+gamma)*k2+alpha)*x(5)*k1*(x(1))^(k2)-alpha*(x(5)*(1-k0)+am-c_);
      (1+gamma+alpha)*(pb+pe)*x(2)+((1+gamma)*k2+alpha)*x(5)*k1*(x(2))^(k2)-alpha*(x(5)*(1-k0)+am-c_);
      (1+gamma+alpha)*(pb+pe)*x(3)+((1+gamma)*k2+alpha)*x(6)*k1*(x(3))^(k2)-alpha*(x(6)*(1-k0)+am-c_);
      b_L-(x(4)*x(1)+(1-x(4))*x(2));
      (1+gamma)*log((pb*x(1)+x(5)*k1*k2*(x(1))^(k2))/((pb+pe)*x(2)+x(5)*k1*k2*(x(2))^(k2)))+alpha*log(x(1)/x(2))+beta*log(x(5)/x(6))-eta;
      tf_L-(x(4)*(1-k0-((gamma*pb*x(1)+(gamma*k2+alpha)*x(5)*k1*x(1)^k2)/(alpha*x(5))))+(1-x(4))*(1-k0-((gamma*(pb+pe)*x(2)+(gamma*k2+alpha)*x(5)*k1*x(2)^k2)/(alpha*x(5)))));
      tf_H-(1-k0-((gamma*(pb+pe)*x(3)+(gamma*k2+alpha)*x(6)*k1*x(3)^k2)/(alpha*x(6))));
      x(5)-A*pi*(N_L*(x(4)*(1-k0-((gamma*pb*x(1)+(gamma*k2+alpha)*x(5)*k1*x(1)^k2)/(alpha*x(5))))+(1-x(4))*(1-k0-((gamma*(pb+pe)*x(2)+(gamma*k2+alpha)*x(5)*k1*x(2)^k2)/(alpha*x(5))))))^(rho-1)*(pi*(N_L*(x(4)*(1-k0-((gamma*pb*x(1)+(gamma*k2+alpha)*x(5)*k1*x(1)^k2)/(alpha*x(5))))+(1-x(4))*(1-k0-((gamma*(pb+pe)*x(2)+(gamma*k2+alpha)*x(5)*k1*x(2)^k2)/(alpha*x(5))))))^rho+(1-pi)*(N_H*(1-k0-((gamma*(pb+pe)*x(3)+(gamma*k2+alpha)*x(6)*k1*x(3)^k2)/(alpha*x(6)))))^rho)^(1/rho-1);
      x(6)-A*(1-pi)*(N_H*(1-k0-((gamma*(pb+pe)*x(3)+(gamma*k2+alpha)*x(6)*k1*x(3)^k2)/(alpha*x(6)))))^(rho-1)*(pi*(N_L*(x(4)*(1-k0-((gamma*pb*x(1)+(gamma*k2+alpha)*x(5)*k1*x(1)^k2)/(alpha*x(5))))+(1-x(4))*(1-k0-((gamma*(pb+pe)*x(2)+(gamma*k2+alpha)*x(5)*k1*x(2)^k2)/(alpha*x(5))))))^rho+(1-pi)*(N_H*(1-k0-((gamma*(pb+pe)*x(3)+(gamma*k2+alpha)*x(6)*k1*x(3)^k2)/(alpha*x(6)))))^rho)^(1/rho-1);
    ];

