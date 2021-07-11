function F = myfun_3_c_p (x,pb,pe,am,c_,k0,k1,k2,N_L,N_H,rho,phi,eta,gamma,beta,A,pi,Q) % 
%known parameter : pb pe am c_ k0 k1 k2 N_L N_H
%assumption paramter : rho 
%caliabration parameter : phi eta gamma beta A pi
%solve parameter: bd_L be_L b_H R wf_L wf_H 
F = [ (1+gamma+(phi*x(1)^(1-eta)))*pb*(1-Q)*x(1)+((1+gamma)*k2+(phi*x(1)^(1-eta)))*x(5)*k1*(x(1))^(k2)-(phi*x(1)^(1-eta))*(x(5)*(1-k0)+(am-(N_L*((x(4)*x(1)+(1-x(4))*x(2)))+N_H*x(3))*pb*Q)-c_) ;
      (1+gamma+(phi*x(2)^(1-eta)))*(pb*(1-Q)+pe)*x(2)+((1+gamma)*k2+(phi*x(2)^(1-eta)))*x(5)*k1*(x(2))^(k2)-(phi*x(2)^(1-eta))*(x(5)*(1-k0)+(am-(N_L*((x(4)*x(1)+(1-x(4))*x(2)))+N_H*x(3))*pb*Q)-c_) ;
      (1+gamma+(phi*x(3)^(1-eta)))*(pb*(1-Q)+pe)*x(3)+((1+gamma)*k2+(phi*x(3)^(1-eta)))*x(6)*k1*(x(3))^(k2)-(phi*x(3)^(1-eta))*(x(6)*(1-k0)+(am-(N_L*((x(4)*x(1)+(1-x(4))*x(2)))+N_H*x(3))*pb*Q)-c_) ;
      (1+gamma)*log(((x(2)/x(1))^(1-eta))*(pb*(1-Q)*x(1)+x(5)*k1*k2*(x(1))^(k2))/((pb*(1-Q)+pe)*x(2)+x(5)*k1*k2*(x(2))^(k2)))+(phi/(1-eta))*(x(1)^(1-eta)-x(2)^(1-eta))+beta*log(x(5)/x(6));
      x(5)-A*(pi)*(N_L*(x(4)*(1-k0-((gamma*pb*(1-Q)*x(1)+(gamma*k2+(phi*x(1)^(1-eta)))*x(5)*k1*(x(1))^(k2))/(x(5)*(phi*x(1)^(1-eta)))))+(1-x(4))*(1-k0-((gamma*(pb*(1-Q)+pe)*x(2)+(gamma*k2+(phi*x(2)^(1-eta)))*x(5)*k1*(x(2))^(k2))/(x(5)*(phi*x(2)^(1-eta)))))))^(rho-1)*(pi*(N_L*(x(4)*(1-k0-((gamma*pb*(1-Q)*x(1)+(gamma*k2+(phi*x(1)^(1-eta)))*x(5)*k1*(x(1))^(k2))/(x(5)*(phi*x(1)^(1-eta)))))+(1-x(4))*(1-k0-((gamma*(pb*(1-Q)+pe)*x(2)+(gamma*k2+(phi*x(2)^(1-eta)))*x(5)*k1*(x(2))^(k2))/(x(5)*(phi*x(2)^(1-eta)))))))^(rho)+(1-pi)*(N_H*(1-k0-((gamma*(pb*(1-Q)+pe)*x(3)+(gamma*k2+(phi*x(3)^(1-eta)))*x(6)*k1*(x(3))^(k2))/(x(6)*(phi*x(3)^(1-eta))))))^rho)^(1/rho-1);
      x(6)-A*(1-pi)*(N_H*(1-k0-((gamma*(pb*(1-Q)+pe)*x(3)+(gamma*k2+(phi*x(3)^(1-eta)))*x(6)*k1*(x(3))^(k2))/(x(6)*(phi*x(3)^(1-eta))))))^(rho-1)*(pi*(N_L*(x(4)*(1-k0-((gamma*pb*(1-Q)*x(1)+(gamma*k2+(phi*x(1)^(1-eta)))*x(5)*k1*(x(1))^(k2))/(x(5)*(phi*x(1)^(1-eta)))))+(1-x(4))*(1-k0-((gamma*(pb*(1-Q)+pe)*x(2)+(gamma*k2+(phi*x(2)^(1-eta)))*x(5)*k1*(x(2))^(k2))/(x(5)*(phi*x(2)^(1-eta)))))))^(rho)+(1-pi)*(N_H*(1-k0-((gamma*(pb*(1-Q)+pe)*x(3)+(gamma*k2+(phi*x(3)^(1-eta)))*x(6)*k1*(x(3))^(k2))/(x(6)*(phi*x(3)^(1-eta))))))^rho)^(1/rho-1);
     ];
