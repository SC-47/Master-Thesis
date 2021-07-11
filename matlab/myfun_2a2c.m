function F = myfun_2a2c(x,alpha_L,alpha_H,am,gamma,k0,k1,k2,wf_L,wf_H,b_L,b_H,N_L,N_H,rho,tf_L,tf_H) 
%known parameter : alpha_L alpha_H am gamma k0 k1 k2 wf_L wf_H b_L b_H tf_L tf_H N_L N_H 
%assumption paramter : rho 
%solve: bd_L be_L R pb pe c_ A pi
F = [ (1+gamma+alpha_L)*x(4)*x(1)+((1+gamma)*k2+alpha_L)*wf_L*k1*(x(1))^(k2)-alpha_L*(wf_L*(1-k0)+am-x(6)) ;
      (1+gamma+alpha_L)*(x(4)+x(5))*x(2)+((1+gamma)*k2+alpha_L)*wf_L*k1*(x(2))^(k2)-alpha_L*(wf_L*(1-k0)+am-x(6)) ;
      (1+gamma+alpha_H)*(x(4)+x(5))*b_H+((1+gamma)*k2+alpha_H)*wf_H*k1*(b_H)^(k2)-alpha_H*(wf_H*(1-k0)+am-x(6)) ;
      b_L-((1-x(3))*x(1)+x(3)*x(2)) ;
      tf_L-((1-x(3))*(1-k0-((gamma*x(4)*x(1)+(gamma*k2+alpha_L)*wf_L*k1*x(1)^k2)/(alpha_L*wf_L)))+x(3)*(1-k0-((gamma*(x(4)+x(5))*x(2)+(gamma*k2+alpha_L)*wf_L*k1*x(2)^k2)/(alpha_L*wf_L)))) ;
      tf_H-(1-k0-((gamma*(x(4)+x(5))*b_H+(gamma*k2+alpha_H)*wf_H*k1*b_H^k2)/(alpha_H*wf_H))) ;
      wf_L-(x(7)*(x(8)*(N_L*((1-x(3))*(1-k0-((gamma*x(4)*x(1)+(gamma*k2+alpha_L)*wf_L*k1*x(1)^k2)/(alpha_L*wf_L)))+x(3)*(1-k0-((gamma*(x(4)+x(5))*x(2)+(gamma*k2+alpha_L)*wf_L*k1*x(2)^k2)/(alpha_L*wf_L)))))^rho+(1-x(8))*(N_H*(1-k0-((gamma*(x(4)+x(5))*b_H+(gamma*k2+alpha_H)*wf_H*k1*b_H^k2)/(alpha_H*wf_H))))^rho)^(1/rho-1)*x(8)*(N_L*((1-x(3))*(1-k0-((gamma*x(4)*x(1)+(gamma*k2+alpha_L)*wf_L*k1*x(1)^k2)/(alpha_L*wf_L)))+x(3)*(1-k0-((gamma*(x(4)+x(5))*x(2)+(gamma*k2+alpha_L)*wf_L*k1*x(2)^k2)/(alpha_L*wf_L)))))^(rho-1)) ;
      wf_H-(x(7)*(x(8)*(N_L*((1-x(3))*(1-k0-((gamma*x(4)*x(1)+(gamma*k2+alpha_L)*wf_L*k1*x(1)^k2)/(alpha_L*wf_L)))+x(3)*(1-k0-((gamma*(x(4)+x(5))*x(2)+(gamma*k2+alpha_L)*wf_L*k1*x(2)^k2)/(alpha_L*wf_L)))))^rho+(1-x(8))*(N_H*(1-k0-((gamma*(x(4)+x(5))*b_H+(gamma*k2+alpha_H)*wf_H*k1*b_H^k2)/(alpha_H*wf_H))))^rho)^(1/rho-1)*(1-x(8))*(N_H*(1-k0-((gamma*(x(4)+x(5))*b_H+(gamma*k2+alpha_H)*wf_H*k1*b_H^k2)/(alpha_H*wf_H))))^(rho-1)) ;
      ];

