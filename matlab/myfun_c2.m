function F = myfun_c2(x,phi,eta,am,gamma,k0,k1,k2,wf_L,wf_H,tf_L,tf_H,b_L,b_H,N_L,N_H,rho) 
%known parameter : phi eta am gamma k0 k1 k2 N_L N_H 
%assumption paramter : rho 
%solve: bd_L be_L R pb pe c_ A pi beta
F = [ (1+gamma+(phi*x(1)^(1-eta)))*x(4)*x(1)+((1+gamma)*k2+(phi*x(1)^(1-eta)))*wf_L*k1*(x(1))^(k2)-(phi*x(1)^(1-eta))*(wf_L*(1-k0)+am-x(6));
      (1+gamma+(phi*x(2)^(1-eta)))*(x(4)+x(5))*x(2)+((1+gamma)*k2+(phi*x(2)^(1-eta)))*wf_L*k1*(x(2))^(k2)-(phi*x(2)^(1-eta))*(wf_L*(1-k0)+am-x(6));
      (1+gamma+(phi*b_H^(1-eta)))*(x(4)+x(5))*b_H+((1+gamma)*k2+(phi*b_H^(1-eta)))*wf_H*k1*(b_H)^(k2)-(phi*b_H^(1-eta))*(wf_H*(1-k0)+am-x(6));     
      b_L-(x(3)*x(1)+(1-x(3))*x(2));
      tf_L-(x(3)*(1-k0-((gamma*x(4)*x(1)+(gamma*k2+(phi*x(1)^(1-eta)))*wf_L*k1*(x(1))^(k2))/(wf_L*(phi*x(1)^(1-eta)))))+(1-x(3))*(1-k0-((gamma*(x(4)+x(5))*x(2)+(gamma*k2+(phi*x(2)^(1-eta)))*wf_L*k1*(x(2))^(k2))/(wf_L*(phi*x(2)^(1-eta))))));
      tf_H-(1-k0-((gamma*(x(4)+x(5))*b_H+(gamma*k2+(phi*b_H^(1-eta)))*wf_H*k1*(b_H)^(k2))/(wf_H*(phi*b_H^(1-eta)))));
      wf_L-x(7)*x(8)*(((N_L*(x(3)*(1-k0-((gamma*x(4)*x(1)+(gamma*k2+(phi*x(1)^(1-eta)))*wf_L*k1*(x(1))^(k2))/(wf_L*(phi*x(1)^(1-eta)))))+(1-x(3))*(1-k0-((gamma*(x(4)+x(5))*x(2)+(gamma*k2+(phi*x(2)^(1-eta)))*wf_L*k1*(x(2))^(k2))/(wf_L*(phi*x(2)^(1-eta))))))))^(rho-1))*(x(8)*((N_L*(x(3)*(1-k0-((gamma*x(4)*x(1)+(gamma*k2+(phi*x(1)^(1-eta)))*wf_L*k1*(x(1))^(k2))/(wf_L*(phi*x(1)^(1-eta)))))+(1-x(3))*(1-k0-((gamma*(x(4)+x(5))*x(2)+(gamma*k2+(phi*x(2)^(1-eta)))*wf_L*k1*(x(2))^(k2))/(wf_L*(phi*x(2)^(1-eta))))))))^(rho)+(1-x(8))*((N_H*(1-k0-((gamma*(x(4)+x(5))*b_H+(gamma*k2+(phi*b_H^(1-eta)))*wf_H*k1*(b_H)^(k2))/(wf_H*(phi*b_H^(1-eta)))))))^(rho))^((1-rho)/rho) ;
      wf_H-x(7)*(1-x(8))*(((N_H*(1-k0-((gamma*(x(4)+x(5))*b_H+(gamma*k2+(phi*b_H^(1-eta)))*wf_H*k1*(b_H)^(k2))/(wf_H*(phi*b_H^(1-eta)))))))^(rho-1))*(x(8)*((N_L*(x(3)*(1-k0-((gamma*x(4)*x(1)+(gamma*k2+(phi*x(1)^(1-eta)))*wf_L*k1*(x(1))^(k2))/(wf_L*(phi*x(1)^(1-eta)))))+(1-x(3))*(1-k0-((gamma*(x(4)+x(5))*x(2)+(gamma*k2+(phi*x(2)^(1-eta)))*wf_L*k1*(x(2))^(k2))/(wf_L*(phi*x(2)^(1-eta))))))))^(rho)+(1-x(8))*((N_H*(1-k0-((gamma*(x(4)+x(5))*b_H+(gamma*k2+(phi*b_H^(1-eta)))*wf_H*k1*(b_H)^(k2))/(wf_H*(phi*b_H^(1-eta)))))))^(rho))^((1-rho)/rho) ;
     %x(9)+((1+gamma)*log(((x(2)/x(1))^(1-eta)*(x(4)*x(1)+wf_L*k1*k2*(x(1))^(k2))/((x(4)+x(5))*x(2)+wf_L*k1*k2*(x(2))^(k2)))+(phi/(1-eta))*(x(1)^(1-eta)-x(2)^(1-eta)))/log(wf_L/wf_H) ;

      ];