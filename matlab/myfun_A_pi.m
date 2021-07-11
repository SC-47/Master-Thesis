function F = myfun_A_pi (x,Z_L,Z_H,wf_L,wf_H,rho) %
%known parameter : Z_L Z_H wf_L wf_H
%assumption paramter : rho 
%calibration parameter : A pi    
F = [ wf_L-x(1)*x(2)*((Z_L)^(rho-1))*(x(2)*(Z_L)^(rho)+(1-x(2))*(Z_H)^(rho))^((1-rho)/rho) ;
      wf_H-x(1)*(1-x(2))*((Z_H)^(rho-1))*(x(2)*(Z_L)^(rho)+(1-x(2))*(Z_H)^(rho))^((1-rho)/rho) ;
    ];
