%Master Thesis R06323052 朝霍 

%calibration aplpha beta gamma pi rho

%price parameter
pb = 0.0982 ; %from calibration 基闽跑计
pe = 0.0760 ; %from calibration 基闽跑计
am = 1      ; %from calibration 基闽跑计 
c_ = 0.6359 ; %from calibration 基闽跑计
%home production
k0 = 0.2021 ; %from home_production.m
k1 = 0.1344 ; %from home_production.m
k2 = 0.5434 ; %from home_production.m
%production
n_L = 5008525 ; % from キАネ▅羆计
n_H =  534096 ; % from キАネ▅羆计
N_L = n_L;
N_H = n_H;
N_L = N_L/(n_L+n_H) ; %夹非てN_L
N_H = N_H/(n_L+n_H) ; %夹非てN_H
%ずネ跑计
b_L  = 2.333  ;  %from キАネ▅羆计(64)
b_H  = 1.691  ;  %from キАネ▅羆计(64)
wf_L = 0.4499 ;  %from calibration 基闽跑计
wf_H = 0.9290 ;  %from calibration 基闽跑计
R    = 0.411  ;  %from PSFD (L->L percentage)
% Assume part

pi  = (wf_L/wf_H)/((wf_L/wf_H)+(N_L*(1-k0-((gamma*pb*b_L+(gamma*k2+alpha)*wf_L*k1*(b_L)^(k2))/(alpha*wf_H))))/(N_H*(1-k0-((gamma*(pb+pe)*b_H+(gamma*k2+alpha)*wf_H*k1*(b_H)^(k2))/(alpha*wf_H))))^(rho-1))   ; % Assumption
rho = 0.401 ; % from Krusell(2000)

