%Master Thesis R06323052 ����� 

%calibration aplpha beta gamma pi rho

%price parameter
pb = 0.0982 ; %from calibration ��������ܼ�
pe = 0.0760 ; %from calibration ��������ܼ�
am = 1      ; %from calibration ��������ܼ� 
c_ = 0.6359 ; %from calibration ��������ܼ�
%home production
k0 = 0.2021 ; %from home_production.m
k1 = 0.1344 ; %from home_production.m
k2 = 0.5434 ; %from home_production.m
%production
n_L = 5008525 ; % from �U�@�N�����ͨ|�`��
n_H =  534096 ; % from �U�@�N�����ͨ|�`��
N_L = n_L;
N_H = n_H;
N_L = N_L/(n_L+n_H) ; %�зǤ�N_L
N_H = N_H/(n_L+n_H) ; %�зǤ�N_H
%�����ܼ�
b_L  = 2.333  ;  %from �U�@�N�����ͨ|�`��(64)
b_H  = 1.691  ;  %from �U�@�N�����ͨ|�`��(64)
wf_L = 0.4499 ;  %from calibration ��������ܼ�
wf_H = 0.9290 ;  %from calibration ��������ܼ�
R    = 0.411  ;  %from PSFD (L->L percentage)
% Assume part

pi  = (wf_L/wf_H)/((wf_L/wf_H)+(N_L*(1-k0-((gamma*pb*b_L+(gamma*k2+alpha)*wf_L*k1*(b_L)^(k2))/(alpha*wf_H))))/(N_H*(1-k0-((gamma*(pb+pe)*b_H+(gamma*k2+alpha)*wf_H*k1*(b_H)^(k2))/(alpha*wf_H))))^(rho-1))   ; % Assumption
rho = 0.401 ; % from Krusell(2000)

