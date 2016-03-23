% CQF Module 2 Exam Part A
mu = [0.04;0.08;0.12;0.15];
sigma = [0.07;0.12;0.18;0.26];
R = [1 0.2 0.5 0.3;0.2 1 0.7 0.4;0.5 0.7 1 0.9;0.3 0.4 0.9 1];
r = 0.03;
Var = R.*(sigma*sigma');
w = (0.1-r)*inv(Var)*(mu-r) / ((mu-r)'*inv(Var)*(mu-r));
sigma_p = sqrt(w'*Var*w);

A = ones(1,4)*inv(Var)*ones(4,1);
B = ones(1,4)*inv(Var)*mu;
C = mu'*inv(Var)*mu;
m = (C - B*r) / (B - A*r);
w_T = (inv(Var)*(mu-r)) / (B - A*r);
sigma_T = sqrt((C-2*r*B+r^2*A)/(B-A*r)^2);
slope = (m-r) / sigma_T;

VaR_norm = w_T'*mu + norminv(0.01,0,1)*sigma_T;
VaR_t = w_T'*mu + tinv(0.01,30)*sigma_T;