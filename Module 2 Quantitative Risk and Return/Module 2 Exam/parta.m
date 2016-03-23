% CQF Module 2 Exam Part A
mu = [0.04;0.08;0.12;0.15];
sigma = [0.07;0.12;0.18;0.26];
R = [1 0.2 0.5 0.3;0.2 1 0.7 0.4;0.5 0.7 1 0.9;0.3 0.4 0.9 1];
r = 0.03;
Var = R.*(sigma*sigma');
w = (0.1-r)*inv(Var)*(mu-r) / ((mu-r)'*inv(Var)*(mu-r));
sigma_p = sqrt(w'*Var*w);