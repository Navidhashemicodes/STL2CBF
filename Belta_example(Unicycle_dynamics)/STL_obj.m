function [J2, Bs] = STL_obj(s0, param, lam_param, beta, q, M, P, T)

n=size(s0,1);

index=0;
Weights{1}=reshape(param(1,index+1:index+M*(n+1)), [n+1, M])';  
index=index+M*(n+1);
Biases{1}= param(1,index+1:index+M)';
index=index+M;
Weights{2}=reshape(param(1,index+1:index+M*P), [M, P])';
index=index+M*P;
Biases{2}=param(1,index+1:index+P)';


eta=(lam_param^2)+1;

exp_sum=0;

Sab{1}=Weights{1}*[s0;0]+Biases{1};
Sa{1}=tanh(Sab{1});
a{1}=Weights{2}*Sa{1}+Biases{2};
s{1}=Dynamics(s0,a{1},q);
exp_sum=exp_sum + exp(-eta*CBF1(s{1}));
Bs(1)=CBF1(s{1});
for i=2:T
    Sab{i}=Weights{1}*[s{i-1};i-1]+Biases{1};
    Sa{i}=tanh(Sab{i});
    a{i}=Weights{2}*Sa{i}+Biases{2};
    s{i}=Dynamics(s{i-1},a{i},q);
    exp_sum=exp_sum + exp(-eta*CBF1(s{i}));
    Bs(i)=CBF1(s{i});
end


sq_sum1=0;
beta_sum1=0;
for i=1:10
    sq_sum1=sq_sum1 + beta(i)^2*CBF2(s{i});
    beta_sum1=beta_sum1+beta(i)^2;
end
w_bar1=sq_sum1/beta_sum1;

sq_sum2=0;
beta_sum2=0;
for i=1:10
    sq_sum2=sq_sum2 + beta(i+10)^2*CBF3(s{i});
    beta_sum2=beta_sum2+beta(i+10)^2;
end
w_bar2=sq_sum2/beta_sum2;

w_bar=(   (beta(21)^2)*w_bar1  +  (beta(22)^2)*w_bar2   )   /    (  beta(21)^2  +  beta(22)^2  );

Bs(T+1)=w_bar;


exp_sum=exp_sum+exp(-eta*w_bar);

J2=(-1/eta)*log(exp_sum);
% if (J2 > 0) 
%     fprintf('J2 value = %d, \n', J2);
%     for jj=1:21
%         fprintf(' B(%d) = %d. \n', jj, Bs(jj));
%     end
% end
end

function z=CBF1(s)
s=s(1:2,1);
Q1=1.5*eye(2);
P1=inv(Q1);
center1=[5;5];
F1=(s-center1)'*P1*(s-center1);
z=1-exp(1-F1);
end

function z=CBF2(s)
s=s(1:2,1);
Q2=1.5*eye(2);
P2=inv(Q2);
center2=[2;8];
F2=(s-center2)'*P2*(s-center2);
z=1-F2;
end

function z=CBF3(s)
s=s(1:2,1);
Q3=1.5*eye(2);
P3=inv(Q3);
center3=[8;2];
F3=(s-center3)'*P3*(s-center3);
z=1-F3;
end