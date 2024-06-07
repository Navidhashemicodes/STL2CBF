function [J2, Bs] = STL_obj(s0, param, lam_param, beta, M, T,A,B, FF)

n=size(s0,1);
P=size(B,2);

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
aa=a{1};
uu=[tan(0.1*tanh(aa(1,1)/10)); tan(0.1*tanh(aa(2,1)/10));   2*tanh(aa(3,1)/10) ];
s{1}=A*s0+B*uu;
exp_sum=exp_sum + exp(-eta*CBF1(s{1},FF));
Bs(1)=CBF1(s{1},FF);
for i=2:T
    Sab{i}=Weights{1}*[s{i-1};i-1]+Biases{1};
    Sa{i}=tanh(Sab{i});
    a{i}=Weights{2}*Sa{i}+Biases{2};
    aa=a{i};
    uu=[tan(0.1*tanh(aa(1,1)/10)); tan(0.1*tanh(aa(2,1)/10));   2*tanh(aa(3,1)/10) ];
    s{i}=A*s{i-1}+B*uu;
    exp_sum=exp_sum + exp(-eta*CBF1(s{i},FF));
    Bs(i)=CBF1(s{i},FF);
end

T_F=floor(0.5*T);

sq_sum1=0;
beta_sum1=0;
for i=1:T_F
    sq_sum1=sq_sum1 + exp(beta(i))*CBF2(s{i},FF);
    beta_sum1=beta_sum1+exp(beta(i));
end
w_bar1=sq_sum1/beta_sum1;

sq_sum2=0;
beta_sum2=0;
for i=1:T_F
    sq_sum2=sq_sum2 + exp(beta(i+T_F))*CBF3(s{i},FF);
    beta_sum2=beta_sum2+exp(beta(i+T_F));
end
w_bar2=sq_sum2/beta_sum2;

w_bar=(   exp(beta(2*T_F+1))*w_bar1  +  exp(beta(2*T_F+2))*w_bar2   )   /    (  exp(beta(2*T_F+1))  +  exp(beta(2*T_F+2))  );

Bs(T+1)=w_bar;

exp_sum=exp_sum+exp(-eta*w_bar);


J2=(-1/eta)*log(exp_sum);


end



function z=CBF1(s,FF)
s=s(1:3,1);
Q1=(FF^2)*1.5*eye(3);
P1=inv(Q1);
center1=FF*[5;5;0];
F1=(s-center1)'*P1*(s-center1);
z=1-exp(1-F1);
end

function z=CBF2(s,FF)
s=s(1:3,1);
Q2=(FF^2)*1.5*eye(3);
P2=inv(Q2);
center2=FF*[2;8;0];
F2=(s-center2)'*P2*(s-center2);
z=1-F2;
end

function z=CBF3(s,FF)
s=s(1:3,1);
Q3=(FF^2)*1.5*eye(3);
P3=inv(Q3);
center3=FF*[8;2;0];
F3=(s-center3)'*P3*(s-center3);
z=1-F3;
end