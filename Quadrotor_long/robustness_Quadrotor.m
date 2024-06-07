function rho = robustness_Quadrotor(s0,param, M, T, A, B,FF)

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



Sab{1}=Weights{1}*[s0;0]+Biases{1};
Sa{1}=tanh(Sab{1});
a{1}=Weights{2}*Sa{1}+Biases{2};
aa=a{1};
uu=[tan(0.1*tanh(aa(1,1)/10)); tan(0.1*tanh(aa(2,1)/10));   2*tanh(aa(3,1)/10) ];
s{1}=A*s0+B*uu;
rho=CBF1(s{1},FF);
for i=2:T
    Sab{i}=Weights{1}*[s{i-1};i-1]+Biases{1};
    Sa{i}=tanh(Sab{i});
    a{i}=Weights{2}*Sa{i}+Biases{2};
    aa=a{i};
    uu=[tan(0.1*tanh(aa(1,1)/10)); tan(0.1*tanh(aa(2,1)/10));   2*tanh(aa(3,1)/10) ];
    s{i}=A*s{i-1}+B*uu;
    rho=min(rho,CBF1(s{i},FF));
end


T_F=floor(0.5*T);

rho_F1=CBF2(s{1},FF);
for i=2:T_F
    rho_F1=max(rho_F1, CBF2(s{i},FF));
end

rho_F2=CBF3(s{1},FF);
for i=2:T_F
    rho_F2=max(rho_F2, CBF3(s{i},FF));
end

rho_F=max(rho_F1, rho_F2);

rho=min(rho, rho_F);


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