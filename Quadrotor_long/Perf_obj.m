function J1 = Perf_obj(s0, param, M, T,A,B,FF)


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


discount=0.9;

Sab{1}=Weights{1}*[s0;0]+Biases{1};
Sa{1}=tanh(Sab{1});
a{1}=Weights{2}*Sa{1}+Biases{2};
aa=a{1};
uu=[tan(0.1*tanh(aa(1,1)/10)); tan(0.1*tanh(aa(2,1)/10));   2*tanh(aa(3,1)/10) ];
s{1}=A*s0+B*uu;
J1=discount^(1)*reward(s{1},FF);
for i=2:T
    Sab{i}=Weights{1}*[s{i-1};i-1]+Biases{1};
    Sa{i}=tanh(Sab{i});
    a{i}=Weights{2}*Sa{i}+Biases{2};
    aa=a{i};
    uu=[tan(0.1*tanh(aa(1,1)/10)); tan(0.1*tanh(aa(2,1)/10));   2*tanh(aa(3,1)/10) ];
    s{i}=A*s{i-1}+B*uu;
    J1=J1+discount^(i)*reward(s{i},FF);
end


end





function z=reward(s,FF)
s=s(1:3,1);
dest=FF*[8;8;-3];
F=(s-dest)'*(s-dest);
z=10*exp(-F/((FF^2)*36));
end