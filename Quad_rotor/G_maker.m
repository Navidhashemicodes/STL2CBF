function G= G_maker(delta_t,A,B)
N=10000;
d_lambda=delta_t/N;
G=zeros(6,6);
A_d=expm(A);
for i=1:N
    lambda=i*d_lambda;
    G=G+(A_d^lambda)*d_lambda;
end
G=G*B;
end