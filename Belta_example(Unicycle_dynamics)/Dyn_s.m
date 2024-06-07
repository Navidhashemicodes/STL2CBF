function d_d_s=Dyn_s(s,a,q)
a(1,1)=0.5*a(1,1);
a(2,1)=0.5*a(2,1);
f1=  (2*sigmoid(a(1,1))/tanh(a(2,1))) * ( cos(s(3,1)+0.5*tanh(a(2,1)))  -  cos(s(3,1)) )  ;
f2=  (2*sigmoid(a(1,1))/tanh(a(2,1))) * ( sin(s(3,1)+0.5*tanh(a(2,1)))  -  sin(s(3,1)) )  ;  

d_d_s=[ q      0     f1 ;...
        0      q     f2 ;...
        0      0      q ];
end