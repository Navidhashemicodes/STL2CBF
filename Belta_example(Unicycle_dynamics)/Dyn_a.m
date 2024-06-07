function d_d_a=Dyn_a(s,a)

a(1,1)=0.5*a(1,1);
a(2,1)=0.5*a(2,1);

f11=exp(-a(1,1))* (  sin(s(3,1)+0.5*tanh(a(2,1))) -  sin(s(3,1))  )    /   (    ((1+exp(-a(1,1)))^2)*tanh(a(2,1))     );

f21=exp(-a(1,1))* (  cos(s(3,1))  -  cos(s(3,1)+0.5*tanh(a(2,1))) )    /   (    ((1+exp(-a(1,1)))^2)*tanh(a(2,1))     );

f31=0;

f12=(  sigmoid(a(1,1))*(1-tanh(a(2,1))^2)*(sin(s(3,1))-sin(s(3,1)+0.5*tanh(a(2,1))))    /   (   tanh(a(2,1))^2   )    )  +...
                 (  0.5*sigmoid(a(1,1))*(1-tanh(a(2,1))^2)*cos(s(3,1)+0.5*tanh(a(2,1)))    /    tanh(a(2,1))      )    ;

f22=(  sigmoid(a(1,1))*(1-tanh(a(2,1))^2)*(cos(s(3,1)+0.5*tanh(a(2,1)))-cos(s(3,1)))    /   (   tanh(a(2,1))^2   )    )  +...
                 (  0.5*sigmoid(a(1,1))*(1-tanh(a(2,1))^2)*sin(s(3,1)+0.5*tanh(a(2,1)))    /    tanh(a(2,1))      )    ;

f32= 0.25 * (1-tanh(a(2,1))^2);

d_d_a=[  f11   f12   ;...
         f21   f22   ;...
         f31   f32   ];
end