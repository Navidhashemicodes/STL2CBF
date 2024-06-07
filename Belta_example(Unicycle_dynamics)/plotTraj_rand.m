function plotTraj_rand(initxx,inityy,initalpha, model_UC, numbers, M, P, Param, T, fignum1, fignum2)

rng(0)

l=30;
u=37.4;
lowest=0.2;

for ii=1:numbers
    q= model_UC(1)+rand*(model_UC(2)-model_UC(1));
    s0 = [initxx(1);inityy(1);initalpha(1)] + rand(3,1).*[initxx(2)-initxx(1);  inityy(2)-inityy(1);  initalpha(2)-initalpha(1)];
    a0 = Actor(s0,0,M,P,Param);
    xx = zeros(3,T);
    xx(:,1)=Dynamics(s0,a0,q);
    for ij=1:T
        a(:,ij)=Actor(xx(:,ij),ij,M,P,Param);
        xx(:,ij+1)=Dynamics(xx(:,ij),a(:,ij),q);
    end
    
    rob = robustness_Belta(s0,Param, q, M, P, T);
    
    if (rob < 0)
        figure(fignum1)
        title('Violating Trajectories')
        color = '-r.' ;
%         violate=[s0 xx];
        plot([s0(1,1) xx(1,:)],[s0(2,1) xx(2,:)], color, 'Linewidth', 1.5);
        
    else
        figure(fignum2)
%         title('Satisfactory Trajectories')
        reward=Perf_obj(s0, Param, q, M, P, T);
        color=[  reward/(l-u)-u/(l-u),   (1-lowest)*reward/(l-u)+(lowest*l-u)/(l-u), reward/(l-u)-u/(l-u)];
        plot([s0(1,1) xx(1,:)],[s0(2,1) xx(2,:)], 'color',color, 'Linewidth', 1.5);
    end
hold on    
end

xlim([-5 15]);
ylim([-5 15]);
% figure(fignum2)
% plot(violate(1,:),violate(2,:), '-r.', 'Linewidth', 1.5);
end