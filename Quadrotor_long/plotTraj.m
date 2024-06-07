function plotTraj(initxx,inityy,delta,M,P,Param,T)
    xxs = initxx(1):delta:initxx(2);
    yys = inityy(1):delta:inityy(2);
    hold on;
    for ii=1:length(xxs)
        initxx = xxs(ii);
        for jj=1:length(yys)
            inityy = yys(jj);
            s0 = [initxx;inityy;pi/2];
            a0 = Actor(s0,0,M,P,Param);
            xx = zeros(3,T);
            xx(:,1)=Dynamics(s0,a0);
            for ij=1:T
                a(:,ij)=Actor(xx(:,ij),ij,M,P,Param);
                xx(:,ij+1)=Dynamics(xx(:,ij),a(:,ij));
            end
            plot([s0(1,1) xx(1,:)],[s0(2,1) xx(2,:)], 'Linewidth', 0.75);
            xlim([-1 10]);
            ylim([-1 10]);
        end
    end
end