function plotTraj(initxx,inityy,initalpha, model_UC, numbers, M, Ma, P, Param, T, fignum1, fignum2)

    xxs = linspace(initxx(1),initxx(2), numbers(1));
    yys = linspace(inityy(1),inityy(2), numbers(2));
    alphas = linspace(initalpha(1),initalpha(2), numbers(3));
    model_UCs = linspace(model_UC(1), model_UC(2), numbers(4));
    
    slope=1/(2*(Ma(1)-Ma(2)));
    intersect = 1- Ma(1)*slope;
    
    for ii=1:length(xxs)
        initxx = xxs(ii);
        for jj=1:length(yys)
            inityy = yys(jj);
            for kk=1:length(alphas)
                initalpha = alphas(kk);
                for ll=1:length(model_UCs)
                    q= model_UCs(ll);
                    s0 = [initxx;inityy;initalpha];
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
                        plot([s0(1,1) xx(1,:)],[s0(2,1) xx(2,:)], color, 'Linewidth', 0.75);
                    
                    else
                        figure(fignum2)
                        title('Satisfactory Trajectories')
                        reward=Perf_obj(s0, Param, q, M, P, T);
                        plot([s0(1,1) xx(1,:)],[s0(2,1) xx(2,:)], 'color',0.9*(slope*reward+intersect)*[0,1,0], 'Linewidth', 1);
                    end

                    hold on
                    xlim([-1 15]);
                    ylim([-1 15]);
                end
            end
        end
    end
end