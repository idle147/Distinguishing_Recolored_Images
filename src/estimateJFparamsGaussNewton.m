
function [p, H] = estimateJFparamsGaussNewton(C,d,Isize,dGamma,maxiter,convergedThresh,makePlots,displayText)

p0 = [Isize(2)/2 Isize(1)/2 1]';
p = p0;
tic
scl = -200;

for iter = 1:maxiter;
    
    %measure error vector, e, EQ(10) in Gloe et al 2010
    eMat = eLCA_JF(C,d,p);
    e = ((eMat(:,1).^2)+(eMat(:,2).^2));
    
    J = eResidJacobian(C,p,d); %jacobian
    delta_p = pinv(J)*e;
    p = p - dGamma*delta_p;
    
    SS(iter) = sum(e.^2); %sum squared of error
    p_save(iter,:) = p; %
    
    if iter > 1
        if (SS(iter-1)-SS(iter)) < convergedThresh
            break
        end
    end
end
eTime = toc;
%%
if makePlots
    H(1) = figure;
    plot3(p_save(:,1),p_save(:,2),p_save(:,3),'-o','LineWidth',2)
    hold on; plot3(p0(1),p0(2),p0(3),'ko','LineWidth',3);
    plot3(p_save(end,1),p_save(end,2),p_save(end,3),'ro','LineWidth',4)
    xlabel('X center'); ylabel('Y Center'); zlabel('Exp. Coef, \alpha')
    grid on
    
    
    %%
    H(2) = figure;
    plot(1:length(SS),SS,'.')
    % set(gca,'YScale','log'); grid on;
    
    grid on
    set(gca,'FontSize',12)
    ylabel('Sum Squared Error');
    xlabel('Iteration');
    %%
    
    H(3) = figure; hold on;
    qMeas = plot_dLCA_vectorField(C,d,scl,'r');
    [Cprime, dJF] = JohnsonFaridLCAmodel(C,p);
    qJF = plot_dLCA_vectorField(C,dJF,scl,'b');
    set(qJF,'LineWidth',1)
    hCenter = plot(p(1),p(2),'k+','LineWidth',2);
    pStr = ['Optical Center, \alpha = ' num2str(p(3))];
    
    xlim([0 max(C(:,1))+64]); ylim([0 max(C(:,2))+64])
    axis ij;
    
    legend([qMeas, qJF, hCenter], 'Measured, Local LCA Displacement','Modeled (Johnson & Farid)',pStr)
    axis equal
    
else
    H = [];
end
%%
if displayText
    disp('Gauss Newton Parameter Estimation: estimateJFparamsGaussNewton.m');
    disp(['    ' num2str(iter) ' iterations, in ' num2str(eTime) ' seconds'])
    dStr{1} = ['alpha = ' num2str(p(3))];
    dStr{2} = ['xdot = ' num2str(p(1))];
    dStr{3} = ['ydot = ' num2str(p(2))];
    disp(dStr)
end
