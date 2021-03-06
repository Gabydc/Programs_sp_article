function temp_loop(dt,totTime,eqs,tol,maxits,output,rhoS,wellRates,Ix,pressureEq,W,f,G,dir)
%% Main simulation
% We solve the equations implicitely. At each time step, the equations are assembled and
% the automatic differentiation framework takes care automatically of the computation of
% the Jacobian.

t = 0; 
step = 0;
nDigits = floor(log10(maxits)) + 1;
while t < totTime,
    t = t + dt;
    step = step + 1;
    fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
            step, convertTo(t - dt, day), convertTo(t, day));

    % Newton loop
    resNorm = 1e99;
    p0  = double(output{1}); % Previous step pressure
    nit = 0;
    while (resNorm > tol) && (nit < maxits)
        % Create equations:
        %eqs = cell([3, 1]);

        eqs{1} = pressureEq(output{1}, p0, dt);
        % Add well contributions in perforated cells:
        eqs{1}(W(1).cells) = eqs{1}(W(1).cells) - wellRates{1}(output{1}, output{2});
        eqs{1}(W(2).cells) = eqs{1}(W(2).cells) - wellRates{2}(output{1}, output{3});
        eqs{1}(W(3).cells) = eqs{1}(W(3).cells) - wellRates{3}(output{1}, output{4});
        eqs{1}(W(4).cells) = eqs{1}(W(4).cells) - wellRates{4}(output{1}, output{5});
        eqs{1}(W(5).cells) = eqs{1}(W(5).cells) - wellRates{5}(output{1}, output{6});

        % Final equation is prescribed bhp at producers (injector)
        eqs{2} = output{2} + 150*barsa;
        eqs{3} = output{3} + 150*barsa;
        eqs{4} = output{4} + 150*barsa;
        eqs{5} = output{5} + 150*barsa;
        eqs{6} = output{6} + 400*barsa;       
        % Sum of wellrates should equal total rate:
        eqs{7} = output{7} - sum(wellRates{1}(output{1}, output{2}))/rhoS;
        eqs{8} = output{8} - sum(wellRates{2}(output{1}, output{3}))/rhoS;
        eqs{9} = output{9} - sum(wellRates{3}(output{1}, output{4}))/rhoS;
        eqs{10} = output{10} - sum(wellRates{4}(output{1}, output{5}))/rhoS;
        eqs{11} = output{11} - sum(wellRates{5}(output{1}, output{6}))/rhoS;
        % Concatenate equations and solve:
        eq  = cat(eqs{:});
        J   = eq.jac{1};  % Jacobian
        res = eq.val;     % residual
       % upd = -(J \ res); % Newton update
        xi=rand(length(res),1);
        tol=10^-7;
        l=ichol(-J);
        e=0;
        iter=200;
        [upd,hc0,hc01,hc02]=ICCG_0(J,res,xi,iter,tol,l,e); 
        
        % Update variables

        
        output{1}.val   = output{1}.val   + upd(Ix{1});
        output{2}.val = output{2}.val + upd(Ix{2});
        output{7}.val  = output{7}.val  + upd(Ix{7});
        output{3}.val = output{3}.val + upd(Ix{3});
        output{8}.val  = output{8}.val  + upd(Ix{8});
        output{4}.val = output{4}.val + upd(Ix{4});
        output{9}.val  = output{9}.val  + upd(Ix{9});
        output{5}.val = output{5}.val + upd(Ix{5});
        output{10}.val  = output{10}.val  + upd(Ix{10});
        output{6}.val = output{6}.val + upd(Ix{6});
        output{11}.val  = output{11}.val  + upd(Ix{11});
        resNorm = norm(res);
        nit     = nit + 1;
        fprintf('  Iteration %*d:  Res = %.4e\n', nDigits, nit, resNorm);
    end

    if nit > maxits,
        error('Newton solves did not converge')
    else
        sol(step+1).time     = t;
        sol(step+1).pressure = double(output{1});
        sol(step+1).bhp      = double(output{2});
        sol(step+1).qS       = double(output{7});
        sol(step+1).bhp2      = double(output{3});
        sol(step+1).qS2       = double(output{8});
        sol(step+1).bhp3      = double(output{4});
        sol(step+1).qS3       = double(output{9});
        sol(step+1).bhp4      = double(output{5});
        sol(step+1).qS4       = double(output{10});
        sol(step+1).bhp5      = double(output{6});
        sol(step+1).qS5       = double(output{11});
        % Plot evolution
        set(0, 'CurrentFigure', f(2));
        subplot(2,1,1), cla, %caxis([120, 205])
        plotCellData(G, convertTo(sol(step+1).pressure, barsa), 'EdgeColor', 'w');

        subplot(2,1,2)
        xaxis = convertTo(sol(step+1).time, day);
        y1 = convertTo(-sol(step+1).qS , meter^3/day);
        y2 = convertTo(-sol(step+1).qS2 , meter^3/day);
        y3 = convertTo(-sol(step+1).qS3 , meter^3/day);
        y4 = convertTo(-sol(step+1).qS4 , meter^3/day);
        y5 = convertTo(sol(step+1).qS5 , meter^3/day);
        h=plot(xaxis, y1,'*', xaxis,y5,'*');
        ylabel('Rate [meter^3/day]','FontSize',14)
    xlabel('Time [days]','FontSize',14)
      legend('P1','I')
        %plot(); hold on
        drawnow
    end
end
file=[dir '/sol/' 'sol'];
 save(file,'sol')
  file='solution';
  B=[dir   file  '.fig'];
 saveas(h,B)
 B=[dir  file  '.jpg'];
 saveas(h,B) 