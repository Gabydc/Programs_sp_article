temporal_5w_sys

while t < totTime,
    pf(1)=sol(step+1).pressure(1);
    nns=1;
    for i=1:nx-1
        nns=nns+1;
        ns=1+nx*i+i;
        pf(nns)=sol(step+1).pressure(ns);
    end
    set(0, 'CurrentFigure', f(2));
    subplot(2,2,2),  %caxis([120, 205])
    plot(convertTo(pf, barsa));
    xlim([1, nx]);
    title('Pressure diagonal [bar]');
    hold on
    subplot(2,2,[3,4]);
    xlim([0, convertTo(totTime,day)]);
    title('Surface volume rate [m^3/day]');
    hold on
    t = t + dt;
    
    step = step + 1;
    fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
        step, convertTo(t - dt, day), convertTo(t, day));
    
    
    % Newton loop
    resNorm = 1e99;
    p0  = double(output{1}); % Previous step pressure
    nit = 0;
    
    if mod(step,5)==0 && step>5;
        counter=counter+1;
    end
    leginf= cell(10,0);
    
    
    while (resNorm > tol) && (nit < maxits)

        % Create equations:
        %eqs = cell([3, 1]);
        miter(step,nit+1)=0;
        eqs{1} = pressureEq(output{1}, p0, dt);
        % Add well contributions in perforated cells:
        eqs{1}(W(1).cells) = eqs{1}(W(1).cells) - wellRates{1}(output{1}, output{2});
        eqs{1}(W(2).cells) = eqs{1}(W(2).cells) - wellRates{2}(output{1}, output{3});
        eqs{1}(W(3).cells) = eqs{1}(W(3).cells) - wellRates{3}(output{1}, output{4});
        eqs{1}(W(4).cells) = eqs{1}(W(4).cells) - wellRates{4}(output{1}, output{5});
        eqs{1}(W(5).cells) = eqs{1}(W(5).cells) - wellRates{5}(output{1}, output{6});
        
        % Final equation is prescribed bhp at producers (injector)
        pw(1)=100;
        pw(2)=100;
        pw(3)=100;
        pw(4)=100;
        pw(5)=600;
        if vwells==1
            if step<=dvec
                if mod(step,1)==0
                    pw(1)=200;
                    pw(5)=500;
                else if mod(step,2)==0
                        pw(2)=200;
                        pw(5)=500;
                    else if mod(step,3)==0
                            pw(3)=200;
                            pw(5)=500;
                        else if mod(step,4)==0
                                pw(4)=200;
                                pw(5)=500;
                            else
                                pw(1)=100;
                                pw(5)=600;
                            end
                        end
                    end
                end
            end
        end
        
        eqs{2} = output{2} - pw(1)*barsa;
        eqs{3} = output{3} - pw(2)*barsa;
        eqs{4} = output{4} - pw(3)*barsa;
        eqs{5} = output{5} - pw(4)*barsa;
        eqs{6} = output{6} - pw(5)*barsa;
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
        resNorm = norm(res);
        xi=zeros(Ix{11},1);
        %    if nit==0
        %         xi=zeros(Ix{11},1);
        % %xi=rand(Ix{11},1);
        % else
        %
        %              xi(1:Ix{11},1)=upd(1:Ix{11},1);
        %
        % end
        if norm(res) < tol

            nit     = nit + 1;
            fprintf('  Iteration %*d:  Res = %.4e\n', nDigits, nit, resNorm);
            figure(700)
             if nit<5
             subplot(2,2,nit) 
            hr=plot(step,resNorm,plotStyle{mod(step,6)+1},'Color', colors(mod(step,6)+1,:));
            hold on
            title( ['NR error, NR Iteration ' num2str(nit)] )
             end
             
        else
            
            % upd = -(J \ res); % Newton update
            l=ichol(J);
            e=0;
            ne=0;

            if  def==0
                [upd,ee,md,hls]=ICCG_te1(-J,res,xi,200,10^-5,l,dir,eig,step,nit,md);

            else
                if step <dvec+1
                    if nit==0
                        dv=dv+1;
                    end
                    [upd,ee,md,hls]=ICCG_te1(-J,res,xi,200,10^-5,l,dir,eig,step,nit,md);
                    zupd(:,dv)=upd;               
                else
                    
                    if pod==1
                        if step==dvec+1
                            [U,S]=defpodf_Dt(zupd,dir);
                        end
                        z=U(:,podv);
                    else z=zupd;
                    end
                    [upd,ee,hls]=DICCG_te1(-J,res,xi,200,10^-8,z,l,eig,step,dir,nit,dvec,md);  
                end
            end
                ne=length(ee);
                miter(step,nit+1)=ne;
                iter=1:ne;
                for i=1:11
                    output{i}.val   = output{i}.val   + upd(Ix{i});
                end
            
            
            file=['miter'];
            filename=[dir file];
            save(filename,file)
            
            
            
            % Update variables
            
            %         for i=1:11
            %         output{i}.val   = output{i}.val   + upd(Ix{i});
            %         end
            
            nit     = nit + 1;
            fprintf('  Iteration %*d:  Res = %.4e\n', nDigits, nit, resNorm);
             figure(700)
             if nit<5
             subplot(2,2,nit) 
            hr=plot(step,resNorm,plotStyle{mod(step,6)+1},'Color', colors(mod(step,6)+1,:));
            hold on
            title( ['NR error, NR Iteration ' num2str(nit)] )
             end
        end
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
        subplot(2,2,1), cla, %caxis([120, 205])
        plotCellData(G, convertTo(sol(step+1).pressure, barsa), 'EdgeColor', 'w');
        
        subplot(2,2,[3,4])
        xaxis = convertTo(sol(step+1).time, day);
        y1 = convertTo(-sol(step+1).qS , meter^3/day);
        y2 = convertTo(-sol(step+1).qS2 , meter^3/day);
        y3 = convertTo(-sol(step+1).qS3 , meter^3/day);
        y4 = convertTo(-sol(step+1).qS4 , meter^3/day);
        y5 = convertTo(-sol(step+1).qS5 , meter^3/day);
        
        hs= plot(xaxis, y1,'x', xaxis,y2,'r*', xaxis,y3,'bo', xaxis,y4,'gp', xaxis,y5,'ms');
        ylabel('Rate [meter^3/day]','FontSize',14)
        xlabel('Time [days]','FontSize',14)
        legend('P1','P2','P3','P4','I1')
        
        drawnow
        
    end
    
end


figure(400)
for i=1:4
    
    if def==1
        subplot(2,2,i)
        hn=plot(1:dv,miter(1:dv,i),'r*');
        hold on
        subplot(2,2,i)
        hn=plot(dv+1:numSteps,miter(dv+1:numSteps,i),'bp');
        hold on
        legend('ICCG','DICCG');
    else
        subplot(2,2,i)
        hn=plot(miter(:,i),'r*');
        legend('ICCG');
    end
    ylabel('Number of iterations','FontSize',16)
    xlabel('Time step ','FontSize',16)
    title(['Iterations N-R' num2str(i) ],'FontSize',16)
end
if ss==1
    file='solution';
    B=[dir  file  '.fig'];
    saveas(hs(1),B)
    B=[dir  file  '.jpg'];
    saveas(hs(1),B)
    
    file='iterations_4NR';
    B=[dir  file  '.fig'];
    saveas(hn,B)
    B=[dir  file  '.jpg'];
    saveas(hn,B)
    
     file='errorNR';
    B=[dir  file  '.fig'];
    saveas(hr,B)
    B=[dir  file  '.jpg'];
    saveas(hr,B)
    
       file='errorLS';
    B=[dir  file  '.fig'];
    saveas(hls,B)
    B=[dir  file  '.jpg'];
    saveas(hls,B)
end