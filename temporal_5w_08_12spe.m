%% Single phase flow simulation using AD
% This example goes through the steps of setting up a single-phase
% simulation with a single horizontal well using the automatic
% differentiation framework.
close all
clear all
clc
% Required modules
try
    require ad-fi
catch
    mrstModule add ad-fi
end
%save graphs
ss=0;
%solve with deflation 1
def=0;
%number of deflation vectors
dvec=10;
%use pod or not
pod=0;
%eigenvectos used
podv=[1:10];
% change the wells
vwells=0;
%compute eigenvalues of the matrix
eig=1;

extra='';

if vwells==1
    extra=[extra 'vw' ];
end
if def==1
    extra=[extra 'dv_' num2str(dvec) ];
    if pod==1
        extra=[extra 'pod' num2str(podv)];
    end
end

% % Setup 10x10x10 grid of 200x200x50 m model.
% nx = 35;    ny = 35;    nz = 1;
% Dx = 35;   Dy = 35;   Dz = 50;
% G = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
% %G = cartGrid([nx, ny, nz]);
% G = computeGeometry(G);

%% Setup rock properties
%
% Assume that the oil compressibility can be approximated as constant in
% the reservoir:
fc=1;
c    = fc*1e-2/barsa;

% % Assume homogeneous/isotropic rock.
% permX = 30*milli*darcy;
% poro  = 0.3;
% rock.perm = repmat(permX, [G.cells.num, 1]);
% rock.poro = repmat(poro , [G.cells.num, 1]);
% sz=nx;
% lsz=round(sz/5);
% per=1;
% for i= 1:2:5
%     rock.perm(1+(i-1)*lsz*nx:i*lsz*nx)  = 30*10^(-per)*milli*darcy();
% end
dir='/home/wagm/cortes/Localdisk/Results/16_10/20/';
%dir='/dev/media/Sphinx/Doctorado_Delft/Results/16_10/20/';
% folder=[ 'size_'   num2str(nx) 'perm_' num2str(per) '_5wells_' 'c_' num2str(fc) 'e-3' extra];
% mkdir([dir], folder)
% dir = [dir folder '/'];
  mrstModule add spe10 coarsegrid;
 nxi=1;
 nyi=1;
 %gridsize
 nx=60;
 ny=220;
 %original size
 nxf=60;
 nyf=220;
 nz=1;
 nzf=1;
 layers = 2;
 fine = [nxf nyf nzf];
 coarse = [nx ny nz];
  name=['SPE10_' num2str(nx) '_' num2str(ny) ];
  mkdir([dir], name )
  dir=[dir name '/'];
%cartDims = [  nx,  ny,   numel(layers)];
    %physDims = [1200, 2200, 2*cartDims(end)] ;   % ft -> m
    %G = cartGrid([nx ny numel(layers)]);
   % G = computeGeometry(G);
    %rock = SPE10_rock(layers);
   [G,rock,hg]=avercells(fine, coarse,layers,3,0);
   savfig(dir,'Grid',hg)
poro  = 0.3;
rock.poro=poro*ones(length(rock.perm),1);
% hp=plotCellData(G, rock.perm);
% break
% Set rock compressibility:
cr = 0e-6/barsa;

%%
% In the case of non-zero rock-compressiblity cr, the input rock porosity
% is taken as reference at a given reference pressure p_r. The grid pore
% volume (pv) becomes a function of pressure given by the differential
% equation
%                 cr = (d pv/d p)/pv
% which results in
%                 pv(p) = pv_r e^( cr(p-p_r) )
% in which pv_r is the reference pore volume (rock.poro x volume) and p_r
% is the reference pressure. We assume the reference pressure is 200 Bar:

pv_r = poreVolume(G, rock);
p_r  = 200*barsa;

% Finally, the pressure-dependent function for pore-volumes becomes:

pv   = @(p) pv_r .* exp( cr .* (p - p_r) );

%% Fluid (oil) properties
% We Assume constant viscosity:
mu   = 5*centi*poise;
% Assume that the oil compressibility can be approximated as constant in
% the reservoir:

% When the compressibility is constant, the fluid density becomes a
% function of pressure given by the differential equation
%                 c = (d rho/d p)/rho
% which results in
%                 rho(p) = rho_r e^( c(p-p_r) )
% in which rho_r is the reference fluid density at reference pressure p_r.
% We assume that rho_r = 850 kg/m^3 at p_r = 200 Bar:

rho_r = 850*kilogram/meter^3;

% Finally define the pressure dependent function for rho:
rho   = @(p) rho_r .* exp( c .* (p - p_r) );

% Computing surface volume rates requires fluid density at surface
% conditions. We assume it is 750 kg/m^3:
rhoS = 750*kilogram/meter^3;

%% Single horizontal well in J direction
% We consider a well with 8 connections in the J direction.

W = [];
nperf = 1;
I = repmat(1, [nperf, 1]);
J = 1;
K = repmat(1, [nperf, 1]);

% Convert IJK-indices to linear index (as used in G)
cellInx = sub2ind(G.cartDims, I, J, K);

W = addWell(W, G, rock, cellInx, 'Name', 'W1', 'InnerProduct', 'ip_tpf', 'sign', -1);

nperf = 1;
I = repmat(nx, [nperf, 1]);
J = 1; %(1 : nperf).';
K = repmat(1, [nperf, 1]);

cellInx = sub2ind(G.cartDims, I, J, K);

W = addWell(W, G, rock, cellInx, 'Name', 'W2', 'InnerProduct', 'ip_tpf', 'sign', -1);

nperf = 1;
I = 1;
J = ny; %(1 : nperf).';
K = repmat(1, [nperf, 1]);

% Convert IJK-indices to linear index (as used in G)
cellInx2 = sub2ind(G.cartDims, I, J, K);

W = addWell(W, G, rock, cellInx2, 'Name', 'W3', 'InnerProduct', 'ip_tpf', 'sign', 1);

nperf = 1;
I = nx;
J = ny;
K = repmat(1, [nperf, 1]);

% Convert IJK-indices to linear index (as used in G)
cellInx = sub2ind(G.cartDims, I, J, K);

W = addWell(W, G, rock, cellInx, 'Name', 'W4', 'InnerProduct', 'ip_tpf', 'sign', -1);

nperf = 1;
I = ceil(nx/2);
J = ceil(ny/2);
K = repmat(1, [nperf, 1]);

% Convert IJK-indices to linear index (as used in G)
cellInx = sub2ind(G.cartDims, I, J, K);

W = addWell(W, G, rock, cellInx, 'Name', 'W5', 'InnerProduct', 'ip_tpf', 'sign', -1);




% Plotting
f = [figure(1), figure(2)];
set(0, 'CurrentFigure', f(1));
clf
hp=plotCellData(G, log(rock.perm(:,1)));

%
  view(0,90)
  colormap jet
  axis equal tight
%plotGrid(G, 'FaceColor', 'g', 'FaceAlpha', .3, 'EdgeColor', 'w');
plotWell(G, W);
axis off;
set(f(1), 'Color', 'w');
camproj perspective;
view(3);
%break


%% Initial conditions
% We assume that the reservoir is initially at equilibrium. This means that
% the following condition must be satisfied:
%          dp/dz = g rho,
% where g is the gavitational accelleration. This relation can be solved
% analytically for p, but alternatively one can solve the above ODE with
% 'initial condtition' p(z_0) = p_r:

gravity reset on
g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,3)) + 5);
press  = ode23(@(z,p) g .* rho(p), [z_0, z_max], p_r);
% We then interpolate onto the grid using cell centers:
p_init = reshape(deval(press, G.cells.centroids(:,3)), [], 1);  clear press

%% Setting up components needed for the simulation
% Since we impose no-flow boundary conditions in this example, we restrict
% connections to interior faces only.
N  = double(G.faces.neighbors);
intInx = all(N ~= 0, 2);
N  = N(intInx, :);
% We will be using the two-point flux approximation. First the one-sided
% transmissibilities are computed, then the harmonic average is taken to
% obtain the two-sided transmissibilites.
hT = computeTrans(G, rock);
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]);
T  = T(intInx);
% In setting up the equations, we need discrete forms of the divergence and
% gradient operators, and we represent these as multiplication by sparse
% matrices. In particular, we construct the 'gradient matrix' C as folows:
n = size(N,1);
C = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[1 -1], n, G.cells.num);
% The discrete gradient and divergence operators are now given by
grad = @(x)-C*x;
div  = @(x)C'*x;

% Additionally, we will need to take the average of cell-based quantities
% in neighboring cells and define the following function
avg = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));

%% Pressure and well equations:
% The pressure equation (without well contributions) is given by:
%
% $\frac{d}{dt}(\phi\rho)+\nabla\cdot(\rho v)=0, \quad v = -\frac{K}{\mu}\nabla(p-g\rho z)$
%
% In discretized form, this leads to
z = G.cells.centroids(:,3); % z-coordinate of grid cells
pressureEq = @(p, p0, dt) (1/dt) .* (pv(p).*rho(p) - pv(p0).*rho(p0)) ...
    - div( avg(rho(p) ./ mu) .* T .* grad(p - g*rho(p).*z) );

% Wellrates are given as Peaceman well-index times pressure drop
% wc = W(1).cells; % perforation grid cells
% WI = W(1).WI;    % well indices
% dz = W(1).dZ;    % perforation depth relative to well reference depth
%
% wc2 = W(2).cells; % perforation grid cells
% WI2 = W(2).WI;    % well indices
% dz2 = W(2).dZ;    % perforation depth relative to well reference depth

% wellRates = {...
%    @(p, bhp) W(1).WI .* (rho(p(W(1).cells)) ./ mu) .* (bhp - p( W(1).cells) + g*dz.*rho(p( W(1).cells)));...
%    @(p, bhp) WI2 .* (rho(p(wc2)) ./ mu) .* (bhp - p(wc2) + g*dz2.*rho(p(wc2)))};
wellRates=cell(5,1);
for i=1:5
    wellRates{i}=@(p, bhp) W(i).WI .* (rho(p(W(i).cells)) ./ mu) .* (bhp - p( W(i).cells) + g*W(i).dZ.*rho(p( W(i).cells)));
end
%% Define ADI variables
% The primary variables are grid-cell pressures, well bhp and surface
% rate.

inputv{1}=p_init;
inputv{2}=p_init( W(1).cells);
inputv{3}=p_init( W(2).cells);
inputv{4}=p_init( W(3).cells);
inputv{5}=p_init( W(4).cells);
inputv{6}=p_init( W(5).cells);
inputv{7}=0;
inputv{8}=0;
inputv{9}=0;
inputv{10}=0;
inputv{11}=0;
%ad_values is a vector containing the resulting initialized adi variables
%[output{1}, output{2}, output{3},output{4}, output{5}, output{6},output{7}, output{8},output{9},output{10},output{11}]=
%[     p_ad,    bhp_ad,   bhp_ad2,  bhp_ad3,   bhp_ad4,   bhp_ad5,    qS_ad,    qS_ad2,   qS_ad3,    qS_ad4,    qS_ad5]
[output{1}, output{2}, output{3},output{4}, output{5}, output{6},output{7}, output{8},output{9},output{10},output{11}]=...
    initVariablesADI(inputv{1},inputv{2},inputv{3},inputv{4},inputv{5},inputv{6},inputv{7},inputv{8},inputv{9},inputv{10},inputv{11});
%return
%[p_ad, bhp_ad, bhp_ad2, qS_ad, qS_ad2] = initVariablesADI(p_init, p_init(wc(1)), p_init(wc2(1)), 0, 0);


%[p_ad,bhp_ad1,bhp_ad2,qS_ad1,qS_ad2]=initVariablesADI(p_init,p_init(wc(1)),p_init(wc(tb)),0,0);
% For convenience, create indices to variables when stacked:
%Ix(1) = pIx, Ix(2)=bhpIx,  Ix(3)=qSIx,Ix(4)=bhpIx2 Ix(5)=qSIx2,
%Ix(6)=bhpIx3, Ix(7)=qSIx3 Ix(8)=bhpIx4, Ix(9)=qSIx4, Ix(10)=bhpIx5,
%Ix(11)= qSIx5 ;
% pIx = 1:G.cells.num;
% bhpIx  = G.cells.num + 1; qSIx  = G.cells.num + 6;
% bhpIx2 = G.cells.num + 2; qSIx2 = G.cells.num + 7;
% bhpIx3 = G.cells.num + 3; qSIx3 = G.cells.num + 8;
% bhpIx4 = G.cells.num + 4; qSIx4 = G.cells.num + 9;
% bhpIx5 = G.cells.num + 5; qSIx5 = G.cells.num + 10;
Ix{1} = 1:G.cells.num;
for i=1:10
    Ix{i+1}=G.cells.num+i;
end

%% Set up simulation parameters
%

numSteps = 52;
totTime  = 156*day;
dt       = totTime / numSteps;
% Tolerance and maximum number of iterations for the Newton solver.
tol      = 1e-5;
maxits   = 10;
miter=zeros(numSteps,maxits);
% Save output in array 'sol'
sol = repmat(struct('time', [], 'pressure', [], 'bhp', [], 'bhp2', [], ...
    'qS', [], 'qS2', []), ...
    [numSteps + 1, 1]);

sol(1).time      = 0;
sol(1).pressure  = double(output{1});
sol(1).bhp       = double(output{2});
sol(1).bhp2      = double(output{3});
sol(1).bhp3      = double(output{4});
sol(1).bhp4      = double(output{5});
sol(1).bhp5      = double(output{6});
sol(1).qS        = double(output{7});
sol(1).qS2       = double(output{8});
sol(1).qS3       = double(output{9});
sol(1).qS4       = double(output{10});
sol(1).qS5       = double(output{11});
% Set up plot
set(0, 'CurrentFigure', f(2));
clf
set(f(2), 'Color', 'w')
subplot(2,2,1),
plotCellData(G, convertTo(p_init, barsa));
title('pressure [bar]', 'EdgeColor', 'w');
%colorbar;
view(0,90);
camproj perspective

subplot(2,2,[3,4]);
xlim([0, convertTo(totTime,day)]);
title('Surface volume rate [m^3/day]');
hold on


%% Main simulation
% We solve the equations implicitely. At each time step, the equations are assembled and
% the automatic differentiation framework takes care automatically of the computation of
% the Jacobian.

t = 0;
step = 0;
counter=0;
dv=0;
dv1=0;
nDigits = floor(log10(maxits)) + 1;
text = [dir 'iterations.txt'];
fileID = fopen(text,'w');
fprintf(fileID,'Time step & Newton iteration & # Iterations solver \\\\ \n');
fclose(fileID);
% fprintf(fileID,' %6s', '%6s ', '%6.2e ', step, nit, iter ' \\');
upd=sol(1).pressure;


while t < totTime,
    pf(1)=sol(step+1).pressure(1);
    nns=1;
    for i=0:ny-2
        if mod(i,3)==0
        nns=nns+1;
        ns=i*nx+nns-1;
       % pause
        pf(nns)=sol(step+1).pressure(ns);
        end
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
    %break
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
        if vwells==1;
            c=0;
        end
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
            for i=1:11
                output{i}.val   = output{i}.val ;
            end
            
            resNorm = norm(res);
            nit     = nit + 1;
            fprintf('  Iteration %*d:  Res = %.4e\n', nDigits, nit, resNorm);
            % pause
        else
            
            % upd = -(J \ res); % Newton update
            l=ichol(J);
            e=0;
            
            
            %         for i=1:11
            %             xi(Ix{i},1)=upd(Ix{i});
            %         end
            
            
            ne=0;
            a{step}=J;
            b{step}=res;
            xs{step}=xi;
            residual{step}=res+J*xi;
            if  def==0
                r0=res+J*xi;
                nr=sqrt(r0'*r0);
                %      nb=sqrt(res'*res)
                nb=sqrt(residual{step}'*residual{step});
                nr/nb;
                upd=-J\res;
                %[upd,ee]=ICCG_t(-J,res,xi,200,10^-5,l,dir,eig,step,nit);
                
                ee=0;
                ne=length(ee);
                miter(step,nit+1)=ne;
                %ee=0;
                fileID = fopen(text,'a');
                %fprintf(fileID,'Time step', 'Newton iteration', '# Iterations solver \\');
                fprintf(fileID,' %6d  &   %6d  &  %6d \\\\ \n', step, nit+1, ne);
                fclose(fileID);
                for i=1:11
                    output{i}.val   = output{i}.val   + upd(Ix{i});
                end
            else
                if step <dvec+1
                    if nit==0
                        dv=dv+1;
                    end
                    [upd,ee]=ICCG_t(-J,res,xi,200,10^-5,l,dir,eig,step,nit);
                    for i=1:11
                        output{i}.val   = output{i}.val   + upd(Ix{i});
                        zout(Ix{i},dv)=output{i}.val;
                    end
                    zupd(:,dv)=upd;
                    
                    
                    
                    ne=length(ee);
                    miter(step,nit+1)=ne;
                    iter=1:ne;
                    %            figure(800)
                    %             if defvect<6
                    %                 subplot(2,3,step)
                    %             else if defvect<11
                    %                     subplot(2,5,step)
                    %
                    %                 else
                    %                     subplot(5,5,step)
                    %                 end
                    %             end
                    %             colors = [1 0 0;0 1 0; 0 1 1; 1 0 1; 0 0 1;1 0 1;0.2 0.5 0.5; 0.5 0.2 1; 1 0.2 0.5; 0.2 1 0.6;];
                    %             hc=semilogy(iter,ee,'o','Color',colors(nit+1,:));
                    %
                    %             hold on
                    %             title( ['Timestep ' num2str(step )] )
                    %             leginf{nit+1} = ['N-R ' num2str(nit+1)];
                    %            legend(leginf{1:nit+1})
                    fileID = fopen(text,'a');
                    %fprintf(fileID,'Time step', 'Newton iteration', '# Iterations solver \\');
                    fprintf(fileID,' %6d  &   %6d  &  %6d \\\\ \n', step, nit+1, ne);
                    fclose(fileID);
                    
                    
                else
                    
                    if pod==1
                        if step==dvec+1
                            [U,S]=defpodf_D(zupd,dir);
                        end
                        z=U(:,podv);
                    else z=zupd;
                    end
                    [upd,ee]=DICCG_te(-J,res,xi,200,10^-5,z,l,eig,step,dir,nit,dvec);
                    for i=1:11
                        output{i}.val   = output{i}.val   + upd(Ix{i});
                    end
                    ne=length(ee);
                    miter(step,nit+1)=ne;
                    iter=1:ne;
                end
            end
            
            
            file=['miter'];
            filename=[dir file];
            save(filename,file)
            
            
            
            %             if mod(step,5)==0
            %                 figure(500)
            %                 if def==1
            %                     title('Iterations DICCG');
            %                 else
            %                     title('Iterations ICCG')
            %                 end
            %                 subplot(2,5,counter)
            %                 colors = [1 0 0;0 1 0; 0 1 1; 1 0 1; 0 0 1;1 0 1;0.2 0.5 0.5; 0.5 0.2 1; 1 0.2 0.5; 0.2 1 0.6;];
            %                 hd=semilogy(iter,ee,'o','Color',colors(nit+1,:));
            %
            %                 hold on
            %                 title( ['Timestep ' num2str(step )] )
            %                 leginf{nit+1} = ['N-R ' num2str(nit+1)];
            %                 legend(leginf{1:nit+1})
            %             end
            %             fileID = fopen(text,'a');
            %             %fprintf(fileID,'Time step', 'Newton iteration', '# Iterations solver \\');
            %             ne=length(ee);
            %
            %             fprintf(fileID,' %6d  &   %6d   & %6d \\\\  \n ', step, nit+1, ne);
            %             fclose(fileID);
            % [upd,hc0,hc01,hc02]=ICCG_t(-J,res,xi,iter,10^-7,l,e,step);
            
            
            %legend(leginf{1:nit+1})
            
            
            % Update variables
            
            %         for i=1:11
            %         output{i}.val   = output{i}.val   + upd(Ix{i});
            %         end
            
            resNorm = norm(res);
            nit     = nit + 1;
            fprintf('  Iteration %*d:  Res = %.4e\n', nDigits, nit, resNorm);
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
    
end