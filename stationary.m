clear all
close all
clc
% Grid and grid information
G = cartGrid([5 5]);
G = computeGeometry(G);
N = G.faces.neighbors;
N = N(all(N ~= 0, 2), :);
nf = size(N,1);
nc = G.cells.num;
rock.perm=ones([G.cells.num, 1]);
hT = computeTrans(G, rock);
cf = G.cells.faces(:,1);
T = 1 ./ accumarray(cf, 1 ./ hT, [G.faces.num, 1]);
T = T(all(N~=0,2),:);

% Operators
C = sparse([(1:nf )'; (1:nf )'], N, ...
ones(nf,1)*[-1 1], nf, nc);
grad = @(x) C*x;
div = @(x) -C'*x;

% Assemble and solve equations
p = initVariablesADI(zeros(nc,1));
q = zeros(nc, 1);
 % source term
q(1) = -1; q(nc) = 1;
 % −> quarter five−spot
eq = -div(grad(p))+q;
 % equation

 % make solution unique
 eq(1) = eq(1) + q(1);
  eq(nc) = eq(nc) + q(nc);
p = eq.jac{1}\eq.val; % solve equation
plotCellData(G,p);
%l=ichol(eq.jac{1});

A=eq.jac{1};
l=diag(length(A));
b=eq.val;
figure
h1=plotCellData(G,p);
figure
 plotCellData(G,A\b);
 [Va,Da] = eigs(inv(l)*A(2:24,2:24)*inv(l'),length(A)-2);
 Da=diag(Da);
 Da=abs(Da);
 conda=condest(inv(l)*A*inv(l'),2)
 condadeff=max(Da)/min(Da)
 figure(40)
subplot(2,1,1)
h21=semilogy(Da,'*');
axis('tight')
title('M_1A')
 [Va,Da] = eigs(A(2:24,2:24),length(A)-2);
 Da=diag(Da);
 Da=abs(Da);
 conda=condest(A,2)
 condadeff=max(Da)/min(Da)
  figure(40)
subplot(2,1,2)
h21=semilogy(Da,'*');
title('A')
axis('tight')
full(A)