clear all

% 方形网格
node1 = [0,0; 1,0; 1,1; 0,1];
elem1 = [2,3,1; 4,1,3];

% L型网格
node2 = [0,0; 0,1; -1,1; -1,0; -1,-1; 0,-1; 1,-1; 1,0];
elem2 = [1,2,3; 4,1,3; 1,4,6; 5,6,4; 6,7,1; 8,1,7];

% 方程数据
pde = mysincosdata;

% 积分精度
option = [];
option.quadorder = 2;       
option.fquadorder = 4; 

% fem
[~,sqrerrL2(1),sqrerrH1(1)] = myfun(elem1,node1,pde,option);
[~,LerrL2(1),LerrH1(1)] = myfun(elem2,node2,pde,option);

for k = 1:5
 [node1,elem1] = uniformrefine(node1,elem1);
 [node2,elem2] = uniformrefine(node2,elem2);
 [~, LerrL2(k+1), LerrH1(k+1)] = myfun(elem2,node2,pde,option);
 [~, sqrerrL2(k+1), sqrerrH1(k+1)] = myfun(elem1,node1,pde,option);
end 

% 绘图
h = [1,1/2,1/4,1/8,1/16,1/32];

figure;
hold on
grid on
plot(h, sqrerrL2, 'r', 'LineWidth',1, 'Marker','+');
plot(h, sqrerrH1, 'b', 'LineWidth',1, 'Marker','*');
ax = gca();
ax.XScale = 'log';
ax.YScale = 'log';
title('方形区域插值误差')
xlabel('h')
ylabel('error')
legend('L2 error', 'H1 error');
legend('show');

figure;
hold on
grid on
plot(h, LerrL2, 'r', 'LineWidth',1, 'Marker','+');
plot(h, LerrH1, 'b', 'LineWidth',1, 'Marker','*');
ax = gca();
ax.XScale = 'log';
ax.YScale = 'log';
title('L形区域插值误差')
xlabel('h')
ylabel('error')
legend('L2 error', 'H1 error');
legend('show');


function [u,errL2,errH1] = myfun(elem,node,pde,option)

[elem2dof,edge,bdDof] = dofP2(elem);
% important constants
N = size(node,1);  NT = size(elem,1); NE = size(edge,1);
Ndof = N + NE;

% \nabla \lambda and |T|
[Dlambda,area] = gradbasis(node,elem); 

% generate sparse pattern
ii = zeros(21*NT,1); 
jj = zeros(21*NT,1); 
index = 0;
for i = 1:6
    for j = i:6
        ii(index+1:index+NT) = double(elem2dof(:,i)); 
        jj(index+1:index+NT) = double(elem2dof(:,j));  
        index = index + NT;
    end
end

[lambda, w] = quadpts(option.quadorder);
nQuad = size(lambda,1);
% compute non-zeros
sA = zeros(21*NT,nQuad);
for p = 1:nQuad
    % Dphi at quadrature points
    Dphip(:,:,6) = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));
    Dphip(:,:,1) = (4*lambda(p,1)-1).*Dlambda(:,:,1);            
    Dphip(:,:,2) = (4*lambda(p,2)-1).*Dlambda(:,:,2);            
    Dphip(:,:,3) = (4*lambda(p,3)-1).*Dlambda(:,:,3);            
    Dphip(:,:,4) = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
    Dphip(:,:,5) = 4*(lambda(p,3)*Dlambda(:,:,1)+lambda(p,1)*Dlambda(:,:,3));
    index = 0;
    for i = 1:6
        for j = i:6
            Aij = 0;            
            Aij = Aij + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2);
            Aij = Aij.*area;
            sA(index+1:index+NT,p) = Aij;
            index = index + NT;
        end
    end
end
sA = sum(sA,2);
% assemble the matrix
diagIdx = (ii == jj);   
upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Ndof,Ndof);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ndof,Ndof);
A = A + AU + AU';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemble the mass matrix M
sM = zeros(21*NT, nQuad);
for p = 1:nQuad
    phip(:,6) = 4*lambda(p,1).*lambda(p,2);
    phip(:,1) = lambda(p,1).*(2*lambda(p,1)-1);
    phip(:,2) = lambda(p,2).*(2*lambda(p,2)-1);
    phip(:,3) = lambda(p,3).*(2*lambda(p,3)-1);
    phip(:,4) = 4*lambda(p,2).*lambda(p,3);
    phip(:,5) = 4*lambda(p,3).*lambda(p,1);
    index = 0;
    for i = 1:6
        for j = i:6
            Mij = 0;
            Mij = Mij + w(p)*phip(:,i).*phip(:,j);
            Mij = Mij.*area;
            sM(index+1:index+NT,p) = Mij;
            index = index + NT;
        end
    end
end
sM = sum(sM,2);

% assemble the mass matrix
M = sparse(ii(diagIdx), jj(diagIdx), sM(diagIdx), Ndof, Ndof);
MU = sparse(ii(upperIdx), jj(upperIdx), sM(upperIdx), Ndof, Ndof);
M = M + MU + MU';

% final system matrix
A = A + 2*M;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% quadrature points in the barycentric coordinate
[lambda,w] = quadpts(option.fquadorder);
nQuad = size(lambda,1);
phi(:,6) = 4*lambda(:,1).*lambda(:,2);
phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
phi(:,4) = 4*lambda(:,2).*lambda(:,3);
phi(:,5) = 4*lambda(:,3).*lambda(:,1);
bt = zeros(NT,6);
for p = 1:nQuad
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);

    fp = pde.f(pxy);   % function handle
    for j = 1:6
        bt(:,j) = bt(:,j) + w(p)*phi(p,j)*fp;
    end
end
bt = bt.*repmat(area,1,6);
b = accumarray(elem2dof(:),bt(:),[Ndof 1]); 

u = zeros(Ndof,1);

% Find Dirichlet boundary dof: fixedDof
fixedDof = []; 
freeDof = [];
isFixedDof = false(Ndof,1); 
    
% isDirichlet(elem2edge(bdFlag(:)==1)) = true;
% isFixedDof(edge(isDirichlet,:)) = true;
% isFixedDof(N + find(isDirichlet')) = true;
% fixedDof = find(isFixedDof);
% freeDof = find(~isFixedDof);    

fixedDof = bdDof;
isFixedDof(fixedDof) = true;
freeDof = find(~isFixedDof);    
   

% Modify the matrix
% Build Dirichlet boundary condition into the matrix AD by enforcing
% |AD(fixedDof,fixedDof)=I, AD(fixedDof,freeDof)=0,
% AD(freeDof,fixedDof)=0|.
if ~isempty(fixedDof)
    bdidx = zeros(Ndof,1); 
    bdidx(fixedDof) = 1;
    Tbd = spdiags(bdidx,0,Ndof,Ndof);
    T = spdiags(1-bdidx,0,Ndof,Ndof);
    AD = T*A*T + Tbd;
else
    AD = A;
end

u(freeDof) = AD(freeDof,freeDof)\b(freeDof);
residual = norm(b - AD*u);

% Dirichlet boundary conditions
idx = (fixedDof > N);                         % index of edge nodes
u(fixedDof(~idx)) = pde.g_D(node(fixedDof(~idx),:)); % bd value at vertex dofs
bdEdgeIdx = fixedDof(idx) - N;
bdEdgeMid = (node(edge(bdEdgeIdx,1),:) + node(edge(bdEdgeIdx,2),:))/2;
u(fixedDof(idx)) = pde.g_D(bdEdgeMid);
b = b - A*u;
b(fixedDof) = u(fixedDof);

errL2 = getL2error(node,elem,pde.exactu,u);
errH1 = getH1error(node,elem,pde.Du,u);
end
