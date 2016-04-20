function [t1N,t2N] = L2ProjectThickness(t1E,t2E,mesh,surfMesh)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

surfh = surfMesh.surfh;
surfX = surfMesh.P;
surfTri = surfMesh.tri;
nele = surfMesh.nele;

[gx,gy,gw]=trigauc([0,1,0],[0,0,1],2); %2nd order for the Mass matrix

dof = 3;
knod = 3; %triangles


%% Assembly stuff
locneq = dof*knod; %local number of equations
%Preallocate stabilization row index matrix that is parallel friendly
rowM = ones((locneq)^2,nele);
colM= rowM;
valM = zeros((locneq)^2,nele);
neq = size(surfX,1)*dof;

rowF1 = zeros(locneq,nele); valF1 = rowF1; %"F"
rowF2 = rowF1; valF2 = rowF1; %"F"
% Test function for right hand side is element indipendent
fi = [1/3,1/3,1/3];
B = zeros(dof,dof*knod);
B(1,1:3:end) = fi;
B(2,2:3:end) = fi;
B(3,3:3:end) = fi;

for iTri = 1:nele
    %local triangle node-numbers 3x1
    iv=surfTri(iTri,:)';
    %local triangle equation numbers 3*3x1
    ieqs = zeros(dof*knod,1);
    ieqs(1:3:end) = iv*3-2;
    ieqs(2:3:end) = iv*3-1;
    ieqs(3:3:end) = iv*3-0;
    
    
    tx = surfX(iv,1);
    ty = surfX(iv,2);
    tz = surfX(iv,3);

    area = surfh(iTri).Area;
    
    %local thickness values
    t1e = t1E(iv); %scalar value
    t2e = t2E(iv); %scalar value
    

    % Right hand side
    f1e = B'*t1e*area; %dof*knod x 1, One point gauss quadrature
    f2e = B'*t2e*area; %dof*knod x 1

    % Mass matrix evaluated at three Gauss points
    Me = zeros(dof*knod); %9x9
    Bf = zeros(dof,dof*knod); %3x9
    for ig = 1:length(gw)
        fi=[1-gx(ig)-gy(ig), gx(ig), gy(ig)];
        Bf(1,1:3:end) = fi;
        Bf(2,2:3:end) = fi;
        Bf(3,3:3:end) = fi;
        Me = Me + (Bf'*Bf) *gw(ig)*area; %9-by-9
    end

    %% Assembly
    len = length(ieqs);
    AX = ieqs(:, ones(1, len));
    AY = AX';
    rowM(:,iTri) = AX(:);
    colM(:,iTri) = AY(:);
    valM(:,iTri) = Me(:);
    
    rowF1(:,iTri) = ieqs(:);
    rowF2(:,iTri) = ieqs(:);
    valF1(:,iTri) = f1e;
    valF2(:,iTri) = f2e;

end
M = sparse(rowM(:), colM(:), valM(:), neq, neq);
% assemble final load vector:
F1 = zeros(neq,1);
F2 = F1;
for iTri=1:nele
    F1(rowF1(:,iTri)) = F1(rowF1(:,iTri)) + valF1(:,iTri);
    F2(rowF1(:,iTri)) = F2(rowF2(:,iTri)) + valF2(:,iTri);
end

t1N = M\F1;
t2N = M\F2;

end

