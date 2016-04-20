function [SU,SF] = Project3DToSurface(surfh,H,F,U)
%Project3DToSurface 
%   [SU,FU] = Project3DToSurface(surfh,H,F,U)
% hej

nele = length([surfh.iel]);
surfh(nele+1:length(surfh)) = [];
nele=size(surfh,2); %number of elements

X = [H.XC, H.YC, H.ZC];
nodes = H.Connectivity;
tri = H.SurfaceP1Triangulation;

SU = zeros(6*nele,3); SF = SU;
for iTri = 1:length(surfh)
    iel = surfh(iTri).iel;
    iv=nodes(iel,:)';
    
    %local Hex equation numbers 24x1
    ieqs = zeros(24,1);
    ieqs(1:3:end) = iv*3-2;
    ieqs(2:3:end) = iv*3-1;
    ieqs(3:3:end) = iv*3-0;
    
    Uloc = U(ieqs); %Local solution
    Floc = F(ieqs); %Local load
    
    % Triangle node coordinates
    TX = surfh(iTri).Xe;

    ux = Uloc(1:3:end);uy = Uloc(2:3:end);uz = Uloc(3:3:end);
    fx = Floc(1:3:end);fy = Floc(2:3:end);fz = Floc(3:3:end);
    
    
    
    [fi, ~, ~, ~, ~] = H.baseHexP1(iel, TX);
    Utri = [sum(fi.*ux(:,ones(3,1)),1)', sum(fi.*uy(:,ones(3,1)),1)', sum(fi.*uz(:,ones(3,1)),1)'];
    Ftri = [sum(fi.*fx(:,ones(3,1)),1)', sum(fi.*fy(:,ones(3,1)),1)', sum(fi.*fz(:,ones(3,1)),1)'];
    
    % P2 triangulation
    ivt = tri(iTri,:);
    
    SU(ivt,:) = Utri;
    SF(ivt,:) = Ftri;

    
end
SU = SU(1:max(tri(:)),:);
SF = SF(1:max(tri(:)),:);

end

