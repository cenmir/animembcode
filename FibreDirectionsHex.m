function m = FibreDirectionsHex(U,mesh,m)
%ElementDirections 
%   Detailed explanation goes here
surfh = mesh.SurfaceP1;
nodes = mesh.Connectivity;
[gx,gy,gw]=trigauc([0,1,0],[0,0,1],m.GO); %Triangle element integration
mu = m.mu;
lambda0 = m.lambda0;
nTri = length(surfh);
xnod = mesh.Points(:,1);
ynod = mesh.Points(:,2);
znod = mesh.Points(:,3);
CutEle = unique([surfh(:).iel])';

triEle(nTri).s1 = [0,0,0];
triEle(nTri).s2 = [0,0,0];

%% ComputeDirections per triangle
for iTri = 1:nTri
    %local Hex1 element number
    iel=surfh(iTri).iel;
    %local Hex1 node-numbers 8x1
    iv=nodes(iel,:)';
    %local Hex equation numbers 24x1
    ieqs = zeros(24,1);
    ieqs(1:3:end) = iv*3-2;
    ieqs(2:3:end) = iv*3-1;
    ieqs(3:3:end) = iv*3-0;
    % Triangle 3 node coordinates
    tx = surfh(iTri).xp;
    ty = surfh(iTri).yp;
    tz = surfh(iTri).zp;
%     TX = []
%     TX = surfh(iTri).Xe;
    
%     tf1 = m.tf1(iTri);
%     tf2 = m.tf2(iTri);
%     t = m.tb+tf1+tf2;
    
    % elemet area
    area = surfh(iTri).Area;
    % Face normal
    n = surfh(iTri).faceNormal;
    P=eye(3)-n'*n;
    
    u = U(ieqs);
    ux = u(1:3:end); uy = u(2:3:end); uz = u(3:3:end);
        
    epsSele = zeros(3);
    for ig = 1:length(gw)
        r = gx(ig); s = gy(ig);
        igw = gw(ig);
        % P1 triangular basis function
        fi = [(1-r-s),r,s];       
        % Gauss point for mother element
        xig = [fi*tx, fi*ty, fi*tz];        
        [fi, fix, fiy, fiz, ~] = mesh.baseHexP1(iel, xig);
        gradFi = [fix';fiy';fiz'];
        gradFiS=P*gradFi;
%         GradS=[gradFiS(1,:)*ux, gradFiS(2,:)*ux, gradFiS(3,:)*ux;
%                gradFiS(1,:)*uy, gradFiS(2,:)*uy, gradFiS(3,:)*uy;
%                gradFiS(1,:)*uz, gradFiS(2,:)*uz, gradFiS(3,:)*uz]  
        GradS = (gradFiS*[ux,uy,uz]).';
        
        %Use epsS
        epsS=(GradS+GradS')/2;
        epsS = P*epsS*P;

        %% epsSele
        epsSele = epsSele + (epsS*area*igw ); % Assemble to element tensor
    end
    epsSele = epsSele/area;
    
    
    
    epsS = epsSele;
    sig = 2*mu*epsS+lambda0*trace(epsS)*eye(3);
    sigP = P*sig*P;
%     sigApprox = sort(eig(sigP));
    
    [V,D]=eig(sigP);
    e = diag(D);
    % Sort by largest eigenvalue, such that the first eigenvector
    % represents the largest principal stress.
    [~,ind] = sort(e,'descend');
%     [~,ind] = sort(e);
    V = V(:,ind);
    s1 = V(:,1)/norm(V(:,1));
    s2 = V(:,2)/norm(V(:,2));
    
    triEle(iTri).s1 = s1;
    triEle(iTri).s2 = s2;
    
    triEle(iTri).epsS = epsS;
    
    
    
    
    % Check out kron, yo
    
    %% Viz
%     if iel == 10
%     xfigure(22)
%     xlabel('X'); ylabel('Y'); zlabel('Z')
%     hp = patch(tx,ty,tz,'b'); hold on; axis equal
% %     hp.FaceColor = b1;
%     [V,D]=eig(sigP)
%     sigP
%     
%     s1
% %     s2
% %     b1
% %     b2
% %     epsS
%     drawnow
%     
%     r = gx'; s = gy';
%     fi = [(1-r-s),r,s];
%     Xig = [fi*tx, fi*ty, fi*tz];
%     Xm = mean(Xig,1);
%     
%     sc = 0.1;
%     quiver3(Xm(1),Xm(2),Xm(3),s1(1),s1(2),s1(3),sc,'color','r')
%     quiver3(Xm(1),Xm(2),Xm(3),s2(1),s2(2),s2(3),sc,'color','c')
%     quiver3(Xm(1),Xm(2),Xm(3),n(1),n(2),n(3),sc,'color','k')
% %     cosTheta = s1.'*s2/(norm(s1)*norm(s2))
% %     thetad = acos(cosTheta)*180/pi
%     a= 1;
%     end
end

%% Average Directions per Hex element
map = [surfh.iel].';
hex = unique(map);

HexEle = m.HexEle;
for iel = hex.'
    iTri = find(map==iel);
    HexEle(iel).iTri = iTri;
    %% Compute orientation
    s1 = [triEle(iTri).s1];
    s2 = [triEle(iTri).s2];
    if dot(s1(:,1),s1(:,2)) < 0
        s1(:,2) = -s1(:,2);
    end
    ms1 = mean(s1.');
    ms1 = ms1/norm(ms1);
    if dot(s2(:,1),s2(:,2)) < 0
        s2(:,2) = -s2(:,2);
    end
    ms2 = mean(s2.');
    ms2 = ms2/norm(ms2);
    
    HexEle(iel).s1 = ms1; 
    HexEle(iel).s2 = ms2;
    
    S1 = ms1.'*ms1;
    S2 = ms2.'*ms2;
    
    %% Hex element mean strain tensor
    nTri = length(iTri);
    epsSh = zeros(3);
    for i=1:nTri
        epsSh = epsSh + 1/nTri * triEle(iTri(i)).epsS;
    end
    
    %% Sigma_h = alpha*Si*Si*epsSh
    alpha1 = m.alpha1; %Young type fibre elastic coefficient
    alpha2 = m.alpha2;
    % TODO: Voigt form
    sigmaS1 = zeros(3,3);
    sigmaS2 = zeros(3,3);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    sigmaS1(i,j) = sigmaS1(i,j) + alpha1*S1(k,l)*S1(i,j)*epsSh(k,l);
                    sigmaS2(i,j) = sigmaS2(i,j) + alpha2*S2(k,l)*S2(i,j)*epsSh(k,l);
                end
            end
        end
    end
    
    %% b1 and b2
    % Double contraction between Sigma_h
    b1 = dot(sigmaS1(:),epsS(:));
    b2 = dot(sigmaS2(:),epsS(:));
    
    HexEle(iel).b1 = b1;
    HexEle(iel).b2 = b2;
    
    
    %% Viz orientations
%     tx = [surfh(iTri).xp];
%     ty = [surfh(iTri).yp];
%     tz = [surfh(iTri).zp];
%     xfigure(22)
%     xlabel('X'); ylabel('Y'); zlabel('Z')
%     hp = patch(tx,ty,tz,'b'); hold on; axis equal
%     sc = 0.2;
% %     quiver3(mean(tx(:,1)),mean(ty(:,1)),mean(tz(:,1)),s2(1,1),s2(2,1),s2(3,1),sc,'color','r')
% %     quiver3(mean(tx(:,2)),mean(ty(:,2)),mean(tz(:,2)),s2(1,2),s2(2,2),s2(3,2),sc,'color','c')
%     
%     quiver3(mean(tx(:)),mean(ty(:)),mean(tz(:)),ms2(1),ms2(2),ms2(3),sc,'color','k')
% %     quiver3(mean(tx(:)),mean(ty(:)),mean(tz(:)),o1(1),o1(2),o1(3),sc,'color','y')
%     drawnow


end

m.HexEle = HexEle;


end



























