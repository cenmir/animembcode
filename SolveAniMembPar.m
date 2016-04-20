function [U,F,S,activeEq] = SolveAniMembPar(mesh,m,iteration)
% SolveAniMemb 
%   TODO: David, Explain this stuff

%% Preallocating
xnod = mesh.XC; ynod=mesh.YC; znod=mesh.ZC;
nodes = mesh.Connectivity;
surfh = mesh.SurfaceP1;
nele=length(surfh); %number of elements



%% Quad Gauss
[qgx,qgy, qgw] = QuadGauss(m.GOGhost); 
% Used for integrating Hex element sides in the stabilization part

%%
[gx,gy,gw]=trigauc([0,1,0],[0,0,1],m.GO); %Triangle element integration
CutEle = unique([surfh(:).iel])';
tb = m.tb; %Base material thickness
mu = m.mu;
lambda0 = m.lambda0;
alpha1 = m.alpha1; %Yung type elasticity coefficient for the fiber material
alpha2 = m.alpha2;

%% Assembly preallocation
ndof = 3;
knod = 8; %Hex element
nneigh = 6;
nsi = length(gw)*nneigh;
nnod = length(xnod);
locneq = ndof*knod; %24
neq = 3*nnod; % number of equations (max)


rowSe = ones((2*locneq)^2,nsi);
rowStab = ones(numel(rowSe),nele); 
colStab = rowStab; 
valStab = 0*rowStab;

rowS = ones((locneq)^2,nele); 
colS = rowS;
valS = zeros((locneq)^2,nele);

rowF = zeros(locneq,nele); valF = rowF; %"F"



% row = zeros(8*neq, 1);col = row;val = row;
% row = zeros(locneq^2*nele, 1);col = row;val = row;
% rowS = zeros((locneq*2)^2*nele*7*length(gw), 1);colS = rowS;valS = rowS;
% up = 0;
% upS = 0;
% F = zeros(neq,1);

% xfigure; hold on;
for iTri=1:nele
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

    fele=zeros(knod*ndof,1);
    % Loop over all triangle Gauss-points
    
    %% Local material Properties
    tf1 = m.HexEle(iel).tf1;
    tf2 = m.HexEle(iel).tf2;
    
    s1 = m.HexEle(iel).s1;
    s2 = m.HexEle(iel).s2;
    S1 = s1'*s1;
    S2 = s2'*s2;
    if iteration == 0
       tf1 = 0;
       tf2 = 0;
    else
        1;
    end
    t = tb+tf1+tf2;
    
    sele = zeros(knod*ndof);
    bw=zeros(ndof,knod*ndof);
    
    
    in = 1; %1:nsi for hex elements
    rowSe = ones((2*locneq)^2,nsi);
    colSe = rowSe;
    valSe = 0*rowSe;
    %% Loop over integration points
    for ig = 1:length(gw)
        r = gx(ig); s = gy(ig);
        igw = gw(ig);
        % P1 triangular basis function
        fi = [(1-r-s),r,s];
        dfidr = [-1,1,0];
        dfids = [-1,0,1];
        
        % Gauss point for mother element
        xig = [fi*tx, fi*ty, fi*tz];
        
        % Compute normal using dxdr and dxds
        dxdr = [dfidr*tx,dfidr*ty,dfidr*tz];
        dxds = [dfids*tx,dfids*ty,dfids*tz];
        n = cross(dxdr,dxds);
        n = n/norm(n);
        % Element normal. Used only to check orientation
        n_ele = surfh(iTri).ElementNormal;
        if sign(dot(n_ele,n)) < 0
            n = -n;
        end
        
        J = [dxdr;dxds;n];
        area = abs(det(J))/2;
        % area = norm(cross(TX(2,:)-TX(1,:),TX(3,:)-TX(1,:)))/2      
        
        [fi, fix, fiy, fiz, ~] = mesh.baseHexP1(iel, xig);
        
        %% Membrane model
        P=eye(3)-n'*n;
        
        %% Create constitutive tensor
        Em = zeros(ndof,ndof,ndof,ndof); % Preallocate some space for the constitutive tensor
        for i=1:3
            for j=1:3
                for k=1:3
                    for l=1:3
                        %Todo: David, Clean this shit up!
                        Ememb = tb/t*( mu*(P(i,k)*P(j,l)+P(i,l)*P(j,k)) + lambda0*P(k,l)*P(i,j) );
                        Efib1 = tf1/t*alpha1*S1(k,l)*S1(i,j);
                        Efib2 = tf2/t*alpha2*S2(k,l)*S2(i,j);
                        Em(i,j,k,l) = Ememb + Efib1 + Efib2;
                    end
                end
            end
        end

        %% Create strain tensor
        epsmat = zeros(2,2,ndof*knod);
        epsmat(1,1,1:3:end) = fix;
        epsmat(1,2,1:3:end) = fiy/2;
        epsmat(1,2,2:3:end) = fix/2;
        epsmat(1,3,1:3:end) = fiz/2;
        epsmat(1,3,3:3:end) = fix/2;
        epsmat(2,1,:) = epsmat(1,2,:);
        epsmat(2,2,2:3:end) = fiy;
        epsmat(2,3,2:3:end) = fiz/2;
        epsmat(2,3,3:3:end) = fiy/2;
        epsmat(3,3,3:3:end) = fiz;
        epsmat(3,1,:) = epsmat(1,3,:);
        epsmat(3,2,:) = epsmat(2,3,:);
        
        sigmat = 0*epsmat;
        for k=1:3
            for l=1:3
                sigmat(1,1,:) = sigmat(1,1,:) + Em(1,1,k,l)*epsmat(k,l,:);
                sigmat(1,2,:) = sigmat(1,2,:) + Em(1,2,k,l)*epsmat(k,l,:);
                sigmat(1,3,:) = sigmat(1,3,:) + Em(1,3,k,l)*epsmat(k,l,:);
                sigmat(2,1,:) = sigmat(2,1,:) + Em(2,1,k,l)*epsmat(k,l,:);
                sigmat(2,2,:) = sigmat(2,2,:) + Em(2,2,k,l)*epsmat(k,l,:);
                sigmat(2,3,:) = sigmat(2,3,:) + Em(2,3,k,l)*epsmat(k,l,:);
                sigmat(3,1,:) = sigmat(3,1,:) + Em(3,1,k,l)*epsmat(k,l,:);
                sigmat(3,2,:) = sigmat(3,2,:) + Em(3,2,k,l)*epsmat(k,l,:);
                sigmat(3,3,:) = sigmat(3,3,:) + Em(3,3,k,l)*epsmat(k,l,:);
            end
        end
        
        %% Double Contratction between sigma and epsilon
%         putvar(sigmat, epsmat)
        for i=1:3
            for j=1:3
                a1 = epsmat(i,j,:);
                a2 = sigmat(i,j,:);
                b1 = a1(:).';
                b2 = a2(:).';
                sele = sele + t * b1.'*b2*igw*area;
            end
        end
        
        %% Stabilization term
        % Gradients of mother element
        B  = zeros(9,3*8);
        B(1,1:3:end)=fix;
        B(2,1:3:end)=fiy;
        B(3,1:3:end)=fiz;
        B(4,2:3:end)=fix;
        B(5,2:3:end)=fiy;
        B(6,2:3:end)=fiz;
        B(7,3:3:end)=fix;
        B(8,3:3:end)=fiy;
        B(9,3:3:end)=fiz;

        nei = mesh.Neighs(iel,:);
        nei = nei(nei>0);
        nei = nei(ismember(nei,CutEle));
        for inei = nei
            ivn=nodes(inei,:);
            ieqsn = zeros(24,1);
            ieqsn(1:3:end) = ivn*3-2;
            ieqsn(2:3:end) = ivn*3-1;
            ieqsn(3:3:end) = ivn*3-0;
            
            % Compute Gauss points
            nod=intersect(iv,ivn);
            p1 = [xnod(nod(1)),ynod(nod(1)),znod(nod(1))];
            p2 = [xnod(nod(2)),ynod(nod(2)),znod(nod(2))];
            p3 = [xnod(nod(3)),ynod(nod(3)),znod(nod(3))];
            p4 = [xnod(nod(4)),ynod(nod(4)),znod(nod(4))];
            p12 = p2-p1;
            p13 = p3-p1;
            
            QX = [p1;p2;p4;p3];
            qx = QX(:,1); qy = QX(:,2); qz = QX(:,3);
            
            ivtot=[ieqs(:);ieqsn(:)];
            stabilTerm = zeros(((3*8)*2)^2,1);
            for iqg = 1:length(qgw)
                xsi = qgx(iqg); eta = qgy(iqg);
                qfi = [(1-xsi)*(1-eta)/4;
                      (1+xsi)*(1-eta)/4;
                      (1+xsi)*(1+eta)/4
                      (1-xsi)*(1+eta)/4]';
               %Gauss point
               xqig = [qfi*qx, qfi*qy, qfi*qz];
               
               [~,fixn,fiyn,fizn,~] = mesh.baseHexP1(inei, [xqig(1),xqig(2),xqig(3)]);
               Bn  = zeros(9,3*8); %
               Bn(1,1:3:end)=fixn;
               Bn(2,1:3:end)=fiyn;
               Bn(3,1:3:end)=fizn;
               Bn(4,2:3:end)=fixn;
               Bn(5,2:3:end)=fiyn;
               Bn(6,2:3:end)=fizn;
               Bn(7,3:3:end)=fixn;
               Bn(8,3:3:end)=fiyn;
               Bn(9,3:3:end)=fizn;
               
               nc=cross(p12,p13);
               
               % we have quadratic faces
               intersectArea = norm(nc);
               
               jump=[B,-Bn];
               
               je=jump'*jump;
               
%                h = area;
               
%                stabilTerm = stabilTerm + gamma1*t*jump(:)*h^3*intersectArea*qgw(iqg);
               stabilTerm = stabilTerm + m.gamma1*t*je(:)*intersectArea*qgw(iqg);

            end
            
            %% Stabilizing the Stiffness Matrix
            lenS = length(ivtot);
            AX = ivtot(:, ones(1, lenS));
            AY = AX';
            
            rowSe(:,in) = AX(:); 
            colSe(:,in) = AY(:); 
            valSe(:,in) = stabilTerm(:); 
            
%             
            
%             lenS = length(ivtot);
%             AXS = ivtot(:, ones(1, lenS));
%             AYS = AXS';
%             nnS = lenS^2;
%             loS = upS + 1;
%             upS = upS + nnS;
%             rowS(loS:upS) = AXS(:);
%             colS(loS:upS) = AYS(:);
%             valS(loS:upS) = stabilTerm;
            in = in +1;
        end
  
        %% Load
        bw(1,1:3:end)=fi;
        bw(2,2:3:end)=fi;
        bw(3,3:3:end)=fi;
        fele = fele+gw(ig)*area*bw'*m.Fload(xig(1),xig(2),xig(3),n(:));
        
    end
    
    %% Assemble stabilization term into element stabilization matrix
    rowStab(:,iTri) = rowSe(:);
    colStab(:,iTri) = colSe(:);
    valStab(:,iTri) = valSe(:);
    
    %% Edge Load
%     if sum(tx == max(surfX(:,1)))>=2
%         force = [1,0,0]';
%         indx = tx == max(tx);
%         edx = tx(indx);
%         edy = ty(indx);
%         edz = tz(indx);
%         edX = [edx,edy,edz];
%         le = norm(edX(1,:)-edX(2,:));
%         for ieg = 1:2
%             [fie, ~, ~, ~, ~] = mesh.baseHexP1(iel, edX(ieg,:) );
%             bwe(1,1:3:end)=fie;
%             bwe(2,2:3:end)=fie;
%             bwe(3,3:3:end)=fie;
%             fele=fele+0.5*le*bwe'*force;
%             
%         end
%         
%     end
    
    %% Assembly
    len = length(ieqs);
    AX = ieqs(:, ones(1, len));
    AY = AX';    
    rowS(:,iTri) = AX(:); 
    colS(:,iTri) = AY(:); 
    valS(:,iTri) = sele(:);
    rowF(:,iTri) = ieqs(:); 
    valF(:,iTri) = fele;
    
%     F(ieqs,1) = F(ieqs,1) + fele;
%     len = length(ieqs);
%     AX = ieqs(:, ones(1, len));
%     AY = AX';
%     nn = len * len;
%     lo = up + 1;
%     up = up + nn;
%     row(lo:up) = AX(:);
%     col(lo:up) = AY(:);
%     val(lo:up) = sele(:);
    
    
    
end
% row2 = [row(1:up);rowS(1:upS)];
% col2 = [col(1:up);colS(1:upS)];
% val2 = [val(1:up);valS(1:upS)];
% S = sparse(row2, col2, val2, neq, neq);
% activeEq = unique([row2(:);col2(:)]);


row = [rowS(:);rowStab(:)];
col = [colS(:);colStab(:)];
val = [valS(:);valStab(:)];
S = sparse(row(:), col(:), val(:), neq, neq);
% assemble final load vector:
F = zeros(neq,1);
for iTri=1:nele
    F(rowF(:,iTri)) = F(rowF(:,iTri)) + valF(:,iTri);
end

activeEq = unique([rowS(:);colS(:)]);

%% Dirichlet Boundary Conditions
% Boundary
% activeNodes = unique(ceil(activeEq/3));
[PN0, ~] = boundaryInds(mesh,{'x',m.x0});
[PN1, ~] = boundaryInds(mesh,{'x',m.x1});
PrescEqs = [3*PN0-2;3*PN0-1;3*PN0-0;  3*PN1-2;3*PN1-1;3*PN1-0];
FreeEqs = setdiff(activeEq,PrescEqs);
Sreduced=S(FreeEqs,FreeEqs); Freduced = F(FreeEqs);

%% Solve the system
U=zeros(neq,1);
Ureduced=full(Sreduced\Freduced);
U(FreeEqs)=Ureduced;

end

