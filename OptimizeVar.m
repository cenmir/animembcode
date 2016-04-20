function [m] = OptimizeVar(S,U,F,mesh,m,k)
%OptimizeVar Summary of this function goes here
%   TODO: David, put some info here on hat this does...
% for iel = 1:nele
% Guess a Lambda
% Compute B1 and B2
% Local Bisect method to get Lambda
% update t1 and t2
% Check Volume constraint
%
nodes = mesh.Connectivity;
surfh = mesh.SurfaceP1;
[gx,gy,gw]=trigauc([0,1,0],[0,0,1],m.GO); %Triangle element integration
mu = m.mu;
lambda0 = m.lambda0;
nele=length(surfh); % number of triangular elements
for iTri = 1:nele
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
    for ig = 1:length(gw)
        r = gx(ig); s = gy(ig);
        igw = gw(ig);
        % P1 triangular basis function
        fi = [(1-r-s),r,s];
        dfidr = [-1,1,0];
        dfids = [-1,0,1];
        
        dxdr = [dfidr*tx,dfidr*ty,dfidr*tz];
        dxds = [dfids*tx,dfids*ty,dfids*tz];
        
        % Gauss point for mother element
        xig = [fi*tx, fi*ty, fi*tz];
        n = surfh(iTri).faceNormal;
        J = [dxdr;dxds;n];
        area = abs(det(J))/2;
        
        [fi, fix, fiy, fiz, ~] = mesh.baseHexP1(iel, xig);
        P=eye(3)-n'*n;
        
        phid=P*[fix';fiy';fiz'];
        fixG=phid(1,:);fiyG=phid(2,:);fizG=phid(3,:);
        
        u = U(ieqs);
        ux = u(1:3:end); uy = u(2:3:end); uz = u(3:3:end);
        
%         eps = [fix'*ux, (fix'*uy+fiy'*ux)/2, (fix'*uz+fiz'*ux)/2;
%                (fix'*uy+fiy'*ux)/2, fiy'*uy, (fiy'*uz+fiz'*uy)/2;
%                (fix'*uz+fiz'*ux)/2, (fiy'*uz+fiz'*uy)/2, fiz'*uz];
           
        GradS=[fixG*ux, fiyG*ux, fizG*ux;
               fixG*uy, fiyG*uy, fizG*uy;
               fixG*uz, fiyG*uz, fizG*uz];
        
        %Use epsS
        epsS=(GradS+GradS')/2
%         epsP=P*epsS*P
        divP=fixG*ux+fiyG*uy+fizG*uz;
        div=fix'*ux+fiy'*uy+fiz'*uz;
        sig = 2*mu*eps+lambda0*div*eye(3)
        e=eig(sig)
        sigP=P*(2*mu*epsS+lambda0*divP*eye(3))*P
        [V,D]=eig(sigP)
        
%         trace(eps*eps)
        epsContract = trace(epsS*epsS)
        
        xfigure(22)
        hp = patch(tx,ty,tz,'b'); hold on; axis equal
        txm=mean(tx); tym=mean(ty); tzm=mean(tz);
        
        sc = 0.1
        V = V*sc;
        quiver3(txm,tym,tzm,V(1,1),V(2,1),V(3,1),1,'color','r')
        quiver3(txm,tym,tzm,V(1,2),V(2,2),V(3,2),1,'color','g')
        quiver3(txm,tym,tzm,V(1,3),V(2,3),V(3,3),1,'color','c')
        
        s1 = V(:,1)/norm(V(:,1))
        s2 = V(:,2)/norm(V(:,2))
        cosTheta = s1.'*s2/(norm(s1)*norm(s2))
        thetad = acos(cosTheta)*180/pi
        pause
        
        
    end
end



end

