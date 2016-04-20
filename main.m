clear % NEVER clear all! Performence issues will ensue
clc
close all

% begin by adding the needed path for this session
AddInputPaths()

%% Input parameters for the cylinder
ref = 2
ne=2^(ref+1); %4,8,16,32
L=4;
pressure=1;
x0=0;x1=L;y0=-1.1;y1=1.1;z0=-1.1;z1=1.1;
R=1;
E=100;nu=0.5;
gamma1 = 1;
GO = 2;
GOGhost = 2;
tb=10^-2;

lambda0=E*nu/(1-nu^2);
% lambda0=E*nu/( (1+nu)*(1-2*nu) );
mu=E/(2*(1+nu));
tf1=10^-3;
tf2=10^-3;
alpha1 = 10*E;
alpha2 = 10*E;
CompTol = 0.00001;
model.tf1min = 0;
model.tf1max = 1;
model.tf2min = 0;
model.tf2max = 1;
model.eta = 0.5;
model.maxThickness = tb;
%% Load
he = ((x1-x0)/(2*ne)) * ((y1-y0)/ne) * ((z1-z0)/ne);
H(ref)=he;
Fload = @(x,y,z,n) pressure*he*(x-x0)*(x1-x)*n(:); %Quadratic load 

%% Create Mesh
nxe = 2*ne; nye = ne; nze = ne;
mesh = Hex1Mesh(x0,x1,nxe,y0,y1,nye,z0,z1,nze);
% hv = mesh.vizMesh('ElementNumbers','NodeNumbers');
mesh.Neighbors('Structured');

%% Surface function
xc = mean([x0,x1]); yc = mean([y0,y1]); zc = mean([z0,z1]);
surfaceFunction = @(x,y,z) ((y-yc).^2+(z - zc).^2).^0.5-R;

%% Define model struct
model.L =L;
model.pressure=pressure;
model.x0=x0;model.x1=x1;
model.y0=y0;model.y1=y1;
model.z0=z0;model.z1=z1;
model.xm = mean([x0,x1]);
model.ym = mean([y0,y1]);
model.zm = mean([z0,z1]);
model.R=R;
model.ne=ne;
model.E=E;
model.nu=nu;
model.GO=GO;
model.GOGhost=GOGhost;
model.tb=tb;
model.tf1=tf1;
model.tf2=tf2;
model.lambda0=lambda0;
model.mu=mu;
model.Fload=Fload;
model.surfaceFunction = surfaceFunction;
model.nxe=nxe;
model.nye=nye;
model.nze=nze;
model.gamma1 = gamma1;
model.alpha1=alpha1;
model.alpha2=alpha2;


%% Discrete surface
xnod = mesh.XC; ynod = mesh.YC; znod = mesh.ZC;
phi = surfaceFunction(xnod,ynod,znod);

%% Extract surface
surfh = mesh.CutP1(phi,0);
% TODO: create class of surfhMesh, add viz stuff and interpolation rutine
[tri,surfX] = mesh.TriangulateP1();
surfMesh.P = surfX;
surfMesh.tri = tri;
surfMesh.surfh = surfh;
surfMesh.nele = length(surfh);
% h = median([surfh.Area])^0.5

scale = 10^(ref-1);

%% Surface area
SurfaceArea = 0;
for iTri = 1:length(surfh)
    SurfaceArea = SurfaceArea + surfh(iTri).Area;
end
%% Volume Constraint
model.Vmax = SurfaceArea*model.maxThickness;
%% Viz Surface
% [h1.fig,xf1] = xfigure;
% h1.patch = patch('faces',tri,'vertices',surfX,'FaceColor','b','FaceLighting','flat','EdgeColor','k'); hold on
% axis equal; view(3); h1.light(1) = light;
% % h1.fig.KeyPressFcn = {@GKPF,xf1,h1};
% ele = 1:mesh.nele;
% Visualize.plotMesh(mesh,ele,6); axis tight equal;

%% Shits getting real here...
% Start the optimization of the fiber thickness and orientation by initializing
ntri = length(surfh);
% Two fiberthickness per element.
HexEle(max(unique([surfh.iel]))).tf1 = 0;
for iel = unique([surfh.iel])
    HexEle(iel).tf1 = model.tf1;
    HexEle(iel).tf2 = model.tf2;
    HexEle(iel).s1  = [0,0,0];
    HexEle(iel).s2  = [0,0,0];
end
model.HexEle = HexEle;
% Two fiber directions per element
model.s1 = zeros(ntri,3);
model.s2 = model.s1;
k = 0;

%% Membrane problem
disp('Solving membrane problem...')
tic
[Uk,Fk,Sk] = SolveAniMembPar(mesh,model,k);
toc
%%
[hs.fig,xf1] = xfigure(334);
[SU,SF] = Project3DToSurface(surfh,mesh,Fk,Uk);
SX = surfX+SU*scale;
CData = sqrt(SU(:,1).^2+SU(:,2).^2+SU(:,3).^2);
hs.patch = patch('Faces',tri,'Vertices',SX, 'FaceVertexCData', CData,'FaceColor','Interp');
axis equal; view(3); h1.light(1) = light;
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; hold on
title(['Deformed Surface Fiber dir 1. Iteration: ',num2str(k)])
colorbar

%% Visualize Solution on Hex1Mesh
% CData = sqrt(U(1:3:end).^2+U(2:3:end).^2+U(3:3:end).^2);
% V = [U(1:3:end),U(2:3:end),U(3:3:end)];
% FN = [F(1:3:end),F(2:3:end),F(3:3:end)];
% 
% ele = 1:size(mesh.Connectivity,1); ele = ele(:);
% fele = [6*ele-5;6*ele-4;6*ele-3;6*ele-2;6*ele-1;6*ele-0;];
% Faces = mesh.Faces(fele(:),:);
% Vertices = [mesh.XC, mesh.YC, mesh.ZC];
% 
% % Load on Hex1Mesh
% xfigure(50);clf; axis equal;
% patch('Faces',Faces,'Vertices',Vertices,'FaceColor','none'); hold on;
% quiver3(xnod,ynod,znod,FN(:,1),FN(:,2),FN(:,3),4,'Color','b');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% patch('faces',tri,'vertices',surfX,'FaceColor','c','FaceLighting','flat','EdgeColor','k'); hold on
% % plot3(H.XC(FreeNodes),H.YC(FreeNodes),H.ZC(FreeNodes),'*r')
% view(0,0)
% 
% % Displacement on Hex1Mesh
% xfigure(55);clf; axis equal; 
% patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData', CData,'FaceColor','interp','FaceAlpha',0.4); hold on;
% colormap jet;
% quiver3(xnod,ynod,znod,V(:,1),V(:,2),V(:,3),4,'Color','r');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% view(0,0)
% 
% return

%% Check Compliance
Comp0 = inf;
dComp = inf;

%% Start Optimizing shit yo
%% Yeah Mr White, Yeah Science!
hcf = xfigure(455);
hcf.Name = 'Convergence of the compliance';
hp = plot([0,0],[0,0],'b-*');
while dComp > CompTol
    k = k+1;
    model.k = k;
    Compk = Fk'*Uk
    dComp = abs(Comp0-Compk)
    Comp0 = Compk
    
    
%     model = OptimizeVar(S,U,F,mesh,model,k);
%     [model, element] = ComputeElementInfo(Uk,mesh,model);
    model = ElementDirections(Uk,mesh,model);
    
    model = UpdateThickness(mesh, model);
    
    % Set Lambda0 = 0.0001, Lambda1 = 10^4
    % Loop until Volume constraint satisfied
    %   Recompute Lambda using bisection method
    %   Recompute t1new and t2new
    % 
    % Use t1new and t2new to compute new U, mesh and model
    % Compute the complinece and check stopping criterion
    % 
    
    Data(k).model = model;
    Data(k).U = Uk;
    Data(k).F = Fk;
    Data(k).Comp = Compk;
    
    xcomp = [Data.Comp];
    xcomp = xcomp(:);
    hp.YData = xcomp;
    hp.XData = 1:k;
    
    
    
    %% Solve
    disp('Solving membrane problem...')
    tic
    [Uk1,Fk1,Sk1] = SolveAniMembPar(mesh,model,k);
    toc
    %% Visualize
    
    [hs,Data] = vizSurface(Uk1,Fk1,surfMesh,mesh,scale,Data,k);
    drawnow
    
    %% Update
    Fk = Fk1;
    Uk = Uk1;
    Sk = Sk1;
    
    %% ConvPlot
    
end

%% Plot 
nit = length(Data);
xfigure
xx = 1:nit;
yy = [Data.Comp]'
plot(xx,yy,'b-o')
axis auto