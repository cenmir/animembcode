function [hs,Data] = vizSurface(U,F,surfMesh,mesh,scale,Data,k)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

surfh = surfMesh.surfh;
surfX = surfMesh.P;
surfTri = surfMesh.tri;
model = Data(k).model;


disp('Transfering displacements to triangles...')
tic
[SU,SF] = Project3DToSurface(surfh,mesh,F,U);
toc
nfig = 3;


%% Fiber Direction
[hs.fig,xf1] = xfigure(k*nfig-2);
% 
% hs.patch = patch('Faces',tri,'Vertices',SX1, 'FaceVertexCData', CData,'FaceColor','interp');
hs.patch = patch('Faces',surfTri,'Vertices',surfX, 'FaceColor','w');
axis equal; view(3); h1.light(1) = light;

xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; hold on
title(['Undeformed surface with fiber directions. Iteration: ',num2str(k)])
% colormap jet
% h1.quiver = quiver3(SX1(:,1),SX1(:,2),SX1(:,3),SF(:,1),SF(:,2),SF(:,3),1);
% h1.fig.KeyPressFcn = {@GKPF,xf1,h1};
% title(['\gamma_1=',num2str(gamma1),', GO=',num2str(GO),',ne=',num2str(ne)])

nHex = length([surfh.iel]);
XM = zeros(nHex,3);
S1 = XM; S2 = XM;
hh = zeros(nHex,1);

nTri = length(surfh);
t1tri = zeros(nTri,1);
t2tri = t1tri;
c = 1;
for iel = [surfh.iel]
   iTri = model.HexEle(iel).iTri;
%    nTri = length(iTri);
   t1tri(iTri) = model.HexEle(iel).tf1;
   t2tri(iTri) = model.HexEle(iel).tf2;
   polygonPoints = unique(reshape([surfh(iTri).Xe].',3,[]).','rows');
   Xm = mean(polygonPoints);
   s1 = model.HexEle(iel).s1;
   s2 = model.HexEle(iel).s2;
   h = sum([surfh(iTri).Area])^.5;
   
   XM(c,:) = Xm;
   S1(c,:) = s1;
   S2(c,:) = s2;
   hh(c) = h;
   c = c+1;
end
h = mean(hh);
hs.hq1 = quiver3(XM(:,1),XM(:,2),XM(:,3),S1(:,1),S1(:,2),S1(:,3),h*scale,'Color','r');
hs.hq2 = quiver3(XM(:,1),XM(:,2),XM(:,3),S2(:,1),S2(:,2),S2(:,3),h*scale,'Color','b');


%% Fiber Density
[hs.fig,xf1] = xfigure(k*nfig-1);
CData1 = t1tri;
SX = surfX+SU*scale;
hs.patch = patch('Faces',surfTri,'Vertices',SX, 'CData', CData1,'FaceColor','Flat');
axis equal; view(3); h1.light(1) = light;
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; hold on
titleTxt = ['Deformed Surface Fiber dir 1. Iteration: ',num2str(k)];
fig = gcf;
set(fig, 'Name', titleTxt, 'NumberTitle','off');
fig.Color = 'w';
colorbar


[hs.fig,xf1] = xfigure(k*nfig-0);
CData2 = t2tri;
hs.patch = patch('Faces',surfTri,'Vertices',SX, 'CData', CData2,'FaceColor','Flat');
axis equal; view(3); h1.light(1) = light;
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; hold on
titleTxt = ['Deformed Surface Fiber dir 2. Iteration: ',num2str(k)];
fig = gcf;
set(fig, 'Name', titleTxt, 'NumberTitle','off');
fig.Color = 'w';
colorbar



%% Save Data
Data(k).S1 = S1;
Data(k).S2 = S2;
Data(k).t1E = t1tri;
Data(k).t2E = t2tri;
Data(k).surfX = surfX;
Data(k).scale = scale;
Data(k).SU = SU;
Data(k).surfTri = surfTri;



end

