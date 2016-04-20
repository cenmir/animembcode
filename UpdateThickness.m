function model = UpdateThickness(mesh, model)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
surfh = mesh.SurfaceP1;
nele = length(surfh);
L0 = 0.0001; L1 = 10^4;
tB1 = zeros(nele,1); 
tB2 = tB1;
eta = model.eta;
t1min = model.tf1min;
t1max = model.tf1max;
t2min = model.tf2min;
t2max = model.tf2max;
Vmax = model.Vmax;
HexEle = model.HexEle;
nHex = length(HexEle);

viz = 0;


if viz
    xfigure(45)
    hp = plot([L0,L1],[0,0],'k*');
    grid on
    hold on;
    axis auto
    plot([L0,L1],Vmax*[1,1],'k-')
end


c = 1;
disp('Finding Lambda...')
time1 = tic;
while (L1-L0) > 10^-4
%     disp(['Iteration: ',num2str(c)])
    Lm = (L0+L1)/2;

    %% Update thickness
    for iel = unique([surfh.iel])
        tf1 = HexEle(iel).tf1;
        tf2 = HexEle(iel).tf2;
        b1 = HexEle(iel).b1;
        b2 = HexEle(iel).b2;
        
        tb1 = tf1 * (1/Lm*b1)^eta;
        tb2 = tf2 * (1/Lm*b2)^eta;
        
        if tb1 <= t1min
            HexEle(iel).tf1new = t1min;
        elseif tb1 >= t1max
            HexEle(iel).tf1new  = t1max;
        else
            HexEle(iel).tf1new = tb1;
        end
        
        if tb2 <= t2min
            HexEle(iel).tf2new  = t2min;
        elseif tb2 >= t2max
            HexEle(iel).tf2new  = t2max;
        else
            HexEle(iel).tf2new  = tb2;
        end
        
    end

    %% Check volume constraint
%     disp('Checking volume constraint...')
%     tic
    Vtest = 0;
    for iTri = 1:nele
       A = surfh(iTri).Area;
       iel = surfh(iTri).iel;
       tf1new = HexEle(iel).tf1new;
       tf2new = HexEle(iel).tf2new;
       Vtest = Vtest + (tf1new+tf2new)*A;
    end
%     toc
    
    if viz
        try
            delete(hlo)
            delete(hup)
    %         delete(hmi)
        end
        hlo = plot([L0,L0],[0,Vmax],'b*-');
        hup = plot([L1,L1],[0,Vmax],'r*-');
    end
    
    if Vtest - model.Vmax > 0
        L0 = Lm;
    else
        L1 = Lm;
    end
    
    if viz 
        hmi = plot(Lm,Vtest,'bo');
        title(['Iteration: ',num2str(c),', \Lambda=',num2str(Lm)])
%         pause(0.1)
    end
    c = c+1;
end
disp('Total Lambda time:')
toc(time1)

model.tf1 = tf1new;
model.tf2 = tf2new;

for iel = unique([surfh.iel])
   HexEle(iel).tf1 = HexEle(iel).tf1new;
   HexEle(iel).tf2 = HexEle(iel).tf2new;
end
model.HexEle = HexEle;


end

