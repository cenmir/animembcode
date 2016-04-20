function AddInputPaths()

curDir = pwd;
mainCodeDir = fileparts(curDir);
inputDir = fullfile(mainCodeDir,'input');

%% Add required paths below
HexMeshDir = fullfile(inputDir,'HexMeshDev','Hex1Mesh');
xFigureDir = fullfile(inputDir, 'xFigure');
otherDir = fullfile(inputDir, 'other');

addpath(HexMeshDir,xFigureDir,otherDir)

end

