clear
close all
%%
% run([pwd,'\read_data2.m']);
% run([pwd,'\mkactintiffs.m']);

%%
folder = 'actin_tiffs';
direc = dir([pwd,'/',folder,'/','*','.','tif']); %Load 640 Images

[filenames{1:length(direc),1}] = deal(direc.name);%create a cell array that contains the names of each .tif image file in the current subdirectory
filenames = sortrows(filenames);%sort the filnames 
z=length(filenames);

for i=1:z
    [pathstr, name, ext] = fileparts(filenames{i});
    fileind(i) = str2num(name);
end

[sortind,idx] = sort(fileind);

for i=1:z
    filesort{i} = filenames{idx(i)};
end

for i=1:z
    tempRGB = imread([pwd,'/',folder,'/',filesort{i}]);
    % RGB images are converted to grayscale using the formula
    % gray=(red+green+blue)/3 or gray=0.299red+0.587green+0.114blue if
    % "Weighted RGB to Grayscale Conversion" is checked in
    % Edit>Options>Conversions.
    % https://imagej.nih.gov/ij/docs/menus/image.html
%     stack8bit(:,:,i) = 0.299.*tempRGB(:,:,1)+0.587.*tempRGB(:,:,2)+0.114.*tempRGB(:,:,3);
    stack8bit(:,:,i) = sum(tempRGB,3)/3;
end

stackinvert = 255 - stack8bit; % Invert Image
stack = uint16(stackinvert); % Doesn't actually do anything - I convert to 16bit in ImageJ for display purposes only

%% Calculate Output Metrics
dt = 2; % seconds

frame = 0:1:size(stack,3)-1;
time = frame.*dt;
ID = sum(reshape(stack,[size(stack,1)*size(stack,2) size(stack,3)]))';
strain = 1-ID./ID(1);
strainrate = diff(strain)./dt;

outs = NaN(length(ID),3);
outs(:,1) = ID;
outs(:,2) = strain;
outs(2:end,3) = strainrate;

mets = NaN(1,3);
mets(1) = ID(1);
mets(2) = max(strain);
mets(3) = max(strainrate);

%% Save Metrics
save([pwd,'\SimulationMetrics.mat'],'ID','strain','strainrate','outs','mets','dt','frame','time');


    





