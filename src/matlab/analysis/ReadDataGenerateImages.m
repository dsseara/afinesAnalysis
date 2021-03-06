clear
close all
%%
run('run_read_data.m');     % Create simdata.mat file
run('velInterp.m');         % Interpolate velocity and displacement fields to a grid
run('divVelocity.m');       % Get the divergence of the velocity and displacement fields to get total normal strain and strain rates
run('divStats.m')           % Get statistics, max and 25-75 slope of strain
run('quiverOverlay.m');     % Generate images, quiver overlaid on filaments and on heatmap of velocity divergence


%run('mkactintiffs.m');

% %%
% folder = 'actin_tiffs';
% direc = dir([pwd,'/',folder,'/','*','.','tif']); %Load 640 Images
% 
% [filenames{1:length(direc),1}] = deal(direc.name);%create a cell array that contains the names of each .tif image file in the current subdirectory
% filenames = sortrows(filenames);%sort the filnames 
% z=length(filenames);
% 
% for i=1:z
%     [pathstr, name, ext] = fileparts(filenames{i});
%     fileind(i) = str2num(name);
% end
% 
% [sortind,idx] = sort(fileind);
% 
% for i=1:z
%     filesort{i} = filenames{idx(i)};
% end
% 
% for i=1:z
%     stack(:,:,:,i) = 
%     
% end
%     
% 
% 
% a2 = im2uint8(a);
% % http://matlab.izmiran.ru/help/techdoc/creating_plots/chimag15.html
% I = .2989*rgb_img(:,:,1)+.5870*rgb_img(:,:,2)+.1140*rgb_img(:,:,3);


