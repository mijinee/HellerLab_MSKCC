% converWLJB  - script to convert /Cube/Images (1x320x256) h5 into png
% Saves a png file 

purge;
%% Load a sequence of WL files:

initDirectory=uigetdir('raw');
cd(initDirectory);

[multiBBOpen,path] = uigetfile('*WL*.h5', 'Select the WL files to process','MultiSelect','on');

if isequal(iscellstr(multiBBOpen),0)
    multiBBOpen = cellstr(multiBBOpen);
end
cd(path);

%% Loop for each BB image
cd(path);
nFiles = size(multiBBOpen,2);
for i = 1:nFiles
    tempBBName = multiBBOpen{1,i};
    tempBBData = h5read(tempBBName,'/Cube/Images'); %256 x 320 x uint16
    h1 = imshow(tempBBData,[]);
    h2 = imcontrast(gca);
    waitfor(h2);
    imgNormalizedAdjust = getimage(gcf); % This is a uint16, as without
    % background subtraction, the actual pixel values are huge.
    imgWrite = uint8 (imgNormalizedAdjust / 256);
    % Might as well write this as a png:
    imwrite(imgWrite ,[tempBBName(1:end-3) ' Nice.png']);

end

%% Close all
disp([num2str(i) ' Nice WL files saved.']); 
clear all; close all;
