% correctWL  - script
% Obtains a list of h5 files with /Image/Data/ present, a h5 WL Off file,
% and a /Image/Data containing file with non-uniform excitation.
% Saves a png file and the corrected h5 file
% Correction is: newData = (oldData - WLOff) / (Smoothed non-uniform image)
% If any of the correction files isn't available, press cancel

purge;
%% Load a sequence of WL files:

initDirectory=uigetdir('raw');
cd(initDirectory);

[multiBBOpen,path] = uigetfile('*WL*.h5', 'Select the WL files to process (excluding blank)','MultiSelect','on');


if isequal(iscellstr(multiBBOpen),0)
    multiBBOpen = cellstr(multiBBOpen);
end
cd(path);
%% Choose blank image with all light sources off.

[darkh5, darkh5path, ~] = uigetfile('*.h5','Choose blank WL image all light sources off');

if isa(darkh5,'double') == 1 % incase the user presses cancel, no blank is available
    imageDark = uint16(zeros(512,320));
   
else
    imageDark = h5read(strcat(darkh5path,darkh5),'/Cube/Images');
end

clear darkh5;


%% Choose image for non-uniform illumination 
if (darkh5path == 0)
    darkh5path = path;
end

cd(darkh5path);
[normalizedh5, normh5path, ~] = uigetfile('*.h5','Choose the NonUniform 32bit WL image generated with createWL');


if isa(normalizedh5,'double') == 1 % incase the user presses cancel, no blank is available
    imageNormalize = double(ones(256,320));
    
else
    imageNormalize = h5read(strcat(normh5path,normalizedh5),'/Image/Data');
end

clear normh5path normalizedh5;
%% Loop for each BB image
cd(path);
nFiles = size(multiBBOpen,2);
for i = 1:nFiles
    tempBBName = multiBBOpen{1,i};
    tempBBData = h5read(tempBBName,'/Cube/Images');
    tempSubtracted = tempBBData - imageDark;
    
    h = fspecial('average',8);
    imageNormSmooth = imfilter(imageNormalize,h,'replicate');
    
    tempNormalized = double(tempSubtracted)./imageNormSmooth;
    h1 = imshow(tempNormalized,[]);
    h2 = imcontrast(gca);
    waitfor(h2);
    imgNormalizedAdjust = getimage(gcf);
    % Might as well write this as a png:
    imwrite(imgNormalizedAdjust,[tempBBName(1:end-3) ' Nice.png']);
    %%
    %newBBName = ['Corrected - ' tempBBName];
    newBBName = [tempBBName(1:end-3) ' Nice.h5'];
    newBBWrite = fullfile(path,newBBName);
    h5create(newBBWrite,'/Image/Data',[256 320],'Datatype','double');
    h5write(newBBWrite,'/Image/Data',imgNormalizedAdjust);
      
    
end

%% Close all
disp([num2str(i) ' Nice WL files saved.']); 
purge;
