% correctBB  - script
% Obtains a list of h5 files with /Image/Data/ present, a h5 Laser Off file,
% and a /Image/Data containing file with non-uniform excitation.
% Saves the corrected h5 file
% Correction is: newData = (oldData - Laser Off) / (non-uniform image)
% If any of the correction files isn't available, press cancel


purge;
%% Load a sequence of BB files:
initDirectory=uigetdir('raw');
cd(initDirectory);

[multiBBOpen,path] = uigetfile('*BB*.h5', 'Select the BB files to process (excluding blank)','MultiSelect','on');


if isequal(iscellstr(multiBBOpen),0)
    multiBBOpen = cellstr(multiBBOpen);
end

%% Choose blank image with all light sources off.
[darkh5, darkh5path, ~] = uigetfile('*.h5','Choose blank image all light sources off');

if isa(darkh5,'double') == 1 % incase the user presses cancel, no blank is available
    imageDark = uint16(zeros(256,320));
   
else
    imageDark = h5read(strcat(darkh5path,darkh5),'/Image/Data');
    %imageDark = h5read(strcat(darkh5path,darkh5),'/Cube/Images');
end

clear darkh5;

%% Choose image for non-uniform illumination 
cd(darkh5path);
[normalizedh5, normh5path, ~] = uigetfile('*.h5','Choose the NonUniform 32bit BB image generated with createBB');

if isa(normalizedh5,'double') == 1 % incase the user presses cancel, no blank is available
    imageNormalize = double(ones(256,320));
    
else
    imageNormalize = h5read(strcat(normh5path,normalizedh5),'/Image/Data');
end

clear normalizedh5 normh5path;
cd(path);

%% Loop for each BB image
cd(initDirectory);

nFiles = size(multiBBOpen,2);
for i = 1:nFiles
    tempBBName = multiBBOpen{1,i};
    tempBBData = h5read(tempBBName,'/Image/Data');
    %tempBBData = h5read(tempBBName,'/Cube/Images');
    tempSubtracted = tempBBData - imageDark;
    tempNormalized = double(tempSubtracted)./imageNormalize;
    figure('name',tempBBName); imshow(tempNormalized,[],'Border','tight');
    %% imshow(image,[]) uses min(I(:)) and max(I(:)) for its low and high values
%     lowI = min(tempNormalized(:));
%     highI = max(tempNormalized(:));
    h1 = imshow(tempNormalized,[]);
    h2 = imcontrast(gca);
    waitfor(h2);
    imgNormalizedAdjust = getimage(gcf);
    close(gcf);
    imgName = [tempBBName(1:end-3) ' Nice.png'];
    imwrite(imgNormalizedAdjust,imgName);
    imgBBName = [tempBBName(1:end-3) ' Nice.h5'];
    imgBBName = fullfile(path,imgBBName);
    h5create(imgBBName,'/Image/Data',[256 320],'Datatype','double');
    h5write(imgBBName,'/Image/Data',tempNormalized);
end

disp(['Written: ' num2str(i) ' Files.']);
purge;
