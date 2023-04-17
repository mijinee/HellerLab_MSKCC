% Z:\Software\Matlab\mouseImager
% correctBBmouseImager
% This is specifically for taking the average projection of a short BB
% movie (to compensate for mouse breathing), and applying corrections to it

%% Obtain h5 file names, remember current path
purge;
clc;
oldPath = cd;

initDirectory=uigetdir('raw');
cd(initDirectory);

[multiBBOpen,path] = uigetfile('*BB*.h5', 'Select the BB movies to process: ','MultiSelect','on');

if isequal(iscellstr(multiBBOpen),0)
    multiBBOpen = cellstr(multiBBOpen);
end

%% Obtain blank image, required for mouse data, set status flags if some are missing
[blankh5, blank5path, ~] = uigetfile('*.h5','Please select the Blank Image for this emission channel: ');
dataBlankImage = double((h5read(strcat(blank5path,blankh5),'/Image/Data'))); % double


%% Obtain non-uniform correction files, set status flags if some are missing
% disp('Please select the Dark Image for this emission channel: ');
[nucih5, nucih5path, ~] = uigetfile('*.h5','Please select the Non-Uniform Correction Image the mouse imager: ');
dataNucImage = (h5read(strcat(nucih5path,nucih5),'/Image/Data')); % double


%% Correct h5 file
nFiles = size(multiBBOpen,2);

for ii = 1:nFiles
    tempMovieName = multiBBOpen{1,ii};
    tempMovieData = double(h5read(tempMovieName,'/Cube/Images')); % double
    
    % tempNumFrames = size(tempMovieData,3);
    
    tempDim1 = size(tempMovieData,1);
    tempDim2 = size(tempMovieData,2);
    
    % Projecting the average of the stack
    tempMovieDataProj = mean(tempMovieData,3);
    % imshow(tempMovieDataProj,[]); % Image check
    
    % Correcting the tempMovieDataProj
    tempMovieDataProjSub = bsxfun(@minus,tempMovieDataProj,dataBlankImage);
    
    % Compensating for non-uniform correction
    tempMovieDataProjSubNorm = bsxfun(@rdivide,tempMovieDataProjSub,dataNucImage);
    
    % saving files 
 
    tempMovieDataOutput = (tempMovieDataProjSubNorm);

    tempWriteName = [tempMovieName(1:end-3) ' - Blank + NUC.h5'];
    tempWriteFullName = fullfile(oldPath,tempWriteName);

    h5create(tempWriteFullName,'/Image/Data',[tempDim1 tempDim2],'Datatype','double');
    h5write(tempWriteFullName,'/Image/Data',tempMovieDataOutput);

    disp('Subtracted, non uniform correction applied, and saved.');

    clear ii tempMovieName tempMovieData tempNumFrames tempDim1 tempDim2 ;
    clear tempMovieDataProj tempMovieDataProjSub tempMovieDataProjSubNorm;
    clear tempWriteName tempWriteFullName; 
end


%% Close variables
clear all;
close all;

disp('All files corrected and written.');
