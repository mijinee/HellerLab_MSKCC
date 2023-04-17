% correctBBmovie.m
% Prakrit V. Jena, 2017 - Made in USA
% Code will correct h5 movies taken by PhySpec and correct using the
% dark (laser off) image, and the non-uniform correction image 
% (created using createBB)

% This version is for movies acquired using the sequencer on the HS
% Specifically - data sets are: [60 x 1040 x 1392] for visible 

%% Obtain h5 file names, remember current path
purge;
clc;
oldPath = cd;

initDirectory=uigetdir('raw');
cd(initDirectory);

[multiBBOpen,path] = uigetfile('*.h5', 'Select the BB movies to process: ','MultiSelect','on');

if isequal(iscellstr(multiBBOpen),0)
    multiBBOpen = cellstr(multiBBOpen);
end

%% Obtain dark image correction files, set status flags if some are missing
[darkh5, darkh5path, ~] = uigetfile('*.h5','Please select the Dark Image for this emission channel: ');

if isequal(darkh5,0)
   disp('No Dark Image present for this channel!')
   darkImageCheck = false;
else
   disp('User selected Dark Image:');
   disp(fullfile(darkh5path, darkh5));
   darkImageCheck = true;
end

%% Obtain non-uniform correction files, set status flags if some are missing
% disp('Please select the Dark Image for this emission channel: ');
[nucih5, nucih5path, ~] = uigetfile('*.h5','Please select the Non-Uniform Correction Image for this emission channel: ');

if isequal(nucih5,0)
   disp('No Non-Uniform Correction Image present for this channel!')
   nucImageCheck = false;
else
   disp('User selected Non-Uniform Correction Image:');
   disp(fullfile(nucih5path, nucih5));
   nucImageCheck = true;
end

%% Verify data format of files
% Now different versions of Physpec/movie vs z stack might have different
% data forms. I shall make cases as needed.

% ii = 1;
% tempMovieName = multiBBOpen{1,ii};
% tempMovieData = h5read(tempMovieName,'/Cube/Images');
% In practice, test for existence; this can be case 1
% [dim1,dim2,dim3] = size(tempMovieData);

% Make a loop for each movie
%% Correct h5 file
nFiles = size(multiBBOpen,2);

for ii = 1:nFiles
    tempMovieName = multiBBOpen{1,ii};
    tempMovieData = double(h5read(tempMovieName,'/Cube/Images')); % double
    tempNumFrames = size(tempMovieData,3);
    tempDim1 = size(tempMovieData,1);
    tempDim2 = size(tempMovieData,2);
    
    % Test for dark image file to use
    
    if darkImageCheck == 1
         dataDarkImage = double((h5read(strcat(darkh5path,darkh5),'/Image/Data')));
    else dataDarkImage = zeros(tempDim1,tempDim2);
    end
    
    % Test for non uniform correction file to use

    if nucImageCheck == 1
         dataNucImage = (h5read(strcat(nucih5path,nucih5),'/Image/Data')); % double
    else dataNucImage = ones(tempDim1,tempDim2);
    end
        
    dataDarkImageStack = repmat(dataDarkImage,[1 1 tempNumFrames]);
    dataNucImageStack = repmat(dataNucImage, [1 1 tempNumFrames]);
   
    % Correcting the tempMovieData stack
   
    tempMovieDataSubtracted = bsxfun(@minus,tempMovieData,dataDarkImageStack);
    tempMovieDataSubtractedNuc = bsxfun(@rdivide,tempMovieDataSubtracted,dataNucImageStack);
    
    tempOutputData = uint16(tempMovieDataSubtractedNuc);
    tempWriteName = [tempMovieName(1:end-3) ' - Dark + NUC.h5'];
    tempWriteFullName = fullfile(oldPath,tempWriteName);
    
    h5create(tempWriteFullName,'/Image/Data',[tempDim1 tempDim2 tempNumFrames],'Datatype','uint16');
    h5write(tempWriteFullName,'/Image/Data',tempOutputData);
    disp(tempWriteFullName);
    disp('Subtracted, non uniform correction applied, and saved.');

    clear ii tempMovieName tempMovieData tempNumFrames tempDim1 tempDim2 dataDarkImage;
    clear dataNucImage dataDarkImageStack dataNucImageStack tempMovieDataSubtracted;
    clear tempMovieDataSubtractedNuc tempWriteName tempWriteFullName tempOutputData; 
end


%% Close variables
clear all;
close all;

disp('All files corrected and written.');


