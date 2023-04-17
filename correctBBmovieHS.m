% correctBBmovie.m ver 1.0
% Prakrit V. Jena, 2018 - Made in USA
% Code will correct h5 movies taken by PhySpec and correct using the
% dark (laser off) image, and the non-uniform correction image 
% (created using createBB)
%
% Output file is 64-bit
%
% All movies being corrected must be of the same type
% i.e. NIR movies with dark image and laser on image (same dimensions)
% or visible images with corresponding correction files

% This version is for movies acquired using the sequencer on the HS
% Specifically - data sets are: [N x X x Y] 
%
% There are no dependencies on other functions and this function is
% stand-alone
%
% Expected updates for ver 1.1
% Output to be 16-bit, same as input

%% Obtain h5 file names, remember current path
clear all; close all; clc;
clc;
oldPath = cd;

initDirectory=uigetdir('raw');
cd(initDirectory);

[multiBBOpen,path] = uigetfile('*.h5', 'Select the BB movies to process: ','MultiSelect','on');

% If the user selects only one file, the code below converts the string
% input to a cell

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

[nucih5, nucih5path, ~] = uigetfile('*.h5','Please select the Non-Uniform Correction Image for this emission channel: ');

if isequal(nucih5,0)
   disp('No Non-Uniform Correction Image present for this channel!')
   nucImageCheck = false;
else
   disp('User selected Non-Uniform Correction Image:');
   disp(fullfile(nucih5path, nucih5));
   nucImageCheck = true;
end


%% Correct h5 files
nFiles = size(multiBBOpen,2);

for ii = 1:nFiles
    tempMovieName = multiBBOpen{1,ii};
    tempMovieData = double(h5read(tempMovieName,'/Cube/Images')); % double
    tempNumFrames = size(tempMovieData,3);
    tempDim1 = size(tempMovieData,1);
    tempDim2 = size(tempMovieData,2);
    
    % Test for dark image file to use
    
    if darkImageCheck == 1
         dataDarkImage = h5read(strcat(darkh5path,darkh5),'/Cube/Images');
        
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
    % Correction algorithm here is:
    % Data = Movie (x by y) x N
    % DarkImage is duplicated to an N stack
    % Non-Uniform Correction is duplicated to an N staack
    % The reason is, vector calculation on stacks is faster than looping
    % through each frame
    
    tempMovieDataSubtracted = bsxfun(@minus,tempMovieData,double(dataDarkImageStack));
    tempMovieDataSubtractedNuc = bsxfun(@rdivide,tempMovieDataSubtracted,dataNucImageStack);
    
   
    %% Case-dependent file naming
    if (darkImageCheck == 1) && (nucImageCheck == 1)
        disp('Dark + NUC');
        tempWriteName = [tempMovieName(1:end-3) ' - Dark + NUC.h5'];
    elseif (darkImageCheck == 1) && (nucImageCheck == 0)
            disp('Dark');
            tempWriteName = [tempMovieName(1:end-3) ' - Dark.h5'];
    elseif (darkImageCheck == 0) && (nucImageCheck == 1)
            disp('NUC');
            tempWriteName = [tempMovieName(1:end-3) ' - NUC.h5'];
    else disp('NO corrections applied');
            tempWriteName = [tempMovieName(1:end-3) ' - NO correction.h5'];
    end
    
    
    %%
    tempWriteFullName = fullfile(oldPath,tempWriteName);
    
    h5create(tempWriteFullName,'/Cube/Images',[tempDim1 tempDim2 tempNumFrames],'Datatype','double');
    h5write(tempWriteFullName,'/Cube/Images',tempMovieDataSubtractedNuc);
    disp(tempWriteFullName);
    disp('Subtracted, non uniform correction applied, and saved.');

    clear ii tempMovieName tempMovieData tempNumFrames tempDim1 tempDim2 dataDarkImage;
    clear dataNucImage dataDarkImageStack dataNucImageStack tempMovieDataSubtracted;
    clear tempMovieDataSubtractedNuc tempWriteName tempWriteFullName; 
end


%% Close variables
clear all;
close all;

disp('All files corrected and written.');


