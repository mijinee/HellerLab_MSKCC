% Core function for analyzing spectra acquired from microscope 3
% 2019-11-21 | Prakrit Jena | prakrit.jena@gmail.com

function content = m3ple()
% function m3ple
% Required files:
%   Data: separate [wavelength intensity] files
%   Blank: Optional - 1 or N [wavelength intensity] files
%   Correction Factors: MAT file, double (512 x 2) - [W I]
%
% Function dependencies:
%   proInput
%   proLoad
%   proCombine
%   proApply
%   anaPlot
%   anaPeakRange
%   anaSubset
%   anaFitSpectra
%   anaFitPar
%       anaFitParameters
%       fn_LS_voigtfun
%   anaOutput
%   plotSmooth


%% Code Block 1: Processing - Obtain location of data, blank and calibration files

% Clearing MATLAB environment

clearvars; close all; clc

% proInput: Function to obtain file locations for the data, blank and calibration
[filesData,filesBlank,filesCalibration] = proInput();

% proLoad: Function to load the data locations obtained by proInput
[contentData,contentBlank,contentCalibrationFull] = proLoad(filesData,filesBlank,filesCalibration);

% proCombine: Function to combine the separate data and blank files
[contentDataFull,contentBlankFull] = proCombine(contentData,contentBlank);

% proApply: Function to apply the pixel-to-wavelength calibration and
% detector efficiency correction factors

content = proApply(contentDataFull,contentBlankFull,contentCalibrationFull);

% User Input: Name file to be saved
content.filesBlank = filesBlank;
content.filesCalibration = filesCalibration;
content.filesData = filesData;

disp('Data processing complete');

% Code to save content
prompt = 'Enter a name for the combined analysis mat file:';
name = input(prompt,'s');
content.saveName = name;
saveName = [name '.mat'];

cd(content.filesData.path);
save(saveName,'content');

disp('Data processing complete and MAT-file saved');

% Clean up
clearvars -except content; clc;
%% Code Block 2: Plotting - User decision on whether to make plots or not
answer = questdlg('Make plots of the processed data?','Initial Analysis','Yes','No','Yes');
% answer = char, 'Yes' or 'No'

switch answer
    case 'Yes'
        disp('Plotting 2 Blanks and 6 Data graphs');
        
        cd(content.filesData.path);
        
        anaPlot(content.BlankRaw,'Blank Raw');
        anaPlot(content.BlankCalibrated,'Blank Calibrated');
        
        anaPlot(content.SpectraRaw,'Raw');
        anaPlot(content.SpectraCalibrated,'Calibrated');
        anaPlot(content.SpectraBlanked,'Blanked');
        anaPlot(content.SpectraSmoothed,'Smoothed');
        anaPlot(content.SpectraNormalized,'Normalized');
        anaPlot(content.SpectraPositive,'Positive');
        
        disp('The data is stored in "content" as SpectraRaw...SpectraBlanked... etc');
    case 'No'
        disp('The data is stored in "content" as SpectraRaw...SpectraBlanked... etc');
end

clearvars -except content;

%% Code Block 3: Fitting - User decision on whether to fit data or not

answer = questdlg('Would you like to fit user-selected peaks with Lorentzian functions?' ...
    ,'Peak-Fitting Procedure','Yes','No','Yes');

switch answer
    case 'Yes'
        disp('Select wavelength ranges containing a single peak');
        
        % anaPeakRange plots spectra (intensity of each file vs wavelength), for
        % the user to click and select wavelength ranges containing the peak to be
        % fit. Each pair of clicks selectes one wavelength range. When all
        % peak-fitting ranges have been selected, press 'Enter'
        content.fitUserSelectedRange = anaPeakRange(content.SpectraBlanked);
        clearvars -except content;
        
        % anaSubset extracts spectra corresponding to the user-selected wavelength
        % range, and formats the data for curve-fitting
        content.fitUserSelectedData = anaSubset(content);
        
        % Prepare space for fitting results using anaFitSpectra
        
        fitUserSelectedDataFitParameters = cell(size(content.fitUserSelectedData));
        fitUserSelectedDataFitParameters(:,1) = content.fitUserSelectedData(:,1);
        fitUserSelectedDataFitParameters(1,:) = content.fitUserSelectedData(1,:);
        fitUserSelectedDataFitCurves = fitUserSelectedDataFitParameters;
        
        % Fit all data via loop using anaFitSpectra
        nRanges = size(content.fitUserSelectedData,1)-1;
        
        for ii = 1:nRanges
            for j = 2:7
                tempData = content.fitUserSelectedData{ii+1,j};
                tempOutputResults = anaFitSpectra(tempData);
                % Replace index values in Row 1 with file names
                tempFileNames = tempData(1,2:end);
                tempOutputResults.Fits(2:end,1) = tempFileNames';
                fitUserSelectedDataFitParameters{ii+1,j} = tempOutputResults.Fits;
                fitUserSelectedDataFitCurves{ii+1,j} = tempOutputResults.Smooth;
            end
        end
        clear j ii;
        % Saving the fit results
        content.fitUserSelectedDataFitParameters = fitUserSelectedDataFitParameters;
        content.fitUserSelectedDataFitCurves = fitUserSelectedDataFitCurves;
        clearvars -except content;
        disp('Fitting complete.');
        
        % Save combined fit results as well
        dataFit = content.fitUserSelectedDataFitParameters;
        content.fitUserSelectedDataFitParametersCombined = anaOutput(dataFit);
        tempSaveName = content.saveName;
        tempSaveNameAppend = [tempSaveName ' with Fits'];
        content.saveName = tempSaveNameAppend;
        save(tempSaveNameAppend,'content')
        
        % Write individual CSV files for the fit results
        dataFits = content.fitUserSelectedDataFitParametersCombined;
        
        % Loop for Raw,Calibrated,Blanked,Smoothed,Normalized,Positive
        for j = 2:7
            fileNameType = dataFits{1,j};
            fileNameFull = ['Fit Results for ' fileNameType '.csv'];
            fileContent = dataFits{2,j};
            writetable(cell2table(fileContent), fileNameFull, ...
                'writevariablenames', false, 'quotestrings', true);
        end
        
        % Plotting - User secision on whether to plot the fit results or not
        answer = questdlg('Make plots of the fitted data?',...
            'Data points + curve fit','Yes','No','Yes');
        % answer = char, 'Yes' or 'No'
        
        switch answer
            case 'Yes'
                disp('Plotting fits of all ranges and all data processing variations');
                nRanges = size(content.fitUserSelectedData,1)-1;
                
                for ii = 2:nRanges+1
                    for j = 2:7
                        tempData = content.fitUserSelectedData{ii,j};
                        tempName = content.fitUserSelectedData{1,j};
                        tempRange = content.fitUserSelectedData{ii,1};
                        tempSmooth = content.fitUserSelectedDataFitCurves{ii,j};
                        plotSmooth(tempData,tempName,tempRange,tempSmooth);
                    end
                end
                clear ii j;
                disp('The data is stored in "content" as fitUserSelectedDataFitParametersCombined');
            case 'No'
                disp('The data is stored in "content" as fitUserSelectedDataFitParametersCombined');
        end
        
        clearvars -except content;
        
    case 'No'
        tempFileName = content.saveName;
        disp(['Processed data has been saved in: ' tempFileName '.mat']);
        
end

%%

end         
% ----------------------End of primary function m3ple----------------------


% ---------------------Local functions called by m3ple---------------------

% proInput
function [filesData,filesBlank,filesCalibration] = proInput()
%   function [filesData,filesBlank,filesCalibration] = proInput()
%   m3plePRO --> proInput
%
%   Function called by m3ple
%   proInput is a function where the user inputs the
%   (i)   data files
%   (ii)  blanks (0, 1, oer N (1 per well))
%   (iii) calibration (Ocean Optics - wavelength-dependent detection)
%
%   proInput returns:
%   filesData = 1x1 struct with
%   filesData.names | filesData.path | with data location
%   filesBlank.names | filesBlank.path| with blank location
%   filesCalibration.names | filesCalibrationpath | with calibration
%% Code Block: Initialization

disp('Running: m3ple-->proInput');
initialPath = cd;

%% Code Block: Acquire data files (0, 1 or N)

[tempFileOpen,tempFilePath] = uigetfile({'*.txt';'*.csv';'*.xlsx'},...
    'Select the spectra files to process, excluding the blank:','MultiSelect','on');

% Code Block: Check if 0, 1 or N files were selected
if isa(tempFileOpen,'double')
    disp('No files selected. Exiting m3ple.');
    return;
elseif isa(tempFileOpen,'char')
    disp('1 file selected.')
    tempFileOpen = cellstr(tempFileOpen);
else
    tempFileOpenNumber = size(tempFileOpen,2);
    disp([num2str(tempFileOpenNumber) ' files selected.']);
end

% Code Block: When 1 or N files are selected, function continues
filesData.names = tempFileOpen;
filesData.path = tempFilePath;

clear tempFileOpen tempFilePath;
%% Code Block: Acquire blank file or files (0, 1 or N)

[tempFileOpen,tempFilePath] = uigetfile({'*.txt';'*.csv';'*.xlsx'},...
    'Select the blank files to use (None, 1 or 1 per well):','MultiSelect','on');

% Code Block: Check if 0, 1 or N files were selected
if isa(tempFileOpen,'double')
    disp('No blank file selected.');
elseif isa(tempFileOpen,'char')
    disp('1 file selected.')
    tempFileOpen = cellstr(tempFileOpen);
else
    tempFileOpenNumber = size(tempFileOpen,2);
    disp([num2str(tempFileOpenNumber) ' files selected.']);
end

% Code Block: When 1 or N files are selected, function continues
filesBlank.names = tempFileOpen;
filesBlank.path = tempFilePath;

clear tempFileOpen tempFilePath
%% Code Block: Acquire ocean-optics intensity calibration file (0 or 1)
cd('\\skimcs.mskcc.org\hellerlab\Software\Matlab\m3ple\')
[tempFileOpen,tempFilePath] = uigetfile({'*.mat'}, ...
    'Select the Ocean-Optics calibration file (None or 1):','MultiSelect','off');

if isa(tempFileOpen,'double')
    disp('No calibration file selected.');
elseif isa(tempFileOpen,'char')
    disp('Calibration file selected.')
    tempFileOpen = cellstr(tempFileOpen);
end

filesCalibration.names = tempFileOpen;
filesCalibration.path = tempFilePath;

clear tempFileOpen tempFilePath
%% Code Block: Return to initial folder
cd(initialPath);
end

% proLoad
function [contentData,contentBlank,contentCalibrationOutput] = proLoad(filesData,filesBlank,filesCalibration)
%   [contentData,contentBlank,contentCalibrationOutput] = proLoad(filesData,filesBlank,filesCalibration)
%   m3ple --> proLoad

%   Function called by m3mple
%   proLoad is a function where the incoming data, blank and
%   calibration file names and paths are loaded into the workspace
%   checked and compared for self-consistency, then applied
%
%   Inputs:
%   filesData = 1x1 struct with the fields:
%   filesData.names = 1 x 'number of files' cell
%   filesData.names{1,1} = 'A2.txt' type of file name
%   filesData.path = string of the file folder path
%   filesBlank and filesCalibration include the same structure and content
%
%
%   proLoad returns:
%   contentData = 2 x 'number of files' cell
%   contentData: row 1 = 'A2.txt' type file names
%   contentData: row 2 = '512 x 2' double containing [W I] for each file
%   contentBlank = 2 x 'number of blanks' cell
%   contentBlank: row 1 = 'Names of Blanks'
%   contentBlank: row 2 = '512 x 2' double containing [W I] for each blank
%   contentCalibration
%% Code Block: Initialization

disp('Running: m3ple-->proCorrections');
initialPath = cd;

%% Code Block: Load Data

tempPath = filesData.path;
cd(tempPath);
tempNumFiles = size(filesData.names,2);
tempDataCell = cell(2,tempNumFiles);

for ii = 1:tempNumFiles
    tempFileName = filesData.names{1,ii};
    tempFileData = dlmread(tempFileName);
    tempDataCell{1,ii} = tempFileName;
    tempDataCell{2,ii} = tempFileData;
end

contentData = tempDataCell;
clear tempPath tempFileName tempFileData tempNumFiles tempDataCell ii;
cd(initialPath);
%% Code Block: Load Blank

% Check if a blank was selected or not

TF = isequal(filesBlank.names,0); % If no blank, then the entry has value 0

if (TF == 1)
    % Generate a 'blank' file with zero data
    disp('No blank files were selected.');
    tempDataFirst = contentData{2,1};
    tempDataFirst(:,2) = 0;
    tempDataCell = cell(2,1);
    tempDataCell{1,1} = 'Zero Blank';
    tempDataCell{2,1} = tempDataFirst;
    contentBlank = tempDataCell;
    clear tempDataFirst;
else
    % Either 1 or N blank files selected
    tempPath = filesBlank.path;
    cd(tempPath);
    tempNumFiles = size(filesBlank.names,2);
    % disp([num2str(tempNumFiles) ' blank files was/were selected.']);
    tempDataCell = cell(2,tempNumFiles);
    
    for ii = 1:tempNumFiles
        tempFileName = filesBlank.names{1,ii};
        tempFileData = dlmread(tempFileName);
        tempDataCell{1,ii} = tempFileName;
        tempDataCell{2,ii} = tempFileData;
    end
    contentBlank = tempDataCell;
    
    clear tempPath tempFileName tempFileData tempNumFiles tempDataCell ii TF;
end

%% Code Block: Load Calibration Data
% _old_v1 allowed the analysis to run without a calibration file
% current version REQUIRES a calibration file of the form:
% 2019-11-07 M3 Well Plate Correction Factors.mat, with content
% correctionFactors = 512 x 2 double, [wavelength CF] values
% Check if a blank was selected or not

TF = isequal(filesCalibration.names,0); % If no blank, then the entry has value 0

if (TF == 1)
    % Generate a 'calibration' file with zero data
    disp('No calibration file was selected.');
    disp('A Calibration File is REQUIRED.');
    return;
    %     tempDataFirst = contentData{2,1};
    %     tempDataFirst(:,2) = 1;
    %     contentCalibration = tempDataFirst;
    %     clear tempDataFirst;
    %
else
    % Load calibration file
    tempPath = filesCalibration.path;
    cd(tempPath);
    load(filesCalibration.names{1,1}); % correctionFactors
    contentCalibration = correctionFactors;
end

contentCalibrationHeader = {'Wavelength','CF'};
contentCalibrationOutput = [contentCalibrationHeader; num2cell(contentCalibration)];

%%
clear TF tempPath correctionFactors;

cd(initialPath);

end

% proCombine
function [contentDataFull,contentBlankFull] = proCombine(contentData,contentBlank)
%   function [contentDataFull,contentBlankFull] = proCombine(contentData,contentBlank)
%   m3ple --> proCombine
%
%   Function called by m3ple
%   proCombine combines the separate data and blank into one cell

%   Inputs:
%   contentData = 2 x 'number of files' cell
%   contentData: row 1 = 'A2.txt' type file names
%   contentData: row 2 = '512 x 2' double containing [W I] for each file
%   contentBlank = 2 x 'number of blanks' cell
%   contentBlank: row 1 = 'Names of Blanks'
%   contentBlank: row 2 = '512 x 2' double containing [W I] for each blank
%
%   Output:
%
%   contentDataFull = 513 x 'number of files + 1' cell
%   Header Row = 'Wavelength', 'file names'
%   Column 1 = Wavelength values
%
%   contentBlankFull = 513 x 'number of blanks + 1' cell
%   Header Row = 'Wavelength', 'file names'
%   Column 1 = Wavelength values
%
%% Code Block:
disp('Running: m3ple-->proCombine');

%% Code Block: Initialize data, combine

tempNumFiles = size(contentData,2);

tempFirstData = contentData{2,1};
tempNumRows = size(tempFirstData,1);

tempWavelength = tempFirstData(:,1);
tempCombinedData = zeros(tempNumRows,tempNumFiles);

for ii = 1:tempNumFiles
    tempCombinedData(:,ii) = contentData{2,ii}(:,2);
end

tempFullSpectra = [tempWavelength tempCombinedData];
tempFullSpectra = flipud(tempFullSpectra);

tempDataHeader = cell(1,tempNumFiles+1);
tempDataHeader{1,1} = 'Wavelength';
tempDataHeader(1,2:tempNumFiles+1) = contentData(1,1:tempNumFiles);

contentDataFull = [tempDataHeader; num2cell(tempFullSpectra)];
clear tempNumFiles tempFirstData tempNumRows tempWavelength; 
clear tempCombinedData ii tempFullSpectra tempDataHeader;
%% Code Block: Initialize blanks, combine

tempNumFiles = size(contentData,2);
tempFirstData = contentData{2,1};
tempWavelength = tempFirstData(:,1);
tempNumRows = size(tempFirstData,1);
tempNumBlanks = size(contentBlank,2);
tempCombinedBlank = zeros(tempNumRows,tempNumFiles);
%%
if (tempNumBlanks == 1)
    disp('One blank being applied to all files.')
    tempBlankData = contentBlank{2,1}(:,2);
    tempCombinedBlank = repmat(tempBlankData,1,tempNumFiles);
    % TODO: Code for repeating single file name many times
    tempBlankHeader = cell(1,tempNumFiles+1);
    tempBlankHeader{1,1} = 'Wavelength';
    tempBlankHeader(1,2:tempNumFiles+1) = contentBlank(1,1);
    tempFullBlank = [tempWavelength tempCombinedBlank];
    tempFullBlank = flipud(tempFullBlank);
else
    for ii = 1:tempNumFiles
        tempCombinedBlank(:,ii) = contentBlank{2,ii}(:,2);
    end
    
    tempFullBlank = [tempWavelength tempCombinedBlank];
    tempFullBlank = flipud(tempFullBlank);
    tempBlankHeader = cell(1,tempNumFiles+1);
    tempBlankHeader{1,1} = 'Wavelength';
    tempBlankHeader(1,2:tempNumFiles+1) = contentBlank(1,1:tempNumFiles);
end

contentBlankFull = [tempBlankHeader; num2cell(tempFullBlank)];

clear tempNumFiles tempFirstData tempWavelength tempNumRows;
clear tempCombinedBlank tempBlankData ii tempFullBlank tempBlankHeader;

end

% proApply
function [content] = proApply(contentDataFull,contentBlankFull,contentCalibrationFull)
%   function [content] = proApply(contentDataFull,contentBlankFull,contentCalibrationFull)
%   m3ple --> proApply
%
%   Function called by m3ple
%   proApply applies the calibration file to populate the calibrated
%   wavelength and corrected intensity into the data and blanks
%
%   Inputs:
%
%   contentDataFull = 513 x ( 'number of files' + 1) cell
%       Header row = 'Wavelength', 'File Names'
%       Column 1 = Wavelength values
%       Column 2 - N = Intensity values
%
%   contentBlankFull = 513 x 1 + (1 OR  'number of files') cell
%                      depending on whether 0, 1 or N blanks were provided
%       Header row = 'Wavelength', 'File Names'
%       Column 1 = Wavelength values
%       Column 2 - N = Intensity values
%
%   contentCalibrationFull = 513 x 2 cell
%       Header row = 'Wavelength', ' CF'
%       Column 1 = calibrated wavelength (512 values)
%       Column 2 = intensity correction factors (512 values)
%
%   Output:
%
%   contentDataFull = 513 x 'number of files + 1' cell
%   Header Row = 'Wavelength', 'file names'
%   Column 1 = Wavelength values
%
%   contentBlankFull = 513 x 'number of blanks + 1' cell
%   Header Row = 'Wavelength', 'file names'
%   Column 1 = Wavelength values
%% Code Block: Generate subset of data matching

%   Calibrating the wavelength will change the [pixel <--> wavelength] map
%   The wavelength used will be 320 pixels long, starting from the index value
%   closest to 900 nm and going to 320 indices down from that

% Matching wavelength range ~ 900 to 1400 nm
disp('Running: m3ple-->proApply');
tempCalibrationWavelengthFull = cell2mat(contentCalibrationFull(2:end,1));
[~,wavelengthStartPixel] = min(abs(900-tempCalibrationWavelengthFull));

% Above gives wavelengthStartPixel, the row with a value closest to 900
wavelengthEndPixel = wavelengthStartPixel + 319;



%% Code Block: Separate wavelength and intensity
contentDataHeader = contentDataFull(1,:);
contentBlankHeader = contentBlankFull(1,:);

tempContentData = cell2mat(contentDataFull(2:end,:));
tempContentSpectra = tempContentData(:,2:end);

tempBlankData = cell2mat(contentBlankFull(2:end,:));
tempBlankSpectra = tempBlankData(:,2:end);

contentCalibration = cell2mat(contentCalibrationFull(2:end,:));
tempCalibratedWavelength = contentCalibration(:,1);
tempCorrectionFactors = contentCalibration(:,2);

% clear tempContentData tempBlankData;
%% Code Block: Reduce data into the right wavelength range
% Wavelength in the Correction Factor file is replacing whatever wavelength
% the individual txt files had

tempCalibratedWavelength = tempCalibratedWavelength(wavelengthStartPixel:wavelengthEndPixel,1);
tempCorrectionFactors = tempCorrectionFactors(wavelengthStartPixel:wavelengthEndPixel,:);

tempContentSpectra = tempContentSpectra(wavelengthStartPixel:wavelengthEndPixel,:);
tempBlankSpectra = tempBlankSpectra(wavelengthStartPixel:wavelengthEndPixel,:);

tempContentWavelength = tempCalibratedWavelength;
% When no correction factor is chosen, CF length = date length = 512x1
% When calibration file is selected, CF length = 267
% %% Code Block: Choose the correct length of the right wavelength range
% tempLengthCF = size(tempCorrectionFactors,1);
% if tempLengthCF == 512
%    tempCorrectionFactors = tempCorrectionFactors(wavelengthStartPixel:wavelengthEndPixel,1);
% end
%
% clear wavelengthEndPixel wavelengthStartPixel;
%% Code Block: Process the spectra in different ways

% contentDataRaw = input data, in the selected wavelength range
% contentBlankRaw = input blank, in the selected wavelength range
contentDataRaw = [contentDataHeader; num2cell([tempContentWavelength tempContentSpectra])];
contentBlankRaw = [contentBlankHeader; num2cell([tempContentWavelength tempBlankSpectra])];

% contentDataCalibrated = raw data calibrated with correction factors
% contentBlankCalibrated = raw blank calibrated with correction factors



tempContentSpectraCalibrated = bsxfun(@times,tempContentSpectra,tempCorrectionFactors);
tempBlankSpectraCalibrated = bsxfun(@times,tempBlankSpectra,tempCorrectionFactors);

contentDataCalibrated = [contentDataHeader; num2cell([tempContentWavelength tempContentSpectraCalibrated])];
contentBlankCalibrated = [contentBlankHeader; num2cell([tempContentWavelength tempBlankSpectraCalibrated])];

% contentDataBlanked = calibrated data with blanks subtracted

tempContentSpectraBlanked = bsxfun(@minus,tempContentSpectraCalibrated,tempBlankSpectraCalibrated);
contentDataBlanked = [contentDataHeader; num2cell([tempContentWavelength tempContentSpectraBlanked])];

% contentDataSmoothed = blanked data smoothed with 'sgolay', default window
% Needs to be done in a loop

tempContentSpectraSmoothed = tempContentSpectraBlanked;
for ii = 1:1:size(tempContentSpectraBlanked,2)
    tempContentSpectraSmoothed(:,ii) = smooth(tempContentWavelength,tempContentSpectraBlanked(:,ii),'sgolay');
end
contentDataSmoothed = [contentDataHeader; num2cell([tempContentWavelength tempContentSpectraSmoothed])];
clear ii;

% contentDataNormalized = blanked data normalized

tempContentSpectraNormalized = normalize(tempContentSpectraBlanked,1,'range');
contentDataNormalized = [contentDataHeader; num2cell([tempContentWavelength tempContentSpectraNormalized])];

% contentDataPositive = blanked data with a global offset (add one minimum
% value to all columns to make the data positive)

tempMinValue = -1 * min(tempContentSpectraBlanked(:));
tempContentSpectraPositive = bsxfun(@plus,tempContentSpectraBlanked,tempMinValue);

contentDataPositive = [contentDataHeader; num2cell([tempContentWavelength tempContentSpectraPositive])];

%% Code Block: Combine data into a structure and return

% Input from user
content.InputData = contentDataFull;
content.InputBlank = contentBlankFull;
content.InputCalibration = contentCalibrationFull;

% Processed
content.BlankRaw = contentBlankRaw;
content.BlankCalibrated = contentBlankCalibrated;

content.SpectraRaw = contentDataRaw;
content.SpectraCalibrated = contentDataCalibrated;
content.SpectraBlanked = contentDataBlanked;
content.SpectraSmoothed = contentDataSmoothed;
content.SpectraNormalized = contentDataNormalized;
content.SpectraPositive = contentDataPositive;



end

% anaPlot
function anaPlot(plotSpectra,plotName)
%   function anaPlot(plotSpectra,plotName)
%   m3ple --> anaPlot
%
%   Function called by m3ple
%   anaPlot plots the data in plotSpectra and labels the plot and file name
%   with plotName
%
%   Input:
%
%   plotName = cell, character identifying the data 'Calibrated'
%
%   plotSpectra = cell, 321 x (number of files + 1)
%   Header row = 'Wavelength', 'File Names'
%   Column 1 = wavelength
%   Column 2...N = Intensity

%% Code Block: Extract data from plotSpectra

spectraNames = plotSpectra(1,2:end);

wavelengthData = cell2mat(plotSpectra(2:end,1));
spectraData = cell2mat(plotSpectra(2:end,2:end));

%% Code Block: Plotting the data above
f1 = figure;
plot(wavelengthData,spectraData);

xlabel('Wavelength (nm)','fontsize',12,'fontweight','b','color','k');
ylabel('Intensity (a.u.)','fontsize',12,'fontweight','b','color','k');

%% Code Block: Axis labels, round to the nearest multiple of 100

xminReal = wavelengthData(1,1);
xmaxReal = wavelengthData(end,1);

xmin = round(xminReal,-2);
xmax = round(xmaxReal,-2);

xlim([xmin xmax]);

%% Code Block: Labels

title(plotName,'fontsize',12,'fontweight','b','color','k','Interpreter','none');
legend(spectraNames,'FontSize',7,'Interpreter','none');

%% Save plot with a specific name
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 8];
saveName = ['Plot - Spectra - ' plotName];
print(saveName,'-dpng','-r0');
%%
saveas(f1,saveName,'png');
close(f1);
end

% anaPeakRange
function selectedRangeFits = anaPeakRange(contentSpectraBlanked)
% function selectedRangeFits = anaPeakRange(contentSpectraBlanked)
% m3ple --> anaPeakRange
%
% Function called by m3ple
% anaPeakRange plots spectra (intensity of each file vs wavelength), for
% the user to click and select wavelength ranges containing the peak to be
% fit. Each pair of clicks selectes one wavelength range. When all
% peak-fitting ranges have been selected, press 'Enter'
%
% Input:
%
% contentSpectraBlanked = cell, 321 x (number of files + 1)
% Row 1 = header
% Column 1 = wavelength
% Column 2...N = intensity
%
% Output:
%
% selectedRangeFits = N x 3 double, N = number of ranges
% Each row: range index, left (starting) index, right (ending) index

%% Code Block: Prepare incoming data

tempWavelength = cell2mat(contentSpectraBlanked(2:end,1));
tempSpectraBlanked = cell2mat(contentSpectraBlanked(2:end,2:end));

%% Code Block: Initial graph for the user to select a wavelength range from

f1 = figure('Name','Left-click to select points defining ranges, Enter to finish',...
    'NumberTitle','off');
plot(tempWavelength,tempSpectraBlanked);
xlabel('Wavelength (nm)','fontsize',12,'fontweight','b','color','k');
ylabel('Intensity (a.u.)','fontsize',12,'fontweight','b','color','k');

% xmin = 900;
% xmax = 1400;

xmin = round(tempWavelength(1,1),-2);
xmax = round(tempWavelength(end,1),-2);
xlim([xmin xmax]);

title('Blanked Calibrated Spectra for Peak Finding','fontsize',12,'fontweight',...
    'b','color','k','Interpreter','none');

clear xlabel ylabel xmin xmax xlim;
%% Code Block: User selection of wavelength ranges

[x,~] = ginput();
numClicks = size(x,1);

%% Code Block: Find wavelength-index closest to user clicks
tempIndex = zeros(numClicks,1);

for ii = 1:numClicks
    [~, tempIndex(ii,1)] = min(abs(tempWavelength - x(ii,1)));
end
clear ii;
%% Code Block: Plot figure with clicked ranges

tempXselection = tempWavelength(tempIndex,1);
tempYselection = tempSpectraBlanked(tempIndex,1);

hold on;
scatter(tempXselection,tempYselection,50,'k','MarkerFaceColor','r');
pngName = 'Plot - Blanked Calibrated Spectra for Peak Finding';
saveas(f1,pngName,'png');

%% Code Block: Transform clicked ranges into start and end

tempClickedData = zeros(numClicks,1);
tempClickedData(:,1) = tempIndex;

tempStartClicks = tempClickedData(1:2:numClicks,1);
tempEndClicks = tempClickedData(2:2:numClicks,1);

numRanges = size(tempStartClicks,1);
saveClicks = [(1:1:numRanges)' tempStartClicks(:,1) tempEndClicks(:,1)];


%% Code Block: Return user selected range
selectedRangeFits = saveClicks;
close(f1);
end

% anaSubset
function outputFittingRange = anaSubset(content)
% function outputFittingRange = anaSubset(content)
% m3ple --> anaSubset
%
% Function called by m3ple
% anaSubset extracts spectra corresponding to the user-selected wavelength
% range, and formats the data for curve-fitting
%
% Input:
%   content is the full and only workspace variable
%
% Output:
%   outputFittingRange = (Number of Ranges + 1) x 7 cell
%   Header Row = Range, Raw, Calibrated, Blanked, Smoothed, Normalized,
%   Positive
%   Column 1 = Wavelength Range Index
%   Cell contents are [Header; ]Wavelength I1 I2...IN]

%% Code Block: Prepare fitting range data holders

% numFiles = size(content.SpectraRaw,2)-1;
selectedRangeFits = content.fitUserSelectedRange;
numRanges = size(selectedRangeFits,1);

tempCellFittingRanges = cell(numRanges+1,7);

tempCellFittingRanges(1,:) = {'Range','Raw','Calibrated','Blanked',...
    'Smoothed','Normalized','Positive'};

%% Code Block: Run loop for number of ranges and fill in data

for ii = 1:numRanges
    tempCellFittingRanges{ii+1,1} = ii;
    tempCellFittingRanges{ii+1,2} = extractRange(content.SpectraRaw,selectedRangeFits(ii,:));
    tempCellFittingRanges{ii+1,3} = extractRange(content.SpectraCalibrated,selectedRangeFits(ii,:));
    tempCellFittingRanges{ii+1,4} = extractRange(content.SpectraBlanked,selectedRangeFits(ii,:));
    tempCellFittingRanges{ii+1,5} = extractRange(content.SpectraSmoothed,selectedRangeFits(ii,:));
    tempCellFittingRanges{ii+1,6} = extractRange(content.SpectraNormalized,selectedRangeFits(ii,:));
    tempCellFittingRanges{ii+1,7} = extractRange(content.SpectraPositive,selectedRangeFits(ii,:));
end
clear ii;

%% Code Block: outputFittingRange
outputFittingRange = tempCellFittingRanges;

    function extractedRange = extractRange(spectralData, spectralRange)
        % function extractedRange = extractRange(spectralData, spectralRange)
        % m3ple --> anaSubset --> extractedRange
        
        spectralDataHeader = spectralData(1,:);
        spectralDataContentFull = cell2mat(spectralData(2:end,:));
        spectralDataContentSubset = spectralDataContentFull(spectralRange(1,2):spectralRange(1,3),:);
        extractedRange = [spectralDataHeader; num2cell(spectralDataContentSubset)];
    end

end

% anaFitSpectra
function outputResults = anaFitSpectra(tempData)
%   function outputResults = anaFitSpectra(tempData)
%   m3ple --> anaFitSpectra
%
%   Function called by m3ple
%   anaFitSpectra fits the data in tempData, and outputs the results of the
%   fitting (parameters + fitted curve)
%
%   Input:
%
%   tempData = cell, [W I1 Iw] and a header row
%
%   Subfunctions:

%   anaFitPar
%       anaFitParameters
%       fn_LS_voigtfun
%
%% Code Block: Initialize variables
tempSpectra = cell2mat(tempData(2:end,2:end));
tempWavelength = cell2mat(tempData(2:end,1));
tempDataHeader = tempData(1,:);

%% Code Block: Using first data to get size of the Smooth Output
numSpectra = size(tempSpectra,2);
tempDataSizeSmooth = [tempWavelength tempSpectra(:,1)];
[~,x_nm_smooth,~] = anaFitPar(tempDataSizeSmooth);
numSmoothX = size(x_nm_smooth,1);

%% Code Block: Prepare output size
outputSmooth = zeros(numSmoothX,numSpectra);
x_smooth_shared = x_nm_smooth;

outputy0 = zeros(size(numSpectra,1));
outputXc = zeros(size(numSpectra,1));
outputArea = zeros(size(numSpectra,1));
outputPeak = zeros(size(numSpectra,1));
outputFWHM = zeros(size(numSpectra,1));
outputR2 = zeros(size(numSpectra,1));
outputIndex = zeros(size(numSpectra,1));

clear ii tempData x_nm_smooth;

%% Code Block: Fitting
parfor ii = 1:numSpectra
    tempData = [tempWavelength tempSpectra(:,ii)];
    [output,~,F_smooth] = anaFitPar(tempData);
    
    outputIndex(ii,1) = ii;
    outputy0(ii,1) = output(1);
    outputXc(ii,1) = output(2);
    outputArea(ii,1) = output(3);
    outputPeak(ii,1) = output(4);
    outputFWHM(ii,1) = output(5);
    outputR2(ii,1) = output(6);
    
    outputSmooth(:,ii) = F_smooth;
    
end

outputFits = [outputIndex outputy0 outputXc outputArea outputPeak outputFWHM outputR2];
outputWavelengthSmooth = [x_smooth_shared outputSmooth];
clear ii tempData outputIndex outputy0 outputXc outputArea outputFWHM outputR2 F_smooth;


%% Code Block: Headers to cell
outputFitHeader = {'File','Y-Offset','Wavelength','Area','Peak','FWHM','R2'};
outputFitCell = [outputFitHeader; num2cell(outputFits)];

outputSmoothCell = [tempDataHeader; num2cell(outputWavelengthSmooth)];
outputResults.Fits = outputFitCell;
outputResults.Smooth = outputSmoothCell;
end

% anaFitPar
function [output,x_nm_smooth,F_smooth] = anaFitPar(tempData)
% function [output,x_nm_smooth,F_smooth] = anaFitPar(tempData)
% anaFitSpectra --> anaFitPar
% anaFitPar --> anaFitParameters
%           --> fn_LS_voigtfun
%
% Function anaFitPar fits one peak with the Lorentzian function
%
% Input:
% tempData = (wavelength,intensity) double
%
% Output:
% output = double, 1 x 6 = [y_final; x_final_nm; area_final; y_peak; gamma_final; R ];
% x_nm_smooth = wavelength values for fitted curve
% F_smooth = intensity values for fitted curve

% Code Block:
x_nm = tempData(:,1); % simply extracting the wavelength column
intensity = tempData(:,2); %column of nRows x 1, intensity
x_eV = 1239842 ./x_nm; % Converting wavelength (nm) into intensity (eV)

% Initial Fit Parameters
[y0,x0_nm,area_0,gamma_0,y_peak] = anaFitParameters(x_nm,intensity);

x0_eV = 1239842/x0_nm; % peak position in eV
x_initial = [y0; x0_eV; area_0; gamma_0]; % set of guess parameters

clear y0 area_0 gamma_0 x0_eV;
%% How to decide if this data set will be fit or not i.e. it'll crash.
F = fn_LS_voigtfun(x_initial,x_eV); % This is either finite or not.
TF = isnan(F);
testTF = sum(TF);

if testTF > 0
    output = [0;0;0;0;0];
    x_nm_smooth = x_nm(1,1):0.5:x_nm(end,1); % 1x161
    sizeSmooth = size(x_nm_smooth,2);
    x_nm_smooth = zeros(sizeSmooth,1);
    F_smooth = x_nm_smooth;
    
else
    
    %x_ub = x0_nm + 8; x_lb = x0_nm - 8;
    %x_ub_eV = 1239842/x_ub; x_lb_eV = 1239842/x_lb;
    % Range in y0 is from -infinity to half the y-value in the data
    %y_lb = -inf; y_ub = max(intensity)/2;
    %gamma_lb = 10;
    %gamma_ub = x_eV(1)-x_eV(end);
    % Define bounds - solution = [y_final x_final_eV area_final gamma_final
    %lb = [y_lb x_ub_eV 0 gamma_lb];   %upper x_nm is lower x_eV
    %ub = [y_ub x_lb_eV inf gamma_ub];
    lb = [];
    ub = [];
    %options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','Display','final-detailed');
    options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','Display','off');
    
    [solution,resnorm] = lsqcurvefit(@fn_LS_voigtfun,x_initial,x_eV,intensity,lb,ub,options);

    
    %% Converting solution into usable format;
    area_final = solution(3); % area of fit
    gamma_final = solution(4); % width of fit
    x_final_eV = solution(2); % peak in eV
    x_final_nm = 1239842 /x_final_eV; % peak in nm
    y_final = solution(1); % y_offset
    ss_tot=sum((intensity(:,1)-mean(intensity(:,1))).^2);
    R = 1-resnorm/ss_tot; %R square calculation
    
    %% Code Block:
    output = [y_final; x_final_nm; area_final; y_peak; gamma_final; R ];
    % Result from fit
    %%
    x_nm_smooth = x_nm(1,1):0.5:x_nm(end,1);
    x_eV_smooth = (1239842./x_nm_smooth)';
    
    F_smooth = fn_LS_voigtfun(solution,x_eV_smooth); %For smooth voigt fit
    
    x_nm_smooth = x_nm_smooth(:); 
    F_smooth = F_smooth(:);      
    
end

    function [y0,xC_nm,area,gamma,y_peak] = anaFitParameters(tempx_nm,tempIntensity)
        % function [y0,xC_nm,area,gamma,y_max] = anaFitParameters(tempx_nm,tempIntensity)
        % Input wavelength(nm) and intensity data, receive fit parameters as
        %y0, xC_nm,area,gamma
        
        %% Data to extract the parameters from
        x_nm_afp = tempx_nm;              % Wavelength (nm)
        intensity_afp = tempIntensity;    % Intensity
        
        %% Parameter extraction
        
        y0 = min(intensity_afp);
        [ymax,index_peak] = max(intensity_afp(:,1));
        xC_nm = x_nm_afp(index_peak);
        
        yhalf = (y0 + ymax)/2;
        [~,r_ind] = min(abs(intensity_afp(index_peak:end,1)-yhalf));
        %r_index = index_peak+r_ind-1; only searching on right side of peak
        [~,l_ind] = min(abs(intensity_afp(1:index_peak,1)-yhalf));
        %l_index = l_ind; only searching on left side of peak
        gamma = x_nm_afp(index_peak+r_ind-1,1)-x_nm_afp(l_ind,1);
        
        Intensity_baseline = bsxfun(@minus,intensity_afp,y0);
        % area = trapz(x_nm,Intensity_baseline);
        y_peak = ymax;
        %%
        try
            area = trapz(x_nm_afp,Intensity_baseline);
        catch
            area = 0;
        end
        %%
        % gamma = 20;
        
    end

    function F = fn_LS_voigtfun(x,xdata_eV)
        % function F = fn_LS_voigtfun(x,xdata_eV)
        % anaFitPar --> fn_LS_voigtfun
        % Functional form of Lorentzian fit
        % temp = (area*gamma./((E-E_0).^2+gamma.^2/4)/(2*pi))+Min;
        % x(1) = y0;
        % x(2) = x0_eV
        % x(3) = area_0;
        % x(4) = gamma_0;
        % x is the initial parameters, xdata_eV is the wavelength in energy
        

        F = (x(3)*x(4)./((xdata_eV-x(2)).^2+x(4).^2/4)/(2*pi))+x(1);
        %
        
        % Looks like this functional form:
        % http://mathworld.wolfram.com/LorentzianFunction.html
        
        % F is the function mapping xdata_eV with parameters x to intensity

    end

end

% anaOutput
function dataOutputCombined = anaOutput(dataFit)
% function dataOutputCombined = anaOutput(dataFit)
% m3ple --> anaOutput
%
% Function called by m3ple
% anaOutput takes the fit results for each wavelength range and combines
% them into a single cell

%% Code Block: Just run for Blanked, for now

nRanges = size(dataFit,1)-1;

for j = 2:7
    outputCellBlank = cell(1,7);
    outputCell = cell(1,7);
    for ii = 2:nRanges+1
        tempCellData = dataFit{ii,j};
        outputCell = [outputCell; tempCellData; outputCellBlank];
    end
    dataFit{nRanges+2,j} = outputCell;
    clear outputCell;
end

%% Code Block: Prepare data to return:

dataOutputHeader = dataFit(1,:);
dataOutputHeader{1,1} = 'Data - ';
dataOutputContent = dataFit(end,:);

dataOutputCombined = [dataOutputHeader; dataOutputContent];

end

% plotSmooth
function plotSmooth(tempData,tempName,tempRange,tempSmooth)
% function plotSmooth(tempData,tempName,tempRange,tempSmooth)

% tempData = cell, [W I1 I2] with header row
% tempName = 'char' string, eg 'Blanked'
% tempRange = double, eg '1', corresponding to wavelength range
%

%% Code Block: Format data for plotting

dataLegend = tempData(1,2:end);
dataSpectra = cell2mat(tempData(2:end,2:end));
dataWavelength = cell2mat(tempData(2:end,1));

dataSpectraSmooth = cell2mat(tempSmooth(2:end,2:end));
dataWavelengthSmooth = cell2mat(tempSmooth(2:end,1));

%% Need to figure out color scale here
nSpectra = size(dataSpectra,2);
c = colormap(parula(nSpectra));

%% Code Block: Line and scatter plot

f1 = figure;
hold on;
for ii = 1:nSpectra
    scatter(dataWavelength,dataSpectra(:,ii),30,c(ii,:),'filled');
end
legend(dataLegend,'Location','northeastoutside','AutoUpdate','off','Interpreter','none');
plot(dataWavelengthSmooth,dataSpectraSmooth,'LineWidth',1,'Color','k','LineStyle','--');

%% Code Block: Axis, title etc addition
xlabel('Wavelength (nm)');
ylabel('Intensity (a.u.)');

titleText = ['Data with Fit of ' tempName ' for Range ' num2str(tempRange)];
title(titleText);

xmin = min(dataWavelength);
xmax = max(dataWavelength);
xlim([xmin xmax]);

%% Code to Save:
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 8];
print(['Plot - ' titleText],'-dpng','-r300');
close(f1);
close all;
end
