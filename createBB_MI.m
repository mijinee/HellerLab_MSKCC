% Version for Mouse Imager (MI)
% Compensate for non uniform excitation by using a (1) BB image with laser off,
% and a (2) BB image with laser on.
% Output is an h5 file, run as script.
purge;
oldPath = cd;
%% Load a blank file with BB on and BB off.
initDirectory=uigetdir('raw');
cd(initDirectory);

[bbOn, bbOnPath, ~] = uigetfile('*.h5','Choose BB image with laser ON.');
[bbOff, bbOffPath, ~] = uigetfile('*.h5','Choose BB image with laser OFF.');

imageBBOn = h5read(strcat(bbOnPath,bbOn),'/Image/Data');
%imageBBOn = h5read(strcat(bbOnPath,bbOn),'/Cube/Images');
imageBBOff = h5read(strcat(bbOffPath,bbOff),'/Image/Data');
%imageBBOff = h5read(strcat(bbOffPath,bbOff),'/Cube/Images');

%% Gaussian smoothing?    
% H=fspecial('gaussian',[3 3],0.5); %Apply Gaussian lowpass filter to smooth data
H=fspecial('average',8); %Apply Gaussian lowpass filter to smooth data
imageSmooth=imfilter(imageBBOn,H,'replicate');
figure('name','Smoothed BB-On Image'); imshow(imageSmooth,[]); colorbar;

%% Processing

imageSubtract = imageSmooth - imageBBOff;
figure('name','BB-On minus BB-Off Image'); imshow(imageSubtract,[]); colorbar;
imageScaled = double(imageSubtract)./(max(max(double(imageSubtract))));
figure('name','After Normalization'); imshow(imageScaled,[]); colorbar;

%% Write as hdf
%nameImage = 'NonUniform 32bit.h5';
%copyfile(bbOff,nameImage);
writeName = fullfile(bbOnPath,'Non Uniform 32bit for BB.h5');
h5create(writeName,'/Image/Data',[256 320],'Datatype','double');
h5write(writeName,'/Image/Data',imageScaled);
disp('Correction file for Non-Uniform BB Excitation written.');
cd(oldPath);



