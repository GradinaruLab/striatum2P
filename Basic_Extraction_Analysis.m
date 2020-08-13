%% Basic Feature Extraction --------------------------------------------------------

% Variables calculated here include:
% peakLocsperMin which gives event rate per minute for each neuron 
% peakMags which gives average magnitude for each neuron

%% Obtain Files

% Run this to obtain a .mat file containing cell activity (using example
% .mat files provided in Example_Cell_Traces

[fileName, filePath] = uigetfile('*.mat');
addpath(filePath)
load(fileName)    

%% Correction for neuropil

% This removes neuropil signal according to the criteria described in
% Suite2P (Pachitariu et al., 2017, bioRxiv)

FCor1 = F - 0.7*Fneu;

numberofcell = size(FCor1);

%% Selection of cells classified using suite2p

% Suite2P allows users to inspect traces obtained from source extraction
% and remove traces that contain noise above a given criterion, the code
% below obtains manually inspected cells.

FCor2 = single.empty;
for n = 1:numberofcell(1,1)

    
    FCorW = FCor1(n,:)*iscell(n,1);
    if sum(FCorW) > 0
        FCor2  = [FCor2;FCorW];
    end
        
    
end
  
%% Obtaining df_f for each cell

% We obtain df_f using the definition provided in Badimon et al., 2020

numberofcellv2 = size(FCor2);
FCor3 = single(zeros(numberofcellv2));

for n = 1:numberofcellv2(1,2)
    
    FCor3(:,n) = (FCor2(:,n) - median(FCor2(:,:),2))./median(FCor2(:,:),2);
        
end

%% Applying median filter

% Traces are denoised using a median filter

FCor4 = single(zeros(numberofcellv2));
for n = 1:numberofcellv2(1,1)

    
    FCor4(n,:) = medfilt1(double(FCor3(n,:)),10);
            
    
end

%% Finding peaks and calculating events per min

% The following code caculates the main variables reported in Badimon et
% al., 2020, including the average magnitude and event rate of calcium
% transients

% Location of peak in trace
peakLocs = single.empty;

% Magnitude of peak (df/f)
peakMags = single.empty;

for n = 1:numberofcellv2(1,1)
    
    [peakLoc, peakMag] = peakfinder(double(FCor4(n,:)));
    
    peakLocX = single.empty;
    for i = 1:length(peakLoc)
    
    if peakMag(i) >= 2*std(FCor4(n,:))
       peakLocX = [peakLocX,peakMag(i)];
    end
    
    end
    
    peakLocmean = length(peakLocX);
    peakMagmean = mean(peakLocX);
    peakLocs = [peakLocs,peakLocmean];
    peakMags = [peakMags,peakMagmean];
        
end


% Peaks per min 

% Example traces provided are imaged at 4 Hz and last 10 mins
% We obtain Event Rate per minute below:

peakLocsperMin = peakLocs/(length(F)/60/4);

%% Results:

disp(['Mean Event Rate per minute: ' num2str(mean(peakLocsperMin)) ])
disp(['Mean Magnitude (df/f) per minute: ' num2str(mean(peakMags))])
