%% Batch Process: Spatial Coorelation Plots -------------------------------------

%% Output Format

% The output of the script Cell_Array contains cells with 4 columns. The
% first and second column represent the cell IDs being compared, e.g Cell 1
% vs Cell 2. The third column is the euclidean distance between cells and
% the fouth column is the correlation coefficient between them. Plots in Badimon 
% can be made by a scatter plot between column 3 and column 4. 

% Written by Aditya Nair (adi.nair@caltech.edu)

%% Input files

% All .mat files in the selected folder are analysed
% Direct to folder Example_Cell_Traces

folder = uigetdir();
fileList = dir(fullfile(folder, '*.mat'));

addpath(fileList.folder)

%%
Cell_Array = {};
for p = 1:length(fileList)
tic    
load(fileList(p).name)    

%% Pipeline Basic Extraction ----

%addpath 'E:\Matlab_2015a\bin'

%Correction for neuropil
FCor1 = F - 0.7*Fneu;

numberofcell = size(FCor1);

%Selection of cells classified by suite2p
FCor2 = single.empty;
for n = 1:numberofcell(1,1)

    
    FCorW = FCor1(n,:)*iscell(n,1);
    if sum(FCorW) > 0
        FCor2  = [FCor2;FCorW];
    end
        
    
end
  
%Obtaining df_f
numberofcellv2 = size(FCor2);
FCor3 = single(zeros(numberofcellv2));

for n = 1:numberofcellv2(1,2)
    
    FCor3(:,n) = (FCor2(:,n) - median(FCor2(:,:),2))./median(FCor2(:,:),2);
        
end

%Applying median filter
FCor4 = single(zeros(numberofcellv2));
for n = 1:numberofcellv2(1,1)

    
    FCor4(n,:) = medfilt1(double(FCor3(n,:)),10);
            
    
end

%%  Obtain Spatial Corr

Celllocs = single.empty;
for n = 1:numberofcell(1,1)

    
    CelllocsW = stat{1,n}.med(1)*iscell(n,1);
    if sum(CelllocsW) > 0
        Celllocs  = [Celllocs;stat{1,n}.med];
    end
        
    
end

Correlationclusall = corr(FCor4(:,:).');    

Cell_Loc_Dist = pdist(Celllocs);
Cell_Loc_Dist_Matrix = squareform(Cell_Loc_Dist);



C1 = triu(Cell_Loc_Dist_Matrix,1);
C2 = triu(Correlationclusall,1);

Cell_Loc_Corr_Sp =  single.empty;
for i = 1:numberofcellv2(1)
       
    for j = 1:numberofcellv2(1)
        if C1(i,j) == 0 && C2(i,j) == 0
        else
            
            Cell_Loc_Corr_Sp = [Cell_Loc_Corr_Sp;[[i,j],C1(i,j),C2(i,j)]];
           
        
        end
    end

end

Cell_Array{p} = Cell_Loc_Corr_Sp;


disp(['Completed Analyzing File ' fileList(p).name])

end
toc