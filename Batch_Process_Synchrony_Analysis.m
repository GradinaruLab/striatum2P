%% Synchrony Calculation

%% Output Format

% The output of the script Cell_Array contains cells with 5 columns. The
% first and second column represent the cell IDs being compared, e.g Cell 1
% vs Cell 2. The third column is the euclidean distance between cells and
% the fouth column is the correlation coefficient between them. 
% The fith column is synchrony calculated as per Badimon et al., 2020.

% Written by Aditya Nair (adi.nair@caltech.edu)

%% Obtain Files

% All .mat files in the selected folder are analysed
% Direct to folder Example_Cell_Traces

folder = uigetdir();
fileList = dir(fullfile(folder, '*.mat'));

addpath(fileList.folder)

%% Find threshold events

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

%if numberofcellv2(1,2) == 7181

FCor4_10min = FCor4(:,1:(numberofcellv2(1,2)/3));
FCor4_20min = FCor4(:,((numberofcellv2(1,2)/3) + 1):(numberofcellv2(1,2)/3*2));
FCor4_30min = FCor4(:,((numberofcellv2(1,2)/3*2) + 1):(numberofcellv2(1,2)));

%else
%    FCor4_10min = FCor4(:,1:(numberofcellv2(1,2)/3+ 1));
%    FCor4_20min = FCor4(:,((numberofcellv2(1,2)/3)):(numberofcellv2(1,2)/3*2));
%    FCor4_30min = FCor4(:,((numberofcellv2(1,2)/3*2)):(numberofcellv2(1,2)));

%end


requiredvar= FCor4_30min;

B1_events = zeros(length(requiredvar(:,1)),length(requiredvar(1,:)));
for i= 1:length(B1_events(:,1))
    
    for j = 1:length(B1_events(1,:))
      x = requiredvar(i,:);
      thresh = 2*std(x);
      idxl = x>=thresh;
      idxl(1) = 0;
      idx = find(idxl);
      yest = x(idx-1)<thresh; 
      final = idx(yest);
      x_req =  zeros(1,length(B1_events(1,:)));
      for k = 1:length(final)
          
         x_req(final(k)) = 1;
         x_req(final(k)-1) = 1;
         x_req(final(k)+1) = 1;
          
      end
      B1_events(i,:) = x_req(1:length(B1_events(i,:)));

    end
   
end

%% Perform synchrony measurement

synchronytest = zeros(length(requiredvar(:,1)),length(requiredvar(:,1)));
for m = 1:length(requiredvar(:,1))

i = B1_events(m,:);
for q = 1:length(requiredvar(:,1)) 

j = B1_events(q,:);
numberact = 0;  
for n = 1:length(requiredvar(1,:))
  if i(n) == 1 && j(n) == 1 
     numberact = numberact+1; 
  end   
end
synchronytest(m,q) = numberact/nnz(i) ;

end

end

synchronfinal = zeros(length(requiredvar(:,1)),length(requiredvar(:,1)));

for m = 1:length(requiredvar(:,1))
for q = 1:length(requiredvar(:,1)) 
synchronfinal(m,q) = (synchronytest(m,q) +synchronytest(q,m))/2 ;
end
end

%% Concatenate

%Select non zero cells

Celllocs = single.empty;
for n = 1:numberofcell(1,1)

    
    CelllocsW = stat{1,n}.med(1)*iscell(n,1);
    if sum(CelllocsW) > 0
        Celllocs  = [Celllocs;stat{1,n}.med];
    end
        
    
end

Correlationclusall = corr(requiredvar(:,:).');    

Cell_Loc_Dist = pdist(Celllocs);
Cell_Loc_Dist_Matrix = squareform(Cell_Loc_Dist);



C1 = triu(Cell_Loc_Dist_Matrix,1);
C2 = triu(Correlationclusall,1);
C3 = triu(synchronfinal,1);
Cell_Loc_Sync =  single.empty;
for i = 1:numberofcellv2(1)
       
    for j = 1:numberofcellv2(1)
        if C1(i,j) == 0 && C2(i,j) == 0
        else
            
            Cell_Loc_Sync = [Cell_Loc_Sync;[[i,j],C1(i,j),C2(i,j),C3(i,j)]];
           
        
        end
    end

end

Cell_Array{p} = Cell_Loc_Sync;

disp(['Completed Analyzing File ' fileList(p).name])
end
toc