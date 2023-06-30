%% Auditory brainstem responses - thresholds
% Plots figure panels 1G and S3B.

% INFO

% 2nd section plots ABR plots sequentially, for manual threshold
% determination. Uncomment to run.

% The subsequent sections use metadata based on manual trehsold 
% determination to visualise thresholds in figure panels S3B and 1G.

clear all
close all

group = 'F'; 
savevar = 0; % save variables for plotting startle/ABR correlations yes or no
dosave = 0; %to save figures
% clicks SD
        dpA = [-27 nan 1 24 39 56]; % list with start of window for respective peaks
        dpB = [-10 nan 23 38 55 76]; % list with ends of window for respective peaks
        
        dpA_ms = dpA/25;
        dpB_ms = dpB/25;
        


%user input -----
plotmarkers = 0; % to plot traces without (0) or with peak markers (1)
path = 'D:\PC-VV001-restoration\manuscript\Repository\ABR_data\';



% Built input matrix (Animals, Time & Paradigm)
    % for all IDs
    % for time 1 2 3
    clear InDet
    anmID = [1902 1903 1904 1709 1710 1711 1712];
    InDet1(1,:) = [1 2 4 8 16 20]; % Line#1 freq
    InDet2(1,:) = [3 3 3 3 3 1]; % Line#2 define paradigm (only right)
    InDet3(1,:) = ones(1, 6); % Line#3 for Time
    InDet4(1,:) = ones(1, 6)*1902; % Line#4 for animal
    

    InDet = [InDet1;InDet2;InDet3;InDet4];
    InDet = [];

    for n = 1:length(anmID)
        if anmID(n) == 1902 || anmID(n) == 1903 || anmID(n) == 1904
            for Time = 1:3
                  InDetAdd = [InDet1;InDet2;ones(1, 6)*Time;ones(1, 6)*anmID(n)];
                  InDet = [InDet InDetAdd];
            end
            
        elseif anmID(n) == 1709 || anmID(n) == 1710 || anmID(n) == 1711 || anmID(n) == 1712
            
            for Time = 1:2
                  InDetAdd = [InDet1;InDet2;ones(1, 6)*Time;ones(1, 6)*anmID(n)];
                  InDet = [InDet InDetAdd];
            end
            
        end
    end
    
% Add column ID in last row
InDet(5,:) = 1:length(InDet);

% Select columns of InDet in a random order
Seq = randperm(length(InDet(1,:)));

% InDet after randomisation
InDetSeq = ["Freq "; "Prdgm"; "Time "; "ID   "; "Column"];
for n  = 1:length(Seq)
InDetSeq = [InDetSeq InDet(:,Seq(n))];
end

% Save randomised file sequence as .mat file
t = datetime;
filetime = [num2str(t.Hour) '-' mat2str(t.Minute)];

filename = ['D:\PC-VV001-restoration\Ferrets\Ephys\final_scripts\ABR_metadata\ABR_RandSequence_' num2str(date) '-' filetime '.mat'];
save( filename,'InDetSeq');



paradigm = {'clickR';'clickL';'narrowR';'narrowL'};
channelcode = {'Left';'Right'};

lastpk = 6; % last peak to plot (example: lastpk = 4 plots peaks 1,2,3,4)


plotall = 1;
savetraces =0; %to save plots of ABR traces
savethreshold = 0; %to save plots with ABR thresholds pre and post NOE
SaveThresholdChange = 0;

% Index for correct loading of raw data from .txt files

    PigletPRE = [ 22 18406 3093]; % [clickR1 clickL1 narrowR narrowL(at21477)] Marks first 90dB line for each protocol,. e.g. clickR.
    Piglet5dPOST = [ 18928  544 21999]; % [clickR1 clickL1 narrowR narrowL(at3615)]
    PigletPOST2 = [22 22 21]; % [clickR1 clickL1 narrowR]

    KangaPRE = [22 21988 3093]; % [ clickR1 clickL1  narrowR narrowL(at25059)] 
    Kanga5dPOST = [18406 22 21477]; % [ clickR1 clickL1 narrowR narrowL(at3093)] 
    KangaPOST2 = [22 22 21]; % [ clickR1 clickL1 narrowR ] 

    HeffalumpPRE = [22 18406 3093 ]; % [clickR1 clickL1 NarrowR (NarrowL(unsure) at 21478)]
    Heffalump5dPOST = [21490 3106 24561]; % [clickR clickL NarrowR NarrowL at 6177]
    HeffalumpPOST2 = [22 22 21]; % [clickR1 clickL narrowR ]

    F1709_PRE = [544 22 21]; % [clickR clickL narrowR ]
    F1709_5dPOST = [22 22 21];

    F1710_PRE = [22 22 21]; % [clickR clickL narrowR ]
    F1710_5dPOST = [22 22 21];

    F1711_PRE = [22 22 21]; % [clickR clickL narrowR ]
    F1711_5dPOST = [22 22 21];

    F1712_PRE = [22 22 21]; % [clickR clickL narrowR ]
    F1712_5dPOST = [22 22 21];

    Index{1904} = [PigletPRE;Piglet5dPOST;PigletPOST2];
    Index{1903} = [KangaPRE;Kanga5dPOST;KangaPOST2];
    Index{1902} = [HeffalumpPRE;Heffalump5dPOST;HeffalumpPOST2];

    Index{1709} = [F1709_PRE; F1709_5dPOST];
    Index{1710} = [F1710_PRE; F1710_5dPOST];
    Index{1711} = [F1711_PRE; F1711_5dPOST];
    Index{1712} = [F1712_PRE; F1712_5dPOST];



%% Plot ABRs sequentially

% % USER INPUT
% 
% % ABR plot titles on/off.
% titleON = 0;
% 
% % define input matrix for this run (either randomised, InDetSeq, or not
% % randomised,InDet.
% 
% CurrentInput = str2double(InDetSeq(:,2:end)); % randomised sequence
% % CurrentInput = InDet;
% tic 
% for trgt = 1:length(CurrentInput)
%    
% % file details
% % for Nmouse = 1:size(IDs,1)
%     Nmouse = trgt; %:size(InDet,2)
%     
%     
%     if CurrentInput(4,Nmouse) == 1904; fname = 'Piglet';
%     elseif CurrentInput(4,Nmouse) == 1903; fname = 'Kanga';
%     elseif CurrentInput(4,Nmouse) == 1902; fname = 'Heffalump'; 
%     elseif CurrentInput(4,Nmouse) == 1709; fname = '1709';
%     elseif CurrentInput(4,Nmouse) == 1710; fname = '1710';
%     elseif CurrentInput(4,Nmouse) == 1711; fname = '1711';
%     elseif CurrentInput(4,Nmouse) == 1712; fname = '1712';
%     end
% 
% %    for time = 1:3 %1:2 [clickR clickL]
%    for time = CurrentInput(3,trgt) % specifcy time
% %         for    prdgm = prdgm_input %1:4
%         for prdgm = CurrentInput(2,trgt) % specifiy paradigm
%             for ch = 2 %1 for left electrode, 2 for right electrode
%                 clear rawtable raw
% 
%                 if prdgm == 3 || prdgm == 4
%                     intlist = [9 8 7 6 5 4]; % ints for tones
%                 elseif prdgm == 1 || prdgm == 2
%                     intlist = [9 8 7 6 5 4 ]; % ints for clicks [9 8 7 6 5 4 ]
%                 end
%                 
%                 %load file
%                 if CurrentInput(4,Nmouse) == 1904 || CurrentInput(4,Nmouse) == 1903 || CurrentInput(4,Nmouse) == 1902
%                     path = 'D:\PC-VV001-restoration\Ferrets\ABRdata\190X_ABRs_pre_1p_2p_3p\';
%                     if time == 1
% %                          rawtable1 = readtable([path fname '_ABRs_preNOE.txt'],'VariableNamingRule','preserve'); % edit 120822: find out how to edit this so strings are preserved after import. It seems that after matlba version 2020, strings are converted to NaNs by default.
%                     
%                            rawtable = readcell([path fname '_ABRs_preNOE.txt'],'VariableNamingRule','preserve'); % edit 120822: find out how to edit this so strings are preserved after import. It seems that after matlba version 2020, strings are converted to NaNs by default.
% 
%                     elseif time == 2
%                          rawtable = readcell([path fname '_ABRs_5DpNOE.txt']);
%                     elseif time == 3
%                         try
%                         rawtable = readcell([path 'txtfiles_pNOE2\' num2str(CurrentInput(4,Nmouse)) '_pNOE2_' paradigm{prdgm} '.txt']);
%                         catch
%                         rawtable = readcell([path 'txtfiles_pNOE2\' num2str(CurrentInput(4,Nmouse)) '_pNOE3_' paradigm{prdgm} '.txt']);
%                         end
%                     end
%                 else
%                     path = 'D:\PC-VV001-restoration\Ferrets\ABRdata\17XX_ABRs_pre_1p\';
%                      if time == 1
%                         rawtable = readcell([path num2str(CurrentInput(4,Nmouse)) '_' paradigm{prdgm} '_BL.txt']);
%                      elseif time == 2
%                         rawtable = readcell([path num2str(CurrentInput(4,Nmouse)) '_' paradigm{prdgm} '_POST1w.txt']);
%                      elseif time == 3
%                          continue
%                      end
%                  end
%                     
%                 %convert to array
% %                 raw = table2array(rawtable);
% %                 raw = table2array(rawtable(:,1)); % use only first column, edit 120822
%                 
%                 % if already array, keep only first, column
%                 raw = rawtable(:,1);
% 
%                 % step to next intensity 
%                 % cklickstep = length(22:531);
%                 % cklickstep = length(22:2872);
%                 if prdgm == 3 || prdgm == 4
%                 prtcstep = length(3093:18383)-1;
%                 elseif prdgm == 1 || prdgm == 2
%                 prtcstep = length(22:3071)-1;
%                 prtcstep = length(544:3593)-1;
% 
%                 end
%                 % define input
%                 % input = raw( Index(prdgm):Index(prdgm)+cklickstep*5);
%                 input = raw( Index{CurrentInput(4,Nmouse)}(time,prdgm):Index{CurrentInput(4,Nmouse)}(time,prdgm)+prtcstep );
%                 
%                 
%                 
%                 % find indexes for level
%                 meta_index{9} = Index{CurrentInput(4,Nmouse)}(time,prdgm) + find(strcmp(input, 'Level = 90 dB'));
%                 meta_index{8} = Index{CurrentInput(4,Nmouse)}(time,prdgm) + find(strcmp(input, 'Level = 80 dB'));
%                 meta_index{7} = Index{CurrentInput(4,Nmouse)}(time,prdgm) + find(strcmp(input, 'Level = 70 dB'));
%                 meta_index{6} = Index{CurrentInput(4,Nmouse)}(time,prdgm) + find(strcmp(input, 'Level = 60 dB'));
%                 meta_index{5} = Index{CurrentInput(4,Nmouse)}(time,prdgm) + find(strcmp(input, 'Level = 50 dB'));
%                 meta_index{4} = Index{CurrentInput(4,Nmouse)}(time,prdgm) + find(strcmp(input, 'Level = 40 dB'));
%                 % meta_index{3} = Index(prdgm) + find(strcmp(input, 'Level = 30 dB')); % for tones
%                 % meta_index{2} = Index(prdgm) + find(strcmp(input, 'Level = 20 dB')); % for tones
% 
%                 
%     % find indexes for intesnities and choose depending on no of intensity
%     % find index for level and frequency
%     clear index
%     if prdgm == 3 || prdgm == 4
%         for int = intlist
%             fcount = 1;
%             for freq = [1 2 4 8 16] % 1 2 4 8 16K sound frequency
%             index{freq}{int} = meta_index{int}(fcount:fcount+1)+1;
%             fcount = fcount+2; %go to value corresponding to next higher freq in idx vector
%             end
%         end
%     elseif prdgm == 1 || prdgm == 2
%         for int = intlist
%             fcount = 1;
%             for freq = [20] % 1 2 4 8 16K sound frequency, or 20 for clicks
%             index{freq}{int} = meta_index{int}(fcount:fcount+1);
%             fcount = fcount+2; %go to value corresponding to next higher freq in idx vector
%             end
%         end
%     end
% 
% 
% %             for freq = freqlist
%             for freq = CurrentInput(1,trgt)
% 
% 
%              % Define plot input matrix
% 
%             if prdgm == 1 || prdgm == 2
%             dpoints = length(24:266);
%             elseif prdgm == 3 || prdgm == 4
%             dpoints = length(15335:15578)-1;
%             end
%                 
%                 
%                             clear metaL metaR MetaL MetaR 
%                 % extract ABR trace for each sound level
%                 for int = intlist
%                     for t = 1:dpoints % number of datapoints per ABRS trace
%                         if prdgm == 1 || prdgm == 2 % for tones
%                             metaL{freq}(int,:) = raw( (index{freq}{int}(1)) : (index{freq}{int}(1)+dpoints) )';
%                             metaR{freq}(int,:) = raw( (index{freq}{int}(2)) : (index{freq}{int}(2)+dpoints) )';
%                         elseif prdgm == 3 || prdgm == 4 % for clicks
%                             metaL{freq}(int,:) = raw( (index{freq}{int}(1)) : (index{freq}{int}(1)+dpoints) )';
%                             metaR{freq}(int,:) = raw( (index{freq}{int}(2)) : (index{freq}{int}(2)+dpoints) )';
%                         end
%                     end
%                 % convert to numerical 
%                     for x = 1:length(metaL{freq}(int,:))
% %                         MetaL{freq}(int,x) = str2num( cell2mat(metaL{freq}(int,x)) );
% %                         MetaR{freq}(int,x) = str2num( cell2mat(metaR{freq}(int,x)) );
%                         % Edit 130822 
%                         MetaL{freq}(int,x) = cell2mat(metaL{freq}(int,x));
%                         MetaR{freq}(int,x) = cell2mat(metaR{freq}(int,x));
%                     end
%                 end
%                 
%                 
%                 
% 
%                 colourmap = [0 0.7 0.7];
%               
% 
%                 %% define input values
%                 % to do: compute per freq
% 
%                 % define if R or L 
%                 clear input
%                 if ch == 2 
%                 input{freq} = MetaR{freq};
%                 elseif ch == 1
%                     input{freq} = MetaL{freq};
%                 end
% 
%                 %% find peaks & throughs
%                 % between datapoints 25-75 
% 
%                 %find local maxima magnitude and locations
%                 clear pks locs trs...
%                     all_ons_locs all_ons_pks ons_locs ons_pks idx_prepk idx_postpk thr_locs thr_mag...
%                     loc_mxpk var
% 
% 
% 
% 
%                 %%  
%                     figure;hold all;
% 
%                     intcl = [ 0 3 6 8 12 15 18 21 24]; % indicators for shifting response window per intensity, relative to highest intesnity
%                     c=0; % intesnity counter
%                 for int = intlist 
%                     c = c+1;
%                     intc = intcl(c);
% %                     for peak = [2 1 3:lastpk]
% 
% 
% 
%               
% 
%                 % --- --- --- PLOT ABR traces --- --- ---
%                 
%                 % figures of ABR traces with markers
%                 if plotall ==1
% 
%                     thresh = 'n.a.';
%                 hold all
%                 if titleON == 1
%                     if time ==1
%                     title([ 'Ch' num2str(ch) ' PRE ' num2str(paradigm{prdgm}) ' - ' num2str(CurrentInput(4,Nmouse)) ' - treshold ' num2str(thresh) '0 dB - ' num2str(freq) 'K'])
%                     elseif time ==2
%                     title([ 'Ch' num2str(ch) ' POST ' num2str(paradigm{prdgm}) ' - ' num2str(CurrentInput(4,Nmouse)) ' - treshold ' num2str(thresh) '0 dB - ' num2str(freq) 'K'])
%                     elseif time ==3
%                     title([ 'Ch' num2str(ch) ' POST2 ' num2str(paradigm{prdgm}) ' - ' num2str(CurrentInput(4,Nmouse)) ' - treshold ' num2str(thresh) '0 dB - ' num2str(freq) 'K'])
%                     end
%                 end
% 
% %                 for peak = 1:lastpk
%                     
% %                
%                 
%                 if plotmarkers == 0
%                     for int = [ 9 8 7 6 5 4 3 2]
% %                         p(int) = plot(input{freq}(int,:)+int/200000,'-','linewidth',2,'Color',colourmap*(int/10) ); %ABR trace
%                         p(int) = plot(input{freq}(int,:)+int/200000,'-','linewidth',2,'Color',colourmap*0.5 ); %ABR trace
% 
%                     % save input to plot ABR traces separately
%                     ABRtrace{Nmouse}{time}{prdgm}{freq}(int,:) = input{freq}(int,:);
% 
%                          plot([1 length(input{freq}(int,:))],[var(int)*3 var(int)*3]+int/200000,'--','linewidth',0.5,'Color','k' );% variance line
%                          plot([1 length(input{freq}(int,:))],[0 0]+int/200000,'-','linewidth',0.5,'Color','k' );% variance line
% 
%                         % plot settings
%                     set(gca,'linewidth',1,'Fontsize',10)
%                     xlabel('ms')
%                     ylabel('Intensity (dB)')
%                     end
% 
%                 end
% 
% 
%                 yticks([2:9]/200000)
%                 yticklabels({'20';'30';'40';'50';'60';'70';'80';'90'})
%                 xticks([25:25:250])
%                 xticklabels({'1';'2';'3';'4';'5';'6';'7';'8';'9';'10'})
%                 
%                 
% 
% 
%               
% 
% %                 end
% 
%                 end
% 
% 
% %                 end
% 
%             end
%             end
%         end
%     end
%    end
%    pause
% %    close all
% end
% 
% 
% toc

%% Load table with threshold metadata

clear all

path = 'D:\PC-VV001-restoration\manuscript\Repository\ABR_data\Metadata\';
savepath = 'D:\PC-VV001-restoration\Ferrets\Behaviour\Figures\Time123\ABRs\';

% INFO: Thresholds_030223_17-19.xlsx are thresholds scored manually on
% 030223 using the ABR-sequence generated on the same day at 17:19
% (correspinding file ABR_RandSequence_03-Feb-2023-17-19.mat).

% load threshold data
T = readtable([ path 'Thresholds_030223_17-19.xlsx']);
Tarray = table2array(T)';

% load ABR_sequence metadata 
load([ path 'ABR_RandSequence_03-Feb-2023-17-19.mat'])


% Add selected thresholds to InDetSeq matrix.
% Remember to use entire Tarray eventually.

Thresholds_str = [InDetSeq(:,2:end); Tarray];
% Convert string matrix to doubles.
Thresholds = str2double(Thresholds_str);

% Plot all tresholds across time (1 2 30

IdxT1 = find(Thresholds(3,:)==1);
IdxT2 = find(Thresholds(3,:)==2);
IdxT3 = find(Thresholds(3,:)==3);



IdxT1BBN = find(Thresholds(3,:)==1 & Thresholds(1,:)==20);
IdxT2BBN = find(Thresholds(3,:)==2 & Thresholds(1,:)==20);
IdxT3BBN = find(Thresholds(3,:)==3 & Thresholds(1,:)==20);


%% ABR values, produce table

% IDX for each freq
for t = 1:3
    for freq = [1 2 4 8 16 20]
        IdxTSA{t}{freq} = [];
        for anm = [1709 1710 1711 1712 1902 1903 1904]
        % sort the idx values by animal number
        IdxTSA{t}{freq} =[IdxTSA{t}{freq}; find(Thresholds(3,:)== t & Thresholds(1,:)==freq & Thresholds(4,:)==anm)];
        end
    end
end
Thresholds(6,IdxTSA{1}{freq})
Thresholds';

for freq = [1 2 4 8 16 20]
threshtab = [Thresholds([6],IdxTSA{1}{freq})' Thresholds([6],IdxTSA{2}{freq})'  [Thresholds([6],IdxTSA{3}{freq})'; [nan nan nan nan]'] ];
% Construct table
Thresh = threshtab;
Thresh_table = array2table(Thresh); % convert to table
Thresh_table.Properties.VariableNames(1:3) = {'BL','Post1','Post2'}; % name the columns

% writetable(Thresh_table,['D:\PC-VV001-restoration\manuscript\Tables\ABR_RandSequence_03-Feb-2023-17-19_thresholds_' num2str(freq) 'K.xlsx']) % convert table to .xlsx file
end


%% Figure 1G: Plot mean thresholds across time for each frequency


close all

% IDX for each freq
for t = 1:3
    for freq = [1 2 4 8 16 20]
        IdxTS{t}{freq} = find(Thresholds(3,:)== t & Thresholds(1,:)==freq);
    end
end

% colourmap
cmap = [.5 .5 1; 1 .5 .5; 0.9 0.6 0.3];

figure
hold all
count = 0;
sctr = [0 0.2 0.4]; % x-axis scatter

% stats
for freq = [1 2 4 8 16 20]
[p2(freq), h2(freq)] = ranksum(Thresholds(6,IdxTS{1}{freq}), Thresholds(6,IdxTS{2}{freq}))
[p3(freq), h3(freq)] = ranksum(Thresholds(6,IdxTS{1}{freq}), Thresholds(6,IdxTS{3}{freq}))
end

for freq = [1 2 4 8 16 20]
[p2(freq), h2(freq)] = signrank(Thresholds(6,IdxTS{1}{freq}), Thresholds(6,IdxTS{2}{freq})  )
[p3(freq), h3(freq)] = signrank(Thresholds(6,IdxTS{1}{freq}), [Thresholds(6,IdxTS{3}{freq}) nan nan nan nan] )
end


for time = [1 2 3]
    
    % shift on x axis
    tinput = [0 .1 .2];
    count = count+1;
    t = tinput(count);
    
b1(time) = plot([1 2 3 4 5 6]+t,[nanmean(Thresholds(6,IdxTS{time}{1})') nanmean(Thresholds(6,IdxTS{time}{2})') nanmean(Thresholds(6,IdxTS{time}{4})')...
     nanmean(Thresholds(6,IdxTS{time}{8})') nanmean(Thresholds(6,IdxTS{time}{16})') nanmean(Thresholds(6,IdxTS{time}{20})')],'-')

b = errorbar([1 2 3 4 5 6]+t,[nanmean(Thresholds(6,IdxTS{time}{1})') nanmean(Thresholds(6,IdxTS{time}{2})') nanmean(Thresholds(6,IdxTS{time}{4})')...
    nanmean(Thresholds(6,IdxTS{time}{8})') nanmean(Thresholds(6,IdxTS{time}{16})') nanmean(Thresholds(6,IdxTS{time}{20})')],...
    [std(Thresholds(6,IdxTS{time}{1})')/sqrt(length(Thresholds(6,IdxTS{time}{1})))...
     std(Thresholds(6,IdxTS{time}{2})')/sqrt(length(Thresholds(6,IdxTS{time}{2})))...
     std(Thresholds(6,IdxTS{time}{4})')/sqrt(length(Thresholds(6,IdxTS{time}{4})))...
     std(Thresholds(6,IdxTS{time}{8})')/sqrt(length(Thresholds(6,IdxTS{time}{8})))...
     std(Thresholds(6,IdxTS{time}{16})')/sqrt(length(Thresholds(6,IdxTS{time}{16})))...
     std(Thresholds(6,IdxTS{time}{20})')/sqrt(length(Thresholds(6,IdxTS{time}{20})))],'k');
% main plot settings
 set(b,{'linew'},{2.5},'Color',cmap(count,:),'CapSize',3.5)
% legend plot settings
set(b1(time),{'linew'},{2.5},'Color',cmap(count,:))



end

    

xlim([0 7]);
ylim([30 105])
 xticks([1 2 3 4 5 6])
 xticklabels({'1 kHz NBN','2 kHz NBN','4 kHz NBN','8 kHz NBN*','16 kHz NBN','Clicks'})
 xtickangle(15)

 yticks([40 50 60 70 80 90]+2)
 yticklabels({'40','50','60','70','80','90'})
 ylabel('ABR threshold (dB SPL)')
 xlabel('Stimulus type')

legend([b1(1) b1(2) b1(3)],'BL','Post 1','Post 2','location','southwest')
set(gca, 'linewidth',1.5,'fontsize',15)
legend boxoff




%% produce SPPS matrix 
% including RTs & with 1,2,3... as IDs
Anmlvalues = [1902 1903 1904 1709 1710 1711 1712 ];

%  HR{time}{freq,counter};
 ID =[]; FREQ =[]; TIME = []; THRESH = [];
 
 IDc = 0; % ID counter
 for time = 1:3
     for freq = [1 2 4 8 16 20]
         for id = 1:7
              ID = [ID;id];
              FREQ = [FREQ;freq];
              try
              THRESH = [THRESH;Thresholds(6,IdxTS{time}{freq}(id))];
              catch
              THRESH = [THRESH;nan];
              end
              TIME=[TIME;time];
         end
     end
 end                            
                    

% Construct table
HRBLm = [ID FREQ TIME THRESH]; % all
HRBLm_table = array2table(HRBLm); % convert to table
HRBLm_table.Properties.VariableNames(1:4) = {'ID','FREQ','TIME','THRESH'}; % name the columns


%% Figure S3B: Plot mean thresholds across time for each frequency
% - As bar plots with error bar, grouped by stimulus on x axis.


% close all

% IDX for each freq
for t = 1:3
    for freq = [1 2 4 8 16 20]
        IdxTS{t}{freq} = find(Thresholds(3,:)== t & Thresholds(1,:)==freq);
    end
end

% colourmap
cmap = [.5 .5 1; 1 .5 .5; 0.9 0.6 0.3];

figure
hold all
count = 0;
sctr = [0 0.2 0.4]; % x-axis scatter

% stats

% RANKSUM
for freq = [1 2 4 8 16 20]
[p2(freq), h2(freq)] = ranksum(Thresholds(6,IdxTS{1}{freq}), Thresholds(6,IdxTS{2}{freq}))
[p3(freq), h3(freq)] = ranksum(Thresholds(6,IdxTS{1}{freq}), Thresholds(6,IdxTS{3}{freq}))
end
% SIGNRANK (paired)
for freq = [1 2 4 8 16 20]
[p2(freq), h2(freq)] = signrank(Thresholds(6,IdxTS{1}{freq}), Thresholds(6,IdxTS{2}{freq})  )
[p3(freq), h3(freq)] = signrank(Thresholds(6,IdxTS{1}{freq}), [Thresholds(6,IdxTS{3}{freq}) nan nan nan nan] )
end

% Produce table with SIGNRANK p values
threshtab1 = []
threshtab2 = []

for freq = [1 2 4 8 16 20]
        threshtab1 = [threshtab1 p2(freq)]; % BL vs Post1 p-values
    try
        threshtab2 = [threshtab2 p3(freq)]; % BL vs Post2 p-values
    catch
        threshtab2 = [threshtab2 nan];
    end
end

% Construct table
Thresh = [threshtab1; threshtab2] ;
Thresh_table = array2table(Thresh); % convert to table
Thresh_table.Properties.VariableNames(1:6) = {'1','2','4','8','16','20'}; % name the columns






for freq = [1 2 4 8 16 20]
    
    % shift on x axis
    tinput = [0 5 10 15 20 25];
    count = count+1;
    t = tinput(count);
    
% errorbar([1 2 3],[mean(Thresholds(6,IdxTS{1}{freq}))' mean(Thresholds(6,IdxTS{2}{freq})')...
%     mean(Thresholds(6,IdxTS{3}{freq})')],[std(Thresholds(6,IdxTS{1}{freq}))' std(Thresholds(6,IdxTS{2}{freq})')...
%     std(Thresholds(6,IdxTS{3}{freq})')],'o-','Color',cmap(count,:),'markerfacecolor',cmap(count,:),'linewidth',2.5)

b = errorbar([1 2 3]+t,[mean(Thresholds(6,IdxTS{1}{freq})') mean(Thresholds(6,IdxTS{2}{freq})') mean(Thresholds(6,IdxTS{3}{freq})') ],...
    [std(Thresholds(6,IdxTS{1}{freq})')/sqrt(length(Thresholds(6,IdxTS{1}{freq})))...
     std(Thresholds(6,IdxTS{2}{freq})')/sqrt(length(Thresholds(6,IdxTS{2}{freq})))...
     std(Thresholds(6,IdxTS{3}{freq})')/sqrt(length(Thresholds(6,IdxTS{3}{freq}))) ],'.k');
set(b,{'linew'},{1.5})

% save mean thresholds and standart errors
ThresholdsTable(freq,:) = [mean(Thresholds(6,IdxTS{1}{freq})') mean(Thresholds(6,IdxTS{2}{freq})') mean(Thresholds(6,IdxTS{3}{freq})') ]
TreshErrTable(freq,:) = [std(Thresholds(6,IdxTS{1}{freq})')/sqrt(length(Thresholds(6,IdxTS{1}{freq})))...
     std(Thresholds(6,IdxTS{2}{freq})')/sqrt(length(Thresholds(6,IdxTS{2}{freq})))...
     std(Thresholds(6,IdxTS{3}{freq})')/sqrt(length(Thresholds(6,IdxTS{3}{freq}))) ]


% Plot BARS
b1 = bar([1]+t,[mean(Thresholds(6,IdxTS{1}{freq})') ],1,'FaceColor',cmap(1,:),'LineWidth',1.5)
b2 = bar([2]+t,[mean(Thresholds(6,IdxTS{2}{freq})') ],1,'FaceColor',cmap(2,:),'LineWidth',1.5)
b3 = bar([3]+t,[mean(Thresholds(6,IdxTS{3}{freq})') ],1,'FaceColor',cmap(3,:),'LineWidth',1.5)



end

% 1K, n.s.
% 2K, BL vs P2, p<0.001
 plot([1+tinput(2) 3+tinput(2)], [94 94],'-k','linewidth',1) % line
 plot(mean([1+tinput(2) 3+tinput(2)])-0.4, 95, '*k','Markersize',4,'linewidth',1.1) % marker1
 plot(mean([1+tinput(2) 3+tinput(2)]), 95, '*k','Markersize',4,'linewidth',1.1) % marker1
 plot(mean([1+tinput(2) 3+tinput(2)])+0.4, 95, '*k','Markersize',4,'linewidth',1.1) % marker1

% 4K, BL vs P2, p <0.001
 plot([1+tinput(3) 3+tinput(3)], [94 94],'-k','linewidth',1) % line
 plot(mean([1+tinput(3) 3+tinput(3)])-0.4, 95, '*k','Markersize',4,'linewidth',1.1) % marker1
 plot(mean([1+tinput(3) 3+tinput(3)]), 95, '*k','Markersize',4,'linewidth',1.1) % marker1
 plot(mean([1+tinput(3) 3+tinput(3)])+0.4, 95, '*k','Markersize',4,'linewidth',1.1) % marker1

% 8K, BL vs P1, P<0.05, vs P2, p<0.001
 plot([1+tinput(4) 2+tinput(4)], [92 92],'-k','linewidth',1) % line
 plot(mean([1+tinput(4) 2+tinput(4)]), 93, '*k','Markersize',4,'linewidth',1.1) % marker1
 
 plot([1+tinput(4) 3+tinput(4)], [94 94],'-k','linewidth',1) % line
 plot(mean([1+tinput(4) 3+tinput(4)])-0.4, 95, '*k','Markersize',4,'linewidth',1.1) % marker1
 plot(mean([1+tinput(4) 3+tinput(4)]), 95, '*k','Markersize',4,'linewidth',1.1) % marker1
 plot(mean([1+tinput(4) 3+tinput(4)])+0.4, 95, '*k','Markersize',4,'linewidth',1.1) % marker1

% 16K, BL vs P1, p<0.05, vs P2, p<0.05
 plot([1+tinput(5) 2+tinput(5)], [92 92],'-k','linewidth',1) % line
 plot(mean([1+tinput(5) 2+tinput(5)]), 93, '*k','Markersize',4,'linewidth',1.1) % marker1
 
 plot([1+tinput(5) 3+tinput(5)], [94 94],'-k','linewidth',1) % line
 plot(mean([1+tinput(5) 3+tinput(5)]), 95, '*k','Markersize',4,'linewidth',1.1) % marker1
 
 
% 20K, BL vs P2, p<0.05
 plot([1+tinput(6) 3+tinput(6)], [94 94],'-k','linewidth',1) % line
 plot(mean([1+tinput(6) 3+tinput(6)]), 95, '*k','Markersize',4,'linewidth',1.1) % marker1

xlim([0 29]);
ylim([35 110])
 xticks([0 5 10 15 20 25]+2)
 xticklabels({'1kHz','2kHz','4kHz','8kHz*','16kHz','Clicks'})
 yticks([40 50 60 70 80 90]+2)
 yticklabels({'40','50','60','70','80','90'})
 ylabel('ABR threshold (dB SPL)')
% legend('1 kHz','2 kHz','4 kHz','8 kHz','16 kHz','BBN','location','southeast')
lgd = legend([b1 b2 b3],'BL','Post1','Post2')
set(lgd, 'box','off','Orientation','horizontal')
set(gca, 'linewidth',1.2,'fontsize',15)


%% Threshold decrease: per animal based on dB differences.
%  Save to file.



clear IdxTSA RatThresh ChanThresh AvChanThresh
% IDX for each freq and anm
anm = [1902 1903 1904];

for ID = 1:3
    freqc = 0; %initiate freq counter
    for freq = [1 2 4 8 16 20]
        freqc = freqc+1;
        for t = 1:3
            IdxTSA{ID}{freq}(t) = find(Thresholds(3,:)== t & Thresholds(1,:)==freq & Thresholds(4,:)==anm(ID));
            ThreshTSA{ID}{freq}(t) = Thresholds(6,IdxTSA{ID}{freq}(t));
        end
        % calculate change in ABR threshold
        RatThresh{ID}{1}(freqc) = ThreshTSA{ID}{freq}(2) - ThreshTSA{ID}{freq}(1); % Post1/BL *100
        RatThresh{ID}{2}(freqc) = ThreshTSA{ID}{freq}(3) - ThreshTSA{ID}{freq}(1); % Post2/BL *100
        ChanThresh{ID}{1}(freqc) = -RatThresh{ID}{1}(freqc); % Threshold increase depicted as negative value ('impairment')
        ChanThresh{ID}{2}(freqc) = -RatThresh{ID}{2}(freqc); % Threshold increase depicted as negative value ('impairment')
       % calculate average change in threshold across freqs
       AvChanThresh{ID}(1) = mean(ChanThresh{ID}{1}); 
       AvChanThresh{ID}(2) = mean(ChanThresh{ID}{2}); 

    end    
end




%  save('ABRthresholdChange','AvChanThresh')



%% Threshold decrease: per animal (selected freqs)
clear IdxTSA RatThresh ChanThresh1
% IDX for each freq and anm
anm = [1902 1903 1904];

for ID = 1:3
    for freq = [1 4 8 16 20]
        for t = 1:3
            IdxTSA{ID}{freq}(t) = find(Thresholds(3,:)== t & Thresholds(1,:)==freq & Thresholds(4,:)==anm(ID));
            ThreshTSA{ID}{freq}(t) = Thresholds(6,IdxTSA{ID}{freq}(t));
        end
        % calculate change in ABR threshold
        RatThresh{ID}{freq}(1) = ThreshTSA{ID}{freq}(2) / ThreshTSA{ID}{freq}(1)*100; % Post1/BL *100
        RatThresh{ID}{freq}(2) = ThreshTSA{ID}{freq}(3) / ThreshTSA{ID}{freq}(1)*100; % Post2/BL *100
        ChanThresh1{ID}{freq} = RatThresh{ID}{freq}-100; % Threshold chnage in % relative to BL
    end    
end

clear IdxTA ThreshTA RatThresh ChanThresh1
for ID = 1:3
        for t = 1:3
            IdxTA{ID}{t} = find(Thresholds(3,:)== t & Thresholds(4,:)==anm(ID) & Thresholds(1,:)~=2);
            ThreshTA{ID}(t) = mean( Thresholds(6,IdxTA{ID}{t}) );
        end
        % calculate change in ABR threshold
        RatThresh{ID}(1) = ThreshTA{ID}(2) / ThreshTA{ID}(1)*100; % Post1/BL *100
        RatThresh{ID}(2) = ThreshTA{ID}(3) / ThreshTA{ID}(1)*100; % Post2/BL *100
         ChanThresh1{ID} = RatThresh{ID}-100; % Threshold chnage in % relative to BL

end







