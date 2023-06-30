%% Auditory brainstem responses

% Plots figure panels 1F, S3A.
% Calculates total ABR magnitude.


clear all
close all

group = 'F';
savevar = 0; % save variables for plotting startle/ABR correlations yes or no
dosave = 0; %to save figures
% clicks SD
        dpA = [-27 nan 1 24 39 56]; % list with start of window for respective peaks
        dpB = [-10 nan 23 38 55 76]; % list with ends of window for respective peaks
        
        dpA_ms = dpA/25
        dpB_ms = dpB/25
        

IDs = [1902;1903;1904;1709;1710;1711;1712];

% make sure the data path is correct
path = 'D:\PC-VV001-restoration\manuscript\Repository\ABR_data';
% savepath = 'D:\PC-VV001-restoration\Ferrets\Behaviour\Figures\Time123\ABRs\';


% user input
freqlist = [20]; % 1,2,4,8,16K or 20 for clicks
prdgm_input = 1  %1; %3 % [1 2 3 4] for clickR, L, and narrowR, L
% specifiy prdgm below [ 1 2 3 4]

paradigm = {'clickR';'clickL';'narrowR';'narrowL'};
channelcode = {'Left';'Right'};

lastpk = 6; % last peak to plot (example: lastpk = 4 plots peaks 1,2,3,4)



plotall = 1;
savetraces =0; %to save plots of ABR traces
savethreshold = 0; %to save plots with ABR thresholds pre and post NOE
SaveThresholdChange = 0;



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



% ---------------
   
% file details
for Nmouse = 1:size(IDs,1)
    
    if IDs(Nmouse) == 1904; fname = 'Piglet';
    elseif IDs(Nmouse) == 1903; fname = 'Kanga';
    elseif IDs(Nmouse) == 1902; fname = 'Heffalump'; 
    elseif IDs(Nmouse) == 1709; fname = '1709';
    elseif IDs(Nmouse) == 1710; fname = '1710';
    elseif IDs(Nmouse) == 1711; fname = '1711';
    elseif IDs(Nmouse) == 1712; fname = '1712';
    end

%     close all
   for time = 1:3 %1:2 [clickR clickL]
        for    prdgm = prdgm_input %1:4
            for ch = 2 %1 for left electrode, 2 for right electrode
                clear rawtable raw

                if prdgm == 3 || prdgm == 4
                    intlist = [9 8 7 6 5 4]; % ints for tones
                elseif prdgm == 1 || prdgm == 2
                    intlist = [9 8 7 6 5 4 ]; % ints for clicks [9 8 7 6 5 4 ]
                end
                
                %load file
                if IDs(Nmouse) == 1904 || IDs(Nmouse) == 1903 || IDs(Nmouse) == 1902
                    path = 'D:\PC-VV001-restoration\Ferrets\ABRdata\190X_ABRs_pre_1p_2p_3p\';
                    if time == 1
%                          rawtable1 = readtable([path fname '_ABRs_preNOE.txt'],'VariableNamingRule','preserve'); % edit 120822: find out how to edit this so strings are preserved after import. It seems that after matlba version 2020, strings are converted to NaNs by default.
                    
                           rawtable = readcell([path fname '_ABRs_preNOE.txt'],'VariableNamingRule','preserve'); % edit 120822: find out how to edit this so strings are preserved after import. It seems that after matlba version 2020, strings are converted to NaNs by default.

                    elseif time == 2
                         rawtable = readcell([path fname '_ABRs_5DpNOE.txt']);
                    elseif time == 3
                        try
                        rawtable = readcell([path 'txtfiles_pNOE2\' num2str(IDs(Nmouse)) '_pNOE2_' paradigm{prdgm} '.txt']);
                        catch
                        rawtable = readcell([path 'txtfiles_pNOE2\' num2str(IDs(Nmouse)) '_pNOE3_' paradigm{prdgm} '.txt']);
                        end
                    end
                else
                    path = 'D:\PC-VV001-restoration\Ferrets\ABRdata\17XX_ABRs_pre_1p\';
                     if time == 1
                        rawtable = readcell([path num2str(IDs(Nmouse)) '_' paradigm{prdgm} '_BL.txt']);
                     elseif time == 2
                        rawtable = readcell([path num2str(IDs(Nmouse)) '_' paradigm{prdgm} '_POST1w.txt']);
                     elseif time == 3
                         continue
                     end
                 end
                    
                %convert to array
%                 raw = table2array(rawtable);
%                 raw = table2array(rawtable(:,1)); % use only first column, edit 120822
                
                % if already array, keep only first, column
                raw = rawtable(:,1);

                % step to next intensity 
                % cklickstep = length(22:531);
                % cklickstep = length(22:2872);
                if prdgm == 3 || prdgm == 4
                prtcstep = length(3093:18383)-1;
                elseif prdgm == 1 || prdgm == 2
                prtcstep = length(22:3071)-1;
                prtcstep = length(544:3593)-1;

                end
                % define input
                % input = raw( Index(prdgm):Index(prdgm)+cklickstep*5);
                input = raw( Index{IDs(Nmouse)}(time,prdgm):Index{IDs(Nmouse)}(time,prdgm)+prtcstep );
                
                
                
                % find indexes for level
                meta_index{9} = Index{IDs(Nmouse)}(time,prdgm) + find(strcmp(input, 'Level = 90 dB'));
                meta_index{8} = Index{IDs(Nmouse)}(time,prdgm) + find(strcmp(input, 'Level = 80 dB'));
                meta_index{7} = Index{IDs(Nmouse)}(time,prdgm) + find(strcmp(input, 'Level = 70 dB'));
                meta_index{6} = Index{IDs(Nmouse)}(time,prdgm) + find(strcmp(input, 'Level = 60 dB'));
                meta_index{5} = Index{IDs(Nmouse)}(time,prdgm) + find(strcmp(input, 'Level = 50 dB'));
                meta_index{4} = Index{IDs(Nmouse)}(time,prdgm) + find(strcmp(input, 'Level = 40 dB'));
                % meta_index{3} = Index(prdgm) + find(strcmp(input, 'Level = 30 dB')); % for tones
                % meta_index{2} = Index(prdgm) + find(strcmp(input, 'Level = 20 dB')); % for tones

                
    % find indexes for intesnities and choose depending on no of intensity
    % find index for level and frequency
    clear index
    if prdgm == 3 || prdgm == 4
        for int = intlist
            fcount = 1;
            for freq = [1 2 4 8 16] % 1 2 4 8 16K sound frequency
            index{freq}{int} = meta_index{int}(fcount:fcount+1)+1;
            fcount = fcount+2; %go to value corresponding to next higher freq in idx vector
            end
        end
    elseif prdgm == 1 || prdgm == 2
        for int = intlist
            fcount = 1;
            for freq = [20] % 1 2 4 8 16K sound frequency, or 20 for clicks
            index{freq}{int} = meta_index{int}(fcount:fcount+1);
            fcount = fcount+2; %go to value corresponding to next higher freq in idx vector
            end
        end
    end


            for freq = freqlist

             % Define plot input matrix

            if prdgm == 1 || prdgm == 2
            dpoints = length(24:266);
            elseif prdgm == 3 || prdgm == 4
            dpoints = length(15335:15578)-1;
            end
                
                
                            clear metaL metaR MetaL MetaR 
                % extract ABR trace for each sound level
                for int = intlist
                    for t = 1:dpoints % number of datapoints per ABRS trace
                        if prdgm == 1 || prdgm == 2 % for tones
                            metaL{freq}(int,:) = raw( (index{freq}{int}(1)) : (index{freq}{int}(1)+dpoints) )';
                            metaR{freq}(int,:) = raw( (index{freq}{int}(2)) : (index{freq}{int}(2)+dpoints) )';
                        elseif prdgm == 3 || prdgm == 4 % for clicks
                            metaL{freq}(int,:) = raw( (index{freq}{int}(1)) : (index{freq}{int}(1)+dpoints) )';
                            metaR{freq}(int,:) = raw( (index{freq}{int}(2)) : (index{freq}{int}(2)+dpoints) )';
                        end
                    end
                % convert to numerical 
                    for x = 1:length(metaL{freq}(int,:))
%                         MetaL{freq}(int,x) = str2num( cell2mat(metaL{freq}(int,x)) );
%                         MetaR{freq}(int,x) = str2num( cell2mat(metaR{freq}(int,x)) );
                        % Edit 130822 
                        MetaL{freq}(int,x) = cell2mat(metaL{freq}(int,x));
                        MetaR{freq}(int,x) = cell2mat(metaR{freq}(int,x));
                    end
                end
                
                
%         

                colourmap = [0 0.7 0.7];

%                 figure
%                 plot(MetaL{freq}(9,:),'-')


                %% define input values
                % to do: compute per freq

                % define if R or L 
                clear input
                if ch == 2 
                input{freq} = MetaR{freq};
                elseif ch == 1
                    input{freq} = MetaL{freq};
                end

                %% find peaks & throughs
                % between datapoints 25-75 

                %find local maxima magnitude and locations
                clear pks locs trs...
                    all_ons_locs all_ons_pks ons_locs ons_pks idx_prepk idx_postpk thr_locs thr_mag...
                    loc_mxpk var




                %%  
                    figure;hold all;

                    intcl = [ 0 3 6 8 12 15 18 21 24]; % indicators for shifting response window per intensity, relative to highest intesnity
                    c=0; % intesnity counter
                for int = intlist 
                    c = c+1;
                    intc = intcl(c);
                    % find magnitude and location of all local peaks
                    [pks{int},locs{int}] = findpeaks(input{freq}(int,:));

                    % find magnitude and location of all local throughs
                    [trs_meta{int},locs_trs{int}] = findpeaks(-input{freq}(int,:));
                    trs{int} = -trs_meta{int}; %reverse negative values again

                    % DEFINE VARIRABLITY VECTOR
                %     var(int) = std(input{freq}(int,1:40)); % fixed: baseline of the trace
                    var(int) = std(input{freq}(36+intc:96+intc)); % moving window: response window



                    %---
                    % FIRST PEAK
                    % select peaks in expected response window (FIRST pos. peak)
                %     metapkidx = find(locs{int}<=65 & locs{int}>= 40 ); % ixd of peak
                %     locations in window used for threshold detection (18/08/21)

                %     metapkidx = find(locs{int}<= 66+intc & locs{int}>= 36+intc ); % eff. 70 and 40; ixd of peak locations in window % old
                    for peak = [2 1 3:lastpk]

                        if peak == 2 % find second peak in given window
                           metapkidx = find(locs{int}<= 60+intc & locs{int}>= 39+intc ); % previously 45, 060921
                        else % find other peaks relative to second peak ( here loc_mxpk{2} )
                            try
                            metapkidx = find(locs{int}<= loc_mxpk{2}{int}+dpB(peak)+intc & locs{int}>= loc_mxpk{2}{int}+dpA(peak)+intc ); %loc_mxpk refers to peak 2 here
                            end
                        end
          

                    % predefine meta matrixes
                    all_ons_locs{int} = [];
                    all_ons_pks{int} = [];

                    for n = 1:length(metapkidx) 
                    all_ons_locs{int}(n) = locs{int}(metapkidx(n)); % peak locations in window
                    all_ons_pks{int}(n) = pks{int}(metapkidx(n)); % peak magnitude in window
                    end  



                    % --- exclude 'early' peaks, rel. to correspnding high int. responses ---
                    if int == 9 % for highest intensity
                    ons_locs{int} = all_ons_locs{int}; % peak locations in window
                    ons_pks{int} = all_ons_pks{int}; % peak magnitude in window
                    end

                    if int < 9 % for lower intensities
                        for n = 1:length(metapkidx) 
                            if isempty(loc_mxpk{peak}{int+1}) % check if higher int has a response
                                continue                %...if not, ignore current peak
                            else
                               if all_ons_locs{int}(n) < loc_mxpk{peak}{int+1}-1 % also ignore peak if it is earlier than higher intensity peak,
                                                                           % granting a 'buffer' of 3 dp
                                 ons_locs{int}(n) = nan; % peak locations in window
                                 ons_pks{int}(n) = nan; % peak magnitude in window
                                   continue
                               else ons_locs{int}(n) = locs{int}(metapkidx(n)); % peak locations in window
                                 ons_pks{int}(n) = pks{int}(metapkidx(n)); % peak magnitude in window
                               end
                            end
                        end
                    end

                    % --- find peak ---
                    % Find magnitude and location of peak closest to corresponding peak for
                    % the pre condtion (for highest intesntity, 90 dB) or corresponding to the next higher intesnity
                    % within the same condtion for intesnity below 90 dB. 

                    % Highest intesnity (90dB)
                    if time == 2 && int == 9 || time == 3 && int == 9
                        % select peak closest to PRE NOE peak
                        try
                        pkmetdiff = abs( ons_locs{int}-pklat{peak}{Nmouse}{1}{freq}(ch,int) );
                        mxpk{peak}{int} = ons_pks{int} (find(pkmetdiff == min(pkmetdiff)));
                        loc_mxpk{peak}{int} = ons_locs{int}(find( ons_pks{int} == mxpk{peak}{int} ));
                           % if there is no corresponding peak detected in the pre
                           % condtion,accept detected peak by default
                           if isnan(pklat{peak}{Nmouse}{1}{freq}(ch,int)); mxpk{peak}{int}=ons_pks{int}; loc_mxpk{peak}{int}=ons_locs{int};end        
                           % if there are two peaks with egual distance to corresp. high
                           % intensity peak, pick the later one
                           if length(mxpk{peak}{int})>1;  mxpk{peak}{int} =  mxpk{peak}{int}(2); loc_mxpk{peak}{int}=loc_mxpk{peak}{int}(2);
                           end
                        catch
                            display('WARNING: PRE NOE day possibly not loaded')
                        mxpk{peak}{int} = [];
                        loc_mxpk{peak}{int} = [];
                        end
                    elseif time == 1 && int == 9 % This is where mxpk and loc_mxpk are defined for the first time.
                        try
                    %find magnitude and location of highest peak (within 25-75 dp)
                        mxpk{peak}{int} = max(ons_pks{int}); % maximum peak for this intensity
                        loc_mxpk{peak}{int} = ons_locs{int}(find( ons_pks{int} == mxpk{peak}{int} )); %location of maximum peak
                        catch
                        mxpk{peak}{int} = [];
                        loc_mxpk{peak}{int} = [];
                        end
                    end

                    % Lower intensities
                    if int < 9
                        try
                        pkmetdiff = abs(ons_locs{int}-loc_mxpk{peak}{int+1});
                        mxpk{peak}{int} = ons_pks{int} (find(pkmetdiff == min(pkmetdiff)));
                        loc_mxpk{peak}{int} = ons_locs{int}(find( ons_pks{int} == mxpk{peak}{int} ));
                           % if there are two peaks with egual distance to corresp. high
                           % intensity peak, pick the later one
                           if length(mxpk{peak}{int})>1;  mxpk{peak}{int} =  mxpk{peak}{int}(2); loc_mxpk{peak}{int}=loc_mxpk{peak}{int}(2);
                           end
                        catch
                    mxpk{peak}{int} = [];
                    loc_mxpk{peak}{int} = [];
                        end
                    end




                    % find neighbrouing throughs
                    if isempty(mxpk{peak}{int}) == 1 || mxpk{peak}{int} == 0 || isnan( mxpk{peak}{int})
                       thr_locs{peak}(int,:) = [ 0 0 ]; % if there is no defined peak, leave throughs data empty
                       thr_mag{peak}(int,:) = [0 0 ];
                    else
                        idx_prepk = find (locs_trs{int} < loc_mxpk{peak}{int});
                        idx_postpk = find (locs_trs{int} > loc_mxpk{peak}{int});
                          % select closest throughs
                          if isempty(idx_prepk); trh_locsA = nan; thr_magA = nan; else trh_locsA = locs_trs{int}(idx_prepk(end)); thr_magA = trs{int}(idx_prepk(end));  end
                          if isempty(idx_postpk); trh_locsB = nan; thr_magB = nan; else trh_locsB = locs_trs{int}(idx_postpk(1)); thr_magB = trs{int}(idx_postpk(1)); end

                        thr_locs{peak}(int,:) = [trh_locsA trh_locsB ]; 
                        thr_mag{peak}(int,:) = [thr_magA thr_magB];
                    end

                    end


               
              
                % figures of ABR traces with markers
                if plotall ==1

                    thresh = 'n.a.'

                hold all
                if time ==1
                title([ 'Ch' num2str(ch) ' PRE ' num2str(paradigm{prdgm}) ' - ' num2str(IDs(Nmouse,:)) ' - treshold ' num2str(thresh) '0 dB - ' num2str(freq) 'K'])
                elseif time ==2
                title([ 'Ch' num2str(ch) ' POST ' num2str(paradigm{prdgm}) ' - ' num2str(IDs(Nmouse,:)) ' - treshold ' num2str(thresh) '0 dB - ' num2str(freq) 'K'])
                elseif time ==3
                title([ 'Ch' num2str(ch) ' POST2 ' num2str(paradigm{prdgm}) ' - ' num2str(IDs(Nmouse,:)) ' - treshold ' num2str(thresh) '0 dB - ' num2str(freq) 'K'])
                end

                for peak = 1:lastpk
                for int = [ 9 8 7 6 5 4 3 2]
                    p(int) = plot(input{freq}(int,:)+int/200000,'-','linewidth',2,'Color',colourmap*(int/10) ); %ABR trace
                %     p(int) = plot(detrend(input{freq}(int,:))+int/200000,'-','linewidth',1.5,'Color','red' ); %ABR trace
                
                % save input to plot ABR traces separately
                ABRtrace{Nmouse}{time}{prdgm}{freq}(int,:) = input{freq}(int,:);

                    try
                    p1(int) = plot(loc_mxpk{peak}{int},mxpk{peak}{int}+int/200000,'*','linewidth',1.5,'Color','r','MarkerSize',5 ); %primary peak
                    p2(int) = plot(thr_locs{peak}(int,:),thr_mag{peak}(int,:)+int/200000,'*','linewidth',1.5,'Color','b','MarkerSize',5 ); %primary throughs
                    catch;end
                    
                % save peak marker info for plotting ABR traces separately
                ABRtrace_loc_mxpk{Nmouse}{time}{prdgm}{freq}{int}{peak} = loc_mxpk{peak}{int};
                ABRtrace_mxpk{Nmouse}{time}{prdgm}{freq}{int}{peak} = mxpk{peak}{int};
                ABRtrace_loc_thr{Nmouse}{time}{prdgm}{freq}{int}{peak} = thr_locs{peak}(int,:);
                ABRtrace_thrmag{Nmouse}{time}{prdgm}{freq}{int}{peak} = thr_mag{peak}(int,:);

                    %     pvar1(int) = plot([1 length(MetaR(int,:))],[var(int) var(int)]+int/200000,'--','linewidth',0.5,'Color','k' )% variance line
                     plot([1 length(input{freq}(int,:))],[var(int)*3 var(int)*3]+int/200000,'--','linewidth',0.5,'Color','k' );% variance line
                     plot([1 length(input{freq}(int,:))],[0 0]+int/200000,'-','linewidth',0.5,'Color','k' );% variance line

                    % plot settings
                set(gca,'linewidth',1,'Fontsize',10)
                xlabel('ms')
                ylabel('Intensity (dB)')
                end


                yticks([2:9]/200000)
                yticklabels({'20';'30';'40';'50';'60';'70';'80';'90'})
                xticks([25:25:250])
                xticklabels({'1';'2';'3';'4';'5';'6';'7';'8';'9';'10'})


%                 if savetraces == 1
%                     if time == 1
%                          print('-r750','-dtiff',[savepath 'traces\Ch' num2str(ch) ' PRE '...
%                           num2str(paradigm{prdgm}) ' - ' num2str(IDs(Nmouse,:)) ],'-painters');
%                     else
%                          print('-r750','-dtiff',[savepath 'traces\Ch' num2str(ch) ' POST '...
%                           num2str(paradigm{prdgm}) ' - ' num2str(IDs(Nmouse,:)) ],'-painters');
%                     end
%                 end




                %% Peak latency 
                for int = intlist 
                    if isempty(loc_mxpk{peak}{int}); pklat{peak}{Nmouse}{time}{freq}(ch,int) = nan;
                    else
                    pklat{peak}{Nmouse}{time}{freq}(ch,int) = loc_mxpk{peak}{int};
                    end
                end

                display NOW PKLAT
                try
                pklat{1}{Nmouse}{1}{freq}(ch,int)
                catch
                    display EMPTY
                end

                % for int = [ 9 8 7 6 5 4 3 2]
                %     if isempty(loc_mxpk2{peak}{int}); pklat2{peak}{Nmouse}{time}{freq}(ch,int) = nan;
                %     else
                %     pklat2{peak}{Nmouse}{time}{freq}(ch,int) = loc_mxpk2{peak}{int};
                %     end
                % end
                %% Calculate peak magnitude
                
       
                for int = [ 9 8 7 6 5 4 3 2]
                    if isempty(mxpk{peak}{int}-thr_mag{peak}(int,2)); pkmag{peak}{Nmouse}{time}{freq}(ch,int) = 0;
                    else
                    pkmag{peak}{Nmouse}{time}{freq}(ch,int) = mxpk{peak}{int}-thr_mag{peak}(int,2);
                    end
                end

                %% Calculate RMS for each intensity

                % range 30 - 200 seems to include the whole response (Afour 90 dB)
                % start w range 40 - 100, but shift it by 4 for every intensity
                if prdgm == 1 || prdgm == 3
                    rmsintc = 0;
                    for int = intlist
                        rmsintc = rmsintc +4;
                    rmslvl_ClR{peak}{Nmouse}{time}{freq}(ch, int) = rms(input{freq}(int,36+rmsintc:96+rmsintc));
                    end
                elseif prdgm == 2 || prdgm == 4
                            rmsintc = 0;
                    for int = intlist
                            rmsintc = rmsintc +4;
                    rmslvl_ClL{peak}{Nmouse}{time}{freq}(ch, int) = rms(input{freq}(int,36+rmsintc:96+rmsintc));
                    end
                end

                end

                end


                end

               

                %based on magnitude of detected peaks
                thresh = 9;
                for peak = 1:5
                    for int = intlist
                        crit = 3*var(int); % threshold criterium based on variance of response window

                        if isempty(mxpk{peak}{int}) % if no peak detected, keep prev intensity as threshold
                        else if int < thresh... % if peak detected and current intesnity is lower then current treshold...
                                && mxpk{peak}{int}-thr_mag{peak}(int,1)>= crit || mxpk{peak}{int}-thr_mag{peak}(int,2)>= crit % ...and if peak is high enough
                                thresh = int; % ...take current intensity as threshold
                             end    
                        end       
                    end
                end

                threshold{Nmouse}{time}{freq}(ch) = thresh;

                hold all
                if time ==1
                title([ 'Ch' num2str(ch) ' PRE ' num2str(paradigm{prdgm}) ' - ' num2str(IDs(Nmouse,:))  num2str(freq) 'K'])
                else
                title([ 'Ch' num2str(ch) ' POST ' num2str(paradigm{prdgm}) ' - ' num2str(IDs(Nmouse,:)) num2str(freq) 'K'])
                end

            end
            end
        end
    end
end


%% Figure 1F: Plot example ABR trace (for all intensities)
% animal 1, click R, 90dB
close all

try
    figure;hold all
for int = 4:9;
prdgm = 1; % clickR
time = 1; % BL
Nmouse = 1;

% define input
ABRplot_input = ABRtrace{Nmouse}{time}{prdgm}{freq}(int,:)*10^6; % *10^6 to convert into uV
% define input for markers
for peak = 1:6
ABR_locmxpk{peak} = ABRtrace_loc_mxpk{Nmouse}{time}{prdgm}{freq}{int}{peak};
ABR_mxpk{peak} = ABRtrace_mxpk{Nmouse}{time}{prdgm}{freq}{int}{peak}*10^6;
ABR_locthr{peak} = ABRtrace_loc_thr{Nmouse}{time}{prdgm}{freq}{int}{peak};
ABR_mxthr{peak} = ABRtrace_thrmag{Nmouse}{time}{prdgm}{freq}{int}{peak}*10^6;
end

% plot
% figure; hold all
% plot(ABRplot_input+int/1,'-','linewidth',3,'Color',colourmap*(int/10) ); %ABR trace
plot(ABRplot_input+int/1,'-','linewidth',3,'Color',[.8 .8 .8]*(int/10) ); %ABR trace
plot(ABRplot_input+int/1,'-','linewidth',3,'Color',[.6 .6 .6]*(int/10) ); %ABR trace

% plot(ABRplot_input,'-','linewidth',3,'Color',[.6 .6 .6]); %ABR trace

if int == 9
for peak = 1:6
plot(ABR_locmxpk{peak},ABR_mxpk{peak}+int/1,'O','linewidth',2,'Color',[0 1 0],'MarkerSize',4,'markerfacecolor',[0 1 0]); %primary peak
plot(ABR_locthr{peak},ABR_mxthr{peak}+int/1,'O','linewidth',2,'Color',[ 0 .5 0],'MarkerSize',4,'markerfacecolor',[ 0 .5 0]); %primary throughs

end
end

% plot settings

xlim([5 6*25])
set(gca,'linewidth',1.5,'Fontsize',14)
xlabel('Time (ms)')
ylabel('Stimulus intesnity (dB SPL)')
% yticks([2:9]/200000)
% yticklabels({'20';'30';'40';'50';'60';'70';'80';'90'})
xticks([25:25:250])
xticklabels({'1';'2';'3';'4';'5';'6';'7';'8';'9';'10'})
yticks([4:1:9])
yticklabels({'40';'50';'60';'70';'80';'90'})


end
end





%% sort into SPSS matrix
% close all
% peak magnitude I-V
% pkmag{peak}{Nmouse}{time}{freq}(ch,int)
PKMAG=[];peakID=[];TIME=[];ID=[];COND=[];PKLAT=[];

for t = 1:3
    for peak = 1:5
            for Nmouse = 1:size(IDs,1)
                clear input
                try
                input1 = pkmag{peak}{Nmouse}{t}{freq}(ch,9);
                input2 = pklat{peak}{Nmouse}{t}{freq}(ch,9);
                catch
                    input1 = nan; input2 = nan;
                end

              PKMAG = [PKMAG; input1];
              PKLAT = [PKLAT; input2/25];% it is 25dP per ms
              peakID = [peakID; peak];
              TIME = [TIME; t];
              ID = [ID; Nmouse ];

            end
    end
end

% Construct SPSS table
VpeakData2 = [ID PKMAG PKLAT peakID TIME];
Vpeak_table2 = array2table(VpeakData2); % convert to table
Vpeak_table2.Properties.VariableNames(1:5) = {'ID','PKMAG','PKLAT','peakID','TIME'}; % name the columns
 



%% Response modulation analysis

%% Plot peak latency, magnitude and RMS - slope & area
% separate plots: 4 plots, slopes post 1 and 2, area post 1 and 2

peak = 2
% dosave= 1

doplot = 0; % 2 to plot fitted slopes and areas under the curve

clear Outarrat2 Outardiff2 Outsloperat2 Outarrat Outardiff Outsloperat Out_slope

% pklat{cond}{peak}{Nmouse}{time}{freq}(ch,int)
% consider excluding peak for 40 dB, as often not present
ylabels = {'Peak latency','Peak magnitude (?V)','RMS'};
for c = 1:3
%     close all
        for Nmouse = 1:size(IDs,1) %1:7 %1:size(IDs,1)
            clear y1 y2 y3 nanidx1 nanidx2 ints 

        if c == 1 % Peak latency
            y1 = pklat{peak}{Nmouse}{1}{freq}(ch,intlist(end):9); % time 1
            y2 = pklat{peak}{Nmouse}{2}{freq}(ch,intlist(end):9); % time 2
            try
            y3 = pklat{peak}{Nmouse}{3}{freq}(ch,intlist(end):9); % time 3
            end
        elseif c == 2 % Peak magnitude
            y1 = pkmag{peak}{Nmouse}{1}{freq}(ch,intlist(end):9) %*10^6; % *10^6 to convert to ?V
            y2 = pkmag{peak}{Nmouse}{2}{freq}(ch,intlist(end):9) %*10^6; % *10^6 to convert to ?V
            try
            y3 = pkmag{peak}{Nmouse}{3}{freq}(ch,intlist(end):9) %*10^6; % *10^6 to convert to ?V
            end
        elseif c== 3 % RMS
            y1 = rmslvl_ClR{peak}{Nmouse}{1}{freq}(ch,intlist(end):9);
            y2 = rmslvl_ClR{peak}{Nmouse}{2}{freq}(ch,intlist(end):9);
            try
            y3 = rmslvl_ClR{peak}{Nmouse}{3}{freq}(ch,intlist(end):9);
            end
        end
        
        ints1 = [intlist(end)*10:10:90];
        ints2 = [intlist(end)*10:10:90];
        ints3 = [intlist(end)*10:10:90];

        % exclude NaNs and crop vectors accordingly (make them comparable)
        nanidx1 = find(isnan(y1)==1); %find nans
        nanidx2 = find(isnan(y2)==1); %find nans
        try
        nanidx3 = find(isnan(y3)==1); %find nans
        end
        
        y1([nanidx1])=[];%crop
        y2([nanidx2])=[];%crop
        try
        y3([nanidx3])=[];%crop
        end
        
        ints1([nanidx1])=[];%crop
        ints2([nanidx2])=[];%crop
        try
        ints3([nanidx3])=[];%crop
        end
        
          % CRITERION: if any value is larger than the adjacent higher
          % intensity value
          outidx1 = []; outidx2 = []; outidx3 = [];
        if c == 2 || c == 3
            for n =[length(y1)-2]:-1:1
                if y1(n) > y1(n+1) % outlier if value larger than previous (high intnesity) value
                        try y1(n)=mean([y1(n-1) y1(n+1)]); catch y1(n)=y1(n+1); end 
                        outidx1 = [outidx1 n];
                end
            end
            for n =[length(y2)-2]:-1:1
                if y2(n) > y2(n+1) % outlier if value larger than previous (high intnesity) value
                        try y2(n)=mean([y2(n-1) y2(n+1)]); catch y2(n)=y2(n+1); end 
                        outidx2 = [outidx2 n];
                end
            end
            try
            for n =[length(y3)-2]:-1:1
                if y3(n) > y3(n+1) % outlier if value larger than previous (high intnesity) value
                        try y3(n)=mean([y3(n-1) y3(n+1)]); catch y3(n)=y3(n+1); end 
                        outidx3 = [outidx2 n];
                end
            end
            end
            
        elseif c == 1 % latencies
            for n =[length(y1)-2]:-1:1
                if y1(n) < y1(n+1) % outlier if value SMALLER than previous (high intnesity) value
                        try y1(n)=mean([y1(n-1) y1(n+1)]); catch y1(n)=y1(n+1); end 
                        outidx1 = [outidx1 n];
                end
            end
            for n =[length(y2)-2]:-1:1
                if y2(n) < y2(n+1) % outlier if value SMALLER than previous (high intnesity) value
                        try y2(n)=mean([y2(n-1) y2(n+1)]); catch y2(n)=y2(n+1); end 
                        outidx2 = [outidx2 n];
                end
            end
            try
            for n =[length(y3)-2]:-1:1
                if y3(n) < y3(n+1) % outlier if value SMALLER than previous (high intnesity) value
                        try y3(n)=mean([y3(n-1) y3(n+1)]); catch y3(n)=y3(n+1); end 
                        outidx3 = [outidx2 n];
                end
            end
           end
        end

        % SLOPES
        
        % only pre NOE plot
            clear slope1 slope 2 slope3 intercept1 intercept2
            if doplot == 2;
            figure(1); hold all
            p01 = plot(ints1, y1, '-Ob','color',[.5 .5 1],'linewidth',2);
            % mark outliers & replacements
            po1 = plot(ints1(outidx1), y1(outidx1), 'Og','linewidth',2)
            end

            [fit1] = polyfit(ints1,y1,1); % polyfit. linear
            slope1 = fit1(1);
            intercept1 = fit1(2);
            fy1 = [slope1*min(ints1)+intercept1 slope1*max(ints1)+intercept1];
            if doplot == 2; plot([min(ints1) max(ints1)],fy1,':','color',[.5 .5 1],'linewidth',3); end
           
            % plot settings
            if doplot == 2;
            try
                legend([p01],'Location','northwest','BL')
            end
                set(gca,'linewidth',1.5,'fontsize',14)
                xticks([40:10:90]); 
                xlabel(['Sound level (dB)'])
                ylabel(ylabels{c})
            end
           
        
            % pre NOE 
            if doplot == 2;
            figure(2); hold all
            p1 = plot(ints1, y1, '-Ob','color',[.5 .5 1],'linewidth',2);
            % mark outliers & replacements
            po1 = plot(ints1(outidx1), y1(outidx1), 'Og','linewidth',2)
            end

            [fit1] = polyfit(ints1,y1,1); % polyfit. linear
            slope1 = fit1(1);
            intercept1 = fit1(2);
            fy1 = [slope1*min(ints1)+intercept1 slope1*max(ints1)+intercept1];
            if doplot == 2; plot([min(ints1) max(ints1)],fy1,':','color',[.5 .5 1],'linewidth',3); end

            % post NOE 
%              p2 = plot(ints2, y2, '-O','color',[1 .5 .5],'linewidth',2);
            % mark outliers & replacements
            if doplot == 2; po1 = plot(ints2(outidx2), y1(outidx2), 'Og','linewidth',2); end

            [fit2] = polyfit(ints2,y2,1); % polyfit, linear
            slope2 = fit2(1);
            intercept2 = fit2(2);
            fy2 = [slope2*min(ints2)+intercept2 slope2*max(ints2)+intercept2];
           if doplot == 1; plot([min(ints2) max(ints2)],fy2,':','color',[1 .5 .5],'linewidth',3); end
           
           % plot settings
            if doplot == 2;
            try
                legend([p1 p2],'Location','northwest','BL','Post1')
            end
                set(gca,'linewidth',1.5,'fontsize',14)
                xticks([40:10:90]); 
                xlabel(['Sound level (dB)'])
                ylabel(ylabels{c})
            end
            
          

            % pre NOE (for post NOE 2 plot)
            if doplot == 2;
            figure; hold all
            p1 = plot(ints1, y1, '-Ob','color',[.5 .5 1],'linewidth',2);
            % mark outliers & replacements
            po1 = plot(ints1(outidx1), y1(outidx1), 'Og','linewidth',2)
            end
            % plot linear regression plot again
            if doplot == 2; plot([min(ints1) max(ints1)],fy1,':','color',[.5 .5 1],'linewidth',3); end

            % post NOE 2
            try
            p3 = plot(ints3, y3, '-O','color',[1 .5 .5],'linewidth',2);
            % mark outliers & replacements
            if doplot == 2; po1 = plot(ints3(outidx3), y1(outidx3), 'Og','linewidth',2); end

            [fit3] = polyfit(ints3,y3,1); % polyfit, linear
            slope3 = fit3(1);
            intercept3 = fit3(2);
            fy3 = [slope3*min(ints3)+intercept3 slope3*max(ints3)+intercept3];
            if doplot == 2; plot([min(ints3) max(ints3)],fy3,':','color',[1 .5 .5],'linewidth',3); end
            end
     
            
            % plot settings
            if doplot == 2;
                legend([p1 p3],'BL','Post2','Location','northwest')
                set(gca,'linewidth',1.5,'fontsize',14)
                xticks([40:10:90]); 
                xlabel(['Sound level (dB)'])
                ylabel(ylabels{c})
            end

            % save slope value for pre and post
            Out_slope{c}{ch}{freq}{Nmouse}{1}= slope1;
            Out_slope{c}{ch}{freq}{Nmouse}{2}= slope2;
            try
            Out_slope{c}{ch}{freq}{Nmouse}{3}= slope3;
            end
            
            % EXCLUSION criteria for modulation plots
            % slope ratio post/pre 
            if slope2/slope1 < 0 % EXCLUSION CRITERION: if modulation of either time goes in the opposite direction (i.e. larger response with lower int)
            Outsloperat{c}(Nmouse,ch)=nan;
            else
            Outsloperat{c}(Nmouse,ch)= slope2/slope1;
            end
            % slope ratio post2/pre 
            try
                if slope3/slope1 < 0 % EXCLUSION CRITERION
                Outsloperat2{c}(Nmouse,ch)=nan;
                else
                Outsloperat2{c}(Nmouse,ch)= slope3/slope1;
                end
            end
            
            
            
            
            

        % AREA BELOW PLOT
        
        
            % polygon envelopes
            % x values
            
            ints_env_frame=[90 80 70 60 50 40 30 20];
            ints_env_A = ints_env_frame(1:length(y1)); ints_env_B = ints_env_A(end);
            ints_env1 = [ints_env_A ints_env_B];
            
            ints_env_A = ints_env_frame(1:length(y2)); ints_env_B = ints_env_A(end);
            ints_env2 = [ints_env_A ints_env_B];
            
            try
            ints_env_A = ints_env_frame(1:length(y3)); ints_env_B = ints_env_A(end);
            ints_env3 = [ints_env_A ints_env_B];
            end
            
            % y values
            y1_env = [zeros(1,length(y1)) y1(1)];
            y2_env = [zeros(1,length(y2)) y2(1)];
            try
            y3_env = [zeros(1,length(y3)) y3(1)];
            end

            % plot polygons
            % pre vs post 1
            if doplot == 2;
            figure(3);hold all
%             title('pre vs post1')
            plot([ints1 ints_env1],[y1 y1_env],'O-','color',[.5 .5 1],'linewidth',2)
            plot([ints2 ints_env2],[y2 y2_env],'O-r','color',[1 .5 .5],'linewidth',2)
                % colour in the area
                a = area(ints1,y1);
                a.FaceAlpha = 0.2; a.FaceColor = [.5 .5 1];
                b = area(ints2,y2);
                b.FaceAlpha = 0.2; b.FaceColor = [1 .5 .5];
            end
            
              % plot settings
                set(gca,'linewidth',1.5,'fontsize',14)
                xticks([40:10:90]); 
                xlabel(['Sound level (dB SPL)'])
                ylabel(ylabels{c})
                legend('BL','Post1','Location','northwest')
                
               
                
                
            % pre vs post 2
            try
            if doplot == 2;
                figure;hold all
%                 title('pre vs post2')
                plot([ints1 ints_env1],[y1 y1_env],'O-','color',[.5 .5 1],'linewidth',2)
                try plot([ints3 ints_env3],[y3 y3_env],'O-r','color',[1 .5 .5],'linewidth',2); end
                    % colour in the area
                    a = area(ints1,y1);
                    a.FaceAlpha = 0.2; a.FaceColor = [.5 .5 1];
                        try
                        b = area(ints3,y3);
                        b.FaceAlpha = 0.2; b.FaceColor = [1 .5 .5];
                        end
            end
            end

            % area under plots
            Out_area{c}{1}(Nmouse,ch) = polyarea([y1 y1_env],[ints1 ints_env1]);
            Out_area{c}{2}(Nmouse,ch) = polyarea([y2 y2_env],[ints2 ints_env2]);
                Out_areaB{c}{1}{ch}{Nmouse} = polyarea([y1 y1_env],[ints1 ints_env1]); % only for plots of ind. pre post values
                Out_areaB{c}{2}{ch}{Nmouse} = polyarea([y2 y2_env],[ints2 ints_env2]); % only for plots of ind. pre post values
            try
            Out_area{c}{3}(Nmouse,ch) = polyarea([y3 y3_env],[ints3 ints_env3]);
                Out_areaB{c}{3}{ch}{Nmouse} = polyarea([y3 y3_env],[ints3 ints_env3]); % only for plots of ind. pre post values
            end

            if slope2/slope1 < 0 % EXCLUSION CRITERION - post 1
            Outarrat{c}(Nmouse,ch) = nan;
            Outardiff{c}(Nmouse,ch) = nan;
            else
%             ratio post / pre
            Outarrat{c}(Nmouse,ch) = Out_area{c}{2}(Nmouse,ch)/Out_area{c}{1}(Nmouse,ch);
%             difference post - pre
            Outardiff{c}(Nmouse,ch) = Out_area{c}{2}(Nmouse,ch)-Out_area{c}{1}(Nmouse,ch);
            end
            
            try
                if slope3/slope1 < 0 % EXCLUSION CRITERION - post 2
                Outarrat2{c}(Nmouse,ch) = nan;
                Outardiff2{c}(Nmouse,ch) = nan;
                else
    %             ratio post / pre
                Outarrat2{c}(Nmouse,ch) = Out_area{c}{3}(Nmouse,ch)/Out_area{c}{1}(Nmouse,ch);
    %             difference post - pre
                Outardiff2{c}(Nmouse,ch) = Out_area{c}{3}(Nmouse,ch)-Out_area{c}{1}(Nmouse,ch);
                end
            end
                % plot settings
                set(gca,'linewidth',1.5,'fontsize',14)
                xticks([40:10:90]); 
                xlabel(['Sound level (dB SPL)'])
                ylabel(ylabels{c})
                legend('BL','Post2','Location','northwest')


     end
end


%% Generate SPSS table

% Factors: condtion, pkltslope, pkltarrat, pklatarrdiff
%                    pkmagslope, pkmagarrat, pkmagardiff
%                    RMSslope,RMSarrat, RMSardiff




PKLTslope = [];PKLTarrat = [];PKLTardiff = [];PKmagslope = [];PKMagarrat = [];PKMAGarrdiff = [];
RMSslope = [];RMSarrat = [];RMSdiff = [];TIME = []; ID = [];
for time = 1:2
    % time=2
    if time == 1
        PKLTslope = [PKLTslope; Outsloperat{1}(:,ch)]; % {1} is PkLt
        PKLTarrat = [PKLTarrat; Outarrat{1}(:,ch)];
        PKLTardiff = [PKLTardiff; Outardiff{1}(:,ch)];

        PKmagslope = [PKmagslope; Outsloperat{2}(:,ch)];
        PKMagarrat = [PKMagarrat; Outarrat{2}(:,ch)];
        PKMAGarrdiff = [PKMAGarrdiff; Outardiff{2}(:,ch)];

        RMSslope = [RMSslope; Outsloperat{3}(:,ch)];
        RMSarrat = [RMSarrat; Outarrat{3}(:,ch)];
        RMSdiff = [RMSdiff; Outardiff{3}(:,ch)];

        ID = [ID 1:size(IDs,1)];
        TIME = [TIME; ones(size(IDs,1),1)*time];
        
    elseif time == 2
        try
        PKLTslope = [PKLTslope; Outsloperat2{1}(:,ch)];
        PKLTarrat = [PKLTarrat; Outarrat2{1}(:,ch)];
        PKLTardiff = [PKLTardiff; Outardiff2{1}(:,ch)];

        PKmagslope = [PKmagslope; Outsloperat2{2}(:,ch)];
        PKMagarrat = [PKMagarrat; Outarrat2{2}(:,ch)];
        PKMAGarrdiff = [PKMAGarrdiff; Outardiff2{2}(:,ch)];

        RMSslope = [RMSslope; Outsloperat2{3}(:,ch)];
        RMSarrat = [RMSarrat; Outarrat2{3}(:,ch)];
        RMSdiff = [RMSdiff; Outardiff2{3}(:,ch)];

        ID = [ID 1:length(Outsloperat2{1}(:,ch))]'; % IDs according to no of aniamls included, satrting from ID1
        TIME = [TIME; ones( length(Outsloperat2{1}(:,ch)) ,1)*time];
        catch
            continue
        end
    end
end

ABRdyn = [ID PKLTslope PKLTarrat PKLTardiff PKmagslope PKMagarrat PKMAGarrdiff...
RMSslope RMSarrat RMSdiff TIME];

% Construct SPSS table
ABRdyn = [ID PKLTslope PKLTarrat PKLTardiff PKmagslope PKMagarrat PKMAGarrdiff RMSslope RMSarrat RMSdiff TIME];
ABRdyn_table = array2table(ABRdyn); % convert to table
ABRdyn_table.Properties.VariableNames(1:11) = {'ID' 'PKLTslope', 'PKLTarrat', 'PKLTardiff', 'PKmagslope', 'PKMagarrat', 'PKMAGarrdiff',...
'RMSslope', 'RMSarrat', 'RMSdiff','TIME'}; % name the columns
 

%% ABR dynamics: Contruct matrix with pre and post values
% SPSS matrix for ABR modulation 

% Out_area{c}{cond}{1}(Nmouse,ch)

ID = []; TIME = [];PKLTslope = [];PKLTar = [];PKMGslope = [];PKMGar = [];
RMSslope = [];RMSar = [];

for t = [1 2 3]
    
    TIME = [TIME; ones( size(Out_areaB{1}{t}{ch},2) ,1)*t];
    ID = [ID; [1:size(Out_areaB{1}{t}{ch},2)]' ];
    
    for anm = 1:length(IDs)
        try
        PKLTslope = [PKLTslope; Out_slope{1}{ch}{freq}{anm}{t}];
        PKLTar = [PKLTar; Out_areaB{1}{t}{ch}{anm}];

        PKMGslope = [PKMGslope; Out_slope{2}{ch}{freq}{anm}{t}];
        PKMGar = [PKMGar; Out_areaB{2}{t}{ch}{anm}];

        RMSslope = [RMSslope; Out_slope{3}{ch}{freq}{anm}{t}];
        RMSar = [RMSar; Out_areaB{3}{t}{ch}{anm}];
        end
    end
end

% for t = 1:2
%     TIME = [TIME; ones(size(Out_area{c}{t}(:,ch),1),1)*t];
%     ID = [ID; [1:size(Out_area{c}{t}(:,ch),1)]' ];
% 
%     PKLTslope = [PKLTslope; Out_slope{1}{ch}{freq}(:,t)];
%     PKLTar = [PKLTar; Out_area{1}{t}(:,ch)];
%     
%     PKMGslope = [PKMGslope; Out_slope{2}{ch}{freq}(:,t)];
%     PKMGar = [PKMGar; Out_area{2}{t}(:,ch)];
%     
%     RMSslope = [RMSslope; Out_slope{3}{ch}{freq}(:,t)];
%     RMSar = [RMSar; Out_area{3}{t}(:,ch)];
% 
% end

% Construct SPSS table
ABRdynpp = [ID TIME PKLTslope PKLTar PKMGslope PKMGar  RMSslope RMSar];
ABRdynpp_table = array2table(ABRdynpp); % convert to table
ABRdynpp_table.Properties.VariableNames(1:8) = {'ID','Time','PKLTslope', 'PKLTar', 'PKmagslope', 'PKMagar',...
'RMSslope', 'RMSar'}; % name the columns





%% Tinnitus index

clear IDX1c IDX2c IDX3c
for anm = 1:length(IDs)
IDX1c{anm} = find(ABRdynpp(:,2)==1 & ABRdynpp(:,1)==anm ); % time 1 
IDX2c{anm} = find(ABRdynpp(:,2)==2 & ABRdynpp(:,1)==anm ); % time 2 
IDX3c{anm} = find(ABRdynpp(:,2)==3 & ABRdynpp(:,1)==anm ); % time 3
end

% Use change in RMS as input for tinnitus index

labellist = {'ID','Time','Peak latency (slope)', 'Peak latency (area)', 'Peak magnitude (slope)', 'Peak magnitude (area)',...
'RMS (slope)', 'RMS (area)'};

clear RMSar1 RMSar2 RMSar3 RMSrat1 RMSrat2
for anm = 1:3
    for m = 8
        RMSar1(anm) = ABRdynpp(IDX1c{anm},m);
        RMSar2(anm) = ABRdynpp(IDX2c{anm},m);
        RMSar3(anm) = ABRdynpp(IDX3c{anm},m);
    end
    % RMS change (ratio)
    RMSrat1(anm) = RMSar2(anm)/RMSar1(anm);
    RMSrat2(anm) = RMSar3(anm)/RMSar1(anm);

end

% save variables to file
fname = ['TI_RMSmag_' num2str(freq) 'K_Aug22'];
save(fname,'RMSrat1','RMSrat2')

