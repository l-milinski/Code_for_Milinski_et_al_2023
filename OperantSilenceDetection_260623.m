%% Operant silence detection, 

% Loads and processes data to produce figures Fig. 1E and S2 C,D
% Produces metrics for tinnitus index.



clear all 
close all
conds = 1:3; % 1:3

dosave = 0; 
doprint = 0; % for printing spss tables (.xlcs files)

% animals = {'1901','1902','1903','1904'};
animals = {'1709','1710','1711','1712'};
% animals = {'1710'};

animals = {'1709','1710','1711','1712','1901','1902','1903','1904'};
stims = 1:3; % 1:3 for all 3 stimulus types AM, NBN, silence 

% Make sure data paths are correct for pre NOE data ('prepath'), post NOE
% ('postpath') and post NOE 2 ('postpath2')
prepath = 'D:\PC-VV001-restoration\manuscript\Repository\Operant_silence_data\pre\'; 
postpath = 'D:\PC-VV001-restoration\manuscript\Repository\Operant_silence_data\post\';
postpath2 = 'D:\PC-VV001-restoration\manuscript\Repository\Operant_silence_data\post2\';

% savepath = 'D:\PC-VV001-restoration\Ferrets\Behaviour\Figures\Time123\Silence\';

nancount=0;
Exccount = 0;
% load files
for anm = 1:length(animals)
    for time = conds % pre & post
    
    for sess = 1:70
        loadcheck = 1
        
        % clear metadata
        clear responses nCorr nAll s2
    
        % load session 
          if time == 2 
              fn = [postpath cell2mat(animals(anm)) 'ttest90B0' num2str(sess)]; % name variant 1 19XX ferrets
              try load(fn); catch; loadcheck = 0 ; end
              if loadcheck == 0
                  fn = [postpath cell2mat(animals(anm)) 'ttest70B0' num2str(sess)]; % name variant 2 19XX ferret
                  try load(fn); catch; loadcheck = 0 ; end
              end
              if loadcheck == 0
                  fn = [postpath cell2mat(animals(anm)) 'ttest70B' num2str(sess)]; % name variant 2 19XX ferret
                  try load(fn); catch; loadcheck = 0 ; end
              end
              if loadcheck == 0
                  fn = [postpath cell2mat(animals(anm)) 'ttest90B' num2str(sess)]; % name variant 2 19XX ferret
                  try load(fn); catch; loadcheck = 0 ; end
              end
              if loadcheck == 0
                  fn = [postpath cell2mat(animals(anm)) 'TTs90E0' num2str(sess)];% name variant 2 17XX ferret
                  try load(fn); catch; loadcheck = 0 ; end
              end
              if loadcheck == 0
                  fn = [postpath cell2mat(animals(anm)) 'TTs90E' num2str(sess)];% name variant 2 17XX ferret
                  try load(fn); catch; loadcheck = 0 ; end
              end
          
          elseif time == 1; fn = [prepath cell2mat(animals(anm)) 'ttest70A' num2str(sess)];
            try load(fn); catch;fn = [prepath cell2mat(animals(anm)) 'ttest70A0' num2str(sess)]; %try with addtional '0' in name
               try load(fn);catch; loadcheck = 0 ; end
            end
            % if still no file found, try filenames of 17XX animals
            if loadcheck == 0; fn = [prepath cell2mat(animals(anm)) 'TTs90D' num2str(sess)];
               try load(fn); catch;fn = [prepath cell2mat(animals(anm)) 'TTs90D0' num2str(sess)]; %try with addtional '0' in name
                  try load(fn);catch; loadcheck = 0 ; end
               end
            end
            
          elseif time == 3 
                      fn = [postpath2 cell2mat(animals(anm)) 'ttest70C' num2str(sess)]; % name variant 1 19XX ferrets
                      try load(fn); catch; loadcheck = 0 ; end
          
                      if loadcheck == 0
                      fn = [postpath2 cell2mat(animals(anm)) 'TTEST70C' num2str(sess)]; % name variant 2 19XX ferret
                      try load(fn); catch; loadcheck = 0 ; end
               end
          end
          
     % stop when session can't be loaded or when less than 5 trials
     if exist('s2') == 0; continue; end
     if  length(s2.StimID) < 5; continue; end

% --- for each session: calculate hit rate for BBN, NBN, Silence trials
        responses(1,:) = s2.StimID;
        responses(2,:) = s2.TargetLoc;
        responses(3,:) = s2.RespLoc;
               
        responses(5,:) = s2.RespTime;


        for n=1:length(responses(1,:))
            % if TargetLOc matches RespLoc (i.e. correct response)
            if responses(2,n) == responses(3,n)
            % ... put 1
                responses(4,n) = 1;
                responses(6,n) = responses(5,n); % add correct RT
                responses(7,n) = nan; % add NaN for incorrect RT

            else responses(4,n) = 0;
                 responses(7,n) = responses(5,n); % add incorrect RT
                 responses(6,n) = nan; % add NaN for correct RT

            end
        end
        
        RTfilt = 20; % exclusion criterium in s
        % --- exclude trials with long response times
        idx = find(responses(5,:)> RTfilt); %delays all trials
        % if delay to large, exclude individual trials
        responses(4,idx) = nan;
        responses(6,idx) = nan;
        responses(7,idx) = nan;
        Exccount = Exccount+length(idx);
        % --- --- --- --- --- ---
        
        % find number correct for each Stim 
        for stim =  stims
        nCorr(stim) = sum(responses(4,(find(responses(1,:)==stim))));
        nAll(stim) = length(responses(4,(find(responses(1,:)==stim))));
        end
        
        % find response time for each Stim 
        for stim =  stims
        ResT{anm}{time}(sess,stim) = nanmean(responses(5,(find(responses(1,:)==stim))));
        cResT{anm}{time}(sess,stim) = nanmean(responses(6,(find(responses(1,:)==stim))));
            try
            iResT{anm}{time}(sess,stim) = nanmean(responses(7,(find(responses(1,:)==stim))));
            catch; iResT{anm}{time}(sess,stim) = nan;
            end
        end

        % calculate hitrate % response times
        for stim =  stims
        HR{anm}{time}(sess,stim) = nCorr(stim) / nAll(stim);
        end
        
    end
    
    for sess = 1:70
        
        % replace empty sessions with NaNs (to not bias averages)
        try
        if HR{anm}{time}(sess,:) == [0 0 0]; HR{anm}{time}(sess,:) = [nan nan nan]; 
            nancount = nancount+1;
        end 
        end
    end
    
    
        try
        meanHR{anm}{time} = nanmean( HR{anm}{time} );
        stdHR{anm}{time} = nanstd( HR{anm}{time} );
        err{anm}{time} = nanstd( HR{anm}{time} )/sqrt(length( HR{anm}{time} ));
        catch
                    meanHR{anm}{time} = [nan nan nan];
        end
    end
end

 

% SPSS matrix A

clear Silence
ID = []; HTR = []; TIME = []; SESSION = []; FREQ = []; RT = []; cRT = []; iRT = [];
for id = 1:length(animals)
    for time = conds
        for freq = stims
            try HR{id}{time}(:,freq); %in case there are no data for one a animal (e.g. 1901 post)
                for sess = 1:length( HR{id}{time}(:,freq) )   
                    ID = [ID; id];
                    TIME = [TIME; time];
                    FREQ = [FREQ; freq];
                    SESSION = [SESSION; sess];
                    
                    HTR = [HTR; HR{id}{time}(sess,freq)];
                    RT = [RT; ResT{id}{time}(sess,freq)];
                    cRT = [cRT; cResT{id}{time}(sess,freq)];
                    iRT = [iRT; iResT{id}{time}(sess,freq)];
                    
%                     if HR{id}{time}(sess,freq) ~= 1 && HR{id}{time}(sess,freq) == HTR(end-1); keyboard; end
                end
            catch
            end
        end
    end
end

Silence = [ID HTR TIME SESSION FREQ RT cRT iRT];

  
%   exlcude NaNs
    inIDX = find(~isnan( Silence(:,2) ) ); %find lines with NaNs
    excidx = find(isnan( Silence(:,2) ) ); 
    Silence_nanex = Silence(inIDX,:); % construct final array without NaNs
    
Silence_table = array2table(Silence_nanex); % convert to table
Silence_table.Properties.VariableNames(1:8) = {'ID','HTR','TIME','SESSION','FREQ','RT','cRT','iRT'}; % name the columns



% SPSS matrix B (animal means separated by freq)


clear freq_anm_idx freq_idx HRanmlmean RTianmmean RTcanmmean RTanmmean
for t = 1:3
    for freq = 1:3
        freq_idx{t}{freq} = find(Silence_nanex(:,5) == freq & Silence_nanex(:,3) == t); % idx by freq
        for id = 1:length(animals)
            freq_anm_idx{t}{freq}{id} = find(Silence_nanex(:,5) == freq & Silence_nanex(:,1) == id & Silence_nanex(:,3) == t); % idx by freq & anm
        
        HRanmlmean{t}{freq}{id} = nanmean( Silence_nanex(freq_anm_idx{t}{freq}{id},2) ) ;
        RTianmlmean{t}{freq}{id} = nanmean( Silence_nanex(freq_anm_idx{t}{freq}{id},8) ) ;
        RTcanmlmean{t}{freq}{id} = nanmean( Silence_nanex(freq_anm_idx{t}{freq}{id},7) ) ;
        RTanmlmean{t}{freq}{id} = nanmean( Silence_nanex(freq_anm_idx{t}{freq}{id},6) ) ;


        end
    end
end

ID = []; TIME = []; FREQ = []; HR = []; RTc = []; RTi = [];
                    

for id = 1:8
    for t = 1:3
        for freq = 1:3
        ID = [ID;id];
        TIME = [TIME; t];
        FREQ = [FREQ; freq];
        
        HR = [HR;HRanmlmean{t}{freq}{id}];
        RTc = [RTc;RTcanmlmean{t}{freq}{id}];
        RTi = [RTi;RTianmlmean{t}{freq}{id}];
        
        end
    end
end

SPSSanmbfreqs = [ID TIME FREQ HR RTc RTi ];

SPSSanmbfreqs_table_CompAcrossTime = array2table(SPSSanmbfreqs); % convert to table
SPSSanmbfreqs_table_CompAcrossTime.Properties.VariableNames(1:6) = {'ID','TIME','FREQ','HR','cRT','iRT'}; % name the columns
 


% SPSS matrix C
% silence only

clear freq_anm_idx HRanmlmean RTianmmean RTcanmmean RTanmmean
for t = 1:3
    for freq = 1:3
        freq_idx{t}{freq} = find(Silence_nanex(:,5) == freq & Silence_nanex(:,3) == t); % idx by freq
        for id = 1:length(animals)
            freq_anm_idx{t}{freq}{id} = find(Silence_nanex(:,5) == freq & Silence_nanex(:,1) == id & Silence_nanex(:,3) == t); % idx by freq & anm
        
        HRanmlmean{t}{freq}{id} = nanmean( Silence_nanex(freq_anm_idx{t}{freq}{id},2) ) ;
        RTianmmean{t}{freq}{id} = nanmean( Silence_nanex(freq_anm_idx{t}{freq}{id},8) ) ;
        RTcanmmean{t}{freq}{id} = nanmean( Silence_nanex(freq_anm_idx{t}{freq}{id},7) ) ;
        RTanmmean{t}{freq}{id} = nanmean( Silence_nanex(freq_anm_idx{t}{freq}{id},6) ) ;


        end
    end
end

ID = []; TIME = []; AM_HR = []; NBN_HR = []; SIL_HR = [];...
                    AM_RTc = []; NBN_RTc = []; SIL_RTc = [];...
                    AM_RTi = []; NBN_RTi = []; SIL_RTi = [];

for id = 1:8
    for t = 1:3
        ID = [ID;id];
        TIME = [TIME; t];
        
        AM_HR = [AM_HR;HRanmlmean{t}{1}{id}];
        NBN_HR = [NBN_HR;HRanmlmean{t}{2}{id}];
        SIL_HR = [SIL_HR;HRanmlmean{t}{3}{id}];

        
        AM_RTc = [AM_RTc;RTcanmmean{t}{1}{id}];
        NBN_RTc = [NBN_RTc;RTcanmmean{t}{2}{id}];
        SIL_RTc = [SIL_RTc;RTcanmmean{t}{3}{id}];

        
        AM_RTi = [AM_RTi;RTianmmean{t}{1}{id}];
        NBN_RTi = [NBN_RTi;RTianmmean{t}{2}{id}];
        SIL_RTi = [SIL_RTi;RTianmmean{t}{3}{id}];

        
    end
end

SPSS_acrossT = [ID TIME AM_HR NBN_HR SIL_HR...
                AM_RTc NBN_RTc SIL_RTc...
                AM_RTi NBN_RTi SIL_RTi ];

Silence_table_CompAcrossTime = array2table(SPSS_acrossT); % convert to table
Silence_table_CompAcrossTime.Properties.VariableNames(1:11) = {'ID','TIME','AM_HR','NBN_HR','SIL_HR','AM_RTc','NBN_RTc','SIL_RTc','AM_RTi','NBN_RTi','SIL_RTi'}; % name the columns
 


%% Tinnitus index: change in silence detection per animal

for t = 1:3
        for anm = 6:8
        SdPerfIDX{t}(anm) = find(SPSS_acrossT(:,2) == t & SPSS_acrossT(:,1) == anm); % idx for ID and TIME
        SdPerf{t}(anm) = SPSS_acrossT(SdPerfIDX{t}(anm),5); % SIL_HR per animal and time
        end
end

% calculate change in silence detection (ratio)
for anm = 6:8
    Sdratio1(anm) = SdPerf{2}(anm)/SdPerf{1}(anm);
    Sdratio2(anm) = SdPerf{3}(anm)/SdPerf{1}(anm);
end

% Silence detection input for tinnitus index
Sdratio1(anm)
Sdratio2(anm)

% A decrease in Silence HR is avidence for tinnitus. Therefore, for the TI
% input, calculat3e 1-SDratio, so a decrease results in a positive value.
for anm = 6:8
    TI1(anm) = 1-Sdratio1(anm);
    TI2(anm) = 1-Sdratio2(anm);
end

% save variables to file
fname = ['TI_silence'];
% save(fname,'TI1','TI2')



%% Figure 1E: plot hit rate by stimulus (CI and errobars)
% different order of plotting: NBN, AM, Silence
% close all

clear freq_idx freq_anm_idx HTR grandmeanHR

conds = 1:3
plotci = 2; % 1 for CIs, 2 for errobars, 3 for lines (animal averages)


HTR = Silence_nanex(:,2);

for t = conds
    for freq = 1:3
        freq_idx{t}{freq} = find(Silence_nanex(:,5) == freq & Silence_nanex(:,3) == t); % idx by freq
        for id = 1:length(animals)
            freq_anm_idx{t}{id}{freq} = find(Silence_nanex(:,5) == freq & Silence_nanex(:,1) == id & Silence_nanex(:,3) == t); % idx by freq & anm
        end
    end
end

figure;hold all
Markers = {'+','o','*','x','v','d','^','s','>','<'};
cmap = [.5 .5 1; 1 .5 .5; 0.9 0.6 0.3];
ts = [-0.4 0 0.4]; % time shift on x axis
xpos = [3 1 5]; % to plot in order NBN, AM, Silence
winput = 0.4; % width for box plots

% adapt box plot settings for only BL plots
if conds == 1
winput = 0.6; % width for box plots
ts = [0 0 0]; % time shift on x axis
cmap = [0 0 1; 1 .5 .5; 0.9 0.6 0.3];
end

clear CI95 yCI95 meanb stdb bSEM
for t = conds
    for freq = [1 2 3]
        % HTR: sessions (sorted by t and freq)
        sessHR{t}{freq} = Silence_nanex(freq_idx{t}{freq},2);
        % HTR: session means (by t and freq)
        sessmeanHR{t}{freq} = nanmean(sessHR{t}{freq});
        
        % CI
        % for errorbar plots
               meanb{t}{freq} = sessmeanHR{t}{freq};
               stdb{t}{freq} = nanstd(sessHR{t}{freq});
               bSEM{t}{freq} = nanstd(sessHR{t}{freq})/sqrt( length(sessHR{t}{freq}) );  
               CI95{t}{freq} = tinv([0.025 0.975], length(sessHR{t}{freq})-1);                    % Calculate 95% Probability Intervals Of t-Distribution
               yCI95{t}{freq} = bsxfun(@times, bSEM{t}{freq}, CI95{t}{freq});              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
            
        


        % save measn across animals
        grandmeanHR{t}(1,freq) = nanmean( [mean(HTR(freq_anm_idx{t}{1}{freq}));mean(HTR(freq_anm_idx{t}{2}{freq}));mean(HTR(freq_anm_idx{t}{3}{freq}));...
                         mean(HTR(freq_anm_idx{t}{4}{freq}));mean(HTR(freq_anm_idx{t}{5}{freq}));mean(HTR(freq_anm_idx{t}{6}{freq}));...
                         mean(HTR(freq_anm_idx{t}{7}{freq}));mean(HTR(freq_anm_idx{t}{8}{freq}))] )
        grandmeanHR{t}(2,freq) = nanstd( [mean(HTR(freq_anm_idx{t}{1}{freq}));mean(HTR(freq_anm_idx{t}{2}{freq}));mean(HTR(freq_anm_idx{t}{3}{freq}));...
                         mean(HTR(freq_anm_idx{t}{4}{freq}));mean(HTR(freq_anm_idx{t}{5}{freq}));mean(HTR(freq_anm_idx{t}{6}{freq}));...
                         mean(HTR(freq_anm_idx{t}{7}{freq}));mean(HTR(freq_anm_idx{t}{8}{freq}))] )
    
   
        Markers = {'+','o','*','x','v','d','^','s','>','<'};
        if plotci == 1
            
            for id = 1:7
            %animal means (markers)
             e(t) =  plot(xpos(freq)+ts(t),mean( HTR(freq_anm_idx{t}{id}{freq}) ),strcat(Markers{id}),'color',cmap(t,:)...
                ,'Markersize',7,'linewidth',1.5)
            end
            
            % plot CIs
            errorbar([xpos(freq)+ts(t)],meanb{t}{freq},yCI95{t}{freq}(2),'linewidth',2,'CapSize',25,'color','k') % Plot 95% Confidence as errorbars
            
             % session means
%                    plot(xpos(freq)+ts(t),sessmeanHR{t}{freq},'O','color',cmap(t,:),'Markerfacecolor',cmap(t,:)...
%                             ,'Markersize',10,'linewidth',2)
           

        elseif plotci == 2
            % plot errorbars or box plots
            
            for id = 1:7
            %animal means (markers)
             e(t) =  plot(xpos(freq)+ts(t),mean( HTR(freq_anm_idx{t}{id}{freq}) ),strcat(Markers{id}),'color',cmap(t,:)...
                ,'Markersize',7,'linewidth',1.5)
            end
            
            % plot box plots
            if sum(conds) == 1
             b = boxplot([mean(HTR(freq_anm_idx{t}{1}{freq}));mean(HTR(freq_anm_idx{t}{2}{freq}));mean(HTR(freq_anm_idx{t}{3}{freq}));...
                         mean(HTR(freq_anm_idx{t}{4}{freq}));mean(HTR(freq_anm_idx{t}{5}{freq}));mean(HTR(freq_anm_idx{t}{6}{freq}));...
                         mean(HTR(freq_anm_idx{t}{7}{freq}));mean(HTR(freq_anm_idx{t}{8}{freq}))]...
                         ,'positions',xpos(freq)+ts(t),'Colors','k','Widths',0.8,'Whisker',0,'Symbol','');
            else
            b = boxplot([mean(HTR(freq_anm_idx{t}{1}{freq}));mean(HTR(freq_anm_idx{t}{2}{freq}));mean(HTR(freq_anm_idx{t}{3}{freq}));...
                         mean(HTR(freq_anm_idx{t}{4}{freq}));mean(HTR(freq_anm_idx{t}{5}{freq}));mean(HTR(freq_anm_idx{t}{6}{freq}));...
                         mean(HTR(freq_anm_idx{t}{7}{freq}));mean(HTR(freq_anm_idx{t}{8}{freq}))]...
                         ,'positions',xpos(freq)+ts(t),'Colors','k','Widths',0.4,'Whisker',0,'Symbol','');
            end 
                     
             set(b,{'linew'},{1.6})

            
            
%                    e(t) =  errorbar([xpos(freq)+ts(t)],sessmeanHR{t}{freq},stdb{t}{freq},'linewidth',2,'CapSize',25,'color',cmap(t,:)) % Plot 95% Confidence as errorbars
%                     plot(xpos(freq)+ts(t),sessmeanHR{t}{freq},'O','color',cmap(t,:),'Markerfacecolor',cmap(t,:)...
%                             ,'Markersize',10,'linewidth',2) 
            
        elseif plotci == 3
            % plot individual animal lines
                plot(xpos,[mean( HTR(freq_anm_idx{t}{id}{1}) ) mean( HTR(freq_anm_idx{t}{id}{2}) ) mean( HTR(freq_anm_idx{t}{id}{3}) )],...
                           '-','color',[.5 .5 .5], 'linewidth',1,'markersize',8)
            
            % plot means across animals
                plot([xpos]+ts(t),[grandmeanHR{t}(1,1) grandmeanHR{t}(1,2) grandmeanHR{t}(1,3)] ,'-','Color',cmap(t,:), 'linewidth',3,'markersize',4,'markerfacecolor','k');
        end
    end
end

if conds == 1
% significance markers
plot([xpos(1) xpos(3)], [1.05 1.05],'-k','linewidth',1) % line
plot(mean([xpos(1) xpos(3)])-0.1, 1.065, '*k','Markersize',8,'linewidth',1.1) % marker1
plot(mean([xpos(1) xpos(3)]), 1.065, '*k','Markersize',8,'linewidth',1.1) % marker2
plot(mean([xpos(1) xpos(3)])+0.1, 1.065, '*k','Markersize',8,'linewidth',1.1) % marker3

% plot(mean([xpos(1) xpos(3)])+0.2, 1.058, '*k','Markersize',8,'linewidth',1.1) % marker2

plot([xpos(2) xpos(3)], [1.02 1.02],'-k','linewidth',1) % line
plot(mean([xpos(2) xpos(3)])-0.1, 1.035, '*k','Markersize',8,'linewidth',1.1) % marker
plot(mean([xpos(2) xpos(3)])+0.1, 1.035, '*k','Markersize',8,'linewidth',1.1) % marker 2

% plot(mean([xpos(2) xpos(3)])+0.1, 1.035, '*k','Markersize',8,'linewidth',1.1) % marker2
% plot(mean([xpos(2) xpos(3)])+0.2, 1.038, '*k','Markersize',8,'linewidth',1.1) % marker2
elseif conds == 1:3
% AM, ns, p=0.03
plot([xpos(1)+ts(1) xpos(1)+ts(3)], [1.05 1.05],'-k','linewidth',1) % line
plot(mean([xpos(1)+ts(1) xpos(1)+ts(3)]), 1.065, '*k','Markersize',8,'linewidth',1.1) % marker

% NBN, ns, ns
% S, ns
end

xlim([0 6]);ylim([0.5 1.1])
set(gca,'linewidth',1.5,'fontsize',14)
xticks([1 3 5])
xticklabels({'NBN','AM','Silence'})
xlabel('Stimulus type')

yticks([0.75 0.8 0.9 1])
ylabel('Proportion of correct responses')
box off; grid off
try
legend([e(1) e(2) e(3)],'BL','Post 1','Post 2','location','southwest')
end
legend off




%% Figure S2 C,D: plot RTs by stimulus (CI and errobars)
   % Set RTid = 2 to plot response times for correct trials, and RTid = 3
   % for response times for incorrect trials.


% close all
clear freq_idx freq_anm_idx RT grandmeanRT


% --- define input RT
RTid = 2; % 1 for RT, 2 for cRT, 3 for iRT
% --- --- ---
plotci = 2; % 1 for CIs, 2 for errobars, 3 for lines (animal averages)


% Input
input = [Silence_nanex(:,6) Silence_nanex(:,7) Silence_nanex(:,8)];
RT = input(:,RTid);

conds = 1:3;


for t = conds
    for freq = 1:3
        freq_idx{t}{freq} = find(Silence_nanex(:,5) == freq & Silence_nanex(:,3) == t); % idx by freq
        for id = 1:length(animals)
            freq_anm_idx{t}{id}{freq} = find(Silence_nanex(:,5) == freq & Silence_nanex(:,1) == id & Silence_nanex(:,3) == t); % idx by freq & anm
        end
    end
end

figure;hold all
Markers = {'+','o','*','x','v','d','^','s','>','<'};
cmap = [.5 .5 1; 1 .5 .5; 0.9 0.6 0.3];
ts = [-0.4 0 0.4]; % time shift on x axis
xpos = [1 3 5];
winput = 0.4; % width for box plots

% adapt box plot settings for only BL plots
if conds == 1
winput = 0.6; % width for box plots
ts = [0 0 0]; % time shift on x axis
cmap = [0 0 1; 1 .5 .5; 0.9 0.6 0.3];
end

clear CI95 yCI95 meanb stdb bSEM
for t = conds
    for freq = 1:3
        % RT: sessions (sorted by t and freq)
        sessRT{t}{freq} = RT(freq_idx{t}{freq});
        % RT: session means (by t and freq)
        sessmeanRT{t}{freq} = nanmean(sessRT{t}{freq});
        
        % CI
        % for errorbar plots
               meanb{t}{freq} = sessmeanRT{t}{freq};
               stdb{t}{freq} = nanstd(sessRT{t}{freq});
               bSEM{t}{freq} = nanstd(sessRT{t}{freq})/sqrt( length(sessRT{t}{freq}) );  
               CI95{t}{freq} = tinv([0.025 0.975], length(sessRT{t}{freq})-1);                    % Calculate 95% Probability Intervals Of t-Distribution
               yCI95{t}{freq} = bsxfun(@times, bSEM{t}{freq}, CI95{t}{freq});              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
            
        
        % save means across animals
        grandmeanRT{t}(1,freq) = nanmean( [mean(RT(freq_anm_idx{t}{1}{freq}));mean(RT(freq_anm_idx{t}{2}{freq}));mean(RT(freq_anm_idx{t}{3}{freq}));...
                         mean(RT(freq_anm_idx{t}{4}{freq}));mean(RT(freq_anm_idx{t}{5}{freq}));mean(RT(freq_anm_idx{t}{6}{freq}));...
                         mean(RT(freq_anm_idx{t}{7}{freq}));mean(RT(freq_anm_idx{t}{8}{freq}))] );
        grandmeanRT{t}(2,freq) = nanstd( [mean(RT(freq_anm_idx{t}{1}{freq}));mean(RT(freq_anm_idx{t}{2}{freq}));mean(RT(freq_anm_idx{t}{3}{freq}));...
                         mean(RT(freq_anm_idx{t}{4}{freq}));mean(RT(freq_anm_idx{t}{5}{freq}));mean(RT(freq_anm_idx{t}{6}{freq}));...
                         mean(RT(freq_anm_idx{t}{7}{freq}));mean(RT(freq_anm_idx{t}{8}{freq}))] );
    
   
        Markers = {'+','o','*','x','v','d','^','s','>','<'};
        if plotci == 1
            
            for id = 1:7
            %animal means (markers)
            plot(xpos(freq)+ts(t),nanmean( RT(freq_anm_idx{t}{id}{freq}) ),strcat(Markers{id}),'color',cmap(t,:)...
                ,'Markersize',8,'linewidth',1.5);
            end
            
            % plot CIs
            if t == 1
             e(t) =  errorbar([xpos(freq)+ts(t)],meanb{t}{freq},yCI95{t}{freq}(2),'linewidth',2,'CapSize',35,'color','k') % Plot 95% Confidence as errorbars
            else
             e(t) =  errorbar([xpos(freq)+ts(t)],meanb{t}{freq},yCI95{t}{freq}(2),'linewidth',2,'CapSize',25,'color','k') % Plot 95% Confidence as errorbars
            set(e(t-1),'CapSize',25)
            end
            
             % session means
%                    plot(xpos(freq)+ts(t),sessmeanHR{t}{freq},'O','color',cmap(t,:),'Markerfacecolor',cmap(t,:)...
%                             ,'Markersize',10,'linewidth',2)
           

        elseif plotci == 2
            
            
            % plot errorbars or box plots
            
            for id = 1:7
            %animal means (markers)
             e(t) =  plot(xpos(freq)+ts(t),nanmean( RT(freq_anm_idx{t}{id}{freq}) ),strcat(Markers{id}),'color',cmap(t,:)...
                ,'Markersize',7,'linewidth',1.5)
            end

            % plot box plots
            if sum(conds) == 1
            b = boxplot([nanmean(RT(freq_anm_idx{t}{1}{freq}));nanmean(RT(freq_anm_idx{t}{2}{freq}));nanmean(RT(freq_anm_idx{t}{3}{freq}));...
                         nanmean(RT(freq_anm_idx{t}{4}{freq}));nanmean(RT(freq_anm_idx{t}{5}{freq}));nanmean(RT(freq_anm_idx{t}{6}{freq}));...
                         nanmean(RT(freq_anm_idx{t}{7}{freq}));nanmean(RT(freq_anm_idx{t}{8}{freq}))]...
                         ,'positions',xpos(freq)+ts(t),'Colors','k','Widths',0.8,'Whisker',0,'Symbol','');
            else
            b = boxplot([nanmean(RT(freq_anm_idx{t}{1}{freq}));nanmean(RT(freq_anm_idx{t}{2}{freq}));nanmean(RT(freq_anm_idx{t}{3}{freq}));...
                         nanmean(RT(freq_anm_idx{t}{4}{freq}));nanmean(RT(freq_anm_idx{t}{5}{freq}));nanmean(RT(freq_anm_idx{t}{6}{freq}));...
                         nanmean(RT(freq_anm_idx{t}{7}{freq}));nanmean(RT(freq_anm_idx{t}{8}{freq}))]...
                         ,'positions',xpos(freq)+ts(t),'Colors','k','Widths',0.4,'Whisker',0,'Symbol','');
            end
             set(b,{'linew'},{1.6})
            
            
            % plot errorbars
%             plot(xpos(freq)+ts(t),sessmeanRT{t}{freq},'O','color',cmap(t,:),'Markerfacecolor',cmap(t,:)...
%                             ,'Markersize',10,'linewidth',2) 
%             if t == 1
%             e(t) =  errorbar([xpos(freq)+ts(t)],sessmeanRT{t}{freq},stdb{t}{freq},'linewidth',2,'CapSize',35,'color',cmap(t,:)) % Plot 95% Confidence as errorbars
%             else
%             e(t) =  errorbar([xpos(freq)+ts(t)],sessmeanRT{t}{freq},stdb{t}{freq},'linewidth',2,'CapSize',25,'color',cmap(t,:)) % Plot 95% Confidence as errorbars
%             set(e(t-1),'CapSize',25)
%             end
            
            
            
            
        elseif plotci == 3
            % plot individual animal lines
                plot(xpos,[mean( HTR(freq_anm_idx{t}{id}{1}) ) mean( HTR(freq_anm_idx{t}{id}{2}) ) mean( HTR(freq_anm_idx{t}{id}{3}) )],...
                           '-','color',[.5 .5 .5], 'linewidth',1,'markersize',8)
            
            % plot means across animals
                plot([xpos]+ts(t),[grandmeanRT{t}(1,1) grandmeanRT{t}(1,2) grandmeanRT{t}(1,3)] ,'-','Color',cmap(t,:), 'linewidth',3,'markersize',4,'markerfacecolor','k');
        end
    end
end

xlim([0 6]);ylim([0 7])
set(gca,'linewidth',1.5,'fontsize',14)
xticks([xpos])
xticklabels({'AM','NBN','Silence'})
yticks([1 2 3 4 5 6])
if RTid == 1; ylabel('Response time (s)'); elseif RTid == 2; ylabel(['Response time (s)';' Correct trials  ']);...
                                       elseif RTid == 3; ylabel(['Response time (s)';'Incorrect trials ']); end
box off; grid off

if conds == 1
ylim([0 8])
    if RTid == 2
        % significance markers
        plot([xpos(1) xpos(3)], [4.4 4.4],'-k','linewidth',1) % line
        plot(mean([xpos(1) xpos(3)])-0.1, 4.55, '*k','Markersize',8,'linewidth',1.1) % marker1
        plot(mean([xpos(1) xpos(3)]), 4.55, '*k','Markersize',8,'linewidth',1.1) % marker2
        plot(mean([xpos(1) xpos(3)])+0.1, 4.55, '*k','Markersize',8,'linewidth',1.1) % marker3

        % plot(mean([xpos(1) xpos(3)])+0.2, 1.058, '*k','Markersize',8,'linewidth',1.1) % marker2

        plot([xpos(2) xpos(3)], [4 4 ],'-k','linewidth',1) % line
        plot(mean([xpos(2) xpos(3)])-0.1, 4.15, '*k','Markersize',8,'linewidth',1.1) % marker
        plot(mean([xpos(2) xpos(3)]), 4.15, '*k','Markersize',8,'linewidth',1.1) % marker2
        plot(mean([xpos(2) xpos(3)])+0.1, 4.15, '*k','Markersize',8,'linewidth',1.1) % marker3

        % plot(mean([xpos(2) xpos(3)])+0.2, 1.038, '*k','Markersize',8,'linewidth',1.1) % marker2
    elseif RTid == 3
            % significance markers
        plot([xpos(1) xpos(3)], [7.4 7.4],'-k','linewidth',1) % line
        plot(mean([xpos(1) xpos(3)])-0.1, 7.55, '*k','Markersize',8,'linewidth',1.1) % marker1
        plot(mean([xpos(1) xpos(3)]), 7.55, '*k','Markersize',8,'linewidth',1.1) % marker2
        plot(mean([xpos(1) xpos(3)])+0.1, 7.55, '*k','Markersize',8,'linewidth',1.1) % marker3

        % plot(mean([xpos(1) xpos(3)])+0.2, 1.058, '*k','Markersize',8,'linewidth',1.1) % marker2

        plot([xpos(2) xpos(3)], [7 7 ],'-k','linewidth',1) % line
        plot(mean([xpos(2) xpos(3)])-0.1, 7.15, '*k','Markersize',8,'linewidth',1.1) % marker
        plot(mean([xpos(2) xpos(3)]), 7.15, '*k','Markersize',8,'linewidth',1.1) % marker2
        plot(mean([xpos(2) xpos(3)])+0.1, 7.15, '*k','Markersize',8,'linewidth',1.1) % marker3

        % plot(mean([xpos(2) xpos(3)])+0.2, 1.038, '*k','Markersize',8,'linewidth',1.1) % marker2
    end
elseif sum(conds) > 1 && RTid == 2
    % AM ns, p<0.001
    plot([xpos(1)+ts(1) xpos(1)+ts(3)], [5.65 5.65],'-k','linewidth',1) % line
    plot(mean([xpos(1)+ts(1) xpos(1)+ts(3)])-0.17, 5.8, '*k','Markersize',8,'linewidth',1.1) % marker1
    plot(mean([xpos(1)+ts(1) xpos(1)+ts(3)]), 5.8, '*k','Markersize',8,'linewidth',1.1) % marker2
    plot(mean([xpos(1)+ts(1) xpos(1)+ts(3)])+0.17, 5.8, '*k','Markersize',8,'linewidth',1.1) % marker3

    % NBN ns, p=0.04
    plot([xpos(2)+ts(1) xpos(2)+ts(3)], [5.65 5.65],'-k','linewidth',1) % line
    plot(mean([xpos(2)+ts(1) xpos(2)+ts(3)]), 5.8, '*k','Markersize',8,'linewidth',1.1) % marker1

    % silence p=0.01, ns
    plot([xpos(3)+ts(1) xpos(3)+ts(2)], [5.65 5.65],'-k','linewidth',1) % line
        plot(mean([xpos(3)+ts(1) xpos(3)+ts(2)]), 5.8, '*k','Markersize',8,'linewidth',1.1) % marker1


elseif sum(conds) > 1 && RTid == 3
    % AM ns, p=0.01
    plot([xpos(1)+ts(1) xpos(1)+ts(3)], [5.65 5.65],'-k','linewidth',1) % line
    plot(mean([xpos(1)+ts(1) xpos(1)+ts(3)]), 5.8, '*k','Markersize',8,'linewidth',1.1) % marker1

    % NBN ns, p<0.001
    plot([xpos(2)+ts(1) xpos(2)+ts(3)], [5.65 5.65],'-k','linewidth',1) % line
    plot(mean([xpos(2)+ts(1) xpos(2)+ts(3)])-0.17, 5.8, '*k','Markersize',8,'linewidth',1.1) % marker1
    plot(mean([xpos(2)+ts(1) xpos(2)+ts(3)]), 5.8, '*k','Markersize',8,'linewidth',1.1) % marker2
    plot(mean([xpos(2)+ts(1) xpos(2)+ts(3)])+0.17, 5.8, '*k','Markersize',8,'linewidth',1.1) % marker3


    % silence ns, ns
end






