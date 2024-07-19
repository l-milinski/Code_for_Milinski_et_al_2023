%% Plot figures related to operant gap detection.

% Run sections in order to plot individual figure panels.


%% Loading and processing data

% Ensure to data paths are correct.
% This script requires the sigm_fit function: R P (2023). sigm_fit (https://www.mathworks.com/matlabcentral/fileexchange/42641-sigm_fit), MATLAB Central File Exchange. Retrieved June 26, 2023.

clear all
close all

% --- user input ---
dosave = 0;
doplot = 0;
tbp = 1; % 1 HR, 2 RT, 3 cRT, 4 iRT
% ------------------

freqvalues = [1 4 8 16 30]; % stimulus frequencies (NBN centre freqs, 30 stands for the BBN)
freqvalues_no8 = [1 4 16 30]; % stimulus frequencies for animals without 8kHz centred NBN stimuli 
Anmlvalues = [1709 1710 1711 1712 1901 1902 1903 1904]; % animal IDs



timepath = {'pre','post','post2'}; %'post'
timecounter = 0;
%define frequency to test (1K,4K or 16K)
close all
SPSSdata=[];

% DATAPATH

% path ='D:\PC-VV001-restoration\Ferrets\Behaviour\operant_gap_detection_data\1901to1904\all_ferrets\Gap1K_4K_8K_16K_30K_1901-1904_1709-1712 - revisedNames\'
path ='D:\PC-VV001-restoration\manuscript\Repository\Operant_gap_data\'


for time = 1:length(timepath)
    timecounter = timecounter +1;


counter = 0;
for animals = [Anmlvalues] 
    
animal = animals; %go through the listed animals
counter = counter + 1; % indicator for current animal No.

% specify which freqs to load for each animal
if counter <= 4
    freqinput = freqvalues_no8; %17XX animals have no 8K sessions
elseif counter > 4 
    freqinput = freqvalues;
end

for freq = [freqinput] 


% addpath 'D:\data_linus\psignifit-master'
% addpath D:\PC-VV001-restoration\Ferrets\Behaviour\final_script\functions\psignifit-master

%predefine matrices to collect individual reponse times of all sessions
indRT{time,freq} = [];
indRTc{time,freq} = [];
indRTi{time,freq} = [];

   
    
    day = 0; %reset training day counter
    recfiles = []; % records succesfully loaded files
    largedelay{time}{freq,counter} = 0;
    
for num = [1:100]%training session numbers to include

    
            % copied from Ferret_operant_gap_TI_time1to3_240122
            % Note that here filenames of ferrets 19XX are specified separately for
            % time 1. 

            if time == 1 
            % pre NOE files
                if animal == 1901 || animal==1902 || animal==1903 || animal==1904
                    fileidx = {'KA'; 'KA0'; 'KBB0'; 'KBB'; 'KE0'; 'KE'}; %input individual file indexes here
                else
                    fileidx = {'KA'; 'KA0'; 'KC0'; 'KC'; 'KBB0'; 'KBB'; 'KE0'; 'KE';'KD';'KD0'}; %input individual file indexes here
                end
            elseif time == 2
            fileidx = {'KG0'; 'KG'; 'KF0';  'KF';'KB0';'KB'}; %input individual file indexes here
           % post NOE files
            elseif time == 3
                if animal == 1901 || animal==1902 || animal==1903 || animal==1904
                    fileidx = {'K0'; 'KC'; 'KD0';  'KD'}; %input individual file indexes here
                else
                    continue
                end
            end
    
    
    loadcheck = 0;
    for n = 2:length(fileidx)
    try
    load([path mat2str(animal) 'gap'...
    num2str(freq) char(fileidx(1)) num2str(num) '.mat'])
     recfiles = [recfiles; num];
     loadcheck = 1;
    catch
    try     
    load([path mat2str(animal) 'gap'...
    num2str(freq) char(fileidx(n)) num2str(num) '.mat']) %account for variability in file name
     recfiles = [recfiles; num];
     loadcheck = 1;
    catch      
      continue
    end
    end
    end
    
    
    close all %close figures popping up during file loading
    
    if loadcheck == 0
        continue
    end
    
    

%read out trialscount per session
try
trials(num) = length(s2.Responses);% same value as s2.nTrials ?
% IF FILE LOADED, INCREMENT day
day = day+1;
catch
    continue
end


%collect output in one matrix
responses = [s2.Locations; s2.Responses; s2.gaps];

%calculate correct trials
for n = 1:length(s2.Responses)
    if s2.Locations(n) == s2.Responses(n) % if target locations match response locations -> correct trial -> 1
    responses(4,n) = 1;
    else
    responses(4,n) = 0; % otherwise -> incorrect trial -> 0
    end
end

%add delays to output matrix
responses  = [responses ; s2.RTs];

%calculate performance per session
Correct(day) = sum(responses(4,:));
performance{counter}(day) = Correct(day)/length(s2.Responses);



% ----- Save reaction times.
% ----- Exclude trials based on long RTs.


%calculate cell array with individual delays saved for each gaplength
    gapvals = [0 3 5 10 20 50 100 270];
    
    % filter out extremely long reaction times
    RTfilt = 5; % 3; in sec (RTs tend to be low 3s, below 5 includes most trials)
for g = 1:8
    
    % save average delays (exclude outliers)
    idx = find(responses(3,:)==gapvals(g)); %delays for correct trials
        delays{g} = responses(5,idx);
        exidx =  delays{g} > RTfilt; % define outliers
        delays{g}(exidx) = NaN;            % exclude outliers
        responses(4,exidx) = NaN;    % exclude outliers (hit rate) 

        % save correct delays (exclude outliers)
    corridx = find(responses(3,:)==gapvals(g) & responses(4,:)== 1); %delays for correct trials
        cdelays{g} = responses(5,corridx);
        exidx =  cdelays{g} > RTfilt; % define outliers
        cdelays{g}(exidx) = NaN;            % exclude outliers
    % save incorrect delays (exclude outliers)
    incorridx = find(responses(3,:)==gapvals(g) & responses(4,:)== 0); %delays for incorrect trials
        idelays{g} = responses(5,incorridx);
        exidx = idelays{g} > RTfilt; % define outliers
        idelays{g}(exidx) = NaN;     % exclude outliers (RT)
        
    %save mean delays (mean delay per session for each gap length)
    try
    mean_delays{time}{freq,counter}(day,g) = nanmean(delays{g});
    mean_cdelays{time}{freq,counter}(day,g) = nanmean(cdelays{g});
    mean_idelays{time}{freq,counter}(day,g) = nanmean(idelays{g});
    end
    
    % exclude trials (exlcude outliers)
    incorridx = find(responses(3,:)==gapvals(g) & responses(4,:)== 0); %delays for incorrect trials
        idelays{g} = responses(5,incorridx);
        exidx = idelays{g} > RTfilt; % define outliers
        idelays{g}(exidx) = NaN;            % exclude outliers
    

end

     % -----


%number correct trials for gaps [0 3 10 50 100 270]
gapCorr = [
sum( length( find(responses(3,:)==0  & responses(4,:)==1) ) )...
sum( length( find(responses(3,:)==3  & responses(4,:)==1) ) )...
sum( length( find(responses(3,:)==5  & responses(4,:)==1) ) )...
sum( length( find(responses(3,:)==10 & responses(4,:)==1) ) )...
sum( length( find(responses(3,:)==20 & responses(4,:)==1) ) )...
sum( length( find(responses(3,:)==50 & responses(4,:)==1) ) )...
sum( length( find(responses(3,:)==100 & responses(4,:)==1) ) )...
sum( length( find(responses(3,:)==270 & responses(4,:)==1) ) )...
];
 
%number incorrect trials for gaps [0 3 10 50 100 270]
gapIncorr = [
sum( length( find(responses(3,:)==0  & responses(4,:)==0) ) )...
sum( length( find(responses(3,:)==3  & responses(4,:)==0) ) )...
sum( length( find(responses(3,:)==5  & responses(4,:)==0) ) )...
sum( length( find(responses(3,:)==10 & responses(4,:)==0) ) )...
sum( length( find(responses(3,:)==20 & responses(4,:)==0) ) )...
sum( length( find(responses(3,:)==50 & responses(4,:)==0) ) )...
sum( length( find(responses(3,:)==100 & responses(4,:)==0) ) )...
sum( length( find(responses(3,:)==270 & responses(4,:)==0) ) )...
];




%collect for each session
% ratioCorrect nCorrect nIncorrect
ratioGapCorr{time,freq}{counter}(day,:) = gapCorr./ (gapCorr+gapIncorr);
sessGapCorr{time,freq}{counter}(day,:) = gapCorr;
sessGapIncorr{time,freq}{counter}(day,:) = gapIncorr;

% delete files where no trials were perfromed for AM or NBN
for gap = 1:8
if gapCorr(gap) == 0 && gapIncorr(gap) == 0
    ratioGapCorr{time,freq}{counter}(day,gap) = NaN;
end
end


%collect number correct & number incorrect trials for all sessions (day) and
%animals (counter)
totGAPcorr{time}{counter}(day,:) = gapCorr;
totGAPincorr{time}{counter}(day,:) = gapIncorr;



end




gap = 1:8; % number of different gaplengths
  
%--------------------------------------------------------------------------
% select sessions with enough trials


%collect number correct & number incorrect trials for all sessions (day) and
%animals (counter)
% for freq = [1 4 16 30];
% for animal  = 1:4
    totsessTrials{time}{freq,counter} = [];
    fewTrials{time}{freq,counter} = [];
    fewTrialsidx{time}{freq,counter} = [];
    goodTrialsidx{time}{freq,counter} = [];
    
    

    
    
% --- MARK SESSIONS WITH FEW TRIALS
    for sess = 1:day % size(ratioGapCorr{counter},1)
        totTrials{freq,counter}(sess,:) = totGAPcorr{time}{counter}(sess,:) + totGAPincorr{time}{counter}(sess,:);
        
% After selecting good sessions, only include sessions until number
% of trials passes 1000.

        % check if current session has sufficient number of trials
        count = 0; % set count to 0
        for n = 1: length( totTrials{freq,counter}(sess,:) )
            if totTrials{freq,counter}(sess,n) < 3
                count = count +1; % count +1 if there are less than 3 trials for that stimulus
            end
        end

        
        if count > 4 % if 5 or more gaplengths have few trials, ignore this session
            fewTrials{time}{freq,counter} = [fewTrials{time}{freq,counter}; totTrials{freq,counter}(sess,:)];
            fewTrialsidx{time}{freq,counter} = [fewTrialsidx{time}{freq,counter} sess];
        else % if sufficient number of trials, add session to output file
            totsessTrials{time}{freq,counter} = [totsessTrials{time}{freq,counter}; totTrials{freq,counter}(sess,:)];
            goodTrialsidx{time}{freq,counter} = [goodTrialsidx{time}{freq,counter} sess];
        end
    end
    
    
    % Now include only sessions with max. 1000 cumulative trials
    metasum =0; gscount = 0;
    for n = 1:size(totsessTrials{time}{freq,counter},1) % number of sessions (for current time, freq and animal)
        metasum = metasum + sum(totsessTrials{time}{freq,counter}(n,:)); %add trials count of current session to cumulative total
        if metasum < 1000 % count to the session before cumulative trials pass 1000 trials
        gscount = gscount + 1;
        end
    end
    
    % Uses gscount as number of sessions to include for current time,
    % freq, animal. Ignore sessions after that.
    
  
  % Collect performance in individual sessions
  %  -> only take sessions until 1000 trials are passed
    for n = 1:length(goodTrialsidx{time}{freq,counter})
      nTotal{time}{freq,counter} = totGAPcorr{time}{counter}( goodTrialsidx{time}{freq,counter}(1:gscount),: ) + ...
               totGAPincorr{time}{counter}( goodTrialsidx{time}{freq,counter}(1:gscount),: );
      HR{time}{freq,counter} = totGAPcorr{time}{counter}( goodTrialsidx{time}{freq,counter}(1:gscount),: )./...
                                nTotal{time}{freq,counter};
      allHR_separated{time}{freq,counter} =  mean(HR{time}{freq,counter});
    end
  
  % Collect RTs for selected individual sessions
  % -> only take sessions until 1000 trials are passed

    for n = 1:length(goodTrialsidx{time}{freq,counter})
    meanRTs{time}{freq,counter} = mean_delays{time}{freq,counter}(goodTrialsidx{time}{freq,counter}(1:gscount),:);
    cRTs{time}{freq,counter} = mean_cdelays{time}{freq,counter}(goodTrialsidx{time}{freq,counter}(1:gscount),:) ;
    iRTs{time}{freq,counter} = mean_idelays{time}{freq,counter}(goodTrialsidx{time}{freq,counter}(1:gscount),:);
    end
  
  % Collect perfromance in early & late sessions
  % Divide into early and late sessions (ignore 1 session if uneven
  % number)
  binsize_raw = length(goodTrialsidx{time}{freq,counter})/2;
  binsize = floor(binsize_raw); 
  
  try
  % early sessions
        HR1{time}{freq,counter} = HR{time}{freq,counter}(1:binsize,:);
        % animal averages
        if size(HR1{time}{freq,counter},1) > 1 
        HR1_means{time}{freq,counter} = mean(HR1{time}{freq,counter});
        else % if only one session included, take no average
        HR1_means{time}{freq,counter} = HR1{time}{freq,counter};
        end

  % late sessions
      try
        HR2{time}{freq,counter} = HR{time}{freq,counter}(binsize+2:end,:);
      catch % if there are only two sessions
        HR2{time}{freq,counter} = HR{time}{freq,counter}(binsize+1:end,:);
      end
        % animal averages
        if size(HR2{time}{freq,counter},1) > 1 
        HR2_means{time}{freq,counter} = mean(HR2{time}{freq,counter});
        else % if only one session included, take no average
        HR2_means{time}{freq,counter} = HR2{time}{freq,counter};
        end
  catch
  end
  
minval = 0.4410; 
maxval = 0.8450;

r1_minval = round(minval, 2);
r1_maxval = round(maxval, 2);
  

end
     
  
end
end


 
% produce SPPS matrix A 

% including RTs & with 1,2,3... as IDs
Anmlvalues = [1709 1710 1711 1712 1901 1902 1903 1904];

%  HR{time}{freq,counter};
 ID =[]; SESSION = []; HTR = []; GAP =[]; FREQ =[]; RT = []; cRT = []; iRT = []; TIME = [];
 
 IDc = 0; % ID counter
 for time = 1:3
     for freq = [1 4 8 16 30]
         for id = 1:8
             for gap = 1:8
                     for sess = 1:size(HR{time}{freq,id},1)
                              ID = [ID;id];
                              GAP =[GAP;gap];
                              SESSION = [SESSION;sess];
                              FREQ = [FREQ;freq];
                              HTR = [HTR;HR{time}{freq,id}(sess,gap)];
                              RT = [RT;meanRTs{time}{freq,id}(sess,gap)];
                              cRT = [cRT;cRTs{time}{freq,id}(sess,gap)];
                              iRT = [iRT;iRTs{time}{freq,id}(sess,gap)];
                              TIME=[TIME;time];

                                  if length(FREQ) == length(HTR)
                                  else
                                      display('STOP')
                                      keyboard
                                  end
                     end
             end
         end
     end
 end
 length(RT)
 length([cRT])
 length([iRT])

% Construct table
HRBLm = [ID SESSION HTR GAP FREQ RT cRT iRT TIME]; % all
HRBLm_table = array2table(HRBLm); % convert to table
HRBLm_table.Properties.VariableNames(1:9) = {'ID','SESSION','HTR','GAP','FREQ','RT','cRT','iRT','TIME'}; % name the columns


% Alternative matrixes
RTBL = [ID SESSION RT GAP FREQ TIME]; % RTs
cRTBL = [ID SESSION cRT GAP FREQ TIME]; % cRTs
iRTBL = [ID SESSION iRT GAP FREQ TIME]; % iRTs

      
            
%% FIGURE  1C: Plot by gap type (HTR, cRTs, iRTs) (Confidence Intervals)
% plots animal means as markers & session measn as box plots


if tbp == 1; HRBL = HRBLm; % to be plotted
elseif tbp ==2; HRBL = RTBL;
elseif tbp ==3; HRBL = cRTBL;
elseif tbp ==4; HRBL = iRTBL;
end


% HRBL = [ID SESSION HTR GAP FREQ];
Markers = {'+','o','*','x','v','d','^','s','>','<'};
cmap = [.5 .5 1; 1 .5 .5; 0.9 0.6 0.3];

clear anmmean freq_gap_anmIDX gaptyIDX gap_anmIDX anmmean1 anmmean2...
    sessmeanIDXnongap sessmeanIDXgap sessmean

time =  1:3

freqc=0;
    freqc = freqc+1;
               

    % indexes for gap type
    gaptyIDX{1} = find(HRBL(:,4) == 1); % no gap stimuli
    gaptyIDX{2} = find(HRBL(:,4) ~= 1); % gap stimuli
   for t = time  
   for anm = 1:8
        % find indexes for freq & gap & animal
%         freq_anmIDX{freqc}{anm} = find(HRBL(:,5) == freq & HRBL(:,1) == Anmlvalues(anm));
        % find indexes for gap type & animal
        if tbp == 1
        gap_anmIDX{t}{anm}{1} = find(HRBL(:,4) == 1 & HRBL(:,1) == anm & HRBL(:,9) == t);
        gap_anmIDX{t}{anm}{2} = find(HRBL(:,4) ~= 1 & HRBL(:,1) == anm & HRBL(:,9) == t);
        else
        gap_anmIDX{t}{anm}{1} = find(HRBL(:,4) == 1 & HRBL(:,1) == anm & HRBL(:,6) == t);
        gap_anmIDX{t}{anm}{2} = find(HRBL(:,4) ~= 1 & HRBL(:,1) == anm & HRBL(:,6) == t);
        end
        % calculate animal means
        anmmean1{t}(anm) = nanmean( HRBL(gap_anmIDX{t}{anm}{1},3) );
        anmmean2{t}(anm) = nanmean( HRBL(gap_anmIDX{t}{anm}{2},3) );
   end
   
   
   % index for session means
   if tbp == 1
   sessmeanIDXnongap{t} = find(HRBL(:,4) == 1 & HRBL(:,9) == t);
   sessmeanIDXgap{t} = find(HRBL(:,4) ~= 1 & HRBL(:,9) == t);
   else
   sessmeanIDXnongap{t} = find(HRBL(:,4) == 1 & HRBL(:,6) == t);
   sessmeanIDXgap{t} = find(HRBL(:,4) ~= 1 & HRBL(:,6) == t);
   end
   
   % session means
   sessmean{t}{1} = HRBL(sessmeanIDXnongap{t},3);
   sessmean{t}{2} = HRBL(sessmeanIDXgap{t},3);
   
   end
   
   tb2 = HRBL;
   

% produce CIs

clear yCI95 meanb bSEM CI95
% close all

for t = time
for pos = 1:2 % 1=gap; 2=nongap
   meanb{t}(pos) = nanmean(sessmean{t}{pos});
   bSEM{t}(pos) = nanstd(sessmean{t}{pos})/sqrt( length(sessmean{t}{pos}) );  
   CI95{t}(pos,:) = tinv([0.025 0.975], length(sessmean{t}{pos})-1);                    % Calculate 95% Probability Intervals Of t-Distribution
   yCI95{t}(pos,:) = bsxfun(@times, bSEM{t}(pos), CI95{t}(pos,:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
end
end

% PLOT

xshift = [-0.2 0 0.2]; % shift of x coordinates by time
figure
hold all

for t = time

% PLOT
if sum(time) > 1
% plot ind animal values (lines)
for anm = 1:8
%         plot([1:2]+xshift(t),[anmmean1{t}(anm) anmmean2{t}(anm) ],'-','color',[.7 .7 .7],'Markersize',7,'linewidth',1.5)
end

% plot ind animal values (markers)
for pos = 1:2
   for anm = 1:8
       if pos == 1
        e(t) = plot([pos]+xshift(t),[anmmean1{t}(anm)],[strcat(Markers{anm})],'color',cmap(t,:),'Markersize',7,'linewidth',1.5)
       else
        plot([pos]+xshift(t),[anmmean2{t}(anm)],[strcat(Markers{anm})],'color',cmap(t,:),'Markersize',7,'linewidth',1.5)
       end
   end
end




% plot mean
% plot([1:2]+xshift(t), meanb{t},'linewidth',3,'color','b') % Plot Mean Of All Experiments

% plot CI as errorbar
for pos = 1:2
errorbar([pos]+xshift(t),meanb{t}(pos),yCI95{t}(pos,2),'k','linewidth',2,'CapSize',25) % Plot 95% Confidence as errorbars
end

else % for BL only
    % plot ind animal values (lines)
    for anm = 1:8
            plot([1 2],[anmmean1{t}(anm) anmmean2{t}(anm)],'-','color',[.7 .7 .7],'Markersize',7,'linewidth',1.5)
    end
    % plot ind animal values (markers)
    for pos = 1:2
       for anm = 1:8
           if pos == 1
            e(t) = plot([pos],[anmmean1{t}(anm)],[strcat(Markers{anm})],'color',cmap(t,:),'Markersize',7,'linewidth',1.5);
           else
            plot([pos],[anmmean2{t}(anm)],[strcat(Markers{anm})],'color',cmap(t,:),'Markersize',7,'linewidth',1.5)
           end
       end
    end
    % connect means
    plot([1 2],[meanb{t}(1) meanb{t}(2)],'k','linewidth',3,'color',cmap(t,:)) 
    % plot CI as errorbar
    for pos = 1:2
    errorbar([pos],meanb{t}(pos),yCI95{t}(pos,2),'k','linewidth',2.5,'CapSize',35) % Plot 95% Confidence as errorbars
    end
    

end

end

if tbp == 1
ylabel({'Proportion of','correct responses'})
ylim([0.5 1])
elseif tbp == 2
ylabel({'Response time (s)'})
% ylim([1 2.5])
elseif tbp == 3
ylabel({'Response time (s)','Correct responses'})
ylim([1.3 2.3])
    if sum(time) > 1
        ylim([1 3])
    end
elseif tbp == 4
ylabel({'Response time (s)','Incorrect responses'})
ylim([1.3 2.3])
    if sum(time) > 1
        ylim([1 3])
    end
end

    
% significance markers
if sum(time) > 1
if tbp == 1
    % no gap, p = 0.028, <0.001
    plot([[1]+xshift(1) [1]+xshift(2)], [0.925 0.925],'-k','linewidth',1) % line
    plot(mean([[1]+xshift(1) [1]+xshift(2)]), 0.934, '*k','Markersize',7,'linewidth',1.1) % marker1

    plot([[1]+xshift(1) [1]+xshift(3)], [0.95 0.95],'-k','linewidth',1) % line
    plot(mean([[1]+xshift(1) [1]+xshift(3)])-0.08, 0.96, '*k','Markersize',7,'linewidth',1.1) % marker1
    plot(mean([[1]+xshift(1) [1]+xshift(3)]), 0.96, '*k','Markersize',7,'linewidth',1.1) % marker2
    plot(mean([[1]+xshift(1) [1]+xshift(3)])+0.08, 0.96, '*k','Markersize',7,'linewidth',1.1) % marker3

    % gap p=0.06, >0.001
    plot([[2]+xshift(1) [2]+xshift(3)], [0.95 0.95],'-k','linewidth',1) % line
    plot(mean([[2]+xshift(1) [2]+xshift(3)])-0.08, 0.96, '*k','Markersize',7,'linewidth',1.1) % marker1
    plot(mean([[2]+xshift(1) [2]+xshift(3)]), 0.96, '*k','Markersize',7,'linewidth',1.1) % marker2
    plot(mean([[2]+xshift(1) [2]+xshift(3)])+0.08, 0.96, '*k','Markersize',7,'linewidth',1.1) % marker3

elseif tbp == 3
    % non-gap, p=0.01, <0.001
    plot([[1]+xshift(1) [1]+xshift(2)], [2.65 2.65],'-k','linewidth',1) % line
    plot(mean([[1]+xshift(1) [1]+xshift(2)]), 2.75, '*k','Markersize',7,'linewidth',1.1) % marker1

    plot([[1]+xshift(1) [1]+xshift(3)], [2.7 2.7],'-k','linewidth',1) % line
    plot(mean([[1]+xshift(1) [1]+xshift(3)])-0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker1
    plot(mean([[1]+xshift(1) [1]+xshift(3)]), 2.75, '*k','Markersize',7,'linewidth',1.1) % marker2
    plot(mean([[1]+xshift(1) [1]+xshift(3)])+0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker3

    % gap, ns, p<0.001
    plot([[2]+xshift(1) [2]+xshift(3)], [2.7 2.7],'-k','linewidth',1) % line
    plot(mean([[2]+xshift(1) [2]+xshift(3)])-0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker1
    plot(mean([[2]+xshift(1) [2]+xshift(3)]), 2.75, '*k','Markersize',7,'linewidth',1.1) % marker2
    plot(mean([[2]+xshift(1) [2]+xshift(3)])+0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker3

elseif tbp == 4
    % no gap, ns, p<0.001
    plot([[1]+xshift(1) [1]+xshift(3)], [2.7 2.7],'-k','linewidth',1) % line
    plot(mean([[1]+xshift(1) [1]+xshift(3)])-0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker1
    plot(mean([[1]+xshift(1) [1]+xshift(3)]), 2.75, '*k','Markersize',7,'linewidth',1.1) % marker2
    plot(mean([[1]+xshift(1) [1]+xshift(3)])+0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker3

    % gap ns, p<0.001
    plot([[2]+xshift(1) [2]+xshift(3)], [2.7 2.7],'-k','linewidth',1) % line
    plot(mean([[2]+xshift(1) [2]+xshift(3)])-0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker1
    plot(mean([[2]+xshift(1) [2]+xshift(3)]), 2.75, '*k','Markersize',7,'linewidth',1.1) % marker2
    plot(mean([[2]+xshift(1) [2]+xshift(3)])+0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker3

end

elseif sum(time) == 1
    % no gap, ns, p<0.001
    plot([1 2], [0.93 0.93],'-k','linewidth',1) % line
    plot(mean([[1 2]])-0.05, 0.94, '*k','Markersize',7,'linewidth',1.1) % marker1
    plot(mean([[1 2]])+0.05, 0.94, '*k','Markersize',7,'linewidth',1.1) % marker2
end


%  plot settings
    xticks([1 2])
    xticklabels({'No Gap','Gap'})
    xlabel('Trial type')
%     ylabel({'Proportion of','correct responses'})
    xlim([0 3])
    title([ 'Performance by trial type' ])
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',14)
    grid off; 
    box off
    % LEGEND
%     if sum(time) > 1
%     legend([e(1) e(2) e(3)],'BL','Post 1','Post 2','location','southeast')
%     end
 
 
            
         
            
                       
%%  SPSS matrix for FA and RT non-gap

clear freq_gap_anmIDX
for t = 1:3
freqc=0;
for freq = [1 4 8 16 30]
    freqc = freqc+1;
           for gap = 1
                for anm = 1:8
                    % find indexes
                    if tbp == 1
                    freq_gap_anmIDX{t}{freqc}{anm}{gap} = find(HRBLm(:,5) == freq & HRBLm(:,1) == anm & HRBLm(:,4) == gap & HRBLm(:,9) == t);
                    end
                end
           end
end
end

% HRBLm = [ID SESSION HTR GAP FREQ RT cRT iRT]; % all

FA = []; RTcont = []; cRTcont = []; iRTcont = []; ID =[];  FREQ = []; SESSION = []; TIME = [];
for t = 1:3
freqc = 0;
for gap = 1
    for freq = [1 4 8 16 30]
        freqc = freqc+1;
       for anm = 1:8
            FA =[FA; 1-HRBLm(freq_gap_anmIDX{t}{freqc}{anm}{gap},3)];
            RTcont =[RTcont; HRBLm(freq_gap_anmIDX{t}{freqc}{anm}{gap},6)];
            cRTcont =[cRTcont; HRBLm(freq_gap_anmIDX{t}{freqc}{anm}{gap},7)];
            iRTcont =[iRTcont; HRBLm(freq_gap_anmIDX{t}{freqc}{anm}{gap},8)];
            for sess = 1:length( HRBLm(freq_gap_anmIDX{t}{freqc}{anm}{gap},8) ) % for each session for that freq, anm at gap = 1
            ID = [ID; anm];
            FREQ = [FREQ; freq];
            SESSION = [SESSION; sess];
            TIME = [TIME;t];
            end
       end
    end
end
end

% Construct table
clear RTcontFA RTcontFA_table


RTcontFA = [ID SESSION FREQ FA RTcont cRTcont iRTcont TIME];
RTcontFA_table = array2table(RTcontFA); % convert to table
RTcontFA_table.Properties.VariableNames(1:8) = {'ID','SESSION','FREQ','FA','RTcont','cRTcont','iRTcont','TIME'}; % name the columns






%% FIGs B: by frequency (confidence interval)
clear yCI95 meanb bSEM CI95 stdb
%  close all
% sess means & CI

time = [1 2 3]; % 1 or [1 2 3]

clear freqIDX
% close all
% tbp = 3
% dosave = 2
% dosave = 0

if tbp == 1; HRBL = HRBLm; % to be plotted
elseif tbp ==2; HRBL = RTBL;
elseif tbp ==3; HRBL = cRTBL;
elseif tbp ==4; HRBL = iRTBL;
end


% ---  prepare input matrix

c = 0; % counter
for t = time
for freq = [1 4 8 16 30]
    c = c+1;
    if tbp==1
       freqIDX{t}{freq} = find(HRBL(:,5) == freq & HRBL(:,9) == t);
    else
       freqIDX{t}{freq} = find(HRBL(:,5) == freq & HRBL(:,6) == t);
    end

end
end


% plot animal means
Markers = {'+','o','*','x','v','d','^','s','>','<'};
clear freqanmIDX anmmeanmean sessmean sessmeanIDX

for t = time
    freqc = 0; % reset
for freq = [1 4 8 16 30]
    freqc = freqc+1;
         % find indexes for session means
         if tbp == 1
         sessmeanIDX{t}{freqc} = find(HRBL(:,5) == freq & HRBL(:,9) == t);
         else
         sessmeanIDX{t}{freqc} = find(HRBL(:,5) == freq & HRBL(:,6) == t);
         end
         
           for anm = 1:8
                % find indexes for animal means
                if tbp ==1
                freqanmIDX{t}{freqc}{anm} = find(HRBL(:,5) == freq & HRBL(:,1) == anm & HRBL(:,9) == t);
                else
                freqanmIDX{t}{freqc}{anm} = find(HRBL(:,5) == freq & HRBL(:,1) == anm & HRBL(:,6) == t);
                end
                % calculate animal means
                anmmeanmean{t}(anm,freqc) = nanmean( HRBL(freqanmIDX{t}{freqc}{anm},3) );

           end
         % session means
%          sessmean{freqc} = HRBL(sessmeanIDX{freqc},3);
         sessmean{t}{freqc} = [HRBL(freqanmIDX{t}{freqc}{1},3);HRBL(freqanmIDX{t}{freqc}{2},3);HRBL(freqanmIDX{t}{freqc}{3},3);HRBL(freqanmIDX{t}{freqc}{4},3);...
             HRBL(freqanmIDX{t}{freqc}{5},3);HRBL(freqanmIDX{t}{freqc}{6},3);HRBL(freqanmIDX{t}{freqc}{7},3);HRBL(freqanmIDX{t}{freqc}{8},3)];

end
end

% plot ------

figure
hold all
% dosave = 2
% tbp = 3
for t = time




for freqc = 1:5
   meanb{t}(freqc) = nanmean(sessmean{t}{freqc});
   stdb{t}(freqc) = nanstd(sessmean{t}{freqc});
   bSEM{t}(freqc) = nanstd(sessmean{t}{freqc})/sqrt( length(sessmean{t}{freqc}) );  
   CI95{t}(freqc,:) = tinv([0.025 0.975], length(sessmean{t}{freqc})-1);                    % Calculate 95% Probability Intervals Of t-Distribution
   yCI95{t}(freqc,:) = bsxfun(@times, bSEM{t}(freqc), CI95{t}(freqc,:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
end



xshift = [-0.2 0 0.2]; % shift of x coordinates by time


% PLOT
if sum(time)>1
    % plot ind animal values (markers)
    freqc = 0;
    for freq = [1 4 8 16 30]
    freqc = freqc+1;
       for anm = 1:8
            plot([freqc]+xshift(t),[anmmeanmean{t}(anm,freqc)],[strcat(Markers{anm}) '-'],'color',cmap(t,:),'Markersize',7,'linewidth',1.5)

       end
    end

    % plot CI as errorbar
    for freqc = 1:5
    errorbar([freqc]+xshift(t),meanb{t}(freqc),yCI95{t}(freqc,2),'color','k','linewidth',2,'CapSize',20) % Plot 95% Confidence as errorbars
    end

else % for BL only
    
    % plot ind animal values (markers)
    freqc = 0;
    for freq = [1 4 8 16 30]
    freqc = freqc+1;
       for anm = 1:8
            plot([freqc],[anmmeanmean{t}(anm,freqc)],[strcat(Markers{anm}) '-'],'color',cmap(t,:),'Markersize',7,'linewidth',1.5)

       end
    end
    % Plot line conecting means
    plot([1:5],[meanb{t}(1) meanb{t}(2) meanb{t}(3) meanb{t}(4) meanb{t}(5)],'color',cmap(t,:),'linewidth',2.5) % Plot 95% Confidence as errorbars
    % plot CI as errorbar
    for freqc = 1:5
    errorbar([freqc],meanb{t}(freqc),yCI95{t}(freqc,2),'color','k','linewidth',2,'CapSize',30) % Plot 95% Confidence as errorbars
    end
    
    


end

end




               
xticks([1:5])
yticks([.1:.1:1])
xticklabels({'1 kHz NBN', '4 kHz NBN','8 kHz NBN*', '16 kHz NBN', 'BBN'})
xtickangle(15)
xlabel('Stimulus type')
xlim([0 6])
title([ 'Performance by stimulus type' ])
ylim([0.5 1])
% if freq == 30
%     title([ 'Performance for BBN' ])
% end
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
grid off; 
box off

if tbp == 1
ylabel({'Proportion of','gap trial correct responses'})
% ylim([0.5 1.1])
elseif tbp == 2
ylabel({'Response time (s)'})
% ylim([1 2.5])
elseif tbp == 3
ylabel({'Response time (s)','Correct responses'})
% ylim([1 2.5])
elseif tbp == 4
ylabel({'Response time (s)','Incorrect responses'})
% ylim([1 2.5])
end

if tbp == 3 || tbp == 4 % for RTs
    ylim([1 3])
    if sum(time) == 1   % for BL RTs
    ylim([1.3 2.3])
    yticks([1.2:.2:2.2])

    end
end

if sum(time) == 1 && tbp == 1
    ylim([0.5 1.1])
end


% significance markers

if sum(time) > 1
if tbp == 1
%1K, p= ns, ns
%4K p= ns, <0.001
plot([[2]+xshift(1) [2]+xshift(3)], [0.95 0.95],'-k','linewidth',1) % line
plot(mean([2+xshift(1) 2+xshift(3)])-0.12, 0.961, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([2+xshift(1) 2+xshift(3)]), 0.961, '*k','Markersize',7,'linewidth',1.1) % marker2
plot(mean([2+xshift(1) 2+xshift(3)])+0.12, 0.961, '*k','Markersize',7,'linewidth',1.1) % marker3

%8K, p<0.001, <0.001
plot([3+xshift(1) [3]+xshift(2)], [0.92 0.92],'-k','linewidth',1) % line
plot(mean([3+xshift(1) 3+xshift(2)])-0.12, 0.931, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([3+xshift(1) 3+xshift(2)]), 0.931, '*k','Markersize',7,'linewidth',1.1) % marker2
plot(mean([3+xshift(1) 3+xshift(2)])+0.12, 0.931, '*k','Markersize',7,'linewidth',1.1) % marker3

plot([[3]+xshift(1) [3]+xshift(3)], [0.95 0.95],'-k','linewidth',1) % line
plot(mean([3+xshift(1) 3+xshift(3)])-0.12, 0.961, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([3+xshift(1) 3+xshift(3)]), 0.961, '*k','Markersize',7,'linewidth',1.1) % marker2
plot(mean([3+xshift(1) 3+xshift(3)])+0.12, 0.961, '*k','Markersize',7,'linewidth',1.1) % marker2

%16K, p= ns, <0.001
plot([[4]+xshift(1) [4]+xshift(3)], [0.95 0.95],'-k','linewidth',1) % line
plot(mean([4+xshift(1) 4+xshift(3)])-0.12, 0.961, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([4+xshift(1) 4+xshift(3)]), 0.961, '*k','Markersize',7,'linewidth',1.1) % marker2
plot(mean([4+xshift(1) 4+xshift(3)])+0.12, 0.961, '*k','Markersize',7,'linewidth',1.1) % marker2

%30K p= ns, ns

elseif tbp == 3 % cRT
%1K, p= ns, <0.001
plot([[1]+xshift(1) [1]+xshift(3)], [2.7 2.7],'-k','linewidth',1) % line
plot(mean([1+xshift(1) 1+xshift(3)])-0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([1+xshift(1) 1+xshift(3)])+0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker2
%4K p= ns, <0.001
plot([[2]+xshift(1) [2]+xshift(3)], [2.7 2.7],'-k','linewidth',1) % line
plot(mean([2+xshift(1) 2+xshift(3)])-0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([2+xshift(1) 2+xshift(3)])+0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker2
%8K, p=0.01, <0.001
plot([[3]+xshift(1) [3]+xshift(2)], [2.65 2.65],'-k','linewidth',1) % line
plot([[3]+xshift(1) [3]+xshift(3)], [2.7 2.7],'-k','linewidth',1) % line
plot(mean([3+xshift(1) 3+xshift(3)])-0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([3+xshift(1) 3+xshift(3)])+0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker2
%16K, p= ns,<0.001
plot([[4]+xshift(1) [4]+xshift(3)], [2.7 2.7],'-k','linewidth',1) % line
plot(mean([4+xshift(1) 4+xshift(3)])-0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([4+xshift(1) 4+xshift(3)])+0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker2
%30K p= ns, <0.001
plot([[5]+xshift(1) [5]+xshift(3)], [2.7 2.7],'-k','linewidth',1) % line
plot(mean([5+xshift(1) 5+xshift(3)])-0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([5+xshift(1) 5+xshift(3)])+0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker2

elseif tbp == 4 % iRT
%1K, p= 0.02, <0.001
plot([[1]+xshift(1) [1]+xshift(2)], [2.65 2.65],'-k','linewidth',1) % line
plot([[1]+xshift(1) [1]+xshift(3)], [2.7 2.7],'-k','linewidth',1) % line
plot(mean([1+xshift(1) 1+xshift(3)])-0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([1+xshift(1) 1+xshift(3)])+0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker2
%4K p= ns,<0.001
plot([[2]+xshift(1) [2]+xshift(3)], [2.7 2.7],'-k','linewidth',1) % line
plot(mean([2+xshift(1) 2+xshift(3)])-0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([2+xshift(1) 2+xshift(3)])+0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker2
%8K, p=0.03, <0.001
plot([[3]+xshift(1) [3]+xshift(2)], [2.65 2.65],'-k','linewidth',1) % line
plot([[3]+xshift(1) [3]+xshift(3)], [2.7 2.7],'-k','linewidth',1) % line
plot(mean([3+xshift(1) 3+xshift(3)])-0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([3+xshift(1) 3+xshift(3)])+0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker2
%16K, p=ns, <0.001
plot([[4]+xshift(1) [4]+xshift(3)], [2.7 2.7],'-k','linewidth',1) % line
plot(mean([4+xshift(1) 4+xshift(3)])-0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([4+xshift(1) 4+xshift(3)])+0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker2
%30K p= ns, <0.001
plot([[5]+xshift(1) [5]+xshift(3)], [2.7 2.7],'-k','linewidth',1) % line
plot(mean([5+xshift(1) 5+xshift(3)])-0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([5+xshift(1) 5+xshift(3)])+0.08, 2.75, '*k','Markersize',7,'linewidth',1.1) % marker2

    
end

elseif sum(time) == 1 & tbp == 1  % for BL plot
    % BBN: p<0.001 to all, except 8K (p=0.03)
plot([[1] [5]], [.98 .98],'-k','linewidth',1) % line
plot(mean([1 5])-0.08, .995, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([1 5])+0.08, .995, '*k','Markersize',7,'linewidth',1.1) % marker2
plot([[3] [5]], [.96 .96],'-k','linewidth',1) % line
    % 8K   p<0.001 to all, except BBN (p=0.03)
plot([[3.05] [4]], [.6 .6],'-k','linewidth',1) % line
plot(mean([3.05 4])-0.08, .585, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([3.05 4])+0.08, .585, '*k','Markersize',7,'linewidth',1.1) % marker2
plot([[1 2.95]], [.6 .6],'-k','linewidth',1) % line
plot(mean([1 2.95])-0.08, .585, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([1 2.95])+0.08, .585, '*k','Markersize',7,'linewidth',1.1) % marker2
plot([[3.05] [5]], [.62 .62],'-k','linewidth',1) % line



end



%% SUPPL. FIGURE S2A & S1C: left responses over gaplength

% For FIGURE S2A set time to 1:3, for S1C set time to 1
time = 1:3;


if tbp == 1; HRBL = HRBLm; % to be plotted
elseif tbp ==2; HRBL = RTBL;
elseif tbp ==3; HRBL = cRTBL;
elseif tbp ==4; HRBL = iRTBL;
end

Markers = {'+','o','*','x','v','d','^','s','>','<'};
clear anmmean freq_gap_anmIDX freq_gapIDX sessmean
xshift = [-0.3 0 0.3]; % shift of x coordinates by time

    freqc = 0;
for freq = [1 4 8 16 30]
%     figure(freq); hold all
    freqc = freqc+1;
    for t = 1:3

        % find indexes for freq & gap
        if tbp == 1
        freq_gapIDX{t}{freqc}{gap} = find(HRBL(:,5) == freq & HRBL(:,4) == gap & HRBL(:,9) == t);
        else
        freq_gapIDX{t}{freqc}{gap} = find(HRBL(:,5) == freq & HRBL(:,4) == gap & HRBL(:,6) == t);
        end
                
        for anm = 1:8
                for gap = 1:8
                    if tbp  == 1
                    % find indexes for freq & gap
                    freq_gapIDX{t}{freqc}{gap} = find(HRBL(:,5) == freq & HRBL(:,4) == gap & HRBL(:,9) == t);
                    % find indexes
                    freq_gap_anmIDX{t}{freqc}{anm}{gap} = find(HRBL(:,5) == freq & HRBL(:,1) == anm & HRBL(:,4) == gap & HRBL(:,9) == t);
                    else
                    % find indexes for freq & gap
                    freq_gapIDX{t}{freqc}{gap} = find(HRBL(:,5) == freq & HRBL(:,4) == gap & HRBL(:,6) == t);
                    % find indexes
                    freq_gap_anmIDX{t}{freqc}{anm}{gap} = find(HRBL(:,5) == freq & HRBL(:,1) == anm & HRBL(:,4) == gap & HRBL(:,6) == t);
                    end
                % calculate animal means
                anmmean{t}{freqc}(anm,gap) = nanmean( HRBL(freq_gap_anmIDX{t}{freqc}{anm}{gap},3) );
                end
         end
           
        % calculate session means
        for gap = 1:8
        sessmean{t}{freqc}(:,gap) = HRBL(freq_gapIDX{t}{freqc}{gap},3);
        end
        % session means for FA (for FA plot, see section below)
        sessmeanFA{t}{freqc} = 1-sessmean{t}{freqc}(:,1)
                   
    end
                 
end



% Plot by gap length (across freqs, left responses)

% close all

% dosave = 1
% dosave = 0
% tbp = 3
if tbp == 1; HRBL = HRBLm; % to be plotted
elseif tbp ==2; HRBL = RTBL;
elseif tbp ==3; HRBL = cRTBL;
elseif tbp ==4; HRBL = iRTBL;
end

Markers = {'+','o','*','x','v','d','^','s','>','<'};
clear anmmean freq_gap_anmIDX freq_gapIDX sessionmean

% HRBLm =         [ID SESSION HTR GAP FREQ RT cRT iRT TIME]; % all

    figure; hold all           

        for t = time  
           for anm = 1:8
                for gap = 1:8
                    if tbp  == 1
                    % find indexes for freq & gap
                    freq_gapIDX{t}{gap} = find(HRBL(:,4) == gap & HRBL(:,9) == t);
                    % find indexes
                    freq_gap_anmIDX{t}{anm}{gap} = find(HRBL(:,1) == anm & HRBL(:,4) == gap & HRBL(:,9) == t);
                    else
                    % find indexes for freq & gap
                    freq_gapIDX{t}{gap} = find(HRBL(:,4) == gap & HRBL(:,6) == t);
                    % find indexes
                    freq_gap_anmIDX{t}{anm}{gap} = find(HRBL(:,1) == anm & HRBL(:,4) == gap & HRBL(:,6) == t);
                    end
                % calculate animal means
                anmmean{t}(anm,gap) = nanmean( HRBL(freq_gap_anmIDX{t}{anm}{gap},3) );
                end
           end
    
      
           
           % session means
           for gap = 1:8
           sessionmean{t}(:,gap) = HRBL(freq_gapIDX{t}{gap},3)
           end
           
           
        
          
if tbp == 1 
    
    if sum(time) > 1
                % plot animal means & box plots (for no gap 1-HR)
%                 b = boxplot([1-anmmean{t}(:,1) anmmean{t}(:,2:end)],'positions',[1:8]+xshift(t),'Colors','k','Widths',0.3,'Whisker',0,'Symbol','','medianstyle','line');
%                 set(b,{'linew'},{1.5})
                
                 % plot session means box plots (for no gap 1-HR)
                b = boxplot([1-sessionmean{t}(:,1) sessionmean{t}(:,2:end)],'positions',[1:8]+xshift(t),'Colors','k','Widths',0.2,'Whisker',0,'Symbol','','medianstyle','line');
                set(b,{'linew'},{1.5})

           for anm = 1:8
                for gap = 1:8
                % plot ind animal values
                if gap ~= 1
                e(t) = plot([gap]+xshift(t),[anmmean{t}(anm,gap)],strcat(Markers{anm}),'color',cmap(t,:),'Markersize',6,'linewidth',1.5)
                elseif gap == 1
                e(t) = plot([gap]+xshift(t),1-[anmmean{t}(anm,gap)],strcat(Markers{anm}),'color',cmap(t,:),'Markersize',6,'linewidth',1.5)
                end
                end
           end
           
    else   % if BL only
        
                % plot session means box plots (for no gap 1-HR)
                b = boxplot([1-sessionmean{t}(:,1) sessionmean{t}(:,2:end)],'positions',[1:8],'Colors','k','Widths',0.6,'Whisker',0,'Symbol','','medianstyle','line');
                set(b,{'linew'},{1.5})
                
           for anm = 1:8
                for gap = 1:8
                % plot ind animal values
                if gap ~= 1
                e(t) = plot([gap],[anmmean{t}(anm,gap)],strcat(Markers{anm}),'color',cmap(t,:),'Markersize',6,'linewidth',1.5)
                elseif gap == 1
                e(t) = plot([gap],1-[anmmean{t}(anm,gap)],strcat(Markers{anm}),'color',cmap(t,:),'Markersize',6,'linewidth',1.5)
                end
                end
           end
           
    end
           
elseif tbp == 3 | tbp == 4
    
    if sum(time) > 1
    
                % plot animal means as box plots
%                 b = boxplot([anmmean(:,1:end)],'positions',[1:8],'Colors','k','Widths',0.4,'Whisker',0,'Symbol','','medianstyle','line');
%                 set(b,{'linew'},{1.5})
                
                % plot session means as box plots 
                  b = boxplot([sessionmean{t}(:,1:end)],'positions',[1:8]+xshift,'Colors','k','Widths',0.3,'Whisker',0,'Symbol','','medianstyle','line');
                  set(b,{'linew'},{1.5})
                
                % plot means across sessions
%                 plot([1:8],nanmmean( [anmmean(:,1:end)] ),'--k','Markersize',10,'linewidth',1.5)

           for anm = 1:8
                for gap = 1:8
                % plot ind animal values
                plot([gap]+xshift(t),[anmmean{t}(anm,gap)],strcat(Markers{anm}),'color',cmap(t,:),'Markersize',6,'linewidth',1.5)
                end
           end
           
    else   % if BL only
            
                % plot session means as box plots 
                  b = boxplot([sessionmean{t}(:,1:end)],'positions',[1:8],'Colors','k','Widths',0.6,'Whisker',0,'Symbol','','medianstyle','line');
                  set(b,{'linew'},{1.5})
           
           for anm = 1:8
                for gap = 1:8
                % plot ind animal values
                plot([gap],[anmmean{t}(anm,gap)],strcat(Markers{anm}),'color',cmap(t,:),'Markersize',6,'linewidth',1.5)
                end
           end
           
    end
           
end


        end


    % plot settings
    if tbp == 1
    ylabel({'Proportion of','left responses'})
    ylim([0 1.1])
    elseif tbp == 2
    ylabel({'Response time (s)'})
    ylim([1 2.5])
    elseif tbp == 3
    ylabel({'Response time (s)','Correct responses'})
    ylim([1 2.5])
        if time == 1
                ylim([1.2 2.3])
        end
    elseif tbp == 4
    ylabel({'Response time (s)','Incorrect responses'})
    ylim([1 2.5])
        if time == 1
                ylim([1.2 2.3])
        end
    end
    xticks([1:8])
    xticklabels({'No gap','3','5','10','20','50','100','270'})
    xlabel('Gap length (ms)')
    
    try
    legend([e(1) e(2) e(3)],'BL','Post 1','Post 2','location','southeast')
    end

    
%     ylabel({'Proportion of','correct responses'})
%     ylim([0 1.1])
    xlim([0 9])
    title([ 'Gap detection (all stimuli)' ])
    
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',14)
    grid off;
    box off



%% FIGURE 1D (False Alarm Rate)

% close all

% HRBL = [ID SESSION HTR GAP FREQ];
Markers = {'+','o','*','x','v','d','^','s','>','<'};
clear anmmean freq_gap_anmIDX freq_gapIDX FA
clear yCI95 meanb bSEM CI95 stdb

xshift = [-0.2 0 0.2]; % shift of x coordinates by time

if tbp == 1
figure; hold on
for t = 1:3
freqc=0;
for freq = [1 4 8 16 30]
    freqc = freqc+1;
           % find indexes
                freq_gapIDX{t}{freqc}{1} = find(HRBL(:,5) == freq & HRBL(:,4) == 1 & HRBL(:,9) == t);
           for gap = 1
                for anm = 1:8
                % find indexes
                freq_gap_anmIDX{t}{freqc}{anm}{gap} = find(HRBL(:,5) == freq & HRBL(:,1) == anm & HRBL(:,4) == gap & HRBL(:,9) == t);
                % calculate animal means
                anmmean{t}{freqc}(anm,gap) = nanmean( HRBL(freq_gap_anmIDX{t}{freqc}{anm}{gap},3) );
                FA{t}{freqc}(anm,gap) = 1-anmmean{t}{freqc}(anm,gap);
                % plot ind animal values
                plot([freqc]+xshift(t),[FA{t}{freqc}(anm,gap)],strcat(Markers{anm}),'color',cmap(t,:),'Markersize',7,'linewidth',1.5)
                end
           end
           
           % produce CIs
           % use FA session means to calculate CIs -> sessmeanFA{t}{freqc}


                meanb{t}(freqc) = nanmean(sessmeanFA{t}{freqc});
                stdb{t}(freqc) = nanstd(sessmeanFA{t}{freqc});
                bSEM{t}(freqc) = nanstd(sessmeanFA{t}{freqc})/sqrt( length(sessmeanFA{t}{freqc}) );  
                CI95{t}(freqc,:) = tinv([0.025 0.975], length(sessmeanFA{t}{freqc})-1);                    % Calculate 95% Probability Intervals Of t-Distribution
                yCI95{t}(freqc,:) = bsxfun(@times, bSEM{t}(freqc), CI95{t}(freqc,:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
                
           



% session means & CI errorbars
  errorbar([freqc]+xshift(t),meanb{t}(freqc),yCI95{t}(freqc,2),'color',[0 0 0],'linewidth',2,'CapSize',20) % Plot 95% Confidence as errorbars


end

end

% significance markers
if tbp == 1
    % 1K p=0.01 & <0.001
plot([1+xshift(1) 1+xshift(2)], [0.68 0.68],'-k','linewidth',1) % line post1
plot(mean([1+xshift(1) 1+xshift(2)]), 0.70, '*k','Markersize',7,'linewidth',1.1) % marker1

plot([1+xshift(1) 1+xshift(3)], [0.73 0.73],'-k','linewidth',1) % line post2
plot(mean([1+xshift(1) 1+xshift(3)])-0.12, 0.755, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([1+xshift(1) 1+xshift(3)]), 0.755, '*k','Markersize',7,'linewidth',1.1) % marker2
plot(mean([1+xshift(1) 1+xshift(3)])+0.12, 0.755, '*k','Markersize',7,'linewidth',1.1) % marker3

    % 4K p=0.01 & <0.001
plot([2+xshift(1) 2+xshift(2)], [0.68 0.68],'-k','linewidth',1) % line post1
plot(mean([2+xshift(1) 2+xshift(2)]), 0.70, '*k','Markersize',7,'linewidth',1.1) % marker1


plot([2+xshift(1) 2+xshift(3)], [0.73 0.73],'-k','linewidth',1) % line post2
plot(mean([2+xshift(1) 2+xshift(3)])-0.12, 0.755, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([2+xshift(1) 2+xshift(3)]), 0.755, '*k','Markersize',7,'linewidth',1.1) % marker2
plot(mean([2+xshift(1) 2+xshift(3)])+0.12, 0.755, '*k','Markersize',7,'linewidth',1.1) % marker2

    % 8K p=n.s & <0.001
plot([3+xshift(1) 3+xshift(3)], [0.73 0.73],'-k','linewidth',1) % line post2
plot(mean([3+xshift(1) 3+xshift(3)])-0.12, 0.755, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([3+xshift(1) 3+xshift(3)]), 0.755, '*k','Markersize',7,'linewidth',1.1) % marker2
plot(mean([3+xshift(1) 3+xshift(3)])+0.12, 0.755, '*k','Markersize',7,'linewidth',1.1) % marker3

    % 16K p=0.08 & <0.001
plot([4+xshift(1) 4+xshift(2)], [0.68 0.68],'-k','linewidth',1) % line post1
plot(mean([4+xshift(1) 4+xshift(2)]), 0.70, '*k','Markersize',7,'linewidth',1.1) % marker1


plot([4+xshift(1) 4+xshift(3)], [0.73 0.73],'-k','linewidth',1) % line post2
plot(mean([4+xshift(1) 4+xshift(3)])-0.12, 0.755, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([4+xshift(1) 4+xshift(3)]), 0.755, '*k','Markersize',7,'linewidth',1.1) % marker2
plot(mean([4+xshift(1) 4+xshift(3)])+0.12, 0.755, '*k','Markersize',7,'linewidth',1.1) % marker3

    % BBN p=n.s & <0.001
plot([5+xshift(1) 5+xshift(3)], [0.73 0.73],'-k','linewidth',1) % line post2
plot(mean([5+xshift(1) 5+xshift(3)])-0.12, 0.755, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([5+xshift(1) 5+xshift(3)]), 0.755, '*k','Markersize',7,'linewidth',1.1) % marker2
plot(mean([5+xshift(1) 5+xshift(3)])+0.12, 0.755, '*k','Markersize',7,'linewidth',1.1) % marker3

end
hold off


% plot settings

    xticks([1:5])
    xticklabels({'1 kHz NBN', '4 kHz NBN','8 kHz NBN*', '16 kHz NBN', 'BBN'})
    xtickangle(15)
    xlabel('Stimulus Type')
    ylabel({'False alarm rate'})
    ylim([0 1])
    xlim([0 6])
    title([ 'False alarm rate' ])
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',14)
    grid off; 
    box off



end
  


%% Gap detection thresholds: pre-processing
% for each animal mean

HRBL = HRBLm; % to select HR as input (not RTs)

% close all
clear input thresh slope freq_gap_anmIDX anmmean input_meta
freqc=0;
gaplabels = [3 5 10 20 50 100 270];
for freq = [1 4 8 16 30]
    freqc = freqc+1;
    for t =  [1 2 3]

           for anm = 1:8;
                for gap = 2:8
                % find indexes
                freq_gap_anmIDX{t}{freqc}{anm}{gap} = find(HRBL(:,5) == freq & HRBL(:,1) == anm & HRBL(:,4) == gap & HRBL(:,9) == t);
                % calculate animal means
                anmmean{t}{freqc}(anm,gap) = nanmean( HRBL(freq_gap_anmIDX{t}{freqc}{anm}{gap},3) );
                %sigmoid scaffold matrix (HR given time and freq)
                input_meta(anm,gap) = anmmean{t}{freqc}(anm,gap);
                % plot ind animal values
        %         plot([xpos(gap)+ts(t)],[anmmean{t}{freqc}(anm,gap)],strcat(Markers{anm}),'color',cmap(t,:),'Markersize',7,'linewidth',1.5)
                end
           end

        % crop input to gap 1-7
        input = input_meta(:,2:end);
        
        for anm =  1:8
%             pause
%             close all
%             figure; hold all
        % fit sigmoid parameters
        fixed_params=[0, 1 , NaN , NaN];
        % fitting sigmoid
        try % try fitting - if input matrix, define thresh and slope as 'nan'
            [param] = sigm_fit(1:7,input(anm,:),fixed_params); % param = [min, max, x50, slope]
                if param(4) <= 0 % is sigmoid is flat, calculate again without parameters
                [param] = sigm_fit(1:7,input(anm,:)); % param = [min, max, x50, slope]
                elseif param(4) > 100 % if sigmoid slope is huge, calculate again without parameters
%                 fixed_params=[NaN, NaN , NaN , NaN]; 
                [param] = sigm_fit(1:7,input(anm,:)); 
                end
            
        % find closest gap to x50 -------
        
            % original
    %       nn = nearestNumber(1:7,param(3));
        % edit 120822,
        a= [1:7];
        n= param(3);
        bydiff = unique(abs(a-n),'stable'); % identify differences from n
        gapidx_meta = find(bydiff == min(bydiff))
        % if two gap length are in equal distance to threshold, take the
        % higher one
        gapidx = max(gapidx_meta);
        
        % -------------------------------
            
            thresh{freq}(anm,t) = gaplabels(gapidx);
                %if sigmoid starts above 50%, set threshold to shortest gap
                %length (3ms)
                bl50 = 0;
                for i = 1:length(input(anm,:))
                        testinput = input(anm,:);
                        if testinput(i) <= 0.5
                        bl50 = bl50+1;
                        end
                end
                if bl50 == 0
                thresh{freq}(anm,t) = 3;
                end 

            % define slope
             slope{freq}(anm,t) = param(4);
             
             % set figure title
             title(['slope = ' num2str(slope{freq}(anm,t)) ])

        catch
             thresh{freq}(anm,t) = nan;
             slope{freq}(anm,t) = nan;
             % set figure title (NAN)
             title(['slope = NaN'])
        end
        
%         plot(1:7,nanmean(input(anm,:))','O','color',cmap(t,:),'Markersize',7,'linewidth',1.5)
        plot(1:7,input(anm,:)','O','color',cmap(t,:),'Markersize',2,'linewidth',1.5)
        end
        
        end

    end
% end




% create SPSS matrix: thresholds


TIME = []; ID = []; FREQ = []; THRESH = []; SLOPE = [];
for time = 1:3
    for freq = [1 4 8 16 30]
        for anm = 1:8
            TIME = [TIME; time];
            FREQ = [FREQ; freq];
            ID = [ID;anm];
            THRESH = [THRESH; thresh{freq}(anm,time)];
            SLOPE = [SLOPE; slope{freq}(anm,time)];

        end
    end
end

ThreshSPSS_meta = [ID TIME FREQ THRESH SLOPE];

% define outliers, NaNs & exclude
out = nanmedian(SLOPE)*10; % define
  outidx = find(SLOPE > out); % find outliers
inidx = find(SLOPE < out); % find non-outliers
  % SLOPE(outidx) = nanmedian(SLOPE)*10; % replace
  % SLOPE(outidx) = nan; % exclude

ThreshSPSS = ThreshSPSS_meta(inidx,:) ; % keep only non-outliers

% OR keep all:
% ThreshSPSS = ThreshSPSS_meta ; % keep only non-outliers


ThreshSPSS_table = array2table(ThreshSPSS); % convert to table
ThreshSPSS_table.Properties.VariableNames(1:5) = {'ID','TIME','FREQ','THRESH','SLOPE'}; % name the columns

% writetable(ThreshSPSS_table,'D:\PC-VV001-restoration\Ferrets\Behaviour\Operant_gap_detection\SPSS\June2022\Ferrets_operant_thresholds_time123_Aug22.xlsx') % convert table to .xlsx file

%% FIGURE 1B: Plot thresholds over time (with individual markers for animal means)


close all

dosave = 0;


plotci = 0; % 1 fore CI plots, 0 for mean & std plots



Markers = {'+','o','*','x','v','d','^','s','>','<'};
cmap = [.5 .5 1; 1 .5 .5; 0.9 0.6 0.3];

% Define indeces -------

% Time IDX
for t= 1:3
timeidx{t} = find(ThreshSPSS(:,2) == t);
end

    % For infromation:
    ThreshSPSS_meta = [ID TIME FREQ THRESH SLOPE];
    % Timeidx picks all thresh values for time 1, 2 or 3, respectively.
    % Since thre are 4 stim frequencies (1,4,8,16Khz NBN,30 KhZ BBN) per animal per
    % time, there would normally be 5*8 thresh values for each time. However, since
    % for example animal 5 was only assessed for BL, and the non-implanted
    % animals were not ssessed for stim 8 kHz, and only 3 animals where
    % assessed for time 3, the number of thresh values varies across
    % time point.
    
    % To produce animal means, find 
    % Animal IDX
    for anm = 1:8
    anmidx{anm} = find(ThreshSPSS(:,1) == anm);
    end
    clear time_anm_idx
    % Animal & time IDX
    for t = 1:3
    for anm = 1:8
    time_anm_idx{t}{anm} = find(ThreshSPSS(:,1) == anm & ThreshSPSS(:,2) == t);
    end
    end
    
    
    % Plot -----------------------
    

cmap = [.5 .5 1; 1 .5 .5; 0.9 0.6 0.3];

% for errorbar plots
clear CI95 yCI95
 figure; hold all
for t = 1:3
   meanb{t} = nanmean(ThreshSPSS(timeidx{t},4));
   stdb{t} = nanstd(ThreshSPSS(timeidx{t},4));
   bSEM{t} = nanstd(ThreshSPSS(timeidx{t},4))/sqrt( length(ThreshSPSS(timeidx{t},4)) );  
   CI95{t} = tinv([0.025 0.975], length(ThreshSPSS(timeidx{t},4))-1);                    % Calculate 95% Probability Intervals Of t-Distribution
   yCI95{t} = bsxfun(@times, bSEM{t}, CI95{t});              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
end

if plotci == 1
for t = 1:3
errorbar([t],meanb{t},yCI95{t}(2),'linewidth',2,'CapSize',25,'color',cmap(t,:)) % Plot 95% Confidence as errorbars
plot(t,mean(ThreshSPSS(timeidx{t},4)),'O','color',cmap(t,:),'Markerfacecolor',cmap(t,:)...
        ,'Markersize',10,'linewidth',2)
end
else
for t = 1:3
errorbar([t],nanmean(ThreshSPSS(timeidx{t},4)),nanstd(ThreshSPSS(timeidx{t},4)),'linewidth',2,'CapSize',25,'color',cmap(t,:)) % Plot 95% Confidence as errorbars
plot(t,nanmean(ThreshSPSS(timeidx{t},4)),'O','color',cmap(t,:),'Markerfacecolor',cmap(t,:)...
        ,'Markersize',16,'linewidth',2.5)
end

% % plot individuals

% plot ind animal values (markers)
    clear input
    
for t = 1:3

   for anm = 1:8
       input{t}{anm} = mean( ThreshSPSS(time_anm_idx{t}{anm},4) ); % for each time, calculate animal means across all frequencies for which trehsilds are available

%        plot(t,ThreshSPSS(timeidx{t},4),[strcat(Markers{anm})],'color',cmap(t,:),'Markersize',7,'linewidth',1.5)
%        plot(t+randi([-5 5])/20,input{t}(anm),[strcat(Markers{anm})],'color','k','Markersize',7,'linewidth',1.5)
       plot(t+randi([-5 5])/20,input{t}{anm},[strcat(Markers{anm})],'color',[.3 .3 .3],'Markersize',7,'linewidth',1.5)

   end
end


end

% significance markers
%     % Bl vs Post1 p=0.02,Bl vs Post2 <0.001, Post1 vs Post2 p=0.03
% plot([1 2], [15.7 15.7],'-k','linewidth',1) % line post1
% plot(mean([1 2]), 16.2, '*k','Markersize',7,'linewidth',1.1) % marker1
% 
% plot([1 3], [17.2 17.2],'-k','linewidth',1) % line post2
% plot(mean([1 3])-0.1, 17.7, '*k','Markersize',7,'linewidth',1.1) % marker1
% plot(mean([1 3]), 17.7, '*k','Markersize',7,'linewidth',1.1) % marker2
% plot(mean([1 3])+0.1, 17.7, '*k','Markersize',7,'linewidth',1.1) % marker2
% 
% plot([2 3], [18.7 18.7],'-k','linewidth',1) % line post1
% plot(mean([2 3]), 19.2, '*k','Markersize',7,'linewidth',1.1) % marker1


plot([1 2], [15.7 15.7]+5,'-k','linewidth',1) % line post1
plot(mean([1 2]), 16.2+5, '*k','Markersize',7,'linewidth',1.1) % marker1

plot([1 3], [17.2 17.2]+5,'-k','linewidth',1) % line post2
plot(mean([1 3])-0.1, 17.7+5, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([1 3]), 17.7+5, '*k','Markersize',7,'linewidth',1.1) % marker2
plot(mean([1 3])+0.1, 17.7+5, '*k','Markersize',7,'linewidth',1.1) % marker2

plot([2 3], [18.7 18.7]+5,'-k','linewidth',1) % line post1
plot(mean([2 3]), 19.2+5, '*k','Markersize',7,'linewidth',1.1) % marker1


    


    
    

    xticks([1 2 3])
    xticklabels({'BL', '1 week','6 months'})

    xlabel('Time')
    ylabel({'Threshold';'Gap length (ms)'})
    ylim([0 25])
    xlim([0 4])
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',14)
    grid off;
    box off
    
    

%% FIGURE 1B: Plot thresholds over time (no animal means)

% close all

dosave = 0;


plotci = 0; % 1 fore CI plots, 0 for mean & std plots

for t= 1:3
timeidx{t} = find(ThreshSPSS(:,2) == t);
end

cmap = [.5 .5 1; 1 .5 .5; 0.9 0.6 0.3];

% for errorbar plots
clear CI95 yCI95
 figure; hold all
for t = 1:3
   meanb{t} = nanmean(ThreshSPSS(timeidx{t},4));
   stdb{t} = nanstd(ThreshSPSS(timeidx{t},4));
   bSEM{t} = nanstd(ThreshSPSS(timeidx{t},4))/sqrt( length(ThreshSPSS(timeidx{t},4)) );  
   CI95{t} = tinv([0.025 0.975], length(ThreshSPSS(timeidx{t},4))-1);                    % Calculate 95% Probability Intervals Of t-Distribution
   yCI95{t} = bsxfun(@times, bSEM{t}, CI95{t});              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
end

if plotci == 1
for t = 1:3
errorbar([t],meanb{t},yCI95{t}(2),'linewidth',2,'CapSize',25,'color',cmap(t,:)) % Plot 95% Confidence as errorbars
plot(t,mean(ThreshSPSS(timeidx{t},4)),'O','color',cmap(t,:),'Markerfacecolor',cmap(t,:)...
        ,'Markersize',10,'linewidth',2)
end
else
for t = 1:3
errorbar([t],nanmean(ThreshSPSS(timeidx{t},4)),nanstd(ThreshSPSS(timeidx{t},4)),'linewidth',2,'CapSize',25,'color',cmap(t,:)) % Plot 95% Confidence as errorbars
plot(t,nanmean(ThreshSPSS(timeidx{t},4)),'O','color',cmap(t,:),'Markerfacecolor',cmap(t,:)...
        ,'Markersize',10,'linewidth',2)
end
end

% significance markers
% if tbp == 1
    % Bl vs Post1 p=0.02,Bl vs Post2 <0.001, Post1 vs Post2 p=0.03
plot([1 2], [15.7 15.7],'-k','linewidth',1) % line post1
plot(mean([1 2]), 16.2, '*k','Markersize',7,'linewidth',1.1) % marker1

plot([1 3], [17.2 17.2],'-k','linewidth',1) % line post2
plot(mean([1 3])-0.1, 17.7, '*k','Markersize',7,'linewidth',1.1) % marker1
plot(mean([1 3]), 17.7, '*k','Markersize',7,'linewidth',1.1) % marker2
plot(mean([1 3])+0.1, 17.7, '*k','Markersize',7,'linewidth',1.1) % marker2

plot([2 3], [18.7 18.7],'-k','linewidth',1) % line post1
plot(mean([2 3]), 19.2, '*k','Markersize',7,'linewidth',1.1) % marker1


% end
    
    

    xticks([1 2 3])
    xticklabels({'BL', '1 week','6 months'})

    xlabel('Time')
    ylabel({'Threshold';'Gap length (ms)'})
    ylim([0 25])
    xlim([0 4])
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',14)
    grid off;
    box off
    
 
    
    
    
    
    
    
    
  
    
    
   