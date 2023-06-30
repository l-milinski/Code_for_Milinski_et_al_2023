% Plot vigilance state durations in baseline condtions: sound vs no sound
% Plots Fig. 3D.
         

clear all
close all

% -----------
% USER INPUT:

anmID = 1; % set this as 1,2 or 3 depending on the animal

% NOTE THE FOLLOWING:
% anmID = 1 corresponds to ferret 2.
% anmID = 2 corresponds to ferret 3.
% anmID = 3 corresponds to ferret 1.


% add the correct path to data
path = 'D:\PC-VV001-restoration\manuscript\Repository\Sleep_architecture_data\';
% -----------



dosave = 0; % 1 to save for 1 animal & BL plots (except episode length&duration)
            % 2 to save BL plots with animals 1:3 (only episode length&duration)
doplot = 0; % to plot hypnograms

d=1; % only fro derivation
conds = [1:2];


ders=strvcat('fro','occ');

% pathfig=[path,'hypnogram']; mkdir(pathfig)
mousename={'Heffalump';'Kanga';'Piglet'};
% mousename={'Heffalump'};

for anm = anmID
    
    
     % define input animal
    if mousename{anm}(1) == 'H' % Heffalump
         dates{1}={'080420_1604'}; % undisturbed day (BL)
         dates{2}={'150520_1938'}; % sound stim day (BL)
    elseif mousename{anm}(1) == 'K' % Kanga
         dates{1}={'170620_1700'}; % undisturbed day (BL)
         dates{2}={'220520_1747'}; % sound stim day (BL)
    elseif mousename{anm}(1) == 'P' % Piglet
         dates{1}={'300420_1817'}; % undisturbed day (BL)
         dates{2}={'080520_1915'}; % sound stim day (BL)
    end
    
    
    
for cond = conds
for dayt = 1
%     dayname={'070420_Matt';'080420'};
    day = dates{cond}{dayt};

ders=strvcat('fro','occ');
f=0:0.25:20;
x=1:1:21599; %for 24h recordings
zermat=zeros(1,21599);
x=x/900;
cols='br';
dernames=strvcat('Frontal','Occipital');

h=0.5:1:24;

pathin=[path];
% figure(1)
% title(['hypnogram ' mousename day])

SWA=[];
    
%     figure
    der=ders(d,:);
    fn=[mousename{anm} ,'_',day,'_',der,'_VSspec'];
    clear nr w r w1 nr2 r2
    eval(['load ',pathin,fn,'.mat ma spectr w nr r w1 nr2 r2 mt -mat']);
%     load([pathin,fn,'.mat']) %load variables from VSspec.mat file
    if size(nr,1)==1 nr=nr'; w=w'; r=r'; w1=w1'; nr2=nr2'; r2=r2'; end
%     W=zermat; W([w;w1])=1;
%     N=zermat; N([nr;nr2])=1; 
%     R=zermat; R([r;r3])=1;


 % find start 1h-bin opf current recording
        strth = str2num(day(8:9));
        startm = str2num(day(10:11));
        % difference of start time from 15:00 (the earliest start time)
        % (ZT0=5am, lights on)
        diffh = strth-15;
        emtybins{anm}(cond) = diffh;

    
    W=zermat; W([w])=1;
    N=zermat; N([nr])=1; 
    R=zermat; R([r])=1;
    R2=zermat; R2([r2])=1;
    
    WA=zermat; WA([w1])=1;
    M=zermat; M([mt])=1;
    NA=zermat; NA([nr2])=1;
    
    if day(1:6)  == '111120' % for Heffalum111120 last 2h of rec are flat, replace NA there with 0s
            NA(19670:end) = zeros(1,length(NA(19670:end)));
    end
    
    % for Kanga_161120 treat M as NA and vice versa
    if day(1:6)  == '161120'
            M=zermat; M([nr2])=1;
            display('NA treated as M')
    end

    
    % save total recording time per day
    rectep = sum([length(w) length(nr) length(r) length(r2) length(w1) length(mt)]);
    rectmin{anm}{cond}{dayt} = rectep*4/60;
    
   
    % --- EDIT SCORING ---
    % Exclude REM2 during W (replace it by W)
    R2Wc = 0;
    for i = 2:length(W)
        if R2(i) == 1 && W(i-1) == 1
            R2Wc = R2Wc+1;
            R2(i) = 0;
            W(i) = 1;
        end     
    end
    % Exclude REM during W (replace it by W)
    RWc = 0;
    for i = 2:length(W)
        if R(i) == 1 && W(i-1) == 1
            RWc = RWc+1;
            R(i) = 0;
            W(i) = 1;
        end     
    end
    % FOR STATE DURATIONS
    % Replace NA with N (if last non-NA epoch was N)

    NNAc = 0;
    for i = 2:length(N)
        if NA(i) == 1 && N(i-1) == 1 
            N(i) = 1;
        end     
    end
    % Replace NA with R (if last non-NA epoch was R)
    RNAc = 0;
    for i = 2:length(N)
        if NA(i) == 1 && R(i-1) == 1 
            R(i) = 1;
        end     
    end
    % Replace NA with R2 (if last non-NA epoch was R2)
    RNAc = 0;
    for i = 2:length(N)
        if NA(i) == 1 && R2(i-1) == 1
            R2(i) = 1;
        end     
    end
    
    % --- --- --- --- ---
    
    % --- wake epsiodes ---
    % merge W and WA
    for nep = 1:length(W)
        if W(nep) == 1 || WA(nep) == 1
        Wep(nep) = 1;
        else
        Wep(nep) = 0;
        end
    end
       
       % code to find consecutive 1s is from here:
       % https://uk.mathworks.com/matlabcentral/answers/114852-finding-consecutive-true-values-in-a-vector
          %test
%         Wep = [ 0 0 1 1 0 0 1 1 1 ]
        a=Wep';
        a0 = [a; 0];

        % Find the end of any consecutive 1's in a0
        ii= strfind(a0',[1 0]);

        a1 = cumsum(a);

        % Cumulative sum at end of any consecutive 1's in a0
        i1 = a1(ii);

        % Places the amount to subtract during cumulative-sum 1-element past the
        % consecutive 1's in a, to produce only the cumulative sum of consecutive
        % 1's in a0.  If this is confusing, output a0 after this step.
        a0(ii+1) = -[i1(1);diff(i1)];
        a0;
        
        % save episode lengths
        idx = find(a0 < 0);
        el = a0(idx)*-1;
        % episode start and end
        eidx2 = idx-1; % last epoch of episode
        eidx1 = idx-el; % first epoch of episode
        
        % ---
        
        % find episodes longer than 1 min
        minep_idx = find(el>15);
        mineps_epoch = el(minep_idx);
        mineps = mineps_epoch/15;
        
        % count episodes longer than 1 min
        clear Epnobin_idx
        eph = (60*60)/4; % epochs per h
        eminidx1 = eidx1(minep_idx); % first epoch of episodes longer than 1 min
        for b = 1:24
                Epnobin_idx{b} = find( (eph*b)-eph <= eminidx1 & eminidx1 < eph*b);
                Epnobin_sum(b) = length(Epnobin_idx{b});
        end
        
        % save epidose number (1h bins)
        EpNobin{anm}{cond}{dayt} = Epnobin_sum;
        
        % save episode number
        EpNo{anm}{cond}{dayt} = length(mineps);
        % save epiode length
        EpL{anm}{cond}{dayt} = mineps;
        avEpL{anm}{cond}{dayt} = mean(mineps);
        stdEpL{anm}{cond}{dayt} = std(mineps);
        errEpL{anm}{cond}{dayt} = std(mineps)/sqrt(length(mineps));
        
         % --- brief awakenings (<16s interruptions of sleep) ---
    % merge W and WA and M
    for nep = 1:length(W)
        if W(nep) == 1 || WA(nep) == 1 || M(nep) == 1
        Bep(nep) = 1;
        else
        Bep(nep) = 0;
        end
    end
       
       % code to find consecutive 1s is from here:
       % https://uk.mathworks.com/matlabcentral/answers/114852-finding-consecutive-true-values-in-a-vector
          %test
%         Wep = [ 0 0 1 1 0 0 1 1 1 ]
        a=Bep';
        a0 = [a; 0];

        % Find the end of any consecutive 1's in a0
        ii= strfind(a0',[1 0]);

        a1 = cumsum(a);

        % Cumulative sum at end of any consecutive 1's in a0
        i1 = a1(ii);

        % Places the amount to subtract during cumulative-sum 1-element past the
        % consecutive 1's in a, to produce only the cumulative sum of consecutive
        % 1's in a0.  If this is confusing, output a0 after this step.
        a0(ii+1) = -[i1(1);diff(i1)];
        a0;
        
        % save episode lengths
        idx = find(a0 < 0);
        el = a0(idx)*-1;
        % episode start and end
        eidx2 = idx-1; % last epoch of episode
        eidx1 = idx-el; % first epoch of episode
        
        % ---
        
        
        % find episodes shorter than 16s (4 epochs)
        minep_idx = find(el<4);
        mineps_epoch = el(minep_idx);
        mineps = mineps_epoch*4; % brief aw duration in sec
        
        % count episodes shorter than 16s
        clear Epnobin_idx eminidx1
        eph = (60*60)/4; % epochs per h
        eminidx1 = eidx1(minep_idx); % first epoch of episodes longer than 1 min
        for b = 1:24
                BEpnobin_idx{b} = find( (eph*b)-eph <= eminidx1 & eminidx1 < eph*b);
                BEpnobin_sum(b) = length(BEpnobin_idx{b});
        end
        
        % save epidose number (1h bins)
        BEpNobin{anm}{cond}{dayt} = BEpnobin_sum;
        
        % save episode number
        BEpNo{anm}{cond}{dayt} = length(mineps);
        % save epiode length
        BEpL{anm}{cond}{dayt} = mineps;
        avBEpL{anm}{cond}{dayt} = mean(mineps);
        stdBEpL{anm}{cond}{dayt} = std(mineps);
        errBEpL{anm}{cond}{dayt} = std(mineps)/sqrt(length(mineps));
        
        
        
        % ----- sleep episodes
        
        % merge all sleep epochs (ecverything not W or Wa)
    for nep = 1:length(W)
        if W(nep) ~= 1 & WA(nep) ~= 1
        SEp(nep) = 1;
        else
        SEp(nep) = 0;
        end
    end
       
       % code to find consecutive 1s is from here:
       % https://uk.mathworks.com/matlabcentral/answers/114852-finding-consecutive-true-values-in-a-vector
          %test
%         Wep = [ 0 0 1 1 0 0 1 1 1 ]
        a=SEp';
        a0 = [a; 0];

        % Find the end of any consecutive 1's in a0
        ii= strfind(a0',[1 0]);

        a1 = cumsum(a);

        % Cumulative sum at end of any consecutive 1's in a0
        i1 = a1(ii);

        % Places the amount to subtract during cumulative-sum 1-element past the
        % consecutive 1's in a, to produce only the cumulative sum of consecutive
        % 1's in a0.  If this is confusing, output a0 after this step.
        a0(ii+1) = -[i1(1);diff(i1)];
        a0;
        
        % save episode lengths
        idx = find(a0 < 0);
        el = a0(idx)*-1;
        % episode start and end
        eidx2 = idx-1; % last epoch of episode
        eidx1 = idx-el; % first epoch of episode
        
        % ---
        
        % find episodes longer than 1 min
        minep_idx = find(el>15);
        mineps_epoch = el(minep_idx);
        mineps = mineps_epoch/15;
        
        % count episodes longer than 1 min
        clear Epnobin_idx
        eph = (60*60)/4; % epochs per h
        eminidx1 = eidx1(minep_idx); % first epoch of episodes longer than 1 min
        for b = 1:24
                Epnobin_idx{b} = find( (eph*b)-eph <= eminidx1 & eminidx1 < eph*b);
                Epnobin_sum(b) = length(Epnobin_idx{b});
        end
        
        % save epidose number (1h bins)
        SEpNobin{anm}{cond}{dayt} = Epnobin_sum;
        
        % save episode number
        SEpNo{anm}{cond}{dayt} = length(mineps);
        % save epiode length
        SEpL{anm}{cond}{dayt} = mineps;
        avSEpL{anm}{cond}{dayt} = mean(mineps);
        stdSEpL{anm}{cond}{dayt} = std(mineps);
        errSEpL{anm}{cond}{dayt} = std(mineps)/sqrt(length(mineps));
        
        
        
        
        
    % amount of epochs per 1h bin
    clear STbin
    eph = (60*60)/4; % epochs per h
    currep = 0;
    for b = 1:24
        if b == 1
            STbin(1,b) = sum(Wep(1:eph*b));
            STbin(2,b) = sum(N(1:eph*b));
            STbin(3,b) = sum(R(1:eph*b));
            STbin(4,b) = sum(R2(1:eph*b));
            STbin(5,b) = sum(M(1:eph*b));
            STbin(7,b) = sum(NA(1:eph*b)); % N artefacts


        elseif 1 < b && b < 24
            STbin(1,b) = sum(Wep(currep:eph*b));
            STbin(2,b) = sum(N(currep:eph*b));
            STbin(3,b) = sum(R(currep:eph*b));
            STbin(4,b) = sum(R2(currep:eph*b));
            STbin(5,b) = sum(M(currep:eph*b));
            STbin(7,b) = sum(NA(currep:eph*b));

        elseif b == 24
            % to account for 21599 epochs at end of recording (instead of
            % 21600)
            STbin(1,b) = sum(Wep(currep:(eph*b)-1));
            STbin(2,b) = sum(N(currep:(eph*b)-1));
            STbin(3,b) = sum(R(currep:(eph*b)-1));
            STbin(4,b) = sum(R2(currep:(eph*b)-1));
            STbin(5,b) = sum(M(currep:(eph*b)-1));
            STbin(7,b) = sum(NA(currep:(eph*b)-1));

        end
        currep = currep+eph; % increase epoch counter
            STbin(6,b) = sum([STbin(1,b) STbin(2,b) STbin(3,b) STbin(4,b) STbin(5,b)]);

    end
    
    % EXCLUDE BINS WITH MANY ARTEFACTS
    % If more than 50% of a 1h bin is NA, exlude this bin (replace with
    % NANs)
    for b = 1:24
        if STbin(7,b) > 450
            STbin(1:6,b) = [nan];
        end
    end
    
    % convert into min
    STbinm{anm}{cond}{dayt} = STbin*4/60;
    
    % mark STbinm where it is empty for all states (less than 30 min.)
    emptyIDX = find( sum(STbinm{anm}{cond}{dayt}) < 30 );
    STbinm{anm}{cond}{dayt}(10,emptyIDX) = 1;
    
    
    

     
    art=[w1;nr2;mt];
    
    %Select FFT frequency range
    swa=mean(spectr(:,3:17),2); % 0.5-4 Hz EEG power
    spindle = mean(spectr(:,61:75),2); % 15-20 Hz EEG power
    
    %Select FFT power for each state; replace zeros with NaNs
    swaW=swa; swaW(W==0)=NaN;
    swaN=swa; swaN(N==0)=NaN;
    swaR=swa; swaR(R==0)=NaN;
    swaR2=swa; swaR2(R2==0)=NaN;
    
    sW=spindle; sW(W==0)=NaN;
    sN=spindle; sN(N==0)=NaN;
    sR=spindle; sR(R==0)=NaN;
    sR2=spindle; sR2(R2==0)=NaN;

    
    
    
    % PLOT
    %organise state dependent FFT values in 4 columns
    swa=[swaW swaN swaR swaR2];
    spindle=[sW sN sR sR2];

    if doplot == 1
        figure
        hold all
        xshift = emtybins{anm}(cond); % in hours
        xshift_ep = xshift*60*15;
        
        % bar depicting dark phase from 20pm to 5AM (9hrs)
        plot([5*899 14*899],[299 299],'-','color',[0 0 0],'linewidth',5)
        plot([5*899 5*899],[-10 299],'--k','linewidth',1)
        plot([14*899 14*899],[-10 299],'--k','linewidth',1)

%     plot(swa(1:899*24,:),'O','linewidth',2,'markersize',1) %plot first 6 hours
        plot([1:length(swa(1:899*24,1))]+xshift_ep,swa(1:899*24,1),'rO','linewidth',2,'markersize',.1,'Color',[0 0.75, 0.75]); % WAKE
        plot([1:length(swa(1:899*24,2))]+xshift_ep,swa(1:899*24,2),'bO','linewidth',2,'markersize',.1,'Color',[0.5, 0.5, 0.5]); % NREM
        plot([1:length(swa(1:899*24,3))]+xshift_ep,swa(1:899*24,3),'gO','linewidth',2,'markersize',.1,'Color',[0 0.5 0]) % REM
        plot([1:length(swa(1:899*24,4))]+xshift_ep,swa(1:899*24,4),'cO','linewidth',2,'markersize',.1,'Color',[0.4660, 0.6740, 0.1880]) % REM2

               
%         % attempt to enlarge legend marker size
%         % https://uk.mathworks.com/matlabcentral/answers/517575-how-to-change-the-line-thickness-inside-the-legend-box#answer_425788
%         ax = axes(); 
%         hLeg = copyobj(h, ax); 
%         set(hLeg, 'XData', NaN, 'YData', NaN, 'LineWidth', 1.5)
%         legend(hLeg)

         %%% define plot dimensions %%%
%          xticks([2:2:26])

         set(gca,'XTick',[1:899:899*26])
         set(gca,'XTick',[1:899*2:899*27])

%          set(gca,'XTickLabel',[0:2:24])
         set(gca,'Linewidth',1.5,'fontsize',14)
         xlim([-500 899*27]) %plot first 9 hours
         ylim([0 300])
         
%           xticks([1:2:27])
%           xticklabels({'10';'12';'14';'16';'18';'20';'22';'0';'2';'4';'6';'8';'10'}) % fromm start time 15:00 (ZT10)
         
         xticklabels({'10';'12';'14';'16';'18';'20';'22';'0';'2';'4';'6';'8';'10'}) % fromm start time 15:00 (ZT10)

         
         % legend marker plots (otuside the displayed plot range)
         p1 = plot(1,600,'rO','linewidth',2,'markersize',5,'Color',[0 0.75, 0.75],'markerfacecolor',[0 0.75, 0.75]); % WAKE
         p2 = plot(1,600,'bO','linewidth',2,'markersize',5,'Color',[0.5, 0.5, 0.5],'markerfacecolor',[0.5, 0.5, 0.5]); % NREM
         p3 = plot(1,600,'gO','linewidth',2,'markersize',5,'Color',[0 0.5 0],'markerfacecolor',[0 0.5 0]); % REM
         p4 = plot(1,600,'cO','linewidth',2,'markersize',5,'Color',[0.4660, 0.6740, 0.1880],'markerfacecolor',[0.4660, 0.6740, 0.1880]); % REM2
         % set legend
         if dayt == 1
         legend([p1 p2 p3 p4],'Wake','NREM','REM','REM2')
         end
         

         
         box off
         
         
         
         
    end
    
    hold on
    
    if d==1;
        ylabel('SWA frontal (% of mean)')
    else
        ylabel('SWA occipital (% of mean)')
    end
    
    xlabel('ZT (hours)')
%     axis([0 24 0 800])
%     title([mousename{anm} day]) 
    
  

        
    
end

end

end

%% MAIN FIGURE 1: Total duration for each state
% to do: show as percentager of total recording time, to account of periods
% of flat signal (e.g. Heffalump111120). Use same input as for timecourse
% of states (see section above).

plotallstates = 0; % 1 to plot for all vigilance states, 0 for only total sleep amount


clear STdurnorm meanSTdurnrom stdSTdurnrom indSTdurnorm p1 p2

statn = {'Wake';'NREM';'REM';'REM2';'M'};
cmap = [.5 .5 1; 1 .5 .5; 0.9 0.6 0.3];
cmap = [ .7 .7 1; .5 .5 1]

rectmin{anm}{cond}{dayt}; % total recording time per day

% recording time over 1 day (per cond and state)
for anm = anmID
    for cond = conds
        rect48{anm}{cond} = sum([rectmin{anm}{cond}{1}]);
        % total recording time
    end
end

STbinm{anm}{cond}{dayt}; % lines 1-6 are state durations (in min.) for W, NR, R, R2, M, total time, in 1h bins

close all
for anm =anmID
    for cond = conds
            for state = 1:5
                 for d = 1
                    STdur{anm}{cond}{d}(:,state) = nansum(STbinm{anm}{cond}{d}(state,:)); %daily duration per state (in min.)
                 end
                 % total duration of each vigilance state (over 48 hours)
                 STdur48{anm}{cond}(state) = STdur{anm}{cond}{1}(:,state) ;
                 % as percent of total recording time
                 STdurnorm{anm}{cond}(:,state) = STdur48{anm}{cond}(state)/rect48{anm}{cond}*100; % as % or daily recording time
                 end
            end
end

% SLEEP DURATION
for anm =anmID
    for state = 1
        figure;hold all
%         plot([1 2 3], [STdurnorm{anm}{1}(:,state) STdurnorm{anm}{2}(:,state)  STdurnorm{anm}{3}(:,state) ],'-')
        b1 = bar([1], [100-STdurnorm{anm}{2}(:,state)],'FaceColor',cmap(2,:),'linewidth',1.5);
        b2 = bar([2], [100-STdurnorm{anm}{1}(:,state)],'FaceColor',[.8 .8 .8],'linewidth',1.5);
%         b3 = bar([3], [STdurnorm{anm}{3}(:,state)],'FaceColor',cmap(3,:),'linewidth',1.5);
        
        % display values on top of bars
        Y = [100-STdurnorm{anm}{2}(:,state) 100-STdurnorm{anm}{1}(:,state)];
        Yround = round(Y,1); % round to one decimal
        txt = text(1:2,Yround,num2str(Yround'),'vert','bottom','horiz','center'); 
        set(txt,'fontsize',16)
        
     xlim([0 3])
     ylim([0 100])
      

     set(gca, 'linewidth',1.5,'fontsize',18)
     ylabel('Sleep (% of recording time)')
     xticks([1 2])
     xticklabels({'Sound','No sound'})
%      title('Sleep amount') 
     box off
     
      %      save figure
     savepath = 'D:\PC-VV001-restoration\Ferrets\Ephys\figures\Chapter 4\'
     savepath = 'D:\PC-VV001-restoration\manuscript\AllFigures\Revision_May23\'
        if dosave == 1
            savename = [savepath mousename{anm} '_SleepDuration_NoSoundVSSound' statn{state} '_' der '.tiff'];
        print('-r750','-dtiff',savename);
        end
     
    end
end

if plotallstates == 1
for anm =anmID
    for state = 1:5
        figure;hold all
%         plot([1 2 3], [STdurnorm{anm}{1}(:,state) STdurnorm{anm}{2}(:,state)  STdurnorm{anm}{3}(:,state) ],'-')
        b1 = bar([1], [STdurnorm{anm}{2}(:,state)],'FaceColor',cmap(2,:),'linewidth',1.5);
        b2 = bar([2], [STdurnorm{anm}{1}(:,state)],'FaceColor',[.8 .8 .8],'linewidth',1.5);
%         b3 = bar([3], [STdurnorm{anm}{3}(:,state)],'FaceColor',cmap(3,:),'linewidth',1.5);
        
        % display values on top of bars
        Y = [STdurnorm{anm}{2}(:,state) STdurnorm{anm}{1}(:,state)];
        Yround = round(Y,1); % round to one decimal
        txt = text(1:2,Yround,num2str(Yround'),'vert','bottom','horiz','center'); 
        set(txt,'fontsize',14)
        


     xlim([0 3])
     if state == 1 ||  state == 2
     ylim([0 80])
     else 
     ylim([0 40])
     end

     set(gca, 'linewidth',1.5,'fontsize',14)
     ylabel('% of recording time)')
     xticks([1 2])
%      xticklabels({'BL','Post1','Post2'})
     xticklabels({'Sound','No sound'})

     title(statn{state}) 
     box off
     
   
     
    end
end
end




