%% Sleep pattern: plots viglance state durations 
% Plots figures S4 A-F, 2G.

         

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


dosave = 0; 
doplot = 0;
d=1; % only fro derivation
conds = [1:3];


ders=strvcat('fro','occ');
mousename={'Heffalump';'Kanga';'Piglet'};




for anm = 1:length(mousename)
    
    % define input animal
    if mousename{anm}(1) == 'H' % Heffalump
         dates{1}={'070420_1604';'080420_1604'};
         dates{2}={'240620_1702';'250620_1703'};
         dates{3}={'101120_1734';'111120_1735'};
    elseif mousename{anm}(1) == 'K' % Kanga
         dates{1}={'160620_1700';'170620_1700'};
         dates{2}={'060720_1600';'070720_1600'};
         dates{3}={'301020_1556';'161120_1531'};
    elseif mousename{anm}(1) == 'P' % Piglet
         dates{1}={'290420_1815';'300420_1817'};
         dates{2}={'300620_1618';'010720_1618'};
         dates{3}={'221020_1719';'011120_1720'};
    end
    
for cond = conds
for dayt = 1:2
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
        
        
        % find episodes of different durations
        % find episodes shorter than 1 min
        minep1_idx1 = find( 0 <= el<15);
        mineps_epoch1 = el(minep_idx);
        mineps1 = mineps_epoch/15;
        % find episodes between 1-2 min.
        minep1_idx2 = find(15 <= el < 30);
        mineps_epoch2 = el(minep_idx);
        mineps2 = mineps_epoch/15;
        % find episodes between 2-5 min.
        minep1_idx3 = find(30 <= el < 75);
        mineps_epoch3 = el(minep_idx);
        mineps3 = mineps_epoch/15;
        
        % find episodes 
        clear mineps2
        clear EpNoH EpLH avEpLH stdEpLH errEpLH
        for duridx = [1 2 3 4 5 6 7 8 9 10 11 12 13]
%             duridx = 13
           
            epn1m = 15; % no of epochs in 1 min
            % epochs in 1 min windows, 10-30min, 30-60 min, >60mins
            durep1 = epn1m*[0 1 2 3 4 5 6 7 8 9 10 30 60];
            durep2 = epn1m*[1 2 3 4 5 6 7 8 9 10 30 60 60*24];
            
        minep_idx2 = find( (durep1(duridx) <= el) & (el < durep2(duridx)) );
        mineps_epoch2 = el(minep_idx2);
        mineps2{duridx} = mineps_epoch2/15;
        
        % save episode number
        EpNoH{anm}{cond}{dayt}(duridx) = length( mineps2{duridx});
        
        
        % save epiode length
        EpLH{anm}{cond}{dayt}{duridx} =  mineps2{duridx};
        avEpLH{anm}{cond}{dayt}(duridx) = mean( mineps2{duridx});
        stdEpLH{anm}{cond}{dayt}(duridx) = std( mineps2{duridx});
        errEpLH{anm}{cond}{dayt}(duridx) = std( mineps2{duridx})/sqrt(length( mineps2{duridx}));
        
        end
        
        % save all epsiode lengths
        Elengths{anm}{cond}{dayt} = el;
        
        
        
        % ---
        
        
        
        
        
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
    
    if doplot == 1
    hold on
    
    if d==1;
        ylabel('SWA frontal (% of mean)')
    else
        ylabel('SWA occipital (% of mean)')
    end
    
    xlabel('ZT (hours)')
    
    end

    
    
      


    
end

end

end


        
%% Figure 2G: Plot distribution of episode lengths (bar histogram, cond123)

dosave = 0

        anm  = anmID;

        % All wake epsiodes over 24 hrs
        Elengths{anm}{cond}{dayt};
        
        % All wake epsiodes over 48 hrs
        for cond = 1:3
        Elengths2d{anm}{cond} = [Elengths{anm}{cond}{1}; Elengths{anm}{cond}{2}];
        end
        
        close all
        cmap = [.5 .5 1; 1 .5 .5; 0.9 0.6 0.3];
        figure
        for cond = 1:3
        hold all        
%         histogram( Elengths2d{anm}{cond}/15,3)
        h = histogram( Elengths2d{anm}{cond}/15,[0 1 2 3 4 5 6 7 8 9 10 15 20 25 30 35 40 45 50])
%         h = histogram( Elengths2d{anm}{cond}/15,[0 1 2 3 4 5 6 7 8 9 10 15 20 25 50])

        h.FaceColor = cmap(cond,:);
        h.LineWidth = 1.5
        h.DisplayStyle = 'bar'; %'stairs' 'bar'
        h.EdgeColor = 'k';
        h.EdgeAlpha = 1
        h.FaceAlpha = .6
        
        %   Histogram binning:
        %   The value X(i) is in the kth bin if EDGES(k) <= X(i) < EDGES(k+1). The 
        %   last bin will also include the right edge such that it will contain X(i)
        %   if EDGES(end-1) <= X(i) <= EDGES(end). (I.e. each bin includes
        %   the left edge, the last bin also includes the right edge.)
        xlabel('Wake episode duration (min)')
        ylabel('Number of episodes')
        ylim([0 180])
        xlim([0 30])
        set(gca, 'linewidth',2.5, 'Fontsize',15)
        box off
        
        
       
        
        end
        
        legend('Baseline','Post 1','Post 2','Box','off')
        
        
        
%% Figure 2G (small panels): Plot count of all waking episodes
        
       dosave = 0
       
       for cond = 1:3
           WepNo(cond) = length(Elengths2d{anm}{cond})
       end
        
%        close all
       figure; hold all
       b1 = bar([1],WepNo(1),'FaceColor',cmap(1,:),'linewidth',1.5);
       b2 = bar([2],WepNo(2),'FaceColor',cmap(2,:),'linewidth',1.5);
       b3 = bar([3],WepNo(3),'FaceColor',cmap(3,:),'linewidth',1.5);
       
       
        % display values on top of bars
        Y = [WepNo(1) WepNo(2)  WepNo(3) ];
        Yround = round(Y,1); % round to one decimal
        if mousename{1}(1) == 'P'
        txt = text(1:3,Yround-30,num2str(Yround'),'vert','bottom','horiz','center'); 
        elseif mousename{1}(1) == 'K'
        txt = text(1:3,Yround-45,num2str(Yround'),'vert','bottom','horiz','center'); 
        elseif mousename{1}(1) == 'H'
        txt = text(1:3,Yround-25,num2str(Yround'),'vert','bottom','horiz','center'); 
        end
        set(txt,'fontsize',35)
        
        xlim([0 4])
%         ylim([0 140])
        set(gca, 'linewidth',2.5,'fontsize',30)
        ylabel('No. of wake ep.')
    %      yticks([])
        xticks([])
%         xticklabels({'BL','Post1','Post2'})
        box off
     

   
%% build SPSS matrix - episode duration

% STbinm{anm}{cond}{d}(state,b);

ID = []; TIME = []; DAY = []; STATE = []; BIN = []; Dur = [];

for anm = 1:length(mousename)
    for t = conds
        for day = 1:2
            for state = 1:5
                for bin = 1:24
                    ID = [ID; anm];
                    BIN = [BIN; bin];
                    STATE = [STATE; state];
                    DAY = [DAY; day];
                    TIME = [TIME; t];
                    % exlude bins with no (short) signal
                    if STbinm{anm}{t}{day}(10,bin) == 1
                        Dur = [Dur; nan];
                    elseif STbinm{anm}{t}{day}(10,bin) == 0
                        Dur = [Dur; STbinm{anm}{t}{day}(state,bin)];
                    end

                end
            end
        end
    end
end

STdurSPSS = [ID TIME DAY STATE BIN Dur];
STdur_table = array2table(STdurSPSS); % convert to table
STdur_table.Properties.VariableNames(1:6) = {'ID','TIME','DAY','STATE','BIN','Dur'}; % name the columns


%% Supp Figure S4D,E,F: state duration timecourse, BL only - 1h bins

% Run the previous section to produce input.

statn = {'Wake';'NREM';'REM ';'REM2';' M  '};

cmapST = [0 0.75, 0.75; 0.5, 0.5, 0.5; 0 0.5 0; 0.4660, 0.6740, 0.1880; 0.2 0.2 0.2];


% close all
clear idx
% 1h bins
for anm = anmID
  for state = 1:5
        for t = conds
            for day = 1:2
                for bidx = 1:24
                    % create idx for 2h-bin entries (consisting of hrs 1,2; hrs 2,3;... etc.)
                    idx{anm}{state}{day}{t}{bidx}(1) = find(STdurSPSS(:,1) == anm & STdurSPSS(:,2) == t & STdurSPSS(:,3) == day...
                        & STdurSPSS(:,4) == state & STdurSPSS(:,5) == bidx); 
                end
            end
        end
  end
end


% creat input for anm1 state 1
clear input
for anm = anmID
    for state = 1:5
        for t = conds
            input{anm}{state}{t} = [];
            for bidx = 1:24
                % 'duration' values per 2h-bin, for day 1 and day 2 together(4
                % values)
                add = [STdurSPSS(idx{anm}{state}{1}{t}{bidx},6); STdurSPSS(idx{anm}{state}{2}{t}{bidx},6)]; 
                input{anm}{state}{t} = [input{anm}{state}{t} add ];
            end
            % produce mean values & errors
            avinput{anm}{state}{t} = nanmean(input{anm}{state}{t});
            stdinput{anm}{state}{t} = nanstd(input{anm}{state}{t});
            errinput{anm}{state}{t} = stdinput{anm}{state}{t} / sqrt( size(input{anm}{state}{t},1) );

        end
    end
end

figure; hold all

for anm = anmID
    figure; hold all

    for state = 1:4
        
        % bar depicting dark phase from 20pm to 5AM (9hrs)
        plot([6 15],[75 75],'-','color',[0 0 0],'linewidth',5)
        plot([6 6],[-10 75],'--k','linewidth',.5)
        plot([15 15],[-10 75],'--k','linewidth',.5)
        
        % x-axis shift due to start time after 15:00
        xshift = emtybins{anm}(1);
        errorbar([1:24]+xshift,avinput{anm}{state}{1},errinput{anm}{state}{1},'color',cmapST(state,:),'linewidth',1.2);
        e(state) = plot([1:24]+xshift,avinput{anm}{state}{1},'color',cmapST(state,:),'linewidth',1.5);

        set(gca,'linewidth',1.5,'fontsize',14)
        xlabel('ZT')
        ylabel(['Minutes'])
%         title(statn{state})
        xlim([0 27])
        ylim([-10 80])
        yticks([0 10:10:60])
        xticks([1:2:27])

%          xticklabels({'2';'4';'6';'8';'10';'12';'14';'16';'18';'20';'22';'24'})
        xticklabels({'10';'12';'14';'16';'18';'20';'22';'0';'2';'4';'6';'8';'10'}) % fromm start time 15:00 (ZT10)
        
        

         


    end
end

        legend([e(1) e(2) e(3) e(4)], 'Wake','NREM','REM','REM2')
        

%% Figure S4A,B,C : Total duration for each state (BL only)

clear STdurnorm meanSTdurnrom stdSTdurnrom indSTdurnorm p1 p2

rectmin{anm}{cond}{dayt}; % total recording time per day

% recording time over 2 days (per cond and state)
for anm = anmID
    for cond = conds
        rect48{anm}{cond} = sum([rectmin{anm}{cond}{1} rectmin{anm}{cond}{2}]);
        % total recording time
    end
end

STbinm{anm}{cond}{dayt}; % lines 1-6 are state durations (in min.) for W, NR, R, R2, M, total time, in 1h bins

% close all
for anm =anmID
    for cond = conds
            for state = 1:5
                 for d = 1:2
                    STdur{anm}{cond}{d}(:,state) = sum(STbinm{anm}{cond}{d}(state,:)); %daily duration per state (in min.)
                 end
                 % total duration of each vigilance state (over 48 hours)
                 STdur48{anm}{cond}(state) = STdur{anm}{cond}{1}(:,state) + STdur{anm}{cond}{2}(:,state);
                 % as percent of total recording time
                 STdurnorm{anm}{cond}(:,state) = STdur48{anm}{cond}(state)/rect48{anm}{cond}*100; % as % or daily recording time
                 end
            end
end


% close all
figure;hold all
for anm =anmID
    for state = 1:5
%         plot([1 2 3], [STdurnorm{anm}{1}(:,state) STdurnorm{anm}{2}(:,state)  STdurnorm{anm}{3}(:,state) ],'-')
        b(state) = bar([state], [STdurnorm{anm}{1}(:,state)],'FaceColor',cmapST(state,:),'linewidth',1.5);
%         b2 = bar([2], [STdurnorm{anm}{2}(:,state)],'FaceColor',cmap(2,:),'linewidth',1.5);
%         b3 = bar([3], [STdurnorm{anm}{3}(:,state)],'FaceColor',cmap(3,:),'linewidth',1.5);

    end
end
        
        % display values on top of bars
        Y = [STdurnorm{anm}{1}(:,1) STdurnorm{anm}{1}(:,2) STdurnorm{anm}{1}(:,3) STdurnorm{anm}{1}(:,4) STdurnorm{anm}{1}(:,5)];
        Yround = round(Y,1); % round to one decimal
        txt = text(1:5,Yround,num2str(Yround'),'vert','bottom','horiz','center'); 
        set(txt,'fontsize',14)
        


     xlim([0 6])
     ylim([0 80])
   

     set(gca, 'linewidth',1.5,'fontsize',14)
     ylabel('% of recording time')
     xticks([1 2 3 4 5])
     xticklabels({'Wake','NREM','REM','REM2','M'})
%      title(statn{state}) 
     box off
     
        


             
 