%% EEG evoked potentials: barplots
% Plots Figs. 3E,F, S11, S16

clear all
close all

% to do: move loaded data into reposiutory folder, adapt loading path

% USER INPUT
anminput = 1;

% NOTE THE FOLLOWING:
% anmID = 1 corresponds to ferret 2.
% anmID = 2 corresponds to ferret 3.
% anmID = 3 corresponds to ferret 1.


loadpath = 'D:\PC-VV001-restoration\manuscript\Repository\EEG_evoked_activity_data\'

inputchannels = [1 2];

nint = 4           
if anminput == 3
condinput = [1 2 3]; % for anm 3 Post2 was recorded as well
else
    condinput = [1 2];
end

dosave = 0;

% FIXED INPUT
STORE = 'RawA';
iter = 20;



%%  Fig. S16: Summary figure - EvPot: chnage during wake and NREM

% load data (produced in previous section)
clear input idx_cond_freq
dosave = 0

% STORE = 'RawA';
iter = 20;


for anm = anminput %[1 2]
    if anm == 1; animal = 'Heffalump'; elseif anm == 2; animal = 'Kanga'; elseif anm == 3; animal = 'Piglet'; end
    for ch = inputchannels% [1 2]
            % load table
            if condinput(end) == 3
                if nint == 4
%                 input_table = readtable([loadpath animal '\SPSS\' animal '_EvPotIntAv_' STORE '_CH' num2str(ch) '_Bootsr' num2str(iter) '_Cond' num2str(condinput) '.xlsx']);
                
                input_table = readtable([loadpath animal '\' animal '_EvPotIntAv_' STORE '_CH' num2str(ch) '_Bootsr' num2str(iter) '_Cond' num2str(condinput) '.xlsx']);

                elseif nint == 40 | 50 | 60 | 65
                input_table = readtable([loadpath animal '\' animal '_EvPot_' STORE '_CH' num2str(ch) '_int' num2str(nint) '_Bootsr' num2str(iter) '_Cond' num2str(condinput) '.xlsx']);
                end
            elseif condinput(end) == 2
                if nint == 4
                input_table = readtable([loadpath animal '\' animal '_EvPotIntAv_' STORE '_CH' num2str(ch) '_Bootsr' num2str(iter) '_Cond' num2str(condinput) '.xlsx']);
                elseif nint == 40 | 50 | 60 | 65
                input_table = readtable([loadpath animal '\' animal '_EvPot_' STORE '_CH' num2str(ch) '_int' num2str(nint) '_Bootsr' num2str(iter) '_Cond' num2str(condinput) '.xlsx']);
                end
            end
            % convert into cell array
            input{anm}{ch} = table2array(input_table);
    end
end



% find indeces with window
for anm = anminput
for ch = inputchannels
for freq = [1 4 8 16]
    for wind = [1 2 3]
        for cond = 1:2
            if nint == 4
            idx_cond_freq_wind{anm}{ch}{cond}{freq}{wind} = find(input{anm}{ch}(:,2) == cond & input{anm}{ch}(:,4) == freq & input{anm}{ch}(:,7) == wind );
            else
            idx_cond_freq_wind{anm}{ch}{cond}{freq}{wind} = find(input{anm}{ch}(:,2) == cond & input{anm}{ch}(:,4) == freq & input{anm}{ch}(:,7) == wind );
            end
            end
    end
end
end
end

% for excel tables with specific intensity, EvPot is in column 9
if nint == 4; cidx = 8; else cidx = 8; end

clear p_mean p_std p_err
% produce figure input 
for anm = anminput
for ch = inputchannels
for freq = [1 4 8 16]
    for cond = 1:2
        for wind = [1 2 3]
        p_mean{anm}{ch}{freq}{wind}(cond) = mean(input{anm}{ch}(idx_cond_freq_wind{anm}{ch}{cond}{freq}{wind},cidx));
        p_std{anm}{ch}{freq}{wind}(cond) = std(input{anm}{ch}(idx_cond_freq_wind{anm}{ch}{cond}{freq}{wind},cidx));
        p_err{anm}{ch}{freq}{wind}(cond) = std(input{anm}{ch}(idx_cond_freq_wind{anm}{ch}{cond}{freq}{wind},cidx))/sqrt( length(input{anm}{ch}(idx_cond_freq_wind{anm}{ch}{cond}{freq}{wind},cidx)) );
        end
    end
end
end
end



% find indeces
clear idx_cond_freq_state
for anm = anminput
for ch = inputchannels
for freq = [1 4 8 16]
    for state = [1 2 3 4]
        for wind = [1 2 3]
            for cond = condinput
                if nint == 4
                idx_cond_freq_state{anm}{ch}{cond}{freq}{state}{wind} = find(input{anm}{ch}(:,2) == cond & input{anm}{ch}(:,4) == freq...
                    & input{anm}{ch}(:,3) == state & input{anm}{ch}(:,7) == wind);
                else
                idx_cond_freq_state{anm}{ch}{cond}{freq}{state}{wind} = find(input{anm}{ch}(:,2) == cond & input{anm}{ch}(:,4) == freq...
                    & input{anm}{ch}(:,3) == state & input{anm}{ch}(:,8) == wind);
                end
            end
        end
    end
end
end
end

clear p_mean p_std p_err
% produce figure input 
for anm = anminput
for ch = inputchannels
for freq = [1 4 8 16]
    for state = 1:4
        for wind = [1 2 3]
            for cond = condinput
                p_mean{anm}{ch}{freq}{state}{wind}(cond) = mean(input{anm}{ch}(idx_cond_freq_state{anm}{ch}{cond}{freq}{state}{wind},cidx));
                p_std{anm}{ch}{freq}{state}{wind}(cond) = std(input{anm}{ch}(idx_cond_freq_state{anm}{ch}{cond}{freq}{state}{wind},cidx));
                p_err{anm}{ch}{freq}{state}{wind}(cond) = std(input{anm}{ch}(idx_cond_freq_state{anm}{ch}{cond}{freq}{state}{wind},cidx))/sqrt( length(input{anm}{ch}(idx_cond_freq_state{anm}{ch}{cond}{freq}{state}{wind},cidx)) );
            end
        end
    end
end
end
end


dosave = 0
% close all
% produce figure input: ratio to cond2
clear p_ratio1m p_std_ratio1m p_err_ratio1m
for anm = anminput
for ch = inputchannels
for freq = [1 4 8 16]
    for state = 1:4
        for wind = [1 2 3]
            for cond = 1:2
                p_ratio1m{anm}{ch}{freq}{wind}(state) = p_mean{anm}{ch}{freq}{state}{wind}(2) / p_mean{anm}{ch}{freq}{state}{wind}(1);
                p_std_ratio1m{anm}{ch}{freq}{wind}(state) = p_std{anm}{ch}{freq}{state}{wind}(2) / p_std{anm}{ch}{freq}{state}{wind}(1);
                p_err_ratio1m{anm}{ch}{freq}{wind}(state) = p_err{anm}{ch}{freq}{state}{wind}(2) / p_err{anm}{ch}{freq}{state}{wind}(1);
            end
        end
    end
end
end
end

% produce mean across response windows
clear p_ratio p_std_ratio p_err_ratio
for ch = inputchannels
    for freq = [1 4 8 16]
        for state  = [1 2 3 4]
        p_ratio{anm}{ch}{freq}(state) = nanmean([ p_ratio1m{anm}{ch}{freq}{1}(state) p_ratio1m{anm}{ch}{freq}{2}(state)...
                                                  p_ratio1m{anm}{ch}{freq}{3}(state)]); % nanmean so empty response windows (NaN) are iognored
        p_std_ratio{anm}{ch}{freq}(state) = nanmean([ p_std_ratio1m{anm}{ch}{freq}{1}(state) p_std_ratio1m{anm}{ch}{freq}{2}(state)...
                                                  p_std_ratio1m{anm}{ch}{freq}{3}(state)]); % nanmean so empty response windows (NaN) are iognored
        p_err_ratio{anm}{ch}{freq}(state) = nanmean([ p_err_ratio1m{anm}{ch}{freq}{1}(state) p_err_ratio1m{anm}{ch}{freq}{2}(state)...
                                                  p_err_ratio1m{anm}{ch}{freq}{3}(state)]); % nanmean so empty response windows (NaN) are iognored
        end
    end
end

figure
Wevpot = [p_ratio{anm}{ch}{1}(1) p_ratio{anm}{ch}{4}(1) p_ratio{anm}{ch}{8}(1) p_ratio{anm}{ch}{16}(1)];
NRevpot = [p_ratio{anm}{ch}{1}(2) p_ratio{anm}{ch}{4}(2) p_ratio{anm}{ch}{8}(2) p_ratio{anm}{ch}{16}(2)];
Revpot = [p_ratio{anm}{ch}{1}(3) p_ratio{anm}{ch}{4}(3) p_ratio{anm}{ch}{8}(3) p_ratio{anm}{ch}{16}(3)];
R2evpot = [p_ratio{anm}{ch}{1}(4) p_ratio{anm}{ch}{4}(4) p_ratio{anm}{ch}{8}(4) p_ratio{anm}{ch}{16}(4)];

EvPot = [Wevpot; NRevpot; Revpot; R2evpot ]*100; % concatenate & convert to %

% CREATE HEATMAP
Freqs = {'1 kHz NBN','4 kHz NBN','8 kHz NBN','16 kHz NBN'};
States = {'W','NR','R','R2'};

                    if ch == 1; der = 'Fro'; elseif ch == 2; der = 'Occ'; end

h = heatmap(Freqs,States,EvPot);
h.Title = ['Post1: ' der ' EEG auditory evoked response (% of BL)'];
h.XLabel = 'Stimulus';
h.YLabel = 'Vigilance state';
set(h,'fontsize',12)

% set colourmap limits according to data range
if anm == 3 | anm == 1
caxis([(80), (150)]); % as mean of pos and meanm of neg data as colour range limits
elseif anm == 2
caxis([(50), (100)]); % as mean of pos and meanm of neg data as colour range limits
end





if condinput(end) == 3
    % produce figure input: ratio to cond3
    clear p_ratio2 p_std_ratio2 p_err_ratio2 p_ratio1m p_std_ratio1m p_err_ratio1m
    for anm = anminput
    for ch = inputchannels
    for freq = [1 4 8 16]
        for state = 1:4
            for wind = [1 2 3]
                for cond = 1:2
                    p_ratio1m{anm}{ch}{freq}{wind}(state) = p_mean{anm}{ch}{freq}{state}{wind}(3) / p_mean{anm}{ch}{freq}{state}{wind}(1);
                    p_std_ratio1m{anm}{ch}{freq}{wind}(state) = p_std{anm}{ch}{freq}{state}{wind}(3) / p_std{anm}{ch}{freq}{state}{wind}(1);
                    p_err_ratio1m{anm}{ch}{freq}{wind}(state) = p_err{anm}{ch}{freq}{state}{wind}(3) / p_err{anm}{ch}{freq}{state}{wind}(1);
                end
            end
        end
    end
    end
    end
    % produce mean across response windows
    clear p_ratio p_std_ratio p_err_ratio
    for ch = inputchannels
        for freq = [1 4 8 16]
            for state  = [1 2 3 4]
            p_ratio2{anm}{ch}{freq}(state) = nanmean([ p_ratio1m{anm}{ch}{freq}{1}(state) p_ratio1m{anm}{ch}{freq}{2}(state)...
                                                      p_ratio1m{anm}{ch}{freq}{3}(state)]); % nanmean so empty response windows (NaN) are iognored
            p_std_ratio2{anm}{ch}{freq}(state) = nanmean([ p_std_ratio1m{anm}{ch}{freq}{1}(state) p_std_ratio1m{anm}{ch}{freq}{2}(state)...
                                                      p_std_ratio1m{anm}{ch}{freq}{3}(state)]); % nanmean so empty response windows (NaN) are iognored
            p_err_ratio2{anm}{ch}{freq}(state) = nanmean([ p_err_ratio1m{anm}{ch}{freq}{1}(state) p_err_ratio1m{anm}{ch}{freq}{2}(state)...
                                                      p_err_ratio1m{anm}{ch}{freq}{3}(state)]); % nanmean so empty response windows (NaN) are iognored
            end
        end
    end

    figure
    Wevpot = [p_ratio2{anm}{ch}{1}(1) p_ratio2{anm}{ch}{4}(1) p_ratio2{anm}{ch}{8}(1) p_ratio2{anm}{ch}{16}(1)];
    NRevpot = [p_ratio2{anm}{ch}{1}(2) p_ratio2{anm}{ch}{4}(2) p_ratio2{anm}{ch}{8}(2) p_ratio2{anm}{ch}{16}(2)];
    Revpot = [p_ratio2{anm}{ch}{1}(3) p_ratio2{anm}{ch}{4}(3) p_ratio2{anm}{ch}{8}(3) p_ratio2{anm}{ch}{16}(3)];
    R2evpot = [p_ratio2{anm}{ch}{1}(4) p_ratio2{anm}{ch}{4}(4) p_ratio2{anm}{ch}{8}(4) p_ratio2{anm}{ch}{16}(4)];

    EvPot = [Wevpot; NRevpot; Revpot; R2evpot ]*100; % concatenate & convert to %

    % CREATE HEATMAP
    Freqs = {'1 kHz NBN','4 kHz NBN','8 kHz NBN','16 kHz NBN'};
    States = {'W','NR','R','R2'};
    
    h = heatmap(Freqs,States,EvPot);
    h.Title = ['Post2: ' der ' EEG auditory evoked response (% of BL)'];
    h.XLabel = 'Stimulus';
    h.YLabel = 'Vigilance state';
    set(h,'fontsize',12)

    
    % set colourmap limits according to data range
    if anm == 3
    caxis([(80), (150)]); % as mean of pos and meanm of neg data as colour range limits
    end
    
   

end        


%% Fig. 3E,F: Plot effect of state & cond



plotCI = 0; % 1 to plot CI, else STD


% find indeces
clear idx_state
for anm = anminput
for ch = inputchannels
    for state = [1 2 3 4]
        for cond = condinput
                idx_state{anm}{ch}{state}{cond} = find(input{anm}{ch}(:,3) == state & input{anm}{ch}(:,2) == cond) ;
        end
        end
end
end

% for excel tables with specific intensity, EvPot is in column 9
if nint == 4; cidx = 8; else cidx = 8; end

clear p_mean p_std p_err CI95 yCI95
% produce figure input 
for anm = anminput
for ch = inputchannels
for freq = [1 4 8 16]
    for state = 1:4
            for cond = condinput
                p_mean{anm}{ch}{cond}(state) = mean(input{anm}{ch}(idx_state{anm}{ch}{state}{cond},cidx));
                p_std{anm}{ch}{cond}(state) = std(input{anm}{ch}(idx_state{anm}{ch}{state}{cond},cidx));
                p_err{anm}{ch}{cond}(state) = std(input{anm}{ch}(idx_state{anm}{ch}{state}{cond},cidx))/sqrt( length(input{anm}{ch}(idx_state{anm}{ch}{state}{cond},cidx)) );
               
                N = size(input{anm}{ch}(idx_state{anm}{ch}{state}{cond},cidx),1);                                      % Number of ‘Experiments’ In Data Set
                CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
                yCI95{anm}{ch}{cond}(:,freq) = bsxfun(@times, p_err{anm}{ch}{cond}(state), CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

            
            
            
            end
    end
end
end
end

% plot
% close all
% condinput = [1:3]
cmap = [ 0 0 1; 1 0 0; 1 .5 0]
for anm = anminput
        for ch = inputchannels
            figure; hold all % new figure for each animal and state
            k = 0;

            for cond = condinput
                if animal(1) ~= 'P'
                if cond == 1; k = 0; elseif cond == 2; k = 0.4; end; % shift on x axis
                end
                if animal(1) == 'P' % for Piglet with 3 bars, position them closer together
                if cond == 2; k = 0.3; elseif cond == 3; k = 0.6; end; % shift on x axis
                end
                
                if ch == 1; der = 'Fro'; elseif ch == 2; der = 'Occ'; end % labels for derivations
                
                if plotCI == 1 % plot means and 95% CIs
                errorbar([1 2 3 4]-0.15+k,p_mean{anm}{ch}{cond}, yCI95{anm}{ch}{cond}(2,:),'o','Color',cmap(cond,:),...
                    'linewidth',2.5,'Markersize',12,'MarkerFaceColor',cmap(cond,:),'capsize',20)
                else
                    
                    if animal(1) ~= 'P'
                         p(cond) = errorbar([1 2 3 4]-0.2+k,p_mean{anm}{ch}{cond}, p_err{anm}{ch}{cond},'o','Color','k','linewidth',2,...
                        'Markersize',1,'MarkerFaceColor',cmap(cond,:))

                         p(cond) = bar([1 2 3 4]-0.2+k,p_mean{anm}{ch}{cond},.3,'FaceColor',cmap(cond,:),'LineWidth',1.5)
        %                 b1 = bar([1]+t,[mean(Thresholds(6,IdxTS{1}{freq})') ],1,'FaceColor',cmap(1,:),'LineWidth',1.5)
                    else
                         p(cond) = errorbar([1 2 3 4]-0.3+k,p_mean{anm}{ch}{cond}, p_err{anm}{ch}{cond},'o','Color','k','linewidth',2,...
                        'Markersize',1,'MarkerFaceColor',cmap(cond,:))

                         p(cond) = bar([1 2 3 4]-0.3+k,p_mean{anm}{ch}{cond},.25,'FaceColor',cmap(cond,:),'LineWidth',1.5)
                    end
                end
%                 set(gca, 'YScale', 'log')
                xlim([0 5]); ylim([1 8]*10^-5);
                if animal(1) == 'H'; ylim([0 7]*10^-5); if nint == 4; ylim([0 7]*10^-5); end; end
                if animal(1) == 'K'; ylim([0 8]*10^-5); end

                xticks([1 2 3 4])

                xticklabels(['Wake';'NREM';'REM ';'REM2'])
                
                ylabel([der ' AER (V)' ] ) 
                xlabel('Vigilance state')                 
                set(gca,'fontsize',18,'linewidth',1.5)
                box off
                title([animal(1:2) ' ' num2str(nint) 'dB' ])
                
                % sig markers
                msize = 5;
                madist = 0.09;
                if anm == 3; h1 = 5.8; h2 = 5.9; elseif anm == 2; h1 = 6.6; h2 = 6.75;
                elseif anm == 1; h1 = 5; h2 = 5.1; end
                
                % wake 
                % BL to Post1
                if anm == 3
                plot([1-0.3+0 1], [h1 h1]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([1-0.3+0 1])-madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([1-0.3+0 1]), [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([1-0.3+0 1])+madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                else
                plot([1-0.2+0 1+0.2], [h1 h1]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([1-0.2+0 1+0.2])-madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([1-0.2+0 1+0.2]), [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([1-0.2+0 1+0.2])+madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                end
                
                % BL to Post2
                if anm == 3 & ch == 1
                plot([1-0.3+0 1-0.15+0.44], [6.1 6.1]*10^-5,'-k','linewidth',1) % line post2
                plot(mean([1-0.3+0 1-0.15+0.44])-madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([1-0.3+0 1-0.15+0.44]), [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([1-0.3+0 1-0.15+0.44])+madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                elseif anm == 3 & ch == 2
                plot([1-0.3+0 1-0.15+0.44], [6.1 6.1]*10^-5,'-k','linewidth',1) % line post2
                plot(mean([1-0.3+0 1-0.15+0.44])-madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([1-0.3+0 1-0.15+0.44]), [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([1-0.3+0 1-0.15+0.44])+madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                end
                
               
                
                % NREM
                % BL to Post1
                if anm == 3
                plot([2-0.3+0 2], [h1 h1]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([2-0.3+0 2])-madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([2-0.3+0 2]), [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([2-0.3+0 2])+madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                else
                plot([2-0.2+0 2+0.2], [h1 h1]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([2-0.2+0 2+0.2])-madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([2-0.2+0 2+0.2]), [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([2-0.2+0 2+0.2])+madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                end
                
                % BL to Post2
                if anm == 3
                plot([2-0.3+0 2-0.15+0.44], [6.1 6.1]*10^-5,'-k','linewidth',1) % line post2
                plot(mean([2-0.3+0 2-0.15+0.44])-madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([2-0.3+0 2-0.15+0.44]), [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([2-0.3+0 2-0.15+0.44])+madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                end
                
                
                % REM
                % BL to Post1
                if anm == 3
                plot([3-0.3+0 3], [h1 h1]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([3-0.3+0 3])-madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([3-0.3+0 3]), [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([3-0.3+0 3])+madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                else
                plot([3-0.2+0 3+0.2], [h1 h1]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([3-0.2+0 3+0.2])-madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([3-0.2+0 3+0.2]), [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([3-0.2+0 3+0.2])+madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                end
                
                % BL to Post2
                if anm == 3
                plot([3-0.3+0 3-0.15+0.44], [6.1 6.1]*10^-5,'-k','linewidth',1) % line post2
                plot(mean([3-0.3+0 3-0.15+0.44])-madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([3-0.3+0 3-0.15+0.44]), [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([3-0.3+0 3-0.15+0.44])+madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                end
                
                % REM2
%                 % BL to Post1
                if anm == 3 & ch == 2 % for Piglet OCC (plot markers higher up)
                plot([4-0.3+0 4], [6.8 6.8]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([4-0.3+0 4])-madist, [6.9]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([4-0.3+0 4]), [6.9]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([4-0.3+0 4])+madist, [6.9]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                elseif anm == 3 & ch == 1 % for Piglet OCC (plot markers higher up)
                plot([4-0.3+0 4], [5.8 5.8]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([4-0.3+0 4])-madist, [5.9]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([4-0.3+0 4]), [5.9]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([4-0.3+0 4])+madist, [5.9]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                else
                plot([4-0.2+0 4+0.2], [h1 h1]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([4-0.2+0 4+0.2])-madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([4-0.2+0 4+0.2]), [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([4-0.2+0 4+0.2])+madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                end
                
                
               
                % BL to Post2
                if anm == 3 & ch == 2 % for Piglet OCC (plot markers higher up)
                plot([4-0.3+0 4-0.15+0.44], [7.1 7.1]*10^-5,'-k','linewidth',1) % line post2
                plot(mean([4-0.3+0 4-0.15+0.44])-madist, [7.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([4-0.3+0 4-0.15+0.44]), [7.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([4-0.3+0 4-0.15+0.44])+madist, [7.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                elseif anm == 3 & ch == 1
                plot([4-0.3+0 4-0.15+0.44], [6.1 6.1]*10^-5,'-k','linewidth',1) % line post2
                plot(mean([4-0.3+0 4-0.15+0.44])-madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([4-0.3+0 4-0.15+0.44]), [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([4-0.3+0 4-0.15+0.44])+madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                end
                
              
            end

            if anm == 3 & der(1) == 'F'
%             legend([p(1) p(2) p(3)],'BL','Post1','Post2','location','northeast','FontSize',12)
%             legend('boxoff')
            elseif anm ~= 3 & der(1) == 'F'
%                 legend([p(1) p(2)],'BL','Post1','location','northeast','FontSize',12)
%                 legend('boxoff')
            end
            
            
            
            
        end    
end



%% Fig. S11: Plot effect of frequency & cond
% Bars.


plotCI = 0; % 1 to plot CI, else STD
% dosave = 1
% dosave = 0

% find indeces
clear idx_freq
for anm = anminput
for ch = inputchannels
    for freq = [1 4 8 16]
        for cond = condinput
                idx_freq{anm}{ch}{freq}{cond} = find(input{anm}{ch}(:,2) == cond & input{anm}{ch}(:,4) == freq) ;
                
        end
        end
end
end

% for excel tables with specific intensity, EvPot is in column 9
if nint == 4; cidx = 8; else cidx = 8; end

clear p_mean p_std p_err yCI95 CI95
% produce figure input 
for anm = anminput
for ch = inputchannels
    freqc = 0;
for freq = [1 4 8 16]
    freqc = freqc+1;
            for cond = condinput
                p_mean{anm}{ch}{cond}(freqc) = mean(input{anm}{ch}(idx_freq{anm}{ch}{freq}{cond},cidx));
                p_std{anm}{ch}{cond}(freqc) = std(input{anm}{ch}(idx_freq{anm}{ch}{freq}{cond},cidx));
                p_err{anm}{ch}{cond}(freqc) = std(input{anm}{ch}(idx_freq{anm}{ch}{freq}{cond},cidx))/sqrt( length(input{anm}{ch}(idx_freq{anm}{ch}{freq}{cond},cidx)) );
                N = size(input{anm}{ch}(idx_freq{anm}{ch}{freq}{cond},cidx),1);                                      % Number of ‘Experiments’ In Data Set
                CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
                yCI95{anm}{ch}{cond}(:,freqc) = bsxfun(@times, p_err{anm}{ch}{cond}(freqc), CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

            
            
            end
end
end
end

% plot
% close all
cmap = [ 0 0 1; 1 0 0; 1 .5 0]
for anm = anminput
        k = 0;
        for ch = inputchannels
             figure; hold all % new figure for each animal and state
             k = 0;
            for cond = condinput
                if animal(1) ~= 'P' % for Piglet with 3 bars, position them closer together
                 if cond == 1; k = 0; elseif cond == 2; k = 0.4; end; % shift on x axis
                end
                if animal(1) == 'P' % for Piglet with 3 bars, position them closer together
                if cond == 2; k = 0.3; elseif cond == 3; k = 0.6; end; % shift on x axis
                end
                
                
                for freqc = 1:4
%                 if cond == 2; k = 0.22; elseif cond == 3; k = 0.44; end; % shift on x axis
                if ch == 1; der = 'Fro'; elseif ch == 2; der = 'Occ'; end
                if plotCI == 1 % plot means and 95% CIs
                errorbar([1 2 3 4]-0.15+k,p_mean{anm}{ch}{cond}, yCI95{anm}{ch}{cond}(2,:),'o','Color',cmap(cond,:),...
                    'linewidth',2.5,'Markersize',3,'MarkerFaceColor',cmap(cond,:),'capsize',20)
                else % plot means and errobars
                 if animal(1) ~= 'P'
                         p(cond) = errorbar([freqc]-0.2+k,p_mean{anm}{ch}{cond}(freqc), p_err{anm}{ch}{cond}(freqc),'o','Color','k','linewidth',2,...
                        'Markersize',1,'MarkerFaceColor',cmap(cond,:))

                         p(cond) = bar([freqc]-0.2+k,p_mean{anm}{ch}{cond}(freqc),.3,'FaceColor',cmap(cond,:),'LineWidth',1.5)
        %                 b1 = bar([1]+t,[mean(Thresholds(6,IdxTS{1}{freq})') ],1,'FaceColor',cmap(1,:),'LineWidth',1.5)
                    else
                         p(cond) = errorbar([freqc]-0.3+k,p_mean{anm}{ch}{cond}(freqc), p_err{anm}{ch}{cond}(freqc),'o','Color','k','linewidth',2,...
                        'Markersize',1,'MarkerFaceColor',cmap(cond,:))

                         p(cond) = bar([freqc]-0.3+k,p_mean{anm}{ch}{cond}(freqc),.25,'FaceColor',cmap(cond,:),'LineWidth',1.5)
                    end
                end
                    
%                 set(gca, 'YScale', 'log')
                xlim([0 5]); ylim([1 8]*10^-5);
                if animal(1) == 'H'; ylim([-0 7]*10^-5); if nint == 4; ylim([-0 7]*10^-5); end; end
                if animal(1) == 'K'; ylim([-0 8]*10^-5); end

                xticks([1 2 3 4])

                xticklabels(['1 kHz NBN ';'4 kHz NBN ';'8 kHz NBN ';'16 kHz NBN'])
                xtickangle(15)
                
                ylabel([der ' AER (V)' ] ) 
                xlabel('Stimulus') 
%                 grid on

                set(gca,'fontsize',14,'linewidth',1.5)
                box off
                title([animal(1:2) ' ' num2str(nint) 'dB' ])
                end
                
                 % sig markers
                msize = 5;
                madist = 0.09;
                if anm == 3; h1 = 5.8; h2 = 5.9; elseif anm == 2; h1 = 6.6; h2 = 6.75;
                elseif anm == 1; h1 = 5; h2 = 5.1; end
                
                % wake 
                % BL to Post1
                if anm == 3
                plot([1-0.3+0 1], [h1 h1]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([1-0.3+0 1])-madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([1-0.3+0 1]), [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([1-0.3+0 1])+madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                else
                plot([1-0.2+0 1+0.2], [h1 h1]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([1-0.2+0 1+0.2])-madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([1-0.2+0 1+0.2]), [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([1-0.2+0 1+0.2])+madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                end
                
                % BL to Post2
                if anm == 3 & ch == 1
                plot([1-0.3+0 1-0.15+0.44], [6.1 6.1]*10^-5,'-k','linewidth',1) % line post2
                plot(mean([1-0.3+0 1-0.15+0.44])-madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([1-0.3+0 1-0.15+0.44]), [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([1-0.3+0 1-0.15+0.44])+madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                elseif anm == 3 & ch == 2
                plot([1-0.3+0 1-0.15+0.44], [6.1 6.1]*10^-5,'-k','linewidth',1) % line post2
                plot(mean([1-0.3+0 1-0.15+0.44])-madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([1-0.3+0 1-0.15+0.44]), [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([1-0.3+0 1-0.15+0.44])+madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                end
                
               
                
                % NREM
                % BL to Post1
                if anm == 3
                plot([2-0.3+0 2], [h1 h1]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([2-0.3+0 2])-madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([2-0.3+0 2]), [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([2-0.3+0 2])+madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                else
                plot([2-0.2+0 2+0.2], [h1 h1]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([2-0.2+0 2+0.2])-madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([2-0.2+0 2+0.2]), [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([2-0.2+0 2+0.2])+madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                end
                
                % BL to Post2
                if anm == 3
                plot([2-0.3+0 2-0.15+0.44], [6.1 6.1]*10^-5,'-k','linewidth',1) % line post2
                plot(mean([2-0.3+0 2-0.15+0.44])-madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([2-0.3+0 2-0.15+0.44]), [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([2-0.3+0 2-0.15+0.44])+madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                end
                
                
                % REM
                % BL to Post1
                if anm == 3
                plot([3-0.3+0 3], [h1 h1]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([3-0.3+0 3])-madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([3-0.3+0 3]), [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([3-0.3+0 3])+madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                else
                plot([3-0.2+0 3+0.2], [h1 h1]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([3-0.2+0 3+0.2])-madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([3-0.2+0 3+0.2]), [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([3-0.2+0 3+0.2])+madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                end
                
                % BL to Post2
                if anm == 3
                plot([3-0.3+0 3-0.15+0.44], [6.1 6.1]*10^-5,'-k','linewidth',1) % line post2
                plot(mean([3-0.3+0 3-0.15+0.44])-madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([3-0.3+0 3-0.15+0.44]), [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([3-0.3+0 3-0.15+0.44])+madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                end
                
                % REM2
%                 % BL to Post1
                if anm == 3 & ch == 2 % for Piglet OCC (plot markers higher up)
                plot([4-0.3+0 4], [6.8 6.8]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([4-0.3+0 4])-madist, [6.9]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([4-0.3+0 4]), [6.9]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([4-0.3+0 4])+madist, [6.9]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                elseif anm == 3 & ch == 1 % for Piglet OCC (plot markers higher up)
                plot([4-0.3+0 4], [5.8 5.8]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([4-0.3+0 4])-madist, [5.9]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([4-0.3+0 4]), [5.9]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([4-0.3+0 4])+madist, [5.9]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                else
                plot([4-0.2+0 4+0.2], [h1 h1]*10^-5,'-k','linewidth',1) % line post1
                plot(mean([4-0.2+0 4+0.2])-madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([4-0.2+0 4+0.2]), [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([4-0.2+0 4+0.2])+madist, [h2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                end
                
                
               
                % BL to Post2
                if anm == 3 & ch == 2 % for Piglet OCC (plot markers higher up)
                plot([4-0.3+0 4-0.15+0.44], [7.1 7.1]*10^-5,'-k','linewidth',1) % line post2
                plot(mean([4-0.3+0 4-0.15+0.44])-madist, [7.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([4-0.3+0 4-0.15+0.44]), [7.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([4-0.3+0 4-0.15+0.44])+madist, [7.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                elseif anm == 3 & ch == 1
                plot([4-0.3+0 4-0.15+0.44], [6.1 6.1]*10^-5,'-k','linewidth',1) % line post2
                plot(mean([4-0.3+0 4-0.15+0.44])-madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker1
                plot(mean([4-0.3+0 4-0.15+0.44]), [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker2
                plot(mean([4-0.3+0 4-0.15+0.44])+madist, [6.2]*10^-5, '*k','Markersize',msize,'linewidth',1.1) % marker3
                end
                
              
            end
                
                
                
                 
                
            end
        end    
        

