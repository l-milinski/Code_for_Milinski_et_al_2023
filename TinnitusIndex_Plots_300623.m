%% Calculates Tinnitus index (TI) and chnages in ABR
% Plots Figs. S3D and 2F.


clear all

anmnames = {'Heffalump';'Kanga';'Piglet'};
savefig = 1;

ID = 1; % Heffalump,Kanga, Piglet

% NOTE THE FOLLOWING:
% ID = 1 corresponds to ferret 2.
% ID = 2 corresponds to ferret 3.
% ID = 3 corresponds to ferret 1.


% --- load TI input for silence detection
load D:\PC-VV001-restoration\manuscript\Repository\Tinnitus_Index_metadata\TI_silence

TI_silence1 = TI1(end-2:end);
TI_silence2 = TI2(end-2:end);
clear TI1 TI2

% INFO on TI: The metric here refers to the change in selence detection
% ability (percent correct in silence trials) after NOE. The metric used
% here is M(silence) = 1-( post/ pre silence detection ability ). Thus positive
% values describe the magnitude of a decrease in silence detection ability (and therefore evidence for tinnitus or hearing loss),
% and a negaitve value the magnitude of an increase in silence detection
% ability.

% --- load TI input for operant gap detection
load D:\PC-VV001-restoration\manuscript\Repository\Tinnitus_Index_metadata\TI_OperantGap_thresh&NGapHR

% TI_op1 = [TI1(end-2:end,1) TI1(end-2:end,2)];
% TI_op2 = [TI2(end-2:end,1) TI2(end-2:end,2)];

% use columns 2 and 3 for FA and NormThreshold
TI_op1 = [TI1(end-2:end,2) TI1(end-2:end,3)];
TI_op2 = [TI2(end-2:end,2) TI2(end-2:end,3)];
clear TI1 TI2 TI_op1_merged TI_op1_merged
for anm = 1:3
TI_op1_merged(anm) = sum([TI_op1(anm,1) TI_op1(anm,2)]);
TI_op2_merged(anm) = sum([TI_op2(anm,1) TI_op2(anm,2)]);
end

% INFO on TI: The metric here refers to the change of ability in detecting
% the continuous sound (non-gap) and the change in gap-detection threshold.The metric for continuous sound detection ability used
% here is M(cont) = ( post/ pre percent correct )-1. Thus postive values
% describe the magnitude of an increase in cont. sound detection ability (and therefore evidence for tinnitus),
% and negative values the magnitude of a decrease in cont. sound detection
% ability. Likewise, the metric for a change in gap detection threshold is
% defined as M(thresh) = ( post/ pre percent correct )-1. Thus postive values
% describe an increase on threshold (evidence for tinnitus and hearing
% loss) and negative values a decrease in threshold.




% calculate overall behavioural TI 
% sum of TI input of all measures

% Beh. tinnitus index input 
for anm = 1:3
TI_all_1_meta{anm} = [TI_silence1(anm) TI_op1_merged(anm) ];
TI_all_2_meta{anm} = [TI_silence2(anm) TI_op2_merged(anm) ];
end


% Beh. tinnitus index 
clear TI_all_1 TI_all_2
for anm = 1:3
TI_all_1(anm) = nansum(TI_all_1_meta{anm});
TI_all_2(anm) = nansum(TI_all_2_meta{anm});
end

% INFO on overall behavioural TI:
% The behavioural tinnitus index (TI) describes the evidence for
% tinnitus and hearing loss based on all behavioural tests. TI is the sum
% of the individual metrics for silence detection and operant gap detection.
% TI is defined as TI = M(silence)+M(cont)+M(GPIAS).


%% load TI input for ABRs (based on RMS total response magnitude)
 % input is not yet in the final format
for freq = [ 1 2 4 8 16 20]
load(['D:\PC-VV001-restoration\manuscript\Repository\Tinnitus_Index_metadata\TI_RMSmag_' num2str(freq) 'K_Aug22'])

% a decrease in ABR RMS should be a neg. value, an increase a pos. value
ABRin1{freq} = RMSrat1-1;
ABRin2{freq} = RMSrat2-1;
end
for anm = 1:3
ABRall1(anm) = mean([ABRin1{1}(anm) ABRin1{2}(anm) ABRin1{4}(anm) ABRin1{8}(anm) ABRin1{16}(anm) ABRin1{20}(anm)]);
ABRall2(anm) = mean([ABRin2{1}(anm) ABRin2{2}(anm) ABRin2{4}(anm) ABRin2{8}(anm) ABRin2{16}(anm) ABRin2{20}(anm)]);

ABRall1RMS(anm) = ABRall1(anm);
ABRall2RMS(anm) = ABRall2(anm);
end

% INFO on ABR measure: The metric here refers to the change in audfitory
% brainstem response (ABR) magnitude after relative to before NOE. The
% measure for response magnitude is the root-measn square of the response in a
% predefined response window, including measures across all stimulus
% intensities ('RMS area', see result chapter 2, Fig. XY). % The metric is defined as M(ABR) = ( post/ pre RMS area )-1. Here, postive values
% indicate the magnitude of an response increase whereas negative values
% indicate a reduced response (evidence for hearing loss).


%% load TI input for ABRs (based on scored threshold changes in mean dB differences)

load(['D:\PC-VV001-restoration\manuscript\Repository\Tinnitus_Index_metadata\ABRthresholdChange_03Feb23_dB.mat'])

% a decrease in ABR threshold relative to BL is a neg. value, an increase a pos. value
for anm = 1:3
ABRall1(anm) = AvChanThresh{anm}(1);
ABRall2(anm) = AvChanThresh{anm}(2);

ABRall1THRESH(anm) = ABRall1(anm);
ABRall2THRESH(anm) = ABRall2(anm);

end


%% Fig. S3D:  plot heatmap (TI, ABR thresh & ABR magnitude)


% close all
format short

savefig = 1;

% ID = [3]
figure
% TI only
subplot(1,2,1)
for anm = ID 
    
    cdata = [ TI_all_1(anm) ;...
              TI_all_2(anm) ];
          

    xvalues = {['TI (idx)']};
    yvalues = {'Post1','Post2'};
    
    xvalues = {['TI (idx)']};
    yvalues = {'Post1','Post2'};

%     h = heatmap(xvalues,yvalues,cdata,'Colormap',winter);
    h = heatmap(xvalues,yvalues,cdata,'Colormap',gray);


%     h.Title = 'Summary of changes in behaviour, ABRS and brain activity (fro) after NOE';
    h.FontSize = 14;
    
     % hide colourbar for all but first animal
    if ID > 1
%     h.ColorbarVisible = 'off'
    end
    % set displayed values to 1 decimal point max
    h.CellLabelFormat = '%.1f';


%     caxis([(negm), (posm)]); % as mean of pos and meanm of neg data as colour range limits
    caxis([(0), (2.5)]); % fix colour range limits

end

% ABR only
subplot(1,2,2)
for anm = ID 
    
%      cdata = [  ABRall1(anm) ;...
%               ABRall2(anm)  ]; % previously: *100
          cdata = [  ABRall1THRESH(anm)  ABRall1RMS(anm)*100 ;...
                     ABRall2THRESH(anm)  ABRall2RMS(anm)*100  ];

    xvalues = {'Hearing Sens.','ABR Mag.'};
    yvalues = {'.','..'};

    h = heatmap(xvalues,yvalues,cdata,'Colormap',gray);

%     h.Title = 'Summary of changes in behaviour, ABRS and brain activity (fro) after NOE';
    h.FontSize = 14;
    
    % hide ylabels
    h.YDisplayLabels{1} = ''
    h.YDisplayLabels{2} = ''
    
     % hide colourbar for all but first animal
    if ID > 1
%     h.ColorbarVisible = 'off'
    end
    % set displayed values to 1 decimal point max
    h.CellLabelFormat = '%.1f';

%     caxis([(negm), (posm)]); % as mean of pos and meanm of neg data as colour range limits
    caxis([(-100), (100)]); % fix colour range limits

    Ax.CellLabelColor='none';


end


%% Fig. 2F

% plots TI for post 1 and 2 in one plot and ABR chnages in separate plot.
% Removes xaxis from plot (or both plots) so plots can be manually merged
% into one figure with two separate yaxes for the different scales between
% TI and ABR.


dosave = 0;

% plot for all animals
% ID = [1];

savepath = 'D:\PC-VV001-restoration\Ferrets\Ephys\final_scripts\ABR_metadata\Plots\';
colornames

% close all
% Plot TI
    for anm = ID 
    
    cdata = [ TI_all_1(anm) ;...
              TI_all_2(anm) ];
          figure; hold all
          p1 = plot([1 1],[0 cdata(1)]);
          p2 = plot([2 2],[0 cdata(2)]);
          
          xlim([0 6])
          ylim([0 3])
          % run colornames_view to find colour
          set([p1 p2],'linewidth',15,'Color',[1 0.98 .8])
%           set([p1 p2],'linewidth',15,'Color',colornames)
          


          set([gca],'linewidth',1.5)
          xticks([])
          
          % remove xaxis
         h = gca;
         h.XAxis.Visible = 'off';
          
         

    end
    
    % Plot ABR
    for anm = ID 
    
    cdata = [  ABRall1THRESH(anm) ; ABRall1RMS(anm)*100 ;...
              ABRall2THRESH(anm) ; ABRall2RMS(anm)*100  ]; % previously: *100
          
          figure; hold all
          p1 = plot([.8 .8],[0 cdata(1)]);
          p2 = plot([1.1 1.1],[0 cdata(2)]);
          p3 = plot([1.8 1.8],[0 cdata(3)]);
          p4 = plot([2.1 2.1],[0 cdata(4)]);
          
          %           p1 = bar([.8],[cdata(1)]);

          
          xlim([0 6])
          ylim([-200 0])
          set([p1 p3],'linewidth',15,'Color',[.8 .8 .8])
          set([p2 p4],'linewidth',15,'Color',[.5 .5 .5],"MarkerEdgeColor","k")

          set([gca],'linewidth',1.5)
         
          % remove xaxis
          h = gca;
          h.XAxis.Visible = 'off';
         
         
          

    end
    





