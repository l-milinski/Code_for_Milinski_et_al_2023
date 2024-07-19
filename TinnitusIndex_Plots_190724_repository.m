%% Calculates Tinnitus index (TI) and chnages in ABR
% Plots Figs.3A,B.



clear all

anmnames = {'Heffalump';'Kanga';'Piglet'};
savefig = 1;

ID = 1; % Heffalump,Kanga, Piglet

PlotWithMCont = 0; % if 1 then tinnitus index includes M(cont) metric, if 0 then it doesn't.

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
    for anm = ID
    TI_op1_merged(anm) = sum([TI_op1(anm,1) TI_op1(anm,2)]);
    TI_op2_merged(anm) = sum([TI_op2(anm,1) TI_op2(anm,2)]);
    end

% IF SIMPLIFIED TI_op1 IS TO BE USED (no M(cont))
if PlotWithMCont == 0
    % use columns 3 (in the TI variable) for NormThreshold
    for anm = 1:3
    TI_op1_merged(anm) = sum([ TI_op1(anm,2)]);
    TI_op2_merged(anm) = sum([ TI_op2(anm,2)]);
    end
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
for anm = ID
TI_all_1_meta{anm} = [TI_silence1(anm) TI_op1_merged(anm) ];
TI_all_2_meta{anm} = [TI_silence2(anm) TI_op2_merged(anm) ];
end


% Beh. tinnitus index 
clear TI_all_1 TI_all_2
for anm = ID
TI_all_1(anm) = nansum(TI_all_1_meta{anm});
TI_all_2(anm) = nansum(TI_all_2_meta{anm});
end

% INFO on overall behavioural TI:
% The behavioural tinnitus index (TI) describes the evidence for
% tinnitus and hearing loss based on all behavioural tests. TI is the sum
% of the individual metrics for silence detection and operant gap detection.
% TI is defined as TI = M(silence)+M(cont)+M(GPIAS).

%% Load TI for aboslute values (behavioural metrics)
clear TI1 TI2
load D:\PC-VV001-restoration\manuscript\Repository\Tinnitus_Index_metadata_absolutes\TI_silence_abs

close all
savefig2 = 0

% TI_silence_abs

for anm = 6:8
    cdataSilence(anm,:) = [ TI_silence_abs(anm,:)];
end



clear TI1 TI2
load D:\PC-VV001-restoration\manuscript\Repository\Tinnitus_Index_metadata_absolutes\TI_OperantGap_thresh&NGapHR_abs
% FAabs
% TreshNorm_abs

% plot absolute metrics for each animal in parallel barplots
         figure
         title(['m (Silence abs) -> low value indicates tinnitus']) % low value -> tinnitus
         hold all
          p3 = bar([1 2 3],[TI_silence_abs(8,:)]);
          p1 = bar([5 6 7],[TI_silence_abs(6,:)]);
          p2 = bar([9 10 11],[TI_silence_abs(7,:)]);

         
          p1.FaceColor = [0.1 0.7 0.7];
          p2.FaceColor = [0.1 0.6 0.7];
          p3.FaceColor = [0.1 0.5 0.7];
          
          xticks([ 2 6 10])
          xticklabels({'Case 1','Case 2','Case 3'})
          
           if savefig2 == 1
               savepath = 'D:\PC-VV001-restoration\manuscript\Revision June 2024'
           print('-r750','-dtiff',[savepath '\T_m(silence)_abs.tiff']);
           end
 
 
          
         figure
         title('m (Cont abs) -> high value indicates tinnitus') % high value -> tinnitus
         hold all
         
          p3 = bar([1 2 3],[FAabs(8,:)]);
          p1 = bar([5 6 7],[FAabs(6,:)]);
          p2 = bar([9 10 11],[FAabs(7,:)]);
          
          p1.FaceColor = [0.1 0.7 0.7];
          p2.FaceColor = [0.1 0.6 0.7];
          p3.FaceColor = [0.1 0.5 0.7];
          
          xticks([ 2 6 10])
          xticklabels({'Case 1','Case 2','Case 3'})
          
          
           if savefig2 == 1
               savepath = 'D:\PC-VV001-restoration\manuscript\Revision June 2024'
           print('-r750','-dtiff',[savepath '\T_m(cond)_abs.tiff']);
           end
 
 
          
          
         figure
         hold all
          p3 = bar([1 2 3],[TreshNorm_abs(8,:)]);
          p1 = bar([5 6 7],[TreshNorm_abs(6,:)]);
          p2 = bar([9 10 11],[TreshNorm_abs(7,:)]);

          p1.FaceColor = [0.1 0.7 0.7];
          p2.FaceColor = [0.1 0.6 0.7];
          p3.FaceColor = [0.1 0.5 0.7];
          
          xticks([ 2 6 10])
          xticklabels({'Case 1','Case 2','Case 3'})
          
          set(gca,'linewidth',1.5)
          set(gca,'fontsize',14)
          
          title('m (Thresh norm abs) -> low value indicates tinnitus','Fontsize',8) % low value -> tinnitus

          
           if savefig2 == 1
               savepath = 'D:\PC-VV001-restoration\manuscript\Revision June 2024'
           print('-r750','-dtiff',[savepath '\T_m(thresh norm)_abs.tiff']);
           end
          
          
 
          
          
          % TO DO 180624: save the plots in this section (label and
          % consider sorting by case in the manuscript; save figures; also
          % save the previous tinnitus index figure with and without
          % m(cont).
          % Anm1 -> case 2
          % Anm2 -> case 3
          % Anm3 -> case 1
        
           
           %% (NEW) Fig. 3b: Plot thresholds (absolute values)
           % individual plots for each animal
           
           savefig2 = 1; % to save plots
           
           Case = 1; % define case to plot
           
           % ----
           sleepdur = [85.6 70.6 71.5];
           wakedur = [28.5 14.4 29.4];
           NREMdur = [41.4 51.1 39.9];
           REMdur = [18.3 23.5 8.1];
           REM2dur = [8.4 7.4 16.2];
           % Info: SEM = std/sqrt(n)
           sleepdurSEM = std(sleepdur)/sqrt(length(sleepdur))
           wakedurSEM = std(wakedur)/sqrt(length(sleepdur))
           NREMdurSEM = std(NREMdur)/sqrt(length(sleepdur))
           REMdurSEM = std(REMdur)/sqrt(length(sleepdur))
           REM2durSEM = std(REM2dur)/sqrt(length(sleepdur))

           
           
           
           
         close all
         figure
         hold all
         
          % Set edge colours by time
          cmap = [.5 .5 1; 1 .5 .5; 0.9 0.6 0.3];
          
          if Case == 1
          p31 = bar([1],[TreshNorm_abs(8,1)],'EdgeColor',cmap(1,:),'Linewidth', 3);
          p32 = bar([3],[TreshNorm_abs(8,2)],'EdgeColor',cmap(2,:),'Linewidth', 3);
          p33 = bar([5],[TreshNorm_abs(8,3)],'EdgeColor',cmap(3,:),'Linewidth', 3);
          elseif Case == 2
          p31 = bar([1],[TreshNorm_abs(6,1)],'EdgeColor',cmap(1,:),'Linewidth', 3);
          p32 = bar([3],[TreshNorm_abs(6,2)],'EdgeColor',cmap(2,:),'Linewidth', 3);
          p33 = bar([5],[TreshNorm_abs(6,3)],'EdgeColor',cmap(3,:),'Linewidth', 3);
          elseif Case == 3
          p31 = bar([1],[TreshNorm_abs(7,1)],'EdgeColor',cmap(1,:),'Linewidth', 3);
          p32 = bar([3],[TreshNorm_abs(7,2)],'EdgeColor',cmap(2,:),'Linewidth', 3);
          p33 = bar([5],[TreshNorm_abs(7,3)],'EdgeColor',cmap(3,:),'Linewidth', 3);
          end
          
          
          % Plot balck line on x axis
          plot([0 6],[0 0],'-k','linewidth',1.8)

%           p11 = bar([8],[TreshNorm_abs(6,1)],'EdgeColor',cmap(1,:),'Linewidth', 1.5);
%           p12 = bar([10],[TreshNorm_abs(6,2)],'EdgeColor',cmap(2,:),'Linewidth', 1.5);
%           p13 = bar([12],[TreshNorm_abs(6,3)],'EdgeColor',cmap(3,:),'Linewidth', 1.5);
%           
%           p21 = bar([15],[TreshNorm_abs(7,1)],'EdgeColor',cmap(1,:),'Linewidth', 1.5);
%           p22 = bar([17],[TreshNorm_abs(7,2)],'EdgeColor',cmap(2,:),'Linewidth', 1.5);
%           p23 = bar([19],[TreshNorm_abs(7,3)],'EdgeColor',cmap(3,:),'Linewidth', 1.5);
%           
%           p1 = bar([5 6 7],[TreshNorm_abs(6,:)]);
%           p2 = bar([9 10 11],[TreshNorm_abs(7,:)]);

%           p1.FaceColor = [0.1 0.7 0.7];
%           p2.FaceColor = [0.1 0.6 0.7];
%           p3.FaceColor = [0.1 0.5 0.7];
          
          xticks([ 1 3 5])
          xticklabels({'BL','1 week','6 months'})

%           xticklabels({'Case 1','Case 2','Case 3'})
          ylabel({'Gap detection ability'; ' (% of normalised threshold) '})
          
          set(gca,'linewidth',2.5)
          set(gca,'fontsize',16)

%            title('m (Thresh norm abs) -> low value indicates tinnitus','Fontsize',8) % low value -> tinnitus
          title(['Case ' num2str(Case)])

%           set([p31 p32 p33 p21 p22 p23 p11 p12 p13],'FaceColor',[1 0.98 .8],'BarWidth', 1.6)
          set([p31 p32 p33 ],'FaceColor',[1 0.95 .8],'BarWidth', 1.2)    


        
         
         
      



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


%% TI with (and without) m(cont)

% plots TI for post 1 and 2 in one plot and ABR chnages in separate plot.
% Removes xaxis from plot (or both plots) so plots can be manually merged
% into one figure with two separate yaxes for the different scales between
% TI and ABR.

% TO DO 19/06/24: save the TI plots with and without cond; but ideally
% include the case number in the title

close all


% savefig = 0;

% plot for all animals
% ID = [1];

colornames
% savefig = 0

% close all
% Plot TI
    for anm = ID
    
    cdata = [ TI_all_1(anm) ;...
              TI_all_2(anm) ];
          
           if PlotWithMCont == 1
               fig1 = figure(1);
             title('TI all metrics')
         elseif PlotWithMCont == 0
               fig11 = figure(1);
             title('TI without M(cont)')
           end
          
          hold all
          
          p1 = bar([1 2],[cdata(1) cdata(2)]);
          
%           p1 = plot([1 1],[0 cdata(1)]);
%           p2 = plot([2 2],[0 cdata(2)]);
          
%           xlim([0 6])
          ylim([0 1])
          % run colornames_view to find colour
          set([p1],'FaceColor',[1 0.98 .8])
%           set([p1 p2],'linewidth',15,'Color',colornames)
          
          set([gca],'linewidth',1.5)
          xticks([ 1.5])
          if ID == 3
              xticklabels({'Case1'})
          elseif ID == 2
              xticklabels({'Case3'})
          elseif ID == 1
              xticklabels({'Case2'})
          end
          
          % remove xaxis
%          h = gca;
%          h.XAxis.Visible = 'off';

    
    
                 savepath = 'D:\PC-VV001-restoration\manuscript\Revision June 2024'
     if savefig == 1
                 if PlotWithMCont == 1
                           savename = [savepath '\TI_ANM' num2str(ID)];
                    print(fig1,'-r750','-dtiff',[savename '.tif'],'-painters'); 
                 else
                           savename = [savepath '\TI_noMcont_ANM' num2str(ID)];
                    print(fig11,'-r750','-dtiff',[savename '.tif'],'-painters'); 
                 end
     end
     
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
    
    
 
    






%% Fig. 2A

% plots TI for post 1 and 2 in one plot and ABR chnages in separate plot.
% Removes xaxis from plot (or both plots) so plots can be manually merged
% into one figure with two separate yaxes for the different scales between
% TI and ABR.



close all


% savefig = 0;

% plot for all animals
% ID = [1];

colornames
savefig = 1

% close all
% Plot TI
    for anm = ID
    
    cdata = [ TI_all_1(anm) ;...
              TI_all_2(anm) ];
          
           if PlotWithMCont == 1
               fig1 = figure(1);
             title('TI all metrics')
         elseif PlotWithMCont == 0
               fig11 = figure(1);
             title('TI without M(cont)')
           end
          
          hold all
          
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

    
    
    
    
                 savepath = 'D:\PC-VV001-restoration\manuscript\Revision June 2024'
%      if savefig == 1
%                  if PlotWithMCont == 1
%                            savename = [savepath '\TI_ANM' num2str(ID)];
%                     print(fig1,'-r750','-dtiff',[savename '.tif'],'-painters'); 
%                  else
%                            savename = [savepath '\TI_noMcont_ANM' num2str(ID)];
%                     print(fig11,'-r750','-dtiff',[savename '.tif'],'-painters'); 
%                  end
%      end
     
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
    
    
 
    





