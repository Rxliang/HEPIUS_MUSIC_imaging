% script to analyze multi-depth Doppler data
bLoad = 1;
if bLoad
    load depthSats4
end

% Look at aggregated Vd and PI across pressures for each depth and for Int
% vs Ext
[numPressures, numDepths, numData] = size(statsInt);

bGetVdBaseline = 0;
bPI = 1;

bStatsByPressure = 0;
bStatsByDepth = 1;

bMatchByPressure = 0;
bMatchByDepth = 0;

bMatchIntExt = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GetVdBaseline

if bGetVdBaseline
    vDmean1 = [];
    vDmean2 = [];
for p = 1
    piExtG = [];
    piExtY = [];
    idx = 1;
    vD1 = [];
    vD2 = [];
    for d = 1:numDepths
        depthLabel = depthFolderName{d};
        vD1 = [];
        vD2 = [];
        for k = 1:numData
            tmpDat = statsInt{p,d,k};
            if ~isempty(tmpDat)
                % create stats and group label for anova for each depth
                if tmpDat.PI.snr > 0
                    vD1 = [vD1 tmpDat.Vd.Vd.all];
                end
            end
            tmpDat = statsExt{p,d,k};
            if ~isempty(tmpDat)
                % create stats and group label for anova for each depth
                if tmpDat.PI.snr > 0
                    vD2 = [vD2 tmpDat.Vd.Vd.all];
                end
            end
        end
        vDmean1(d) = mean(vD1);
        vDmean2(d) = mean(vD2);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% StatsByPressure

if bStatsByPressure
for d = 1:numDepths
    piExtG = [];
    piExtY = [];
    idx = 1;
    for p = 1:numPressures
        presLabel = folderNamePres{p};
        for k = 1:numData
            tmpDat = statsInt{p,d,k};
            if ~isempty(tmpDat)
                % create stats and group label for anova for each depth
                if tmpDat.PI.snr > 0
                    if bPI
                        pis = tmpDat.peakTiming.peakTiming.all;
                    else
                        pis = tmpDat.Vd.Vd.all/vDmean1(d);
                    end
                    if isempty(pis)
                        errordlg('Take a look!!');
                    else
                        numPis = length(pis);
                        piExtY(idx:idx+numPis-1) = pis;
                        for a = 1:numPis
                            piExtG{idx+a-1} = presLabel;
                        end
                        idx = idx + numPis;
                    end
                end
            end
        end
    end
    [p,tbl,stats] = anova1(piExtY, piExtG);
    depthName = depthFolderName{d};
    title(['Depth: ' depthName(1:2) ' mm']);
    pause;
    close(gcf);
    close(gcf);
    for p = 1:numPressures
        presLabel = folderNamePres{p};
        for k = 1:numData
            tmpDat = statsExt{p,d,k};
            if ~isempty(tmpDat)
                % create stats and group label for anova for each depth
                if tmpDat.PI.snr > 0
                    if bPI
                        pis = tmpDat.peakTiming.peakTiming.all;
                    else
                        pis = tmpDat.Vd.Vd.all/vDmean2(d);
                    end
                    if isempty(pis)
                        errordlg('Take a look!!');
                    else
                        numPis = length(pis);
                        piExtY(idx:idx+numPis-1) = pis;
                        for a = 1:numPis
                            piExtG{idx+a-1} = presLabel;
                        end
                        idx = idx + numPis;
                    end
                end
            end
        end
    end
    [p,tbl,stats] = anova1(piExtY, piExtG);
    depthName = depthFolderName{d};
    title(['Depth: ' depthName(4:5) ' mm']);
    pause;
    close(gcf);
    close(gcf);
    [c,~,~,gnames] = multcompare(stats);
    pause;
    close(gcf);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% StatsByDepth

if bStatsByDepth
for p = 1:numPressures
    piExtG = [];
    piExtY = [];
    idx = 1;
    for d = 1:numDepths
        depthLabel = depthFolderName{d};
        depthName1 = depthLabel(1:2);
        depthName2 = depthLabel(4:5);
        for k = 1:numData
            tmpDat = statsInt{p,d,k};
            if ~isempty(tmpDat)
                % create stats and group label for anova for each depth
                if tmpDat.PI.snr > 0
                    if bPI
                        pis = tmpDat.PI.PI.all;
                    else
                        pis = tmpDat.Vd.Vd.all/vDmean1(d);
                    end
                    if isempty(pis)
                        errordlg('Take a look!!');
                    else
                        numPis = length(pis);
                        piExtY(idx:idx+numPis-1) = pis;
                        for a = 1:numPis
                            piExtG{idx+a-1} = depthName1;
                        end
                        idx = idx + numPis;
                    end
                end
            end
            tmpDat = statsExt{p,d,k};
            if ~isempty(tmpDat)
                % create stats and group label for anova for each depth
                if tmpDat.PI.snr > 0
                    if bPI
                        pis = tmpDat.PI.PI.all;
                    else
                        pis = tmpDat.Vd.Vd.all/vDmean2(d);
                    end
                    if isempty(pis)
                        errordlg('Take a look!!');
                    else
                        numPis = length(pis);
                        piExtY(idx:idx+numPis-1) = pis;
                        for a = 1:numPis
                            piExtG{idx+a-1} = depthName2;
                        end
                        idx = idx + numPis;
                    end
                end
            end
        end
    end
    [p1,tbl,stats] = anova1(piExtY, piExtG);
    title(['Pressure: ' folderNamePres{p}],'Interpreter', 'none');
    pause;
    close(gcf);
    close(gcf);
%     [c,~,~,gnames] = multcompare(stats);
%     pause;
%     close(gcf);
    plot(stats.means)
    title(['Pressure: ' folderNamePres{p}],'Interpreter', 'none');
    x_tick = 1:1:length(stats.gnames);
    set(gca,'XTick', x_tick, 'XTickLabel',stats.gnames)
    hold on;
    pause;
    %close(gcf);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MatchByPressure

if bMatchByPressure
for d = 1:numDepths
    piExtG = [];
    piExtY = [];
    idx = 1;
    for p = 1:numPressures
        presLabel = folderNamePres{p};
        for k = 1:numData
            tmpDat = statsMatch{p,d,k};
            if ~isempty(tmpDat)
                % create stats and group label for anova for each depth
                
                pis = abs(tmpDat.dPI);
                if isempty(pis)
                    %        errordlg('Take a look!!');
                else
                    numPis = length(pis);
                    piExtY(idx:idx+numPis-1) = pis;
                    for a = 1:numPis
                        piExtG{idx+a-1} = presLabel;
                    end
                    idx = idx + numPis;
                end
                
            end
        end
    end
    [p,tbl,stats] = anova1(piExtY, piExtG);
    title(['Depth: ' depthFolderName{d}]);
    pause;
    close(gcf);
    close(gcf);
    [c,~,~,gnames] = multcompare(stats);
    pause;
    close(gcf);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MatchByDepth

if bMatchByDepth
for p = 1:numPressures
    piExtG = [];
    piExtY = [];
    idx = 1;
    for d = 1:numDepths
        depthLabel = depthFolderName{d};
        for k = 1:numData
            tmpDat = statsMatch{p,d,k};
            if ~isempty(tmpDat)
                % create stats and group label for anova for each depth
                
                pis = abs(tmpDat.dPI);
                if isempty(pis)
                    %errordlg('Take a look!!');
                else
                    numPis = length(pis);
                    piExtY(idx:idx+numPis-1) = pis;
                    for a = 1:numPis
                        piExtG{idx+a-1} = depthLabel;
                    end
                    idx = idx + numPis;
                end
                
            end
        end
    end
    [p1,tbl,stats] = anova1(piExtY, piExtG);
    title(['Pressure: ' folderNamePres{p}]);
    pause;
    close(gcf);
    close(gcf);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MatchIntExt

if bMatchIntExt
    for d = 1:numDepths-1
        piExtY = [];
        idx = 1;
        for p = 1:numPressures
            presLabel = folderNamePres{p};
            for k = 1:numData
                tmpDat = statsInt{p,numDepths,k};
                if ~isempty(tmpDat)
                    % create stats and group label for anova for each depth
                    if tmpDat.PI.snr > 0
                        pis = tmpDat.PI.PI.all;
                        
                        numPis = length(pis);
                        piExtY(idx:idx+numPis-1) = pis;
                        idx = idx + numPis;
                        
                    end
                end
            end
            idxExt = mean(piExtY);
            piExtY = [];
            idx = 1;
            
            for k = 1:numData
                tmpDat = statsInt{p,d,k};
                if ~isempty(tmpDat)
                    % create stats and group label for anova for each depth
                    if tmpDat.PI.snr > 0
                        pis = tmpDat.PI.PI.all;
                        
                        numPis = length(pis);
                        piExtY(idx:idx+numPis-1) = pis;
                        idx = idx + numPis;
                        
                    end
                end
            end
            idx1 = mean(piExtY);
            piExtY = [];
            idx = 1;
            
            balIdx1(p) = (idx1 - idxExt);
            
            for k = 1:numData
                tmpDat = statsExt{p,d,k};
                if ~isempty(tmpDat)
                    % create stats and group label for anova for each depth
                    if tmpDat.PI.snr > 0
                        pis = tmpDat.PI.PI.all;
                        
                        numPis = length(pis);
                        piExtY(idx:idx+numPis-1) = pis;
                        idx = idx + numPis;
                        
                    end
                end
            end
            idx2 = mean(piExtY);
            piExtY = [];
            idx = 1;
            
            balIdx2(p) = (idx2 - idxExt);
        end
        
        pres = [0 12 24 36];
        pres1 = 0:.5:36;
        plot(pres, balIdx1,'o')
        [PP,SS] = polyfit(pres, balIdx1, 2);
        polyData = polyval(PP,pres1);
        
        hold on;
        plot(pres1, polyData,'r')
        [~,minIdx] = min(polyData);
        if minIdx > 3 && minIdx < length(polyData)
            line([pres1(minIdx) pres1(minIdx)], ylim);
        end
        depthName = depthFolderName{d};
        title(['Depth: ' depthName(1:2) ' mm'],'Interpreter', 'none');
        xlabel('P_e [mmHg]')
        ylabel('PI_{ext} - PI_{int}')
        pause;
        close(gcf);
        plot(pres, balIdx2,'o')
        [PP,SS] = polyfit(pres, balIdx1, 2);
        polyData = polyval(PP,pres1);
        hold on;
        plot(pres1, polyData,'r')
        [~,minIdx] = min(polyData);
        if minIdx > 1 && minIdx < length(polyData)
            line([pres1(minIdx) pres1(minIdx)], ylim);
        end
        depthName = depthFolderName{d};
        title(['Depth: ' depthName(4:5) ' mm'],'Interpreter', 'none');
        xlabel('P_e [mmHg]')
        ylabel('PI_{ext} - PI_{int}')
        pause;
        close(gcf);
       
    end
end
