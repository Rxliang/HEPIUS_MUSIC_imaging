function [outputArg1,outputArg2] = spectralDoppler_JP20231215(IData,QData)
persistent counter
persistent ax
persistent x_mm
persistent z_mm
persistent spectAx
persistent veldop
persistent spectData
persistent t

%This assumes that PData 2 is for doppler reconstruction
dopPData = 2; 
PData = evalin('base','PData');
Trans = evalin('base','Trans');
dop = evalin('base','dop');
datcols = 256;
Receive = evalin('base','Receive');
bkgdUSrcv = evalin('base','bkgdUSRcv');
Resource = evalin('base','Resource');
%keyboard

if ~exist('counter','var') || isempty(counter)
    dopFig = findobj('Name','MUSIC Doppler');
    ax = dopFig.Children;
    for i = 1:length(ax.Children)
        set(ax(i).Children,'HitTest','off');
        set(ax(i).Children,'ButtonDownFcn',@clicker)
    end
    set(ax,'ButtonDownFcn',@clicker);
    counter = 1;
    x_mm = Trans.wvl*(PData(dopPData).Origin(1):PData(dopPData).PDelta(1):PData(dopPData).Origin(1)+PData(dopPData).PDelta(1)*(PData(dopPData).Size(2)-1));
    z_mm = Trans.wvl*(PData(dopPData).Origin(3):PData(dopPData).PDelta(3):PData(dopPData).Origin(3)+PData(dopPData).PDelta(3)*(PData(dopPData).Size(1)-1));
    
    f = figure();
    spectAx = axes(f);
    f2vel = Resource.Parameters.speedOfSound/(2*Receive(bkgdUSrcv+1).demodFrequency*10^6)*100;

    veldop = linspace(-dop.PRF/2,dop.PRF/2,size(IData,4))*f2vel;
    t = (1:datcols)*1000/dop.PRF;
    spectData = zeros(length(veldop),datcols);
end

if length(findobj(ax.Children,'Type','Line')) == 2
    objs = findobj(ax.Children,'Type','Line');
    pos = zeros(2,2);
    for i = 1:length(objs)
        pos(i,1) = mean(objs(i).XData);
        pos(i,2) = mean(objs(i).YData);
    end
    [~,Ix] = min(abs(pos(1,1)-x_mm),[],2);
    [~,Iz] = min(abs(pos(1,2)-z_mm),[],2);
    ang = atan((pos(2,1)-pos(1,1))/(pos(2,2)-pos(1,2)));
    
    IQ = squeeze(complex(IData(Iz,Ix,1,:,1),QData(Iz,Ix,1,:,1)));
    dataPos = mod(counter,datcols);
    dataPosNext = dataPos + 1;
    if dataPos == 0
        dataPos = datcols;
    end
    spectData(:,dataPos) = fftshift(fft(IQ));
    spectData(:,dataPosNext) = zeros(length(veldop),1);
    lcspect = 20*log10(abs(spectData));
    lcspect(isinf(lcspect(:)))= 0;
    if counter < datcols
        tmp = lcspect(:,1:counter);
        clims = [prctile(tmp(:),25),prctile(tmp(:),97)];
    else
        clims = [prctile(lcspect(:),25),prctile(lcspect(:),97)];
    end
    imagesc(spectAx,t,veldop/cos(ang),lcspect,clims);
    xlabel(spectAx,'Slow Time (ms)');
    ylabel(spectAx,'Velocity (cm/s)');
    counter = counter + 1;
end

end

function clicker(h,cpt)

    if strcmp(h.Parent.SelectionType,'normal')
        rectpts = [cpt.IntersectionPoint(1:2)-0.25;cpt.IntersectionPoint(1)-0.25,cpt.IntersectionPoint(2)+0.25; cpt.IntersectionPoint(1:2)+0.25;cpt.IntersectionPoint(1)+0.25,cpt.IntersectionPoint(2)-0.25;cpt.IntersectionPoint(1:2)-0.25];
        objs = findobj(h.Children,'Type','Line','Color',[1,0,0,0.35]);
        if length(objs)==1
            objs.XData = rectpts(:,1);
            objs.YData = rectpts(:,2);
        else
            hold(h,'on');
            plot(h,rectpts(:,1),rectpts(:,2),'LineWidth',2,'Color',[1,0,0,0.35],'HitTest','off')
            hold(h,'off');
        end
    elseif strcmp(h.Parent.SelectionType, 'alt')
        rectpts = [cpt.IntersectionPoint(1:2)-0.25;cpt.IntersectionPoint(1)-0.25,cpt.IntersectionPoint(2)+0.25; cpt.IntersectionPoint(1:2)+0.25;cpt.IntersectionPoint(1)+0.25,cpt.IntersectionPoint(2)-0.25;cpt.IntersectionPoint(1:2)-0.25];
        objs = findobj(h.Children,'Type','Line','Color',[0,1,0,0.35]);
        if length(objs)==1
            objs.XData = rectpts(:,1);
            objs.YData = rectpts(:,2);
        else
            hold(h,'on');
            plot(h,rectpts(:,1),rectpts(:,2),'LineWidth',2,'Color',[0,1,0,0.35],'HitTest','off')
            hold(h,'off');
        end
    else
        disp('unrecognized selection type');
    end

%     clickT = get(h, 'selectiontype');
    % 'normal' for left moue button
    % 'alt' for right mouse button
    % 'extend' for middle mouse button
    % 'open' on double click

%     clickpt = get(h, 'currentpoint');
    % Current mouse location, in pixels from the lower left.
    % When the units of the figure are 'normalized', the
    % coordinates will be [0 0] inb lower left, and [1 1] in
    % the upper right.

end