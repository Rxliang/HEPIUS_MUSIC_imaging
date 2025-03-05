% Script to process vdStats
% Investigate relationship of Vd and Vclutter with pressure

% load vdStats
numStudies = length(vdStatsInt);

statsObj.int.p0.vd = [];
statsObj.int.p4.vd = [];
statsObj.int.p8.vd = [];
statsObj.int.p12.vd = [];
statsObj.int.p16.vd = [];
statsObj.int.p20.vd = [];
statsObj.int.p24.vd = [];
statsObj.int.p28.vd = [];
statsObj.int.p32.vd = [];

statsObj.int.p0.pDiff = [];
statsObj.int.p4.pDiff = [];
statsObj.int.p8.pDiff = [];
statsObj.int.p12.pDiff = [];
statsObj.int.p16.pDiff = [];
statsObj.int.p20.pDiff = [];
statsObj.int.p24.pDiff = [];
statsObj.int.p28.pDiff = [];
statsObj.int.p32.pDiff = [];

statsObj.ext.p0.vd = [];
statsObj.ext.p4.vd = [];
statsObj.ext.p8.vd = [];
statsObj.ext.p12.vd = [];
statsObj.ext.p16.vd = [];
statsObj.ext.p20.vd = [];
statsObj.ext.p24.vd = [];
statsObj.ext.p28.vd = [];
statsObj.ext.p32.vd = [];

statsObj.ext.p0.pDiff = [];
statsObj.ext.p4.pDiff = [];
statsObj.ext.p8.pDiff = [];
statsObj.ext.p12.pDiff = [];
statsObj.ext.p16.pDiff = [];
statsObj.ext.p20.pDiff = [];
statsObj.ext.p24.pDiff = [];
statsObj.ext.p28.pDiff = [];
statsObj.ext.p32.pDiff = [];


for a = 1:numStudies
    pres = vdStatsInt{a}.pressure;
 %   if sum(isnan(vdStatsInt{a}.vd)) == 0
        for b = 1:length(pres)
            switch(pres(b))
                case 0
                    if ~isnan(vdStatsInt{a}.vd(b))
                        statsObj.int.p0.vd = [statsObj.int.p0.vd; vdStatsInt{a}.vd(b)];
                        statsObj.int.p0.pDiff = [statsObj.int.p0.pDiff; vdStatsInt{a}.pdiff(b)];
                    end
                case 4
                    if ~isnan(vdStatsInt{a}.vd(b))
                        statsObj.int.p4.vd = [statsObj.int.p4.vd; vdStatsInt{a}.vd(b)];
                        statsObj.int.p4.pDiff = [statsObj.int.p4.pDiff; vdStatsInt{a}.pdiff(b)];
                    end
                case 8
                    if ~isnan(vdStatsInt{a}.vd(b))
                        statsObj.int.p8.vd = [statsObj.int.p8.vd; vdStatsInt{a}.vd(b)];
                        statsObj.int.p8.pDiff = [statsObj.int.p8.pDiff; vdStatsInt{a}.pdiff(b)];
                    end
                case 12
                    if ~isnan(vdStatsInt{a}.vd(b))
                        statsObj.int.p12.vd = [statsObj.int.p12.vd; vdStatsInt{a}.vd(b)];
                        statsObj.int.p12.pDiff = [statsObj.int.p12.pDiff; vdStatsInt{a}.pdiff(b)];
                    end
                case 16
                    if ~isnan(vdStatsInt{a}.vd(b))
                        statsObj.int.p16.vd = [statsObj.int.p16.vd; vdStatsInt{a}.vd(b)];
                        statsObj.int.p16.pDiff = [statsObj.int.p16.pDiff; vdStatsInt{a}.pdiff(b)];
                    end
                case 20
                    if ~isnan(vdStatsInt{a}.vd(b))
                        statsObj.int.p20.vd = [statsObj.int.p20.vd; vdStatsInt{a}.vd(b)];
                        statsObj.int.p20.pDiff = [statsObj.int.p20.pDiff; vdStatsInt{a}.pdiff(b)];
                    end
                case 24
                    if ~isnan(vdStatsInt{a}.vd(b))
                        statsObj.int.p24.vd = [statsObj.int.p24.vd; vdStatsInt{a}.vd(b)];
                        statsObj.int.p24.pDiff = [statsObj.int.p24.pDiff; vdStatsInt{a}.pdiff(b)];
                    end
                case 28
                    if ~isnan(vdStatsInt{a}.vd(b))
                        statsObj.int.p28.vd = [statsObj.int.p28.vd; vdStatsInt{a}.vd(b)];
                        statsObj.int.p28.pDiff = [statsObj.int.p28.pDiff; vdStatsInt{a}.pdiff(b)];
                    end
                case 32
                    if ~isnan(vdStatsInt{a}.vd(b))
                        statsObj.int.p32.vd = [statsObj.int.p32.vd; vdStatsInt{a}.vd(b)];
                        statsObj.int.p32.pDiff = [statsObj.int.p32.pDiff; vdStatsInt{a}.pdiff(b)];
                    end
            end
%        end
    end
    % ext
  %  if sum(isnan(vdStatsExt{a}.vd)) == 0
        for b = 1:length(pres)
            switch(pres(b))
                case 0
                    if ~isnan(vdStatsExt{a}.vd(b))
                        statsObj.ext.p0.vd = [statsObj.ext.p0.vd; vdStatsExt{a}.vd(b)];
                        statsObj.ext.p0.pDiff = [statsObj.ext.p0.pDiff; vdStatsExt{a}.pdiff(b)];
                    end
                case 4
                    if ~isnan(vdStatsExt{a}.vd(b))
                        statsObj.ext.p4.vd = [statsObj.ext.p4.vd; vdStatsExt{a}.vd(b)];
                        statsObj.ext.p4.pDiff = [statsObj.ext.p4.pDiff; vdStatsExt{a}.pdiff(b)];
                    end
                case 8
                    if ~isnan(vdStatsExt{a}.vd(b))
                        statsObj.ext.p8.vd = [statsObj.ext.p8.vd; vdStatsExt{a}.vd(b)];
                        statsObj.ext.p8.pDiff = [statsObj.ext.p8.pDiff; vdStatsExt{a}.pdiff(b)];
                    end
                case 12
                    if ~isnan(vdStatsExt{a}.vd(b))
                        statsObj.ext.p12.vd = [statsObj.ext.p12.vd; vdStatsExt{a}.vd(b)];
                        statsObj.ext.p12.pDiff = [statsObj.ext.p12.pDiff; vdStatsExt{a}.pdiff(b)];
                    end
                case 16
                    if ~isnan(vdStatsExt{a}.vd(b))
                        statsObj.ext.p16.vd = [statsObj.ext.p16.vd; vdStatsExt{a}.vd(b)];
                        statsObj.ext.p16.pDiff = [statsObj.ext.p16.pDiff; vdStatsExt{a}.pdiff(b)];
                    end
                case 20
                    if ~isnan(vdStatsExt{a}.vd(b))
                        statsObj.ext.p20.vd = [statsObj.ext.p20.vd; vdStatsExt{a}.vd(b)];
                        statsObj.ext.p20.pDiff = [statsObj.ext.p20.pDiff; vdStatsExt{a}.pdiff(b)];
                    end
                case 24
                    if ~isnan(vdStatsExt{a}.vd(b))
                        statsObj.ext.p24.vd = [statsObj.ext.p24.vd; vdStatsExt{a}.vd(b)];
                        statsObj.ext.p24.pDiff = [statsObj.ext.p24.pDiff; vdStatsExt{a}.pdiff(b)];
                    end
                case 28
                    if ~isnan(vdStatsExt{a}.vd(b))
                        statsObj.ext.p28.vd = [statsObj.ext.p28.vd; vdStatsExt{a}.vd(b)];
                        statsObj.ext.p28.pDiff = [statsObj.ext.p28.pDiff; vdStatsExt{a}.pdiff(b)];
                    end
                case 32
                    if ~isnan(vdStatsExt{a}.vd(b))
                        statsObj.ext.p32.vd = [statsObj.ext.p32.vd; vdStatsExt{a}.vd(b)];
                        statsObj.ext.p32.pDiff = [statsObj.ext.p32.pDiff; vdStatsExt{a}.pdiff(b)];
                    end
            end
        end
 %   end
    
end

% concatenate for box plot
% int
L = length(statsObj.int.p0.pDiff);
groupVec = ones(L,1)*0;
intBox = statsObj.int.p0.pDiff;

L = length(statsObj.int.p4.pDiff);
groupVec = [groupVec; ones(L,1)*4];
intBox = [intBox; statsObj.int.p4.pDiff];

L = length(statsObj.int.p8.pDiff);
groupVec = [groupVec; ones(L,1)*8];
intBox = [intBox; statsObj.int.p8.pDiff];

L = length(statsObj.int.p12.pDiff);
groupVec = [groupVec; ones(L,1)*12];
intBox = [intBox; statsObj.int.p12.pDiff];

L = length(statsObj.int.p16.pDiff);
groupVec = [groupVec; ones(L,1)*16];
intBox = [intBox; statsObj.int.p16.pDiff];

L = length(statsObj.int.p20.pDiff);
groupVec = [groupVec; ones(L,1)*20];
intBox = [intBox; statsObj.int.p20.pDiff];

L = length(statsObj.int.p24.pDiff);
groupVec = [groupVec; ones(L,1)*24];
intBox = [intBox; statsObj.int.p24.pDiff];

% [p,tbl,stats] = anova1(intBox, groupVec);
% [c,~,~,gnames] = multcompare(stats);

% ext
L = length(statsObj.ext.p0.pDiff);
groupVec1 = ones(L,1)*0;
extBox = statsObj.ext.p0.pDiff;

L = length(statsObj.ext.p4.pDiff);
groupVec1 = [groupVec1; ones(L,1)*4];
extBox = [extBox; statsObj.ext.p4.pDiff];

L = length(statsObj.ext.p8.pDiff);
groupVec1 = [groupVec1; ones(L,1)*8];
extBox = [extBox; statsObj.ext.p8.pDiff];

L = length(statsObj.ext.p12.pDiff);
groupVec1 = [groupVec1; ones(L,1)*12];
extBox = [extBox; statsObj.ext.p12.pDiff];

L = length(statsObj.ext.p16.pDiff);
groupVec1 = [groupVec1; ones(L,1)*16];
extBox = [extBox; statsObj.ext.p16.pDiff];

L = length(statsObj.ext.p20.pDiff);
groupVec1 = [groupVec1; ones(L,1)*20];
extBox = [extBox; statsObj.ext.p20.pDiff];

L = length(statsObj.ext.p24.pDiff);
groupVec1 = [groupVec1; ones(L,1)*24];
extBox = [extBox; statsObj.ext.p24.pDiff];

[p,tbl,stats] = anova1(extBox, groupVec1);
[c,~,~,gnames] = multcompare(stats);