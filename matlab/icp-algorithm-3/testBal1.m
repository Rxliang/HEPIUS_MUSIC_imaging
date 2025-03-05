
for a = 1:length(balSet.y)
    f = [];
    count = 1;
    yy = nanmedian(balSetExt{a}.PIenv') - nanmedian(balSetInt{a}.PIenv');
    P = pressureList{a};
    if length(P) < 2 || length(yy) < 2 || sum(~isnan(yy)) < 2
        P = [];
    end
    if ~isempty(P)
        if length(P) ~= length(yy)
            k = min([length(P) length(yy)]);
            P = P(1:k);
            yy = yy(1:k);
        end
        
        ss = sum(isnan(yy));
        if ss > 0
            fNan = find(~isnan(yy));
            yy = yy(fNan);
            if ~isempty(P)
                P = P(fNan);
            end
            
        end
        
        L = length(yy);
        if yy(2) > yy(1)
            f(count) = 1;
            count = count + 1;
        end
        if yy(L) < yy(L-1)
            f(count) = L;
            count = count + 1;
        end
        if L == 5
            if yy(3) > max(yy(1:2))
                f(count) = 3;
                count = count + 1;
            end
            if yy(3) > max(yy(L-1:L))
                f(count) = 3;
                count = count + 1;
            end
        end
        if L == 6
            [m,I] = min(yy(3:4));
            b = 0;
            if m > max(yy(1:2))
                f(count) = 3;
                count = count + 1;
                f(count) = 4;
                count = count + 1;
                b = 1;
            end
            if m > max(yy(L-1:L))
                f(count) = 3;
                count = count + 1;
                f(count) = 4;
                count = count + 1;
                b = 1;
            end
            if ~b
                [m,I] = max(yy(3:4));
                if m > max(yy(1:2))
                    f(count) = 2+I;
                    count = count + 1;
                end
            end
            if m > max(yy(L-1:L))
                f(count) = 2+I;
                count = count + 1;
            end
        end
%         if L > 6
%             warndlg('More than 6 pressure?')
%         end
        if ~isempty(f)
            f = unique(f);
            f = find(~ismember([1:1:L],f));
        end
        
        
        pVec = [P(1):.1:P(end)];
        if ~isempty(f)
            [Pf, S] = polyfit(P(f), yy(f), 2);
        else
            [Pf, S] = polyfit(P, yy, 2);
        end
        polynomFit = polyval(Pf, pVec);
        [m, I] = min(polynomFit);
        ICP(a) = pVec(I);
    else
        ICP(a) = NaN;
    end
end
