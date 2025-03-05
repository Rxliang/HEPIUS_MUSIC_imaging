% Script to test cross validation of tree ensmble

predictorNames = {'balFact1'; 'balFact2'; 'balFact3'; 'balFact4'; 'balFact5'; 'balFact6'; 'balFact7'; 'balFact8'; 'balFact9'; 'balFact10'; 'balFact11'};

X = balSet.x;
Y = balSet.y';
L = length(Y);

% tree = fitrtree(X, Y, 'MaxNumSplits', 30, 'Prune', 'On', 'MinParentSize', 3, 'PredictorNames', predictorNames);
% yHat = predict(tree, X);
% c1 = corrcoef(Y, yHat)
% 
% return;

for a = 1:L
    yTrain = Y([1:a-1, a+1:L]); % leave one out
    xTrain = X([1:a-1, a+1:L], :); % leave one out
    yTest = Y(a);
    xTest = X(a, :);
    
    % PCA
    %[coeff,score,latent] = pca(xTrain, 'NumComponents', 5);
    
    %t = templateTree('MaxNumSplits',20,'Prune','on', 'MinParentSize', 3);
    %Mdl = fitensemble(xTrain, yTrain, 'LSBoost', numTrees,t, 'Type','regression', 'LearnRate', 0.25);
    %Mdl = TreeBagger(30, score, yTrain, 'Method','regression', 'MaxNumSplits', 20, 'MinLeafSize', 1);
    
    Mdl = fitrtree(xTrain, yTrain, 'MaxNumSplits', 30, 'Prune', 'On', 'MinParentSize', 3);
    %XtestPCA = (xTest - mean(xTrain,1))*coeff;
    %yHat(a) = Mdl.predict(XtestPCA);
    
    %Mdl = fitensemble(xTrain, yTrain, 'LSBoost', 10, 'Tree');
    %Mdl = fitrtree(xTrain, yTrain, 'MaxNumSplits', 100, 'Prune', 'On', 'MinParentSize', 1, 'PredictorNames', predictorNames);
    yHat(a) = predict(Mdl, xTest);
    
end
c2 = corrcoef(Y, yHat)


return;

MdlDeep = fitrtree(X,Y,'CrossVal','on','MergeLeaves','off',...
    'MinParentSize',1,'Surrogate','on');
MdlStump = fitrtree(X,Y,'MaxNumSplits',1,'CrossVal','on','Surrogate','on');

n = size(X,1);
m = floor(log2(n - 1));
lr = [0.1 0.25 0.5 1];
maxNumSplits = 2.^(0:m);
numTrees = 150;
Mdl = cell(numel(maxNumSplits),numel(lr));
rng(1); % For reproducibility
for k = 1:numel(lr);
    for j = 1:numel(maxNumSplits);
        t = templateTree('MaxNumSplits',maxNumSplits(j),'Surrogate','on');
        Mdl{j,k} = fitensemble(X,Y,'LSBoost',numTrees,t,...
            'Type','regression','KFold',5,'LearnRate',lr(k));
    end;
end;

kflAll = @(x)kfoldLoss(x,'Mode','cumulative');
errorCell = cellfun(kflAll,Mdl,'Uniform',false);
error = reshape(cell2mat(errorCell),[numTrees numel(maxNumSplits) numel(lr)]);
errorDeep = kfoldLoss(MdlDeep);
errorStump = kfoldLoss(MdlStump);

mnsPlot = [1 round(numel(maxNumSplits)/2) numel(maxNumSplits)];
figure;
for k = 1:3;
    subplot(2,2,k);
    plot(squeeze(error(:,mnsPlot(k),:)),'LineWidth',2);
    axis tight;
    hold on;
    h = gca;
    plot(h.XLim,[errorDeep errorDeep],'-.b','LineWidth',2);
    plot(h.XLim,[errorStump errorStump],'-.r','LineWidth',2);
    plot(h.XLim,min(min(error(:,mnsPlot(k),:))).*[1 1],'--k');
    h.YLim = [10 50];
    xlabel 'Number of trees';
    ylabel 'Cross-validated MSE';
    title(sprintf('MaxNumSplits = %0.3g', maxNumSplits(mnsPlot(k))));
    hold off;
end;
hL = legend([cellstr(num2str(lr','Learning Rate = %0.2f'));...
        'Deep Tree';'Stump';'Min. MSE']);
hL.Position(1) = 0.6;

[minErr,minErrIdxLin] = min(error(:));
[idxNumTrees,idxMNS,idxLR] = ind2sub(size(error),minErrIdxLin);

fprintf('\nMin. MSE = %0.5f',minErr)
fprintf('\nOptimal Parameter Values:\nNum. Trees = %d',idxNumTrees);
fprintf('\nMaxNumSplits = %d\nLearning Rate = %0.2f\n',...
    maxNumSplits(idxMNS),lr(idxLR))


