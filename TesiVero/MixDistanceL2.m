
function D = MixDistance(X,Y)
%
% Mixed categorical/numerical distance 
%
% INPUT:
% X = matrix of features, nObs x (nCategorical + nNumerical)
%        NOTE: categorical features must be
%                               1) grouped together
%                               2) the first block 
% Y = X (for most applications)
%
% OUTPUT:
% D = matrix of distances (nObsCat+nObsNum) x (nObsCat+nObsNum)

%% Find the number of categorical and numerical features
% The idea is that categorical variables are encoded, so they are
% represented by dummy/binary variables,
% and the sum of the possibile values == 1

nFeatures = size(X,2);
nCat = 0;
for i = 1:nFeatures
    if sum(unique(X(:,i))) == 1
        nCat = nCat + 1;
    end
end
nNum = nFeatures - nCat;

%% Compute distances, separately
DCat = pdist2(X(:,1:9), Y(:,1:9), 'hamming');
DInv = pdist2(X(:,10:11), Y(:,10:11), 'hamming');
DNum = pdist2(X(:,12:end), Y(:,12:end), 'euclidean');
% Compute relative weight based on the number of categorical variables
coeffInv = 1;
wCat = (nCat - 2)/(nCat - 2 + nNum + (2 * coeffInv)); 
wNum = nNum/(nCat - 2 + nNum + (2 * coeffInv)); 
wInv = (2 * coeffInv)/(nCat - 2 + nNum + (2 * coeffInv)); 
D = wCat*DCat + wNum*DNum + DInv*wInv;
end


