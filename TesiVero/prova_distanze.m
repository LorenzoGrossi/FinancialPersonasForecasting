
close all
%clear all
clc
path = '/Users/Lorenzo/Desktop/Tesi/BankClients.mat';

load BankClients.mat
%indici_da_rimuovere = [222, 1007, 2813];

% Rimozione delle righe specificate
%BankClients(indici_da_rimuovere, :) = [];
%BankClients = BankClients(~[222,1007,2813],:);
% Numerical
N = vartype('numeric');
Data = BankClients(:,2:end); % Exclude 1st col = ID
% Removing the gender and the family size
%Data(:, "Gender") = [];
%Data(:, "FamilySize") = [];
NumFeatures = Data(:,N); % subtable of numerical features
XNum = NumFeatures.Variables; % create a matrix 
XNum = rescale(XNum); % normalize in [0, 1]
% Categorical
C = vartype('categorical');
CatFeatures = Data(:,C); % subtable of categorical features
% Encoding (ie create dummy variables)
%grp2idx create a vector of integer values from 1 up to the number K of distinct groups.
XCat = [];
for i = 1:size(CatFeatures, 2)
    XCat_i = dummyvar(grp2idx(CatFeatures{:,i}));
    XCat = [XCat XCat_i(:,1:end-1)]; % To use the dummy variables in a model, we must delete a column
end

X = [XCat XNum];

CatFeatures = table2array(CatFeatures);
NumFeatures = rescale(table2array(NumFeatures));




%% tecnica Zenti

%% Compute distances, separately
DCat = pdist2(X(:,1:11), X(:,1:11), 'hamming');
DNum = pdist2(X(:,12:end), X(:,12:end), 'euclidean');
% Compute relative weight based on the number of categorical variables
%wCat = nCat/(nCat + nNum); 
%D = wCat*DCat + (1 - wCat)*DNum;


%% tecnica nostra

% Compute distances, separately
DCat2 = pdist2(CatFeatures(:,1:end-1), CatFeatures(:,1:end-1), 'hamming');
DInv2 = pdist2(CatFeatures(:,end), CatFeatures(:,end), 'hamming');
DNum2 = pdist2(NumFeatures, NumFeatures, 'euclidean');

DCatDef2 = DCat2 + DInv2;

% Compute relative weight based on the number of categorical variables
coeffInv = 1;
wCat2 = (nCat - 2)/(nCat - 2 + nNum + (2 * coeffInv)); 
wNum2 = nNum/(nCat - 2 + nNum + (2 * coeffInv)); 
wInv2 = (2 * coeffInv)/(nCat - 2 + nNum + (2 * coeffInv)); 
D2 = wCat2*DCat2 + wNum2*DNum2 + DInv2*wInv2;


