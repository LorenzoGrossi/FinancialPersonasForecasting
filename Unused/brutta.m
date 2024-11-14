%% Import

close all;
clc;
clear;
path = '/Users/Lorenzo/Desktop/Tesi/BankClients.mat';

load(path)

%% Visualize dependencies of numeric variables

%{
ID	Numerical ID
Age	Age, in years
Gender	Gender (Female = 1, Male = 0)
Job	1=non occupato, 2=impiegato/operaio, 3=quadro/dirigente, 4=imprenditore/libero professionista, 5=pensionato)
Area	1=Nord, 2=Centro, 3=Sud/Isole
CitySize	1=piccolo centro, 2=cittadina medie dimensioni, 3=grande città >200k
FamilySize	Number of components
Income	Normalized Income (percentiles)
Wealth	Normalized Wealth (percentiles)
Debt	Normalized Debt (percentiles)
FinEdu	Normalized Financial Education (percentiles)
ESG	Normalized ESG propensity (percentiles)
Digital	Normalized Digital propensity (percentiles)
BankFriend	Normalized Bank Friendliness (percentiles)
LifeStyle	Normalized Lifestyle Index (percentiles); higher values = 
Luxury	Normalized Luxury spending (percentiles)
Saving	Normalized Saving propensity (percentiles)
Investments	Investment type: 1 = no investments; 2 = mostly lump sum; 3 = mostly capital accumulation
%}

Datatoplot=table2array(BankClients(:,(7:17)));

Datatoplot=horzcat(table2array(BankClients(:,2)),Datatoplot);

%numerical variables

xnames = {'Age','familysize','income','wealth', 'debt', 'finedu','esg','digital','bankfriend','lifestyle','luxury','savings'};

%%
figure
gplotmatrix(Datatoplot,[],BankClients.Gender, [],[],[],[],'grpbars', xnames)
title('Grouping by Gender')

% gender does not show any particular relation

%%
figure
gplotmatrix(Datatoplot,[],BankClients.Investments, [],[],[],[],'grpbars', xnames)
title('Grouping by Investments')

% SEE WITH DEBT
% investments does not show any particular relation with others var

%%
figure
gplotmatrix(Datatoplot,[],BankClients.Job, [],[],[],[],'grpbars', xnames)
title('Grouping by Job')

%interesting; to be investigated

%%
figure
gplotmatrix(Datatoplot,[],BankClients.Area, [],[],[],[],'grpbars', xnames)
title('Grouping by Area')

% FORSE NORD ALZA I VALORI (LUX,...)?

%%
figure
gplotmatrix(Datatoplot,[],BankClients.CitySize, [],[],[],[],'grpbars', xnames)
title('Grouping by CitySize')

% FORSE GRANDE CITTà ALZA

%% si passa a PCA per vedere gli score maggiori da quali variabili sono composti

% Normalizza i dati
age_da_standardizzare = table2array(BankClients(:, 2));
famsize_da_standardizzare = table2array(BankClients(:, 7));

% Applica rescale solo alla colonna estratta
age_standardizzata = rescale(age_da_standardizzare);
famsize_standardizzata = rescale(famsize_da_standardizzare);

% Aggiorna la colonna estratta con i valori standardizzati
BankClients_st = table2array(BankClients(:,(8:17)));

BankClients_st=horzcat(age_standardizzata,famsize_standardizzata,BankClients_st);

%%

rownames=BankClients.Properties.VariableNames;

rownames=rownames(:,(2:18));

BankClients_ml=table(age_standardizzata, BankClients.Gender, BankClients.Job,...
    BankClients.Area, BankClients.CitySize, famsize_standardizzata, BankClients.Income...
    ,BankClients.Wealth, BankClients.Debt, BankClients.FinEdu, BankClients.ESG...
    ,BankClients.Digital, BankClients.BankFriend, BankClients.LifeStyle, BankClients.Luxury...
    ,BankClients.Saving, BankClients.Investments, 'VariableNames',rownames);

%% Esegui PCA sui dati normalizzati

[coeff, score, latent, ~, explained] = pca(BankClients_st);

% Continua con l'analisi...

% Supponiamo che 'coeff' contenga i coefficienti delle componenti principali
% e che 'BankClients_std' sia il dataset standardizzato

% Imposta il numero di componenti da visualizzare
num_components = 5;

% Visualizza le prime due componenti principali
for i = 1:num_components
    figure;
    bar(coeff(:, i)); % Grafico a barre dei coefficienti della componente principale i
    xlabel('Variabile'); % Etichetta dell'asse x
    ylabel('Coefficiente'); % Etichetta dell'asse y
    title(['Coefficienti della Componente Principale ', num2str(i)]); % Titolo del grafico
    xticklabels(xnames); % Etichette dell'asse x con i nomi delle variabili (assumendo che tu abbia un vettore di nomiVar)
    xtickangle(45); % Angolo di rotazione delle etichette dell'asse x per una migliore leggibilità
    grid on; % Attiva la griglia
end

varianza_spiegata_cumulativa = cumsum(explained);

% Visualizza la varianza spiegata cumulativamente
plot(varianza_spiegata_cumulativa, 'bo-');

% non sembra che ci siano variabili poco significative; 

%% correlazioni

correlazione = corrcoef(BankClients_st);

corr_1 = correlazione-eye(size(correlazione,1));

% Visualizza la matrice di correlazione come una heatmap
figure;
heatmap(correlazione,'XData',xnames,'YData',xnames);
title('Matrice di Correlazione');

max(max((abs(corr_1)))); %0.5666, non ci sono variabili che correlano a livelli da considerare di scartarne alcune

%% controlli di variabili singolarmente, per vedere se sono opportune trasformazioni

figure
histogram(BankClients.Age)
title('AGE')

%%

[h, pval, w]=swtest(BankClients.Age, 0.01); %reject

%%
figure
histogram(BankClients.Income)
title('Income')
%%
[h, pval, w]=swtest(BankClients.Income, 0.1) %reject

%%
figure
histogram(BankClients.FamilySize)
title('FamilySize')
%%
[h, pval, w]=swtest(BankClients.FamilySize, 0.1) %reject
%%
figure
histogram(BankClients.Wealth)
title('Wealth')
%%
[h, pval, w]=swtest(BankClients.Wealth, 0.1) %reject
%% qui cambiare, provare trasformazioni
figure
histogram(BankClients.Debt)
title('Debt')
%%
[h, pval, w]=swtest(BankClients.Debt, 0.1) %reject
%%
figure
histogram(BankClients.FinEdu)
title('FinEdu')

%%
[h, pval, w]=swtest(BankClients.FinEdu, 0.1) %reject
%%
figure
histogram(BankClients.ESG)
title('ESG')

%%
figure
histogram(BankClients.Digital)
title('Digital')

%%
figure
histogram(BankClients.BankFriend)
title('BankFriend')

%%
figure
histogram(BankClients.LifeStyle)
title('lifestyle')
%%  controlla
figure
histogram(BankClients.Luxury)
title('Luxury')

%%
figure
histogram(BankClients.Saving)
title('savings')

%% ML ALGORITHMS

% divido il dataset;

nObs = size(BankClients_ml, 1);
idxPermutation = randperm(nObs);
X = BankClients_ml(idxPermutation,(1:16)); % random permutation
Y = BankClients_ml(idxPermutation,17); % random permutation
nTrain = round(0.60*nObs);
nPureCV = round(0.20*nObs);
nTest = round(0.20*nObs);

Xtrain = X((1:nTrain),:);
XpureCV = X((nTrain+1:nPureCV+nTrain),:);
Xtest = X((nTrain+nPureCV+1:nObs),:);

Ytrain = Y((1:nTrain),:);
YpureCV = Y((nTrain+1:nPureCV+nTrain),:);
Ytest = Y((nTrain+nPureCV+1:nObs),:);

Train = BankClients_ml((1:nTrain),:);
PureCV = BankClients_ml((nTrain+1:nPureCV+nTrain),:);
Test = BankClients_ml((nTrain+nPureCV+1:nObs),:);

%% testo modelli di machine learning

% Definisci i modelli di machine learning che desideri testare
% Ad esempio, supponiamo che tu voglia testare un modello di regressione lineare e un classificatore SVM
%mdl_linear = fitlm(Xtrain, Ytrain); % Modello di regressione lineare
mdl_svm = fitcecoc(Xtrain, Ytrain); % Classificatore SVM

% Valuta le prestazioni dei modelli sui dati di test
%pred_linear = predict(mdl_linear, Xtest);
pred_svm = predict(mdl_svm, Xtest);

pred_svm_num = categorical(pred_svm);

Ytest_num = table2array(Ytest);

% Calcola le metriche di valutazione (ad esempio, per problemi di regressione puoi usare l'errore quadratico medio, per problemi di classificazione l'accuratezza, la precisione, il richiamo, ecc.)
% Esempio per un problema di regressione
%mse_linear = mean((pred_linear - Ytest).^2);

%mse_svm = mean(1-(double(pred_svm) == Ytest).^2);

% Esempio per un problema di classificazione
%accuracy_linear = sum(pred_linear == Ytest) / numel(Ytest);
accuracy_svm = sum(pred_svm_num == Ytest_num) / numel(Ytest);

%%
% Puoi anche utilizzare funzioni MATLAB integrate come 'confusionmat' per problemi di classificazione per ottenere la matrice di confusione e calcolare altre metriche di valutazione

% Ora puoi valutare i modelli utilizzando la cross-validation
% Per esempio, usando k-fold cross-validation
num_folds = 5; % Numero di fold per la cross-validation
cv_linear = crossval(mdl_linear, 'KFold', num_folds);
cv_svm = crossval(mdl_svm, 'KFold', num_folds);

% Calcola le metriche di valutazione della cross-validation (ad esempio, l'errore medio, l'accuratezza media, ecc.)

% Ora puoi analizzare i risultati e confrontare le prestazioni dei diversi modelli

