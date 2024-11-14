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
LifeStyle	Normalized Lifestyle Index (percentiles);
Luxury	Normalized Luxury spending (percentiles)
Saving	Normalized Saving propensity (percentiles)
Investments	Investment type: 1 = no investments; 2 = mostly lump sum; 3 = mostly capital accumulation
%}

Datatoplot=table2array(BankClients(:,(8:17)))*100;

Datatoplot=horzcat(table2array(BankClients(:,2)),Datatoplot);

Datatoplot=horzcat(table2array(BankClients(:,7)),Datatoplot);

%numerical variables

xnames = {'familysize','age','income','wealth', 'debt', 'finedu','esg','digital','bankfriend','lifestyle','luxury','savings'};


%%

rownames=BankClients.Properties.VariableNames;

rownames=rownames(:,(2:18));

BankClients_ml=table(BankClients.Age, BankClients.Gender, BankClients.Job, BankClients.Area, BankClients.CitySize,...
    BankClients.FamilySize,  BankClients.Income*100,...
    BankClients.Wealth*100, BankClients.Debt*100, BankClients.FinEdu*100, BankClients.ESG*100,...
    BankClients.Digital*100, BankClients.BankFriend*100, BankClients.LifeStyle*100, BankClients.Luxury*100,...
    BankClients.Saving*100, BankClients.Investments, 'VariableNames',rownames);

%% Es 1 - Analisi Monovariata Preliminare

% Andiamo a graficare ogni parametro di Van Cauter con ogni covariata per farci un'idea su quali potrebbero essere 
% le covariate maggiormente relazionate con la variabile dipendente e sul tipo di relazione che potrebbe legarle.

% Le covariate "Stato di salute" e "Sesso" sono variabili categoriche, e per questo motivo le rappresentiamo con dei 
% boxplot. Le altre covariate sono variabili continue che possono essere rappresentate con degli scatterplot.

% Prima di fare una regressione lineare --> si fanno tanti scatter plot o boxplot per vedere se le variabili sono tra loro correlate, almeno a coppie.

% inseriamo in 2 vettori i nomi dei parametri e dei regressori (covariate)
nomi_reg = {'Age' 'Gender' 'Job' 'Area' 'CitySize' 'FamilySize'...
    'Income' 'Wealth' 'Debt' 'FinEdu' 'ESG' 'Digital' 'BankFriend'...
    'LifeStyle' 'Luxury' 'Savings' 'Age*Gender' 'Age*Job' 'Age*Area'...
    'Age*CitySize' 'Age*FamilySize' 'Age*Income' 'Age*Wealth' 'Age*Debt'...
    'Age*FinEdu' 'Age*ESG' 'Age*Digital' 'Age*BankFriend' 'Age*LifeStyle'...
    'Age*Luxury' 'Age*Savings' 'Gender*Job' 'Gender*Area' 'Gender*CitySize'...
    'Gender*FamilySize' 'Gender*Income' 'Gender*Wealth' 'Gender*Debt'...
    'Gender*FinEdu' 'Gender*ESG' 'Gender*Digital' 'Gender*BankFriend'...
    'Gender*LifeStyle' 'Gender*Luxury' 'Gender*Savings' 'Job*Area'...
    'Job*CitySize' 'Job*FamilySize' 'Job*Income' 'Job*Wealth' 'Job*Debt'...
    'Job*FinEdu' 'Job*ESG' 'Job*Digital' 'Job*BankFriend' 'Job*LifeStyle'...
    'Job*Luxury' 'Job*Savings' 'Area*CitySize' 'Area*FamilySize' 'Area*Income'...
    'Area*Wealth' 'Area*Debt' 'Area*FinEdu' 'Area*ESG' 'Area*Digital'...
    'Area*BankFriend' 'Area*LifeStyle' 'Area*Luxury' 'Area*Savings'...
    'CitySize*FamilySize' 'CitySize*Income' 'CitySize*Wealth' 'CitySize*Debt'...
    'CitySize*FinEdu' 'CitySize*ESG' 'CitySize*Digital' 'CitySize*BankFriend'...
    'CitySize*LifeStyle' 'CitySize*Luxury' 'CitySize*Savings' 'FamilySize*Income'...
    'FamilySize*Wealth' 'FamilySize*Debt' 'FamilySize*FinEdu' 'FamilySize*ESG'...
    'FamilySize*Digital' 'FamilySize*BankFriend' 'FamilySize*LifeStyle'...
    'FamilySize*Luxury' 'FamilySize*Savings' 'Income*Wealth' 'Income*Debt'...
    'Income*FinEdu''Income*ESG' 'Income*Digital' 'Income*BankFriend' 'Income*LifeStyle'...
    'Income*Luxury' 'Income*Savings' 'Wealth*Debt' 'Wealth*FinEdu' 'Wealth*ESG'...
    'Wealth*Digital' 'Wealth*BankFriend' 'Wealth*LifeStyle' 'Wealth*Luxury' 'Wealth*Savings'...
    'Debt*FinEdu' 'Debt*ESG' 'Debt*Digital' 'Debt*BankFriend' 'Debt*LifeStyle'...
    'Debt*Luxury' 'Debt*Savings' 'FinEdu*ESG' 'FinEdu*Digital' 'FinEdu*BankFriend'...
    'FinEdu*LifeStyle' 'FinEdu*Luxury' 'FinEdu*Savings' 'ESG*Digital' 'ESG*BankFriend'...
    'ESG*LifeStyle' 'ESG*Luxury' 'ESG*Savings' 'Digital*BankFriend' 'Digital*LifeStyle'...
    'Digital*Luxury' 'Digital*Savings' 'BankFriend*LifeStyle' 'BankFriend*Luxury'...
    'BankFriend*Savings' 'LifeStyle*Luxury' 'LifeStyle*Savings' 'Luxury*Savings'};
nomi_par = {'Investments'};

dummy = [2:5 17:20 32:81]; % dove ho le dummy variables 
%%
% divido parametri (y) da covariate (x)
x = BankClients_ml(:,1:16); % regressori
y = table2array(BankClients_ml(:,17)); % parametri

for i=1:width(y)
    figure('Name',nomi_par{i})
    for j=1:width(x)
        subplot(4,4,j)
        if j>1 && j<6
            %boxplot(table2array(x(:,i)),table2array(y(:,j))); hold on;
%             switch j
%                 case 2
%                     xlabel({'0=normale, 1=obeso, 2=diabetico'})
%                 case 3
%                     xlabel({'0=M, 1=F'})
%             end
        else
            scatter(table2array(x(:,j)),y(:,i)), grid on;
            xlabel(nomi_reg{j})
        end
        ylabel(nomi_par{i})
        title(nomi_reg{j})
    end
end



%% Es 2 - Modello di regressione + stima dei parametri + CV e intervalli di confidenza

% Dobbiamo determinare il modello di regressione --> per farlo dobbiamo
% scegliere un sottoinsieme adatto di regressori

% per far ciò dobbiamo implementare a mano un algoritmo, che è la stepwise
% E’ la soluzione più usata quando ci sono parecchi regressori e quindi la ricerca attraverso 
% l’elencazione di tutti i possibili modelli è troppo onerosa per essere effettuata

% Usiamo la stepwise forward:
% Nella procedura forward si parte da un modello semplice e si aggiungono (se opportuno) 
% i regressori più promettenti uno ad uno.

% Si sceglie un criterio per decidere di volta in volta quale è il modello
% migliore, ad esempio l'adjusted R^2.

% L'R^2 tiene conto del numero di parametri --> ci permette di valutare il modello migliore 
% tenendo conto della complessità del modello.
% adjusted R^2 = R^2 - k/(n-k-1)(1- R^2)


% Utilizzo le dummy variables per categoriche

gndr=BankClients_ml.Gender;

gender=dummyvar(gndr);

gender=gender(:, 1:end-1);

jb=BankClients_ml.Job;

job=dummyvar(jb);

job=job(:, 1:end-1);

ar=BankClients_ml.Area;

area=dummyvar(ar);

area=area(:, 1:end-1);

ctsz=BankClients_ml.CitySize;

citysize=dummyvar(ctsz);

citysize=citysize(:, 1:end-1);

%%

%% PROSEGUIRE QUI
% definizione cell array con i regressori candidati
H_cand = {gender [dati(:,2)] [dati(:,3)] [dati(:,4)] [dati(:,5)] [dati(:,6)] [dati(:,7)] [salute.*dati(:,3)] [salute.*dati(:,4)] [salute.*dati(:,5)] [salute.*dati(:,6)] [salute.*dati(:,7)] [dati(:,2).*dati(:,3)] [dati(:,2).*dati(:,4)] [dati(:,2).*dati(:,5)] [dati(:,2).*dati(:,6)] [dati(:,2).*dati(:,7)] [dati(:,3).*dati(:,4)] [dati(:,3).*dati(:,5)] [dati(:,3).*dati(:,6)] [dati(:,3).*dati(:,7)] [salute.*dati(:,2)] [salute.*dati(:,2).*dati(:,3)] [salute.*dati(:,2).*dati(:,4)] [salute.*dati(:,2).*dati(:,5)] [salute.*dati(:,2).*dati(:,6)] [salute.*dati(:,2).*dati(:,7)]};

% Matrice sigma = identità --> non definito il modello dell'errore

ind_corr = cell(1,length(nomi_par));  % preallocazione
X = cell(1,length(nomi_par));  % preallocazione
beta = cell(1,length(nomi_par));  % preallocazione
int_conf = cell(1,length(nomi_par));  % preallocazione
sigma_p = cell(1,length(nomi_par));  % preallocazione
CV_beta = cell(1,length(nomi_par));  % preallocazione
ind_ok = cell(1,length(nomi_par));  % preallocazione
F = cell(1,length(nomi_par));  % preallocazione
Vc = cell(1,length(nomi_par));  % preallocazione
X_end = cell(1,length(nomi_par));  % preallocazione
beta_end = cell(1,length(nomi_par));  % preallocazione
Y_hat_new = cell(1,length(nomi_par));  % preallocazione
reg_modello = cell(1,length(nomi_par));  % preallocazione
modello = cell(1,length(nomi_par));  % preallocazione

disp('Modelli di Regressione Lineare Multipla')

% Ciclo sui 4 parametri (y)

for m = 1:width(y)

    ind_cand = 1:length(nomi_reg);
    ind_corr{m} = 0;

    clear s_beta Y betaciclo R2 R2adj 

    X{m} = ones(length(y),1);

    Y = y(:,m); 
    maxR2_old(m) = 0;  %inizializzo

    for j = 1:length(ind_cand)
        clear R2_adj R2
        for i = 1:length(ind_cand)
            [betaciclo{i}, R2(i), R2_adj(i)] = regressione([X{m} H_cand{ind_cand(i)}],Y);  
        end
        [maxR2, indice] = max(R2_adj);
        ind_max = ind_cand(indice);

        if maxR2 > maxR2_old(m)
           ind_cand(ind_cand == ind_max) = [];  % tolgo dai candidati
           ind_corr{m} = [ind_corr{m}, ind_max];
           maxR2_old(m) = maxR2;
           X{m} = [X{m} H_cand{ind_max}];
        else 
           break
        end
    end

    beta{m} = ((X{m}'*X{m})^-1)*X{m}'*Y;  % parametri del modello

    k = width(X{m})-1; % Numero di regressori (tolgo 1 perchè la prima è la colonna di uni e non la considero)
    n = length(X{m});  % Numero di osservazioni

    % L'incertezza della stima la descrivo attraverso la matrice di covarianza delle stime (non dell'errore), 
    % la sigmap, matrice di dimensioni nXn, che è uguale a (H'*sigma_v^-1*H)^-1
    % con sigma_v = matrice identità * α e H = X;
    % lo scalare α è necessario per avere la precisione delle stime che
    % richiede la conoscenza di Σv.
    % α lo calcoliamo a posteriori come il rapporto tra WRSS(p), cioè la somma dei quadrati dell'errore, e i gradi di libertà
    Y_hat{m} = X{m}*beta{m};
    res = Y-Y_hat{m};

    SSE_full = sum((Y - Y_hat{m}).^2);  % Somma dei quadrati dell'errore (modello completo, cioè prima del test di Fisher)
    coeff_alfa = SSE_full/(n-k-1);
    
    sigma_p{m} = full(inv(X{m}'*coeff_alfa*eye(length(dati))*X{m}));

    % deviazione standard dei beta stimati
    s_beta = sqrt(diag(sigma_p{m}));  % dev standard dei beta (troviamo le varianze sulla diagonale della matrice di covarianza delle stime)

    % spesso si preferisce riportare l'incertezza in termini relativi --> calcoliamo il CV della stima
    CV_beta{m} = (s_beta./beta{m})*100;

  
    % Intervalli di confidenza: β_hat+-t(n-k-1, alfa/2)*s_β
    alfa = 0.05;

    t_crit = tinv(alfa/2, n-k-1);  % dove n-k-1 rappresenta i gradi di libertà e α/2 il livello di significatività

    p_inf = zeros(length(beta{m}),1);  % preallocazione
    p_sup = zeros(length(beta{m}),1);  
    int_conf{m} = zeros(length(beta{m}),2);

    % Intervallo di confidenza per ogni coefficiente beta della regressione multipla
    for i = 1:length(beta{m})
        p_inf(i)=beta{m}(i)+(tinv(alfa/2, n-k-1))*s_beta(i); 
        p_sup(i)=beta{m}(i)-(tinv(alfa/2, n-k-1))*s_beta(i);  
        int_conf{m}(i,:) = [p_inf(i) p_sup(i)];
    end


    % Valutiamo con il test F (a posteriori sul nostro sottoinsieme di regressori) la significatività dei parametri.
    % La statistica F indica la capacità predittiva del modello di regressione determinando se tutti i coefficienti 
    % di regressione del modello sono significativamente diversi da 0.
    % Spesso la procedura è questa:
    % Identificare un modello con molti regressori, vedere le variabili che mostrano un p-value non piccolo a 
    % quindi fare il test F eliminando quelle variabili per decidere se eliminarle o no dal modello completo.
    
    % H0: β_k = 0

    p = k;  % numero di regressori del modello completo
    clear k

    for i = 2:length(ind_corr{m})  % da 2 perché il primo si riferisce a beta_0
        ind_rid = ind_corr{m};
        ind_rid(i) = [];

        X_rid{m} = ones(length(y),1);
        
        for j = 2:length(ind_rid)
            X_rid{m} = [X_rid{m} H_cand{ind_rid(j)}];
        end

        k = width(X_rid{m})-1;  % numero di regressori del modello ridotto 

        [beta_rid, ~, ~] = regressione(X_rid{m},Y); 

        Y_hat_new{m}(:,i) = X_rid{m} * beta_rid;
        SSE_red = sum((Y - Y_hat_new{m}(:,i)).^2);

        % Statistica del test: F = [(SSEreduced – SSEfull)/(p – k)]/[SSEfull/(n – p – 1)]
        F{m}(i) = ((SSE_red - SSE_full)/(p-k))/(SSE_full/(n-p-1));  % statistica
        Vc{m} = finv(0.95,p-k,n-p-1);
        
        % valuto se la statistica del test F sta a destra o sinistra del valore critico, per capire se rientra
        % nella zona di rifiuto o di accettazione dell'ipotesi nulla H0:
        if F{m}(i) > Vc{m}
            % se la statistica cade nella zona di rifiuto --> rifiuto l'H0 e dico che non ho elementi per dire 
            % beta non è significativo e lo tengo.
            ind_ok{m} = [ind_ok{m} ind_corr{m}(i)];
        end
    end
    X_end{m} = ones(size(y,1),1);  % preallocazione

    for i = 1:length(ind_ok{m})
        X_end{m} = [X_end{m} H_cand{ind_ok{m}(i)}];
    end

    beta_end{m} = ((X_end{m}'*X_end{m})^-1)*X_end{m}'*Y;  % parametri finali del modello
    
    % Stampiamo sullo schermo, per ogni parametro di Van Cauter, il rispettivo modello di regressione lineare
    plus = {' + '};    
    for i = 1:length(ind_ok{m})
        reg_modello{m}{i} = nomi_reg{ind_ok{m}(i)};
    end

    modello{m} = num2str(1);

    for i = 1:length(reg_modello{m})
        modello{m} = [modello{m}  plus{1} reg_modello{m}{i}];
    end
    disp([nomi_par{m}, ' ~ ', modello{m}])
    
end



%% Es 3 - Parametrizzazione personale: predizione dei parametri della cinetica del CP usando i nostri dati come regressori

% Il miglior modello di regressione lineare multipla che abbiamo
% determinato è:

myY = cell(1,length(nomi_par));  % preallocazione
myX = cell(1,length(nomi_par));  % preallocazione

Clinical = [0 0];
Sex = 0;
Age = 21;
Altezza = 1.77;
Peso = 73;

BMI = Peso/(Altezza^2);
BSA = 0.20247*Altezza^0.725*Peso^0.425;

miei_dati = [0 Sex Age Altezza Peso BMI BSA];

totX = {[Clinical] [miei_dati(:,2)] [miei_dati(:,3)] [miei_dati(:,4)] [miei_dati(:,5)] [miei_dati(:,6)] [miei_dati(:,7)] [Clinical.*miei_dati(:,3)] [Clinical.*miei_dati(:,4)] [Clinical.*miei_dati(:,5)] [Clinical.*miei_dati(:,6)] [Clinical.*miei_dati(:,7)] [miei_dati(:,2).*miei_dati(:,3)] [miei_dati(:,2).*miei_dati(:,4)] [miei_dati(:,2).*miei_dati(:,5)] [miei_dati(:,2).*miei_dati(:,6)] [miei_dati(:,2).*miei_dati(:,7)] [miei_dati(:,3).*miei_dati(:,4)] [miei_dati(:,3).*miei_dati(:,5)] [miei_dati(:,3).*miei_dati(:,6)] [miei_dati(:,3).*miei_dati(:,7)] [Clinical.*miei_dati(:,2)] [Clinical.*miei_dati(:,2).*miei_dati(:,3)] [Clinical.*miei_dati(:,2).*miei_dati(:,4)] [Clinical.*miei_dati(:,2).*miei_dati(:,5)] [Clinical.*miei_dati(:,2).*miei_dati(:,6)] [Clinical.*miei_dati(:,2).*miei_dati(:,7)]};

for m = 1:4
    myX{m} = 1;
    for i = 1:length(ind_ok{m})
        myX{m} = [myX{m} totX{ind_ok{m}(i)}];
    end
    myY{m} = myX{m}*beta_end{m};
end

V = myY{1};  % Volume di distribuzione
t_corto = myY{2};  % tempo di emivita corto (di distribuzione)
t_lungo = myY{3};  % tempo di emivita lungo (biologica)
Fraction = myY{4};  % Fraction


% 1 Parametrizzazione: Parametri di Von Cauter
disp(' ')
disp(['I-Parametri di Van Cauter: V = ',num2str(V),';  T_corto = ',num2str(t_corto),';  T_lungo = ',num2str(t_lungo),';  Fraction = ',num2str(Fraction),'.']);


% 2 Parametrizzazione: Parametri del modello bicompartimentale

% La cinetica del C-peptide è descritta da un modello bicompartimentale;
% quindi la risposta impulsiva, cioè la concentrazione dato un impulso unitario come bolo, è del tipo A1e^−α1t + A2e^−α2t
% Dobbiamo ricavare i parametri A1, A2, alfa1 e alfa2

alfa1 = 0.693/t_corto; % tempo di emivita di distribuzione = ln2/α1
alfa2 = 0.693/t_lungo; % tempo di emivita biologica = ln2/α2

% Volume di Distribuzione: V = 1/(A1+A2) [Al numeratore avremmo la dose D, ma stiamo considerando la risposta all'impulso unitario]
% Fraction = A1/(A1+A2)

% Risolvendo il sistema di due equazioni in due incognite, otteniamo:
A1 = Fraction/V;
A2 = 1/V - A1;

disp(' ')
disp(['II-Parametri della cinetica a 2 esponenziali: A1 = ',num2str(A1),';  A2 = ',num2str(A2),';  α1 = ',num2str(alfa1),';  α2 = ',num2str(alfa2),'.']);


% 3 Parametrizzazione: modello con clearance

% trovo le clearance unica cosa che manca nella nuova parametrizzazione, Vd gia noto 

Cl_1=V*alfa1;
Cl_2=V*alfa2;

disp(' ')
disp(['III-Parametri della cinetica a 2 esponenziali con la clearance: V= ',num2str(V),';  CL1 = ',num2str(Cl_1),';  CL2 = ',num2str(Cl_2),'.']);

