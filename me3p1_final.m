clear all
close all
clc
profile on

disp('Empirical pset 1, Metrics 3');

%% 0. Importing Data

data = readtable("DataTeam3&4.xlsx");

unemp = table2array(data(:,3));
gdp = table2array(data(:,4));
gtsp = table2array(data(:,5)).*100;
ivbp = table2array(data(:,8));
c1 = table2array(data(:,6));
c2 = table2array(data(:,7));
qtrs = table2array(data(:,2));

% ts1 = timeseries(unemp,qtrs,'name','Unemployment');
% ts2 = timeseries(gdp,qtrs,'name','GDP');
% ts3 = timeseries(gtsp,qtrs,'name','Government Spending');
% ts4 = timeseries(ivbp,qtrs,'name','Blanchard-Perotti IV');
% tsc = tscollection({ts1 ts2 ts3 ts4},'name','okun');

%% Data Exploration

disp('Section 1');
disp('Data Plots')

%sub 1

f11 = figure;
set(f11,'Visible','off');
subplot(2,2,1);
plot(qtrs,unemp);
%plot(ts1,'-xb','Displayname',ts1.Name)
title('Subplot 1: Unemployment series')
xlabel('Quarters')
ylabel('Unemployment percentage')
grid on


subplot(2,2,2);
plot(qtrs,gdp)
%plot(ts2,'-.xm','Displayname',ts2.Name)
title('Subplot 2: GDP series')
xlabel('Quarters')
ylabel('GDP growth percentage')
grid on

subplot(2,2,3);
plot(qtrs,gtsp)
%plot(ts3,'-.xm','Displayname',ts3.Name)
title('Subplot 3: Govt. Spending series')
xlabel('Quarters')
ylabel('Govt. Spending in percentage of GDP')
grid on

subplot(2,2,4);
plot(qtrs,ivbp)
%plot(ts4,'-.xm','Displayname',ts4.Name)
title('Subplot 4: BP series')
xlabel('Quarters')
ylabel('BP in ---')
grid on

% legend('show','Location','NorthWest')
% xlabel('Time (Quarters)')
% ylabel('')
% grid on

% saveas(f11,'Figure 1.1.png');
% disp('Figure 1.1')

%--------------------------------------------------------------------------

%sub 2

f12 = figure;
set(f12,'Visible','off');
subplot(2,2,1);
autocorr(unemp)
title('Subplot 1: Unemployment series')

subplot(2,2,2);
autocorr(gdp)
title('Subplot 2: GDP series')

subplot(2,2,3);
autocorr(gtsp)
title('Subplot 3: Govt. Spending series')

subplot(2,2,4);
autocorr(ivbp)
title('Subplot 4: BP series')

% saveas(f12,'Figure 1.2.png');
% disp('Figure 1.2')

%--------------------------------------------------------------------------

%sub 3

f13 = figure;
set(f13,'Visible','off');
subplot(2,2,1);
histfit(unemp,25,'kernel');
line([mean(unemp), mean(unemp)], ylim, 'LineWidth',1,'Color','r','LineStyle','-.')
line ([mean(unemp)+std(unemp) mean(unemp)+std(unemp) NaN mean(unemp)-std(unemp) mean(unemp)-std(unemp)] , [ylim NaN   ylim],'LineWidth', 0.5, 'Color', 'g','Displayname','St. Dev.')
a=annotation('textbox',...
    [0.30 0.8 0.5 0.04],...
    'String',{'Mean = 6.6401', 'Variance = 16.132', 'Skewness = 2.1696','Kurtosis = 8.3348'},...
    'FitBoxToText','on','LineStyle','none');
a.FontSize=5.5;
t=title('Unemployment series');
t.FontSize=10;
t.FontWeight = 'bold';


subplot(2,2,2);
histfit(gdp,25,'kernel');
line([mean(gdp), mean(gdp)], ylim, 'LineWidth', 1, 'Color','r','LineStyle','-.')
line ([mean(gdp)+std(gdp) mean(gdp)+std(gdp) NaN mean(gdp)-std(gdp) mean(gdp)-std(gdp)] , [ylim NaN   ylim],'LineWidth', 0.5, 'Color', 'g','Displayname','St. Dev.')
a=annotation('textbox',...
    [0.57 0.8 0.5 0.04],...
    'String',{'Mean = 0.8037', 'Variance = 5.1100', 'Skewness = -0.4877','Kurtosis = 6.0247'},...
    'FitBoxToText','on','LineStyle','none');
a.FontSize=5.5;
t=title('GDP series');
t.FontSize=10;
t.FontWeight = 'bold';

subplot(2,2,3);
histfit(gtsp,20,'kernel');
line([mean(gtsp), mean(gtsp)], ylim, 'LineWidth', 1, 'Color', 'r','LineStyle','-.')
line ([mean(gtsp)+std(gtsp) mean(gtsp)+std(gtsp) NaN mean(gtsp)-std(gtsp) mean(gtsp)-std(gtsp)] , [ylim NaN   ylim],'LineWidth', 0.5, 'Color', 'g','Displayname','St. Dev.')
a=annotation('textbox',...
    [0.3 0.3 0.5 0.04],...
    'String',{'Mean = 16.9780', 'Variance = 89.031', 'Skewness = 1.9299','Kurtosis = 10.1116'},...
    'FitBoxToText','on','LineStyle','none');
a.FontSize=5.5;
t=title('Govt. Spending series');
t.FontSize=10;
t.FontWeight = 'bold';

subplot(2,2,4);
histfit(ivbp,25,'kernel')
line([mean(ivbp), mean(ivbp)], ylim, 'LineWidth', 1, 'Color', 'r','LineStyle','-.')
line ([mean(ivbp)+std(ivbp) mean(ivbp)+std(ivbp) NaN mean(ivbp)-std(ivbp) mean(ivbp)-std(ivbp)] , [ylim NaN   ylim],'LineWidth', 0.5, 'Color', 'g','Displayname','St. Dev.')
a=annotation('textbox',...
    [0.57 0.3 0.5 0.04],...
    'String',{'Mean = -5.0201e-07', 'Variance = 1.2948e-04', 'Skewness = 0.1305','Kurtosis = 20.4250'},...
    'FitBoxToText','on','LineStyle','none');
a.FontSize=5.5;
lgd=legend('Distribution','Kernel Density','Mean','St.Dev','Location',[0.27 0.01 0.5 0.04],'Orientation','horizontal');
lgd.FontSize = 6;
t=title('BP series');
t.FontSize=10;
t.FontWeight = 'bold';

saveas(f13,'Figure 1.3.png');
% disp('Figure 1.3')


%% Estimation

disp('Estimation');

one = ones(503,1);
AvarTemp2 = 0;
deltaTemp2 = 0;
AvarTemp31 = 0;
AvarTemp32 = 0;
deltaTemp31 = 0;
deltaTemp32 = 0;
AvarTemp41 = 0;
AvarTemp42 = 0;
deltaTemp41 = 0;
deltaTemp42 = 0;

f = waitbar(0,'Please wait...');
pause(.5)

for h = 0:10
    
    y2 = gdp(5+h:end);
    X2 = [unemp(5:end-h) one(5:end-h) unemp(4:end-h-1) unemp(3:end-h-2) unemp(2:end-h-3) unemp(1:end-h-4) gdp(4:end-h-1) gdp(3:end-h-2) gdp(2:end-h-3) gdp(1:end-h-4) gtsp(4:end-h-1) gtsp(3:end-h-2) gtsp(2:end-h-3) gtsp(1:end-h-4) c1(4:end-h-1) c1(3:end-h-2) c1(2:end-h-3) c1(1:end-h-4) c2(4:end-h-1) c2(3:end-h-2) c2(2:end-h-3) c2(1:end-h-4)];

    b2(:,h+1) = (X2'*X2)\(X2'*y2);
    
    B = diag((y2(:)-X2*b2(:,h+1)),0);
    Avar = inv(X2'*X2)*((X2'*B)*(X2'*B)')*inv(X2'*X2);
    cf21(h+1) = b2(1,h+1) + 1.96*sqrt(Avar(1,1));
    cf22(h+1) = b2(1,h+1) - 1.96*sqrt(Avar(1,1));
    
    delta2(h+1) = b2(1,h+1) + deltaTemp2;
    deltaTemp2 = delta2(h+1);
    
    cfd21(h+1) = delta2(h+1) + 1.96*sqrt(Avar(1,1)+AvarTemp2);     %for delta
    cfd22(h+1) = delta2(h+1) - 1.96*sqrt(Avar(1,1)+AvarTemp2);
    AvarTemp2 = AvarTemp2 + Avar(1,1);
    
    %3
    y31 = gdp(5+h:end);
    X31 = [gtsp(5:end-h) one(5:end-h) unemp(4:end-h-1) unemp(3:end-h-2) unemp(2:end-h-3) unemp(1:end-h-4) gdp(4:end-h-1) gdp(3:end-h-2) gdp(2:end-h-3) gdp(1:end-h-4) gtsp(4:end-h-1) gtsp(3:end-h-2) gtsp(2:end-h-3) gtsp(1:end-h-4) c1(4:end-h-1) c1(3:end-h-2) c1(2:end-h-3) c1(1:end-h-4) c2(4:end-h-1) c2(3:end-h-2) c2(2:end-h-3) c2(1:end-h-4)];
        
    y32 = unemp(5+h:end);
    X32 = [gtsp(5:end-h) one(5:end-h) unemp(4:end-h-1) unemp(3:end-h-2) unemp(2:end-h-3) unemp(1:end-h-4) gdp(4:end-h-1) gdp(3:end-h-2) gdp(2:end-h-3) gdp(1:end-h-4) gtsp(4:end-h-1) gtsp(3:end-h-2) gtsp(2:end-h-3) gtsp(1:end-h-4) c1(4:end-h-1) c1(3:end-h-2) c1(2:end-h-3) c1(1:end-h-4) c2(4:end-h-1) c2(3:end-h-2) c2(2:end-h-3) c2(1:end-h-4)];
        
    b31(:,h+1) = (X31'*X31)\(X31'*y31);
    b32(:,h+1) = (X32'*X32)\(X32'*y32);
    
    B31 = diag((y31(:)-X31*b31(:,h+1)),0);
    Avar31 = inv(X31'*X31)*((X31'*B31)*(X31'*B31)')*inv(X31'*X31);
    cf311(h+1) = b31(1,h+1) + 1.96*sqrt(Avar31(1,1));
    cf312(h+1) = b31(1,h+1) - 1.96*sqrt(Avar31(1,1));
    
    B32 = diag((y32(:)-X32*b32(:,h+1)),0);
    Avar32 = inv(X32'*X32)*((X32'*B32)*(X32'*B32)')*inv(X32'*X32);
    cf321(h+1) = b32(1,h+1) + 1.96*sqrt(Avar32(1,1));
    cf322(h+1) = b32(1,h+1) - 1.96*sqrt(Avar32(1,1));
    
    delta31(h+1) = b31(1,h+1) + deltaTemp31;
    deltaTemp31 = delta31(h+1);
    
    delta32(h+1) = b32(1,h+1) + deltaTemp32;
    deltaTemp32 = delta32(h+1);
    
    cfd311(h+1) = delta31(h+1) + 1.96*sqrt(Avar31(1,1)+AvarTemp31);     %for delta
    cfd312(h+1) = delta31(h+1) - 1.96*sqrt(Avar31(1,1)+AvarTemp31);
    AvarTemp31 = AvarTemp31 + Avar31(1,1);
    
    cfd321(h+1) = delta32(h+1) + 1.96*sqrt(Avar32(1,1)+AvarTemp32);     %for delta
    cfd322(h+1) = delta32(h+1) - 1.96*sqrt(Avar32(1,1)+AvarTemp32);
    AvarTemp32 = AvarTemp32 + Avar32(1,1);
    
    kappa3(h+1) = delta31(h+1)/delta32(h+1);
    
    deltamethod3(:,h+1) = [1/delta32(h+1) -delta31(h+1)/delta32(h+1)^2];
    AvarK312 = inv(X32'*X32)*((X32'*B31)*(X32'*B32)')*inv(X32'*X32);
    AvarK321 = inv(X32'*X32)*((X32'*B32)*(X32'*B31)')*inv(X32'*X32);
    BK3(:,:,h+1) = [AvarTemp31 AvarK312(1,1); AvarK321(1,1) AvarTemp32];
    AvarK3(h+1) = deltamethod3(:,h+1)'*BK3(:,:,h+1)*deltamethod3(:,h+1);
    
    cfk31(h+1) = kappa3(h+1) + 1.96*sqrt(AvarK3(h+1));
    cfk32(h+1) = kappa3(h+1) - 1.96*sqrt(AvarK3(h+1));
    
    
    %4
    Z41 = [ivbp(5:end-h) one(5:end-h) unemp(4:end-h-1) unemp(3:end-h-2) unemp(2:end-h-3) unemp(1:end-h-4) gdp(4:end-h-1) gdp(3:end-h-2) gdp(2:end-h-3) gdp(1:end-h-4) gtsp(4:end-h-1) gtsp(3:end-h-2) gtsp(2:end-h-3) gtsp(1:end-h-4) c1(4:end-h-1) c1(3:end-h-2) c1(2:end-h-3) c1(1:end-h-4) c2(4:end-h-1) c2(3:end-h-2) c2(2:end-h-3) c2(1:end-h-4)];
    X41 = X31;
    y41 = y31;
    
    Z42 = [ivbp(5:end-h) one(5:end-h) unemp(4:end-h-1) unemp(3:end-h-2) unemp(2:end-h-3) unemp(1:end-h-4) gdp(4:end-h-1) gdp(3:end-h-2) gdp(2:end-h-3) gdp(1:end-h-4) gtsp(4:end-h-1) gtsp(3:end-h-2) gtsp(2:end-h-3) gtsp(1:end-h-4) c1(4:end-h-1) c1(3:end-h-2) c1(2:end-h-3) c1(1:end-h-4) c2(4:end-h-1) c2(3:end-h-2) c2(2:end-h-3) c2(1:end-h-4)];    
    X42 = X32;
    y42 = y32;
    
    b41(:,h+1) = inv((X41'*Z41)*inv(Z41'*Z41)*(Z41'*X41))*((X41'*Z41)*inv(Z41'*Z41)*(Z41'*y41));
    b42(:,h+1) = inv((X42'*Z42)*inv(Z42'*Z42)*(Z42'*X42))*((X42'*Z42)*inv(Z42'*Z42)*(Z42'*y42));
        
    B41 = diag((y41(:)-X41*b41(:,h+1)),0);
    Avar41 = inv((X41'*Z41)*inv(Z41'*Z41)*(Z41'*X41))*((X41'*Z41)*inv(Z41'*Z41)*((Z41'*B41)*(Z41'*B41)')*inv(Z41'*Z41)*(Z41'*X41))*inv((X41'*Z41)*inv(Z41'*Z41)*(Z41'*X41));
    cf411(h+1) = b41(1,h+1) + 1.96*sqrt(Avar41(1,1));
    cf412(h+1) = b41(1,h+1) - 1.96*sqrt(Avar41(1,1));
    
    B42 = diag((y42(:)-X42*b42(:,h+1)),0);
    Avar42 = inv((X42'*Z42)*inv(Z42'*Z42)*(Z42'*X42))*((X42'*Z42)*inv(Z42'*Z42)*((Z42'*B42)*(Z42'*B42)')*inv(Z42'*Z42)*(Z42'*X42))*inv((X42'*Z42)*inv(Z42'*Z42)*(Z42'*X42));
    cf421(h+1) = b42(1,h+1) + 1.96*sqrt(Avar42(1,1));
    cf422(h+1) = b42(1,h+1) - 1.96*sqrt(Avar42(1,1));
    
    delta41(h+1) = b41(1,h+1) + deltaTemp41;
    deltaTemp41 = delta41(h+1);
    
    delta42(h+1) = b42(1,h+1) + deltaTemp42;
    deltaTemp42 = delta42(h+1);
    
    cfd411(h+1) = delta41(h+1) + 1.96*sqrt(Avar41(1,1)+AvarTemp41);     %for delta
    cfd412(h+1) = delta41(h+1) - 1.96*sqrt(Avar41(1,1)+AvarTemp41);
    AvarTemp41 = AvarTemp41 + Avar41(1,1);
    
    cfd421(h+1) = delta42(h+1) + 1.96*sqrt(Avar42(1,1)+AvarTemp42);     %for delta
    cfd422(h+1) = delta42(h+1) - 1.96*sqrt(Avar42(1,1)+AvarTemp42);
    AvarTemp42 = AvarTemp42 + Avar42(1,1);
    
    kappa4(h+1) = delta41(h+1)/delta42(h+1);
    
    deltamethod4(:,h+1) = [1/delta32(h+1) -delta31(h+1)/delta32(h+1)^2];
    AvarK412 = inv((X42'*Z42)*inv(Z42'*Z42)*(Z42'*X42))*((X42'*Z42)*inv(Z42'*Z42)*((Z42'*B41)*(Z42'*B42)')*inv(Z42'*Z42)*(Z42'*X42))*inv((X42'*Z42)*inv(Z42'*Z42)*(Z42'*X42));
    AvarK421 = inv((X42'*Z42)*inv(Z42'*Z42)*(Z42'*X42))*((X42'*Z42)*inv(Z42'*Z42)*((Z42'*B42)*(Z42'*B41)')*inv(Z42'*Z42)*(Z42'*X42))*inv((X42'*Z42)*inv(Z42'*Z42)*(Z42'*X42));
    BK4(:,:,h+1) = [AvarTemp41 AvarK412(1,1); AvarK421(1,1) AvarTemp42];
    AvarK4(h+1) = deltamethod4(:,h+1)'*BK4(:,:,h+1)*deltamethod4(:,h+1);
    
    cfk41(h+1) = kappa4(h+1) + 1.96*sqrt(AvarK4(h+1));
    cfk42(h+1) = kappa4(h+1) - 1.96*sqrt(AvarK4(h+1));


    
    %Checking the first stage
    
    X41_mean = ones(499-h,1)*mean(X41(:,1));
    bfs4(:,h+1)= inv(Z41'*Z41)*Z41'*X41(:,1);
    SSreg = (Z41*bfs4(:,h+1)-X41_mean)'*(Z41*bfs4(:,h+1)-X41_mean);
    SStot = (X41(:,1)-X41_mean)'*(X41(:,1)-X41_mean);
    Rfs(h+1) = SSreg/SStot; %This is essentially one for all the
    %regressions, so we are gettig perfect fit in the first stage. Then in
    %the second stage when we substitute zp by its fitted value, it makes
    %no difference compared to OLS.
    
    %z only on controls
    
    W = Z41(:,2:end);
    bzw(:,h+1) = inv(W'*W)*W'*X41(:,1);
    SSregzw = (W*bzw(:,h+1)-X41_mean)'*(W*bzw(:,h+1)-X41_mean);
    SStotzw = (X41(:,1)-X41_mean)'*(X41(:,1)-X41_mean);
    Rzw(h+1) = SSregzw/SStotzw; 
    
    waitbar(h*0.1,f,'Calc co-eff');
    pause(1)
end
close(f)

disp('Done');

%% Plots

disp('Creating Plots')

    f21 = figure;
    set(f21,'Visible','off');
    plot(linspace(0,10,11),b2(1,:),linspace(0,10,11),cf21(:),linspace(0,10,11),cf22(:))
    str = sprintf('beta');
    title(str)
    xlabel('h')
    legend('Beta','Upper bound','Lower Bound','Location','NorthWest')
%     saveas(f211,'Figure 2.1.1.png');

    f22 = figure;
    set(f22,'Visible','off');
    plot(linspace(0,10,11),delta2(:),linspace(0,10,11),cfd21(:),linspace(0,10,11),cfd22(:))
    str = sprintf('delta');
    title(str)
    xlabel('h')
    legend('Delta','Upper bound','Lower Bound','Location','North')
%     saveas(f22,'Figure 2.2.png');

    f31 = figure;
    set(f31,'Visible','off');
    plot(linspace(0,10,11),b31(1,:),linspace(0,10,11),cf311(:),':',linspace(0,10,11),cf312(:),':')
    str = sprintf('beta31');
    title(str)
    xlabel('h')
%     saveas(f31,'Figure 3.1.png');

    f32 = figure;
    set(f32,'Visible','off');
    plot(linspace(0,10,11),b32(1,:),linspace(0,10,11),cf321(:),':',linspace(0,10,11),cf322(:),':')
    str = sprintf('beta32');
    title(str)
    xlabel('h')
%     saveas(f32,'Figure 3.2.png');

    f33 = figure;
    set(f33,'Visible','off');
    plot(linspace(0,10,11),kappa3(:),linspace(0,10,11),cfk31(:),linspace(0,10,11),cfk32(:))
    str = sprintf('kappa');
    title(str)
    xlabel('h')
%     saveas(f33,'Figure 3.3.png');

    f41 = figure;
    set(f41,'Visible','off');
    plot(linspace(0,10,11),b41(1,:),linspace(0,10,11),cf411(:),':',linspace(0,10,11),cf412(:),':')
    str = sprintf('beta41');
    title(str)
    xlabel('h')
%     saveas(f41,'Figure 4.1.png');

    f42 = figure;
    set(f42,'Visible','off');
    plot(linspace(0,10,11),b42(1,:),linspace(0,10,11),cf421(:),':',linspace(0,10,11),cf422(:),':')
    str = sprintf('beta42');
    title(str)
    xlabel('h')
%     saveas(f42,'Figure 4.2.png');

   f43 = figure;
    set(f43,'Visible','off');
    plot(linspace(0,10,11),kappa4(:),linspace(0,10,11),cfk41(:),linspace(0,10,11),cfk42(:))
    str = sprintf('kappa4');
    title(str)
    xlabel('h')
%     saveas(f43,'Figure 4.3.png');

    f44 = figure;
    set(f44,'Visible','off');
    plot(linspace(0,10,11),kappa4(:),linspace(0,10,11),kappa3(:))
    str = sprintf('kappa VS kappa4');
    title(str)
    xlabel('h')
%     saveas(f44,'Figure 4.4.png');

disp('Saving Plots')

%% Appendix

%Save all plots
saveas(f11,'Figure 1.1.png');
saveas(f12,'Figure 1.2.png');
saveas(f13,'Figure 1.3.png');
saveas(f21,'Figure 2.1.1.png');
saveas(f22,'Figure 2.1.2.png');
saveas(f31,'Figure 3.1.png');
saveas(f32,'Figure 3.2.png');
saveas(f33,'Figure 3.3.png');
saveas(f41,'Figure 4.1.png');
saveas(f42,'Figure 4.2.png');
saveas(f43,'Figure 4.3.png');
saveas(f44,'Figure 4.4.png');

profile viewer
%profsave;