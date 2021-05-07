% Written by Angel Rodes -- CIAF-SUERC
% angelrodes@gmail.com
% First version: February 2014
% Last updated May, 2021

%% Clear
clear % clear previous data
close all hidden % close previous graphs

%% Parameters
tm10=1390000; % 10Be half life
tm26=705000; % 26Al half life
l10=log(2)/tm10; % lambda 10Be
l26=log(2)/tm26; % lambda 26Al
Rinit2610=6.75; % 26Al/10Be production rate ratio (e.g. Balco & Rovey http://dx.doi.org/10.2475/10.2008.02)
Rinit2110=4.08; % 21Ne/10Be production rate ratio (e.g. Balco & Shuster http://dx.doi.org/10.1016/j.epsl.2009.02.006)

% Production rates
P10=10.116; % random value
P26=P10*Rinit2610;
P21=P10*Rinit2110;
RDSPrate=0.08; % percentage of error in the individual production rates (0.05 is 5%). Ratios will have an uncertanty of ~RDSPrate*1.41


%% Choose profile data

xls_files = dir(fullfile('*.xls'));
XLS_files = dir(fullfile('*.XLS'));
xls_str = [{'HELP'} {xls_files.name} {XLS_files.name}];
[s,v] = listdlg('Name','Input data','PromptString','Select a .xls file:',...
    'SelectionMode','single',...
    'ListString',xls_str);
profilefile=xls_str(s);
profilefile=profilefile{1};
if s==1
    stringhelp=sprintf('The XLS file must contain ONLY the following information: \n Column 1: sample name \n Column 2: 10Be concentration \n Column 3: 10Be concentration error \n Column 4: 26Al concentration \n Column 5: 26Al concentration error \n Column 6: 21Ne concentration \n Column 7: 21Ne concentration error \n Column 8: depth \n Column 9: depth error \n The first row (column titles) will be ignored. \n \nSelect only the samples containing data for the 3 isotopes.');
    helpdlg(stringhelp,'Help')
    return
end

[num, txt, raw]= xlsread(profilefile);
numsamples=size(num,1)-1;
Samplename=txt(2:numsamples+2,1);
N10=num(1:numsamples+1,1);
d10=num(1:numsamples+1,2);
N26=num(1:numsamples+1,3);
d26=num(1:numsamples+1,4);
N21=num(1:numsamples+1,5);
d21=num(1:numsamples+1,6);
Depth=num(1:numsamples+1,7);
dDepth=num(1:numsamples+1,8);


%% Plot profile

figure(1)
hold on

% Build legend
plot([0 0],[0 0],'Color',[0.8 0.8 0.8]);
plot(0,0,'ok');
plot(0,0,'ob');
plot(0,0,'or');
legend(['Spallation' char(10) 'attenuation'],'[^{10}Be]','[^{26}Al]','[^{21}Ne]','Location','NorthWest','AutoUpdate','off')

% Spallation reference
for C0=logspace(log10(mean(N10)),log10(100*mean(N26)),10)
    z=-linspace(0,max(Depth+dDepth),100);
    plot(C0*exp(z*1.9/160),z,'Color',[0.8 0.8 0.8]);
end

% samples
for sample=1:numsamples+1
    theta=linspace(0,2*pi,100);
    x10=N10(sample)+d10(sample)*cos(theta);
    x26=N26(sample)+d26(sample)*cos(theta);
    x21=N21(sample)+d21(sample)*cos(theta);
    y=-(Depth(sample)+dDepth(sample)*sin(theta));
    plot(x10,y,'-k','LineWidth',2)
    plot(x26,y,'-b','LineWidth',2)
    plot(x21,y,'-r','LineWidth',2)
    text(max(x10),min(y),Samplename(sample));
    text(max(x26),min(y),Samplename(sample));
    text(max(x21),min(y),Samplename(sample));
end

% Beautify
set(gca, 'xscale', 'log')
ylim([-300 25])
xlim([min(N10)/2 max(N21)*2])
xlabel('Concentration [atoms/g]')
ylabel('Depth [cm]')


%% select samples
[samples,v] = listdlg('Name','Input data','PromptString','Select samples:',...
    'SelectionMode','multiple',...
    'ListString',Samplename);
if length(samples)<2
    errordlg('I need at least 2 samples!','Error');
    return
end

%% Linear 3D regression ( Mote Carlo )
Nx=N10(samples);
dx=d10(samples);
Ny=N26(samples);
dy=d26(samples);
Nz=N21(samples);
dz=d21(samples);

% Mote Carlo approach (models)
nmodels=50000; % number of models to run
Rs=0*(1:nmodels);
onesigma=0*(1:nmodels);
Ps=zeros(nmodels,3);
P0s=zeros(nmodels,3);
for model=1:nmodels
    xyz=[random('norm',Nx,dx),random('norm',Ny,dy),random('norm',Nz,dz)]; % random concentrations from gaussian distributions
    if max(max(((xyz-[Nx,Ny,Nz])./[dx,dy,dz])).^2)<1 % check if model is inside 1 sigma limits
        onesigma(model)=1;
    else
        onesigma(model)=0;
    end
    P0 = mean(xyz,1); % position of the line (3D linear regression)
    [~,~,v] = svd(xyz-repmat(P0,size(xyz,1),1)); 
    P = v(:,1)'; % slope (3D linear regression)
    
    P0s(model,:)=P0; % save value
    Ps(model,:)=P; % save value
    Rs(model)=(1-P(2)/P(3)/(P26/P21))/(1-P(1)/P(3)/(P10/P21)); % Calculate R*
end
%% Calculate ages from models

% Formula to calculate ages from R
tref=logspace(3,7.4,1000); % t reference
Rref=(1-exp(-l26*tref))./(1-exp(-l10*tref)); % R* reference
ages=interp1(Rref,tref,Rs); % Interpolate R from models in the reference formula

% Formula to calculate ages from R (including production rate error)
P102=random('norm',P10+0*Ps(:,1),P10*RDSPrate+0*Ps(:,1)); % random (gaussian) 10Be production rates
P262=random('norm',P26+0*Ps(:,2),P26*RDSPrate+0*Ps(:,1)); % random (gaussian) 26Al production rates
P212=random('norm',P21+0*Ps(:,3),P21*RDSPrate+0*Ps(:,1)); % random (gaussian) 21Ne production rates
Rs2=(1-Ps(:,2)./Ps(:,3)./(P262./P212))./(1-Ps(:,1)./Ps(:,3)./(P102./P212)); % R* with random production rate values
ages2=interp1(Rref,tref,Rs2); % Interpolate R from models in the reference formula

%% plot 3D results
figure(33)
hold on
% Plot approaches
t=[-1e10,1e10];
for model=1:10:nmodels
    P0m=P0s(model,:);
    Pm=Ps(model,:);
    if onesigma(model)==1 % plot only the "best" models (otherwise the samples will be hidden)
        plot3(P0s(model,1)+t*Ps(model,1),P0s(model,3)+t*Ps(model,3),P0s(model,2)+t*Ps(model,2),'Color',[0.8 0.8 0.8])
    end
    Rs(model)=(1-P(2)/P(3)/(P26/P21))/(1-P(1)/P(3)/(P10/P21));
    
end

% best fit: P(t) = P0 + t*P
t=linspace(-max(xyz(:)),max(xyz(:)),10000);
P=mean(Ps);
P0=mean(P0s);
X=P0(1)+t*P(1);
Y=P0(2)+t*P(2);
Z=P0(3)+t*P(3);
positive=find((X>0).*(Y>0).*(Z>0));
X=X(positive); Y=Y(positive); Z=Z(positive);
plot3(X,Z,Y,'-b','LineWidth',1) 
% Production rate ratios
X=t*P10;
Y=t*P26;
Z=t*P21;
positive=find((X>0).*(Y>0).*(Z>0));
X=X(positive); Y=Y(positive); Z=Z(positive);
plot3(X,Z,Y,'-k','LineWidth',2)
% legend('slope','Production rate ratios')
% sample elipsoids
xyz=[Nx,Ny,Nz];
plot3(xyz(:,1),xyz(:,3),xyz(:,2),'.k')
[theta phi]=meshgrid(linspace(0,2*pi,30),linspace(0,2*pi,30));
for sample=1:length(samples)
    plot3(Nx(sample)+dx(sample)*cos(theta).*sin(phi),Nz(sample)+dz(sample)*cos(phi),Ny(sample)+dy(sample)*sin(theta).*sin(phi),'k')
    text(Nx(sample),Nz(sample),Ny(sample),['   ' Samplename(samples(sample))])
end

xlabel('[^{10}Be]')
zlabel('[^{26}Al]')
ylabel('[^{21}Ne]')
xlim([0 max(Nx+dx)*1.2])
zlim([0 max(Ny+dy)*1.2])
ylim([0 max(Nz+dz)*1.2])
grid on
% box on

% %% Plot simple ratios/Rinit (projections of the 3D lines)
% figure(2)
% 
% subplot(1,3,1) % Be-Al
% Nx=N10(samples);
% dx=d10(samples);
% Ny=N26(samples);
% dy=d26(samples);
% Ratio=P26/P10;
% hold on
% theta=linspace(0,2*pi,100);
% plot(Nx,Ny,'.k');
% plot(max(Nx)*[0,1.2],max(Nx)*[0,1.2]*Ratio,'-k');
% xindex=1;
% yindex=2;
% for model=find(onesigma==1)
% plot(P0s(model,xindex)+t*Ps(model,xindex),P0s(model,yindex)+t*Ps(model,yindex),'-','Color',[0.8 0.8 0.8])
% end
% Ri=Ps(:,yindex)./Ps(:,xindex);
% p = polyfit(Nx,Ny,1);
% yfit =  p(1) * Nx + p(2);
% yresid = Ny - yfit;
% SSresid = sum(yresid.^2);
% SStotal = (length(Ny)-1) * var(Ny);
% rsq = 1 - SSresid/SStotal;
% plot(Nx,yfit,'--b')
% for sample=1:length(samples)
%     x=(Nx(sample)+dx(sample)*cos(theta));
%     y=(Ny(sample)+dy(sample)*sin(theta));
%     plot(x,y,'-k')
%     text(mean(x),mean(y),Samplename(sample));
% end
% tref=logspace(3,7,1000);
% Rref=exp(-abs(l10)*tref);
% meanR=mean(Ri)/Ratio;
% dmeanR=std(Ri)/Ratio;
% meantb=interp1(Rref,tref,meanR);
% dmeantb=abs(interp1(Rref,tref,meanR+dmeanR)-interp1(Rref,tref,meanR-dmeanR));
% title(['R/Rinit=' num2str(meanR,3) ' +/- ' num2str(dmeanR,3) ' t_{b}=' num2str(meantb/1e6,3)  ' +/- ' num2str(dmeantb/1e6,3) ' Ma'])
% xlabel('N10 [atoms/g]')
% ylabel('N26 [atoms/g]')
% ylim([0 max(y)*1.5])
% t2610=[num2str(meantb/1e6,3)  ' +/- ' num2str(dmeantb/1e6,3) ' Ma, R2=' num2str(rsq)]
% 
% subplot(1,3,2) % Ne-Be
% Nx=N21(samples);
% dx=d21(samples);
% Ny=N10(samples);
% dy=d10(samples);
% Ratio=P10/P21;
% hold on
% theta=linspace(0,2*pi,100);
% plot(Nx,Ny,'.k');
% plot(max(Nx)*[0,1.2],max(Nx)*[0,1.2]*Ratio,'-k');
% xindex=3;
% yindex=1;
% for model=find(onesigma==1)
% plot(P0s(model,xindex)+t*Ps(model,xindex),P0s(model,yindex)+t*Ps(model,yindex),'-','Color',[0.8 0.8 0.8])
% end
% Ri=Ps(:,yindex)./Ps(:,xindex);
% p = polyfit(Nx,Ny,1);
% yfit =  p(1) * Nx + p(2);
% yresid = Ny - yfit;
% SSresid = sum(yresid.^2);
% SStotal = (length(Ny)-1) * var(Ny);
% rsq = 1 - SSresid/SStotal;
% plot(Nx,yfit,'--b')
% for sample=1:length(samples)
%     x=(Nx(sample)+dx(sample)*cos(theta));
%     y=(Ny(sample)+dy(sample)*sin(theta));
%     plot(x,y,'-k')
%     text(mean(x),mean(y),Samplename(sample));
% end
% tref=logspace(3,7,1000);
% Rref=exp(-abs(l10-l26)*tref);
% meanR=mean(Ri)/Ratio;
% dmeanR=std(Ri)/Ratio;
% meantb=interp1(Rref,tref,meanR);
% dmeantb=abs(interp1(Rref,tref,meanR+dmeanR)-interp1(Rref,tref,meanR-dmeanR));
% title(['R/Rinit=' num2str(meanR,3) ' +/- ' num2str(dmeanR,3) ' t_{b}=' num2str(meantb/1e6,3)  ' +/- ' num2str(dmeantb/1e6,3) ' Ma'])
% xlabel('N21 [atoms/g]')
% ylabel('N10 [atoms/g]')
% ylim([0 max(y)*1.5])
% xdatasimple=[mean(Ri) std(Ri)]/Ratio;
% t1021=[num2str(meantb/1e6,3)  ' +/- ' num2str(dmeantb/1e6,3) ' Ma, R2=' num2str(rsq)]
% 
% subplot(1,3,3) % Ne-Al
% Nx=N21(samples);
% dx=d21(samples);
% Ny=N26(samples);
% dy=d26(samples);
% Ratio=P26/P21;
% hold on
% theta=linspace(0,2*pi,100);
% plot(Nx,Ny,'.k');
% plot(max(Nx)*[0,1.2],max(Nx)*[0,1.2]*Ratio,'-k');
% xindex=3;
% yindex=2;
% for model=find(onesigma==1)
% plot(P0s(model,xindex)+t*Ps(model,xindex),P0s(model,yindex)+t*Ps(model,yindex),'-','Color',[0.8 0.8 0.8])
% end
% Ri=Ps(:,yindex)./Ps(:,xindex);
% p = polyfit(Nx,Ny,1);
% yfit =  p(1) * Nx + p(2);
% yresid = Ny - yfit;
% SSresid = sum(yresid.^2);
% SStotal = (length(Ny)-1) * var(Ny);
% rsq = 1 - SSresid/SStotal;
% plot(Nx,yfit,'--b')
% for sample=1:length(samples)
%     x=(Nx(sample)+dx(sample)*cos(theta));
%     y=(Ny(sample)+dy(sample)*sin(theta));
%     plot(x,y,'-k')
%     text(mean(x),mean(y),Samplename(sample));
% end
% tref=logspace(3,7,1000);
% Rref=exp(-abs(l26)*tref);
% meanR=mean(Ri)/Ratio;
% dmeanR=std(Ri)/Ratio;
% meantb=interp1(Rref,tref,meanR);
% dmeantb=abs(interp1(Rref,tref,meanR+dmeanR)-interp1(Rref,tref,meanR-dmeanR));
% title(['R/Rinit=' num2str(meanR,3) ' +/- ' num2str(dmeanR,3) ' t_{b}=' num2str(meantb/1e6,3)  ' +/- ' num2str(dmeantb/1e6,3) ' Ma'])
% xlabel('N21 [atoms/g]')
% ylabel('N26 [atoms/g]')
% ylim([0 max(y)*1.5])
% ydatasimple=[mean(Ri) std(Ri)]/Ratio;
% t2621=[num2str(meantb/1e6,3)  ' +/- ' num2str(dmeantb/1e6,3) ' Ma, R2=' num2str(rsq)]


%% plot mote carlo results (ages)
figure(34)
hold on
% lages=log10(ages);
[n,xout] =hist(ages,20);
% hist(ages,20)
bar(xout,n);
title('Age (a)')
pd=fitdist(ages','normal');
x=linspace(min(ages),max(ages),100);
PDF=pdf(pd,x);
PDF=PDF/max(PDF)*max(n);
plot(x,PDF,'r-','LineWidth',2)
agespdf=[(pd.mu) (pd.sigma)]; % mu sigma
legend('models',['t= ' num2str(round(agespdf(1)/1000)/1000) '+/-' num2str(round(agespdf(2)/1000)/1000) ' Ma'])
t261021=[num2str(round(agespdf(1)/1000)/1000) ' +/- ' num2str(round(agespdf(2)/1000)/1000) ' Ma']

%% plot mote carlo results (ages including production rate errors)
figure(35)
hold on
% lages=log10(ages);
[n,xout] =hist(ages2,20);
% hist(ages,20)
bar(xout,n);
title(['Age (a) including ' num2str(RDSPrate*100) ' % of production rate errors'])
pd=fitdist(ages2,'normal');
x=linspace(min(ages2),max(ages2),100);
PDF=pdf(pd,x);
PDF=PDF/max(PDF)*max(n);
plot(x,PDF,'r-','LineWidth',2)
agespdf2=[(pd.mu) (pd.sigma)]; % mu sigma
legend('models',['t= ' num2str(round(agespdf2(1)/1000)/1000) '+/-' num2str(round(agespdf2(2)/1000)/1000) ' Ma'])
t261021includingallerrors=[num2str(round(agespdf2(1)/1000)/1000) ' +/- ' num2str(round(agespdf2(2)/1000)/1000) ' Ma']


% %% plot mote carlo results (R*s)
% figure(36)
% hold on
% % lages=log10(ages);
% [n,xout] =hist(Rs,20);
% % hist(ages,20)
% bar(xout,n);
% title('R*')
% pd=fitdist(Rs','normal');
% x=linspace(min(Rs),max(Rs),100);
% PDF=pdf(pd,x);
% PDF=PDF/max(PDF)*max(n);
% plot(x,PDF,'r-','LineWidth',2)
% Rspdf=[(pd.mu) (pd.sigma)]; % mu sigma
% legend('models',['R*= ' num2str((Rspdf(1))) '+/-' num2str((Rspdf(2))) ' '])

%% graphical model plot
figure(5)
hold on
     

% data from models
xdat=(Ps(:,1)./Ps(:,3))/(P10/P21); % Slope(10Be/21Ne)/Rinit
ydat=(Ps(:,2)./Ps(:,3))/(P26/P21); % Slope(26Al/21Ne)/Rinit
pointsperaxis=30; % resolution for the density plot
xplot=linspace(quantile(xdat,0.01),quantile(xdat,0.99),pointsperaxis); % axis for the density plot
yplot=linspace(quantile(ydat,0.01),quantile(ydat,0.99),pointsperaxis); % axis for the density plot
[XPLOT,YPLOT] = meshgrid(xplot,yplot);
ZPLOT=0.*XPLOT; % density grid
for i=1:length(xplot)-1
    for j=1:length(yplot)-1
        ZPLOT(j,i)=sum(xdat>xplot(i) & xdat<xplot(i+1) & ydat>yplot(j) & ydat<yplot(j+1)); % calculate densities (x and y defines the corner of the pixel)
    end
end
limits(1)=max((cumsum(sort(ZPLOT(:)))<sum(ZPLOT(:))*(1-0.6827)).*sort(ZPLOT(:))); % calculate one sigma limits
limits(2)=max((cumsum(sort(ZPLOT(:)))<sum(ZPLOT(:))*(1-0.9545)).*sort(ZPLOT(:))); % calculate one sigma limits
ZPLOTs=ZPLOT/max(ZPLOT(:)); % scale to maximum
Zcolors=(1-ZPLOTs*0.7)*256; % 256 greyscale colors
image(xplot,yplot,Zcolors) % plot density
colormap(hot(256))
% contour(XPLOT,YPLOT,ZPLOT,limits,'Color','r') % contour one and two sigma areas
% plot(mean(xdat),mean(ydat),'.k') % plot mean

% constant tb models
for tb=[0.5 1 2 3 4 5 6 7 8 10]*1e6
    expratio=logspace(-5,0,300);
    RBe=(exp(-tb.*l10).*P10.*(1-expratio)+expratio.*P10)./(P21);
    RAl=(exp(-tb.*l26).*P26.*(1-expratio)+expratio.*P26)./(P21);
    RBe=RBe/(P10/P21);
    RAl=RAl/(P26/P21);
    text(RBe(1),RAl(1),[num2str(round(tb/1e6*10)/10) 'Ma'],'HorizontalAlignment','center','VerticalAlignment','bottom','Color','k','Clipping','on')
    plot(RBe,RAl,'-k');
end

% constant texp models
for expratio=[0.2 0.4 0.6 0.8]
    tb=linspace(1,10e6,300);
    RBe=(exp(-tb.*l10).*P10.*(1-expratio)+expratio.*P10)./(P21);
    RAl=(exp(-tb.*l26).*P26.*(1-expratio)+expratio.*P26)./(P21);
    RBe=RBe/(P10/P21);
    RAl=RAl/(P26/P21);
    text(RBe(length(tb)),RAl(length(tb)),[' ' num2str(round(expratio*100)) '%'],'HorizontalAlignment','right','VerticalAlignment','bottom','Color','k','Clipping','on')
    plot(RBe,RAl,'--k');
end

% 0 texp model
expratio=0;
tb=logspace(0,9,300);
RBe=(exp(-tb.*l10).*P10.*(1-expratio)+expratio.*P10)./(P21);
RAl=(exp(-tb.*l26).*P26.*(1-expratio)+expratio.*P26)./(P21);
RBe=RBe/(P10/P21);
RAl=RAl/(P26/P21);
plot(RBe,RAl,'-','Color','k','LineWidth',2);
midline=round(length(tb)*0.65);
text(RBe(midline),RAl(midline),'   Simple exposure+burial','Color','k');

% other basic lines
plot([0 1],[0 1],'-k','LineWidth',2) % infinite burial
plot(1,1,'p','LineWidth',2,'MarkerSize',10, 'MarkerEdgeColor','k','MarkerFaceColor','w') % surface
text(1,1,[' Surface (Prod. rate ratios)   '],'HorizontalAlignment','right','VerticalAlignment','top','Color','k','Clipping','on')      

xlabel('R_{10/21}')
ylabel('R_{26/21}')
title('Axis normalized to production rate ratios')

box('on')
set(gcf,'color','w');
box('on')
set(gcf,'color','w');

%% Check negative slopes (terrible data)
if min(mean(Ps(:,2)./Ps(:,3)./(P262./P212)),mean(Ps(:,1)./Ps(:,3)./(P102./P212)))<0
    errordlg('Negative slopes! Check input data...','Error');
end