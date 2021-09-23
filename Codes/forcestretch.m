function [output] = forcestretch(Filename,n,p,q)

% n is moving mean window
% p and q are to shift the data such that they center around 0 on the y
% axis

%parsing
F1x = h5read(Filename,'/Force HF/Force 1x');
F2x = h5read(Filename,'/Force HF/Force 2x');
trapposition = h5read(Filename,'/Trap position/1X');

timeindex = [0:1:length(F1x)-1];
timeindex = timeindex';
t = timeindex./78125;

%moving average
F1x = movmean(F1x,n);
F2x = movmean(F2x,n);
trapposition = trapposition-trapposition(1);
trapposition = movmean(trapposition,n);
trapposition = trapposition/0.249; %trap position in microns

%fitting
options1 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-Inf,-Inf],...
               'Upper',[Inf,Inf],...
               'Startpoint',[1 1],'algorithm','Trust-Region',...
               'DiffMinChange',(1.0E-8),'DiffMaxChange',(0.1),'MaxFunEvals',(600),'MaxIter',(400),...
               'TolFun',(1.0E-6),'TolX',(1.0E-6));

options2 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-Inf,-Inf],...
               'Upper',[Inf,Inf],...
               'Startpoint',[1 1],'algorithm','Trust-Region',...
               'DiffMinChange',(1.0E-8),'DiffMaxChange',(0.1),'MaxFunEvals',(600),'MaxIter',(400),...
               'TolFun',(1.0E-6),'TolX',(1.0E-6));
           
options3 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-Inf,-Inf],...
               'Upper',[Inf,Inf],...
               'Startpoint',[1 1],'algorithm','Trust-Region',...
               'DiffMinChange',(1.0E-8),'DiffMaxChange',(0.1),'MaxFunEvals',(600),'MaxIter',(400),...
               'TolFun',(1.0E-6),'TolX',(1.0E-6));
           
g1 = fittype('f1*t + b1',...
    'dependent',{'F1x'},'independent',{'t'},...
    'coefficients',{'f1','b1'});

g2 = fittype('f2*t + b2',...
    'dependent',{'F2x'},'independent',{'t'},...
    'coefficients',{'f2','b2'});

g3 = fittype('v*t + x0',...
    'dependent',{'trapposition'},'independent',{'t'},...
    'coefficients',{'v','x0'});

[c1,gof1,output1] = fit(t,F1x,g1,options1);
f1 = c1.f1;
b1 = c1.b1; 
F1xfit = f1*t + b1;

[c2,gof1,output1] = fit(t,F2x,g2,options2);
f2 = c2.f2;
b2 = c2.b2; 
F2xfit = f2*t + b2;

[c3,gof1,output1] = fit(t,trapposition,g3,options3);
v = c3.v;
x0 = c3.x0; 
trapfit = v*t + x0;

c_sys0 = (f2-f1)/(2*-v);

%plotting
figure
yyaxis left
plot(t,F1x+p,'k.','MarkerSize',1)
grid on
hold
plot(t,F1xfit+p,'r','LineWidth',2)
plot(t,F2x-q,'k.','MarkerSize',1)
plot(t,F2xfit-q,'c','LineWidth',2)
xlabel('Time (s)')
ylabel('Force 1x & 2x (pN)')

yyaxis right
plot(t,trap,'g.','MarkerSize',1)
plot(t,trapfit,'g','MarkerSize',2)
grid on

% The output gives the slope and intercept of the force profiles of the two
% beads, the chi0 and the chi_sys0 values. These are further processed to
% obtain the surface tension of the droplet.
output = [f1 b1 f2 b2 v x0 c_sys0];

end