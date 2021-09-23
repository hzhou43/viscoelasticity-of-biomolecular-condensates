function [phasediff,output] = oscillationfit(Filename,initforce,initamp,frequency)
% Input the filename (e.g. 'forcefile.h5'), an initial value for the force
% signal, an initial value for the amplitude and the frequency of
% oscillation.

Force1x = h5read(Filename,'/Force HF/Force 1x');

trapposition = h5read(Filename,'/Trap position/1X');

%The time index is converted to time in seconds with the sampling rate =
%78125 Hz.
timeindex = [0:1:length(Force1x)-1];
timeindex = timeindex';
t = timeindex./78125;

% Both the x component of force signal from trap 1 and the trap 1 position are fit to a sine function.
options1 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-Inf,-Inf,-Inf,0],...
               'Upper',[Inf,Inf,Inf,2*pi],...
               'Startpoint',[initforce initamp frequency pi],'algorithm','Trust-Region',...
               'DiffMinChange',(1.0E-8),'DiffMaxChange',(0.1),'MaxFunEvals',(600),'MaxIter',(400),...
               'TolFun',(1.0E-6),'TolX',(1.0E-6));

options2 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,-Inf,0],...
               'Upper',[Inf,Inf,Inf,2*pi],...
               'Startpoint',[4.2 0.125 freq pi],'algorithm','Trust-Region',...
               'DiffMinChange',(1.0E-8),'DiffMaxChange',(0.1),'MaxFunEvals',(600),'MaxIter',(400),...
               'TolFun',(1.0E-6),'TolX',(1.0E-6));
           
f1 = fittype('a0 + abs(a).*sin(2*pi.*f1.*t + abs(phi1))',...
    'dependent',{'Force1x'},'independent',{'t'},...
    'coefficients',{'a0','a','f1','phi1'});

f2 = fittype('b0 + abs(b).*sin(2*pi.*f2.*t + abs(phi2))',...
    'dependent',{'trapposition'},'independent',{'t'},...
    'coefficients',{'b0','b','f2','phi2'});

[c1,gof1,output1] = fit(t,Force1x,f1,options1);
a0 = c1.a0;
a = c1.a; 
f1 = c1.f1;
phi1 = c1.phi1;
a = abs(a);
phi1 = abs(phi1);
Force1xfit = a0 + a.*sin(2*pi.*f1.*t + phi1);

[c2,gof2,output2] = fit(t,trapposition,f2,options2);
b0 = c2.b0;
b = c2.b; 
f2 = c2.f2;
phi2 = c2.phi2;
b = abs(b);
phi2 = abs(phi2);
trappositionfit = b0 + b.*sin(2*pi.*f2.*t + phi2);

phasediff = rad2deg(phi1-phi2);

% Both the force signal response and the trap oscillation are plotted on
% the same graph.
figure

yyaxis left
plot(t,Force1x,'k.','MarkerSize',1)
hold
plot(t,Force1xfit,'r','LineWidth',2)
xlabel('Time (s)')
ylabel('Force 1x (pN)')
grid on

trapposition = normalize(trapposition/0.249,'center','mean');
yyaxis right
plot(t,trapposition,'k.','MarkerSize',1)
hold
trappositionfit = normalize(trappositionfit/0.249,'center','mean');
plot(t,trappositionfit,'g','LineWidth',2)
xlabel('Time (s)')
ylabel('Trap position (\mum)')
grid on

% The output is further processed to obtain the complex moduli
output = [num2str(a0),' ',num2str(a),' ',num2str(f1),' ',num2str(phi1),' ',num2str(b0),' ',num2str(b),' ',num2str(f2),' ',num2str(phi2),' ',num2str(phasediff)];

end