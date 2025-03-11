%% This script is for MANUAL FITTING and DATA ANALYSIS
% of the FTIR data for use OUTSIDE of the app.
%%
cd ~  % This is here because sometimes MATLAB gets confused 
% finding the Isilon folder so you have to reset the current folder to
% somewhere on the disk first.
spectra_range = [1:131]; 
cd('/Volumes/CHEM-SGR/sgr-ftir.chem.pitt.edu/2025/2025-02-26')
[data1,freq] = LoadSpectra('/Volumes/CHEM-SGR/sgr-ftir.chem.pitt.edu/2025/2025-02-26',...
    'PMNTF2EMIM_20250226_75C_',spectra_range);
freq = freq(:,1);

if freq(2) - freq(1) > 0
    freq = flip(freq);
end
% [data1,freq] = LoadSpectra();

% Subtract to the initial spectrum
sub_data = data1 - data1(:,1);

% INITIALIZE OBJECT
f = FTIRexperiment(sub_data,freq,0,25,1500,30,"PMIM NTF2","2025-02-26","Matt");
f = f.timeAxis('/Volumes/CHEM-SGR/sgr-ftir.chem.pitt.edu/2025/2025-02-26',...
    'PMNTF2EMIM_20250226_75C_',spectra_range);

fprintf("Successfully imported " + size(f.data,2) + " spectra.\n")

clear spectra_range
%% make initial guesses
% have the user select which spectrum to guess from
ii = 115;

% set the fit range
range1 = [2290 2390];

% set starting point using values from the user
center = 2340;
wg = 1.7; 
wl = 1.7;
a1 = 1.75;  % main peak height
a2 = 0.07; % expected Boltzmann factor for bend
a3 = 0.0; % gas lines
c0 = 0.005;
c1 = 0; % baseline slope

%fit function requires fliipped inputs
freq = flip(f.freqAxis);
s = flip(f.data(:,ii));


%get x and y data for the fit
ind1 = find(freq>=range1(1) & freq<range1(2));
x = freq(ind1);
ydata = s(ind1);

%plot the fitted function using user parameters
yfit = co2GasLineFitFunction(x,center,wg,wl,a1,a2,a3,c0,c1);
res = ydata-yfit;
sse = sum(res.^2);

figure(1);clf
plot(x,ydata,'o',x,yfit,x,res-0.1,'r-o')
%app.UIAxes3.Title = (sprintf('Initial guess SSE = %f',sse));
%% do the gas line fit
T = tic; %time the fitting for later display
f = gasLineFit(f,center,wg,wl,...
    a1,a2,a3,c0,...
    c1);
stop = toc(T);

%selecte 4 evenly placed fits to plot
n_spectra = size(f.data,2);
iis = ceil([1 n_spectra/4 n_spectra*3/4 n_spectra]);
figure(2);clf
for ii = iis
    plot(f.fittedSpectra(ii).x,f.fittedSpectra(ii).ydata,'o',...
        f.fittedSpectra(ii).x,f.fittedSpectra(ii).yfit,...
        f.fittedSpectra(ii).x,f.fittedSpectra(ii).res-0.1,'ro')
    hold on
end
hold off

%let the user know how it went
review = "";
tl = 0;
for ii = 1:n_spectra
    if f.fittedSpectra(ii).O.exitflag < 1
        review = [review;'Spectrum '+ii+' did not converge!!! Results might not be trustworthy.'];
        tl = tl+1;
    end
end
if tl==0
    review = "All fits were successful.";
end
review = [review;"Fitting took "+stop+" seconds."];
review
%% plotting the fits
figure(3);clf

% number of spectra to show
n = size(f.data,2);

%find the indicies for the amount of spectra desired
spectraIndicies = zeros(1,n);
interval = ceil(size(f.data,2)/n);
for ii = 1:n
    spectraIndicies(ii) = (ii*interval);
end

for ii = spectraIndicies
    temp = f.fittedSpectra(ii).fobj;
    pf = co2GasLineFitFunction(f.fittedSpectra(ii).x,temp.center,temp.w_g,temp.w_l,...
        temp.a1,temp.a2,0,0,0);
    plot(subplot(2,1,1),f.fittedSpectra(ii).x,pf)
    hold on
end
title('Fitted Spectra')
xlabel('Wavenumbers (cm^{-1})')
ylabel('Absorbance (AU)')
box off
set(gca,'TickDir','out')
hold off

plot(subplot(2,1,2),f.timePts,concOverTime(f),'o-','color','blue');
hold on
title('Concentration Over Time')
xlabel('Time (s)')
ylabel('Concentration (M)')
box off
set(gca,'TickDir','out')
hold off

set(gcf,'Units','normalized')
set(gcf,'Color','w')
set(gcf,'Position',[0.5 0 0.35 1])
%% final conc if applicable

% f.finalSpectrum = LoadSpectra();
% f = f.getFinalConc;

%% fit for diffusion coefficient
%get parameters ready
t = f.timePts;
%         t = t(1:end-3);
%           t = t(1:end-15);
%         t = t-t(1);
y = f.concOverTime;
%         y = y(4:end);
%           y = y(1:end-15);
A = f.radius;
C = f.finalConc;
nmax = 150;
rres = 50;
rlim = 350;
sigma = 704;
dx = 0;
dy = 0;
sp = [45.7 0.1556 480 -2142]; % put guess here
ub = [1e5 1e3 0.5*f.radius 1e5];
lb = [0 0 0 -1e5];

figure(728);clf
plot(t,y)
hold on
plot(t,diffusion_moving_beam(t,sp(1),f.radius,sp(2),nmax,sigma,sp(3),dy,"rlim",rlim,"t0",sp(4)))


%% Actually do the fit

%set up options and type
opts = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',lb,'Upper',ub,'StartPoint',sp,...
    'Display','Iter');

ft = fittype(@(D,C,dx,t0,t) diffusion_moving_beam(t,D,A,C,nmax,sigma,dx,dy,"rlim",rlim,"t0",t0),...
    'independent',{'t'},...
    'dependent','absorbance',...
    'coefficients',{'D','C','dx','t0'},...
    'options',opts);

%set up structure for storing output
out = struct('x',[],'ydata',[],'yfit',[],'res',[],...
    'fobj',[],'G',[],'O',[]);

tic

%do the fit
[fobj,G,O] = fit(t,y',ft);

toc

%get results
yfit = fobj(t);
out.x = t;
out.ydata = y;
out.yfit = yfit;
out.res = y - yfit;
out.fobj = fobj;
out.G = G;
out.O = O;

if out.O.exitflag < 1
    warning('Curve fit did not converge!!! Results might not be trustworthy.');
end

f.diffusionFitResult = out;
%% display fit result
figure(4);clf

plot(f.diffusionFitResult.x,f.diffusionFitResult.ydata,...
    'o','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(f.diffusionFitResult.x,f.diffusionFitResult.yfit,...
    'red','LineWidth',1.5)
residuals = f.diffusionFitResult.yfit - f.diffusionFitResult.ydata(:);
plot(f.diffusionFitResult.x,(residuals*10 - 0.02),'o','MarkerEdgeColor','red')
legend('Data points','Fitted curve','Location','northwest')
xlabel('Time (s)')
ylabel('Concentration (M)')
hold off


% get confidence intervals
ci = confint(f.diffusionFitResult.fobj);

readout = [string(f.diffusionFitResult.fobj.D)]
others = ["95% Confidence Interval is "+ci(1)+" to "+ci(2)+".";...
    "R^2 = "+string(f.diffusionFitResult.G.rsquare)]

f.diffusionFitResult.fobj

%% Update lab notebook with results
f.fitMethod = 'diffusion_moving_beam.m';
cd("/Volumes/CHEM-SGR/sgr-ftir.chem.pitt.edu/2025/"+f.dateString);
save(f.dateString,"f")

obj = labarchivesCallObj('notebook','Matt Lab Notebook',...
    'folder','Experiments',...
    'page','2025-02-19 Diffusion of CO2 in PMIM NTF2 at 35 C');
figure(3)
obj = obj.updateFigureAttachment;
figure(4)
caption = "";
coeffs = coeffnames(f.diffusionFitResult.fobj);
units = ["um^2/s" "M" "um" "s"];
if numel(units) ~= numel(coeffs)
    error("Cannot match all fitting parameters with a unit.")
end
for ii = 1:numel(coeffs)
   std_devs{ii} = (ci(2,ii) - ci(1,ii))/4;
   caption = caption + coeffs{ii} + " = " + f.diffusionFitResult.fobj.(coeffs{ii})...
       + " ± " + std_devs{ii} + " " + units(ii) + ", ";
end
D_std = (ci(2,1) - ci(1,1))/4;
C_std = (ci(2,2) - ci(1,2))/4;
dx_std = (ci(2,3) - ci(1,3))/4;
obj = obj.updateFigureAttachment('caption',caption);

%% Automatically save this script to lab notebook too
