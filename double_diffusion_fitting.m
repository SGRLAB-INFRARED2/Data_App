%% DOUBLE DIFFUSION COEFFICIENT FITTING
% To run this script, must have fully analyzed single diffusion coefficient
% analysis already done by full_diffusion_analysis.m.

% ---- Variables that require user manual input will be surrounded ----
    % like this
% ----

%% Load in the mat file with the experiment data

% ---- IMPORTANT!!! Input the correct data right from the start. This will put
% all the data in the right place. If the date is wrong, it could overwrite
% previous analysis.
date_of_experiment = "2025-02-26";
% ----
year_of_experiment = year(datetime(date_of_experiment));

% ---- here put your path to CHEM-SGR (end with a slash). This will be necessary
% throughout the script.
isilon_path = "/Volumes/CHEM-SGR/";
% ----
data_path = isilon_path + "sgr-ftir.chem.pitt.edu/" + year_of_experiment + "/" + date_of_experiment;

cd(data_path)
load(date_of_experiment + "_single_diffusion_fitting.mat");
load(date_of_experiment + "_single_diffusion_fitting_params.mat")
%% Guesses for uptake curve fitting, by eye
%get parameters ready
t = f.timePts;
y = f.concOverTime;
A = f.radius;
nmax = double_diffusion_fitting_params.nmax;
rlim = double_diffusion_fitting_params.rlim;
sigma = double_diffusion_fitting_params.sigma;
dy = 0;
dx = double_diffusion_fitting_params.dx;

% ---- User input starting values
%     D1  D2  C1    C2     t0
sp = [2   2   0.1225 0.1225 -0]; % put guess here
ub = [1e5 1e5 1e3 1e3 max(t)];
lb = [0 0 0 0 -max(t)];
% ----

T = tic;

figure(728);clf
plot(t,y)
hold on
plot(t,...
    two_step_diffusion(t,sp(1),sp(2),A,sp(3),sp(4),nmax,sigma,dx,dx,dy,rlim,'t0',sp(5)))
plot(t,diffusion_moving_beam(t,sp(1),A,sp(3),nmax,sigma,dx,dy,"rlim",rlim,'t0',sp(5)),'red')
plot(t,diffusion_moving_beam(t,sp(2),A,sp(4),nmax,sigma,dx,dy,"rlim",rlim,'t0',sp(5)),'green')
legend("data","composite fit curve","fit curve 1","fit curve 2","Location",'northwest')

delta_t = toc(T);
fprintf("Time to plot was " + delta_t + " seconds.\n")

%% Actually do the fit

%set up options and type
opts = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',lb,'Upper',ub,'StartPoint',sp,...
    'Display','Iter');

ft = fittype(@(D1,D2,C1,C2,t0,t) two_step_diffusion(t,D1,D2,A,C1,C2,nmax,sigma,dx,dx,dy,rlim,'t0',t0),...
    'independent',{'t'},...
    'dependent','absorbance',...
    'coefficients',{'D1','D2','C1','C2','t0'},...
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
figure(14);clf

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

% readout = [string(uh.diffusionFitResult.fobj.D)]
% others = ["95% Confidence Interval is "+ci(1)+" to "+ci(2)+".";...
%     "R^2 = "+string(uh.diffusionFitResult.G.rsquare)]

f.diffusionFitResult.fobj
%% Save the data

% ---- PUT IN THE CORRECT NOTEBOOK PAGE TITLE ---
lab_notebook = 'Matt Lab Notebook';
folder = 'Experiments';
page_title = '2025-02-26 Diffusion of CO2 in PMNTF2 EMIM at 75 C';
% ----

double_diffusion_fitting_params.sp = sp;
double_diffusion_fitting_params.ub = ub;
double_diffusion_fitting_params.lb = lb;
double_diffusion_fitting_params.dx = dx;
double_diffusion_fitting_params.nmax = nmax;
double_diffusion_fitting_params.rlim = rlim;
double_diffusion_fitting_params.sigma = sigma;

f.fitMethod = 'two_step_diffusion.m';

cd(data_path);
save(date + "_double_diffusion_fitting.mat","f")
save(date + "_double_diffusion_fitting_params.mat",'double_diffusion_fitting_params')

obj = labarchivesCallObj('notebook',lab_notebook,...
    'folder',folder,...
    'page',page_title);
figure(14)
caption = "";
coeffs = coeffnames(f.diffusionFitResult.fobj);
units = ["um^2/s" "um^2/s" "M" "M" "s"];
if numel(units) ~= numel(coeffs)
    error("Cannot match all fitting parameters with a unit.")
end
for ii = 1:numel(coeffs)
   std_devs{ii} = (ci(2,ii) - ci(1,ii))/4;
   caption = caption + coeffs{ii} + " = " + f.diffusionFitResult.fobj.(coeffs{ii})...
       + " ± " + std_devs{ii} + " " + units(ii) + ", ";
end
obj = obj.updateFigureAttachment('caption',caption);