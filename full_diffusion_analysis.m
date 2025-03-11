%% FULL DATA ANALYSIS of diffusion experiment - from start to finish - 

% ---- Variables that require user manual input will be surrounded ----
    % like this
% ----

%%
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
%% Step 1: process the microscope image

%% Load the image from Isilon

% ---- path to the image within CHEM-SGR (end with a slash) ----
image_path = "leica_stereoscope/Matt_Kallol/20250228/";
% ----

% ---- name of the image file (use the .tif one) ----
image_filename = 'PMNTF2EMIM_20250226.tif';
% ----

cd(isilon_path + image_path)
I = imread(image_filename);
I = I(:,:,1:3);
I = rgb2gray(I);
%% Process the image
figure(11);clf
imshow(I)
set(gcf,'Position',[813   332   868   615])
hold on

% ---- When selecting edges, the points you select must be fairly evenly
% spaced. I recommend selecting four points first, one in each cardinal
% direction, then selecting four more points, halfway between each pair of
% cardinal points. These 8 points typically suffice from my testing. ----

% select gel edge
fprintf("Now selecting gel edge...Press enter to continue.\n")
guide = annotation('textbox',[0.2 0 0.1 0.1],'String',"Select gel edge. Press enter to continue.");
guide.FontSize = 18;
guide.Color = 'white';
guide.EdgeColor = 'white';
gel_edge = ginput;
scatter(gel_edge(:,1),gel_edge(:,2),'filled','MarkerFaceColor','green')
gel_center = [mean(gel_edge(:,1)) mean(gel_edge(:,2))];
scatter(gel_center(1),gel_center(2),100,'filled','MarkerFaceColor','red')
delete(findall(figure(11),'type','annotation'))
% select pinhole edge
fprintf("Now selecting pinhole edge...Press enter to continue.\n")
guide = annotation('textbox',[0.2 0 0.1 0.1],'String',"Select pinhole edge. Press enter to continue.");
guide.FontSize = 18;
guide.Color = 'white';
guide.EdgeColor = 'white';
pinhole_edge = ginput;
scatter(pinhole_edge(:,1),pinhole_edge(:,2),'filled','MarkerFaceColor','cyan')
pinhole_center = [mean(pinhole_edge(:,1)) mean(pinhole_edge(:,2))];
scatter(pinhole_center(1),pinhole_center(2),100,'filled','MarkerFaceColor','red')
delete(findall(figure(11),'type','annotation'))

% Calculate values
delete(findall(figure(11),'type','annotation'))
% calculate radius
radii = zeros(size(gel_edge,1),1);
for ii = 1:size(gel_edge,1)
    delta_x = gel_edge(ii,1) - gel_center(1);
    delta_y = gel_edge(ii,2) - gel_center(2);
    radii(ii) = sqrt(delta_x^2 + delta_y^2);
end
radius_pixels = mean(radii);
figure(11)
plot([gel_center(1) gel_center(1)-radius_pixels],[gel_center(2) gel_center(2)],'LineWidth',2,'Color','green')
a = annotation('textbox',[0.1 0.8 0.1 0.1],'String',"r = " + radius_pixels + " pixels");
a.FontSize = 16;
a.Color = 'green';
a.EdgeColor = 'white';

% find pinhole displcement
delta_x = pinhole_center(1) - gel_center(1);
delta_y = pinhole_center(2) - gel_center(2);
displacement_pixels = sqrt(delta_x^2 + delta_y^2);
figure(11)
plot([gel_center(1) pinhole_center(1)],[gel_center(2) pinhole_center(2)],'LineWidth',2','Color','cyan')
b = annotation('textbox',[0.1 0.75 0.1 0.1],'String',"dx = " + displacement_pixels + " pixels");
b.FontSize = 16;
b.Color = 'cyan';
b.EdgeColor = 'white';
%% Get the scale bar

% ---- put in the length of scale bar (in mm) manually ----
%              --------
scale_bar_mm = 0.18678;
%              --------

% ---- change these pixel bounds so that ONLY the scale bar is included.
% not even the number. it will be a VERY small figure, but necessary. ----
vertical_region = 1479:1482;
horizontal_region = 1940:2000;
% ----

figure(2);clf
sub_image = I(vertical_region,horizontal_region);
imshow(sub_image)
set(gcf,'Position',[1   490   681   457])
threshold = 180;
scale_bar = sub_image < threshold;
scale_bar = mean(scale_bar,1);
scale_bar_pixels = numel(scale_bar(scale_bar > 0));
mm_per_pixel = scale_bar_mm/scale_bar_pixels;
%% Calculate the radius
figure(11)
radius_mm = radius_pixels * mm_per_pixel;
a.String = a.String + ", " + radius_mm + " mm";
displacement_mm = displacement_pixels * mm_per_pixel;
b.String = b.String + ", " + displacement_mm + " mm";
%% Save the data
% make a neat structure
image_processing.path = image_path + image_filename;
image_processing.gel_center = gel_center;
image_processing.gel_edge = gel_edge;
image_processing.pinhole_center = pinhole_center;
image_processing.pinhole_edge = pinhole_edge;
image_processing.radius_mm = radius_mm;
image_processing.displacement_mm = displacement_mm;
image_processing.mm_per_pixel = mm_per_pixel;
image_processing.scale_bar_mm = scale_bar_mm;
image_processing.scale_bar_vertical_region = vertical_region;
image_processing.scale_bar_horizontal_region = horizontal_region;

cd(data_path)
save(date_of_experiment + "_image_processing_data.mat",'image_processing')
%% ---- END OF STEP 1 ----

%% Step 2: Fit the FTIR peaks to obtain the uptake curve

%% Load in the spectra
cd ~
% --- the indicies of the spectra you wish to use ----
spectra_range = [1:469]; 
% ----

% --- the spectra file prefix ---
file_prefix = 'PMNTF2EMIM_20250226_75C_';
% ----

% --- experimental parameters ---
volume = 0;  % in microliters
spacer_size = 25;  % in microns
gel_radius = radius_mm*1000;  % as previously calculated, but can be overridden
time_delay = 30;  % between spectra, in seconds
sample_name = "PMNTF2 EMIM";
your_name = "Matt";
% ---

cd(data_path)
[data1,freq] = LoadSpectra(data_path,file_prefix,spectra_range);
freq = freq(:,1);

if freq(2) - freq(1) > 0
    freq = flip(freq);
end

% Subtract to the initial spectrum
sub_data = data1 - data1(:,1);

% INITIALIZE OBJECT
f = FTIRexperiment(sub_data,freq,volume,spacer_size,gel_radius,...
    time_delay,sample_name,date_of_experiment,your_name);
f = f.timeAxis(data_path,file_prefix,spectra_range);

fprintf("Successfully imported " + size(f.data,2) + " spectra.\n")
%% Guesses for FTIR peak fitting, by eye
% ---- Which spectrum will you match to? Usually the last one is good.
trial_spectrum = 469;
% ----

% set the fit range. Usually doesn't need to be changed
range1 = [2290 2390];

% ---- User-input starting point values ----
sp.center = 2340;
sp.wg = 1.7; 
sp.wl = 2;
sp.a1 = 3.6;  % main peak height
sp.a2 = 0.13; % expected Boltzmann factor for bend
sp.a3 = 0.35; % gas lines
sp.c0 = 0.03;
sp.c1 = 0; % baseline slope
% ----

%fit function requires fliipped inputs
freq = flip(f.freqAxis);
s = flip(f.data(:,trial_spectrum));


%get x and y data for the fit
ind1 = find(freq>=range1(1) & freq<range1(2));
x = freq(ind1);
ydata = s(ind1);

%plot the fitted function using user parameters
yfit = co2GasLineFitFunction(x,sp.center,sp.wg,sp.wl,sp.a1,sp.a2,sp.a3,sp.c0,sp.c1);
res = ydata-yfit;
sse = sum(res.^2);

figure(1);clf
plot(x,ydata,'o',x,yfit,x,res-0.1,'r-o')
%% Do the FTIR peak fit
T = tic; %time the fitting for later display
f = gasLineFit(f,sp.center,sp.wg,sp.wl,sp.a1,sp.a2,sp.a3,sp.c0,sp.c1);
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
%% Plotting the uptake curve for viewing
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
%% Save the data
FTIR_peak_fitting_params.sp = sp;
FTIR_peak_fitting_params.trial_spectrum = trial_spectrum;
FTIR_peak_fitting_params.fit_range = range1;
FTIR_peak_fitting_params.spectra_range = spectra_range;
FTIR_peak_fitting_params.file_prefix = file_prefix;
cd(data_path)
save(date_of_experiment + "_FTIR_peak_fitting.mat",'FTIR_peak_fitting_params')
%% ---- END OF STEP 2 ----

%% Step 3: Fit for single diffusion coefficient

%% Guesses for uptake curve fitting, by eye
t = f.timePts;
y = f.concOverTime;
A = f.radius;
nmax = 150;
rres = 50;
rlim = 350;
sigma = 704;
dy = 0;

% ---- User input starting values
dx = displacement_mm*1000;  % from the image analysis. can be overridden
%      D    C   t0
sp = [2 0.245 -0]; % put guess here
ub = [1e5 1e3 max(t)];
lb = [0 0 -max(t)];
% ----

figure(728);clf
plot(t,y)
hold on
plot(t,diffusion_moving_beam(t,sp(1),f.radius,sp(2),nmax,sigma,dx,dy,"rlim",rlim,"t0",sp(3)))
%% Do the uptake curve fitting

%set up options and type
opts = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',lb,'Upper',ub,'StartPoint',sp,...
    'Display','Iter');

ft = fittype(@(D,C,t0,t) diffusion_moving_beam(t,D,A,C,nmax,sigma,dx,dy,"rlim",rlim,"t0",t0),...
    'independent',{'t'},...
    'dependent','absorbance',...
    'coefficients',{'D','C','t0'},...
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
figure(144);clf

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
%% Save the data

% ---- PUT IN THE CORRECT NOTEBOOK PAGE TITLE ---
page_title = '2025-02-26 Diffusion of CO2 in PMNTF2 EMIM at 75 C';
% ----

single_diffusion_fitting_params.sp = sp;
single_diffusion_fitting_params.ub = ub;
single_diffusion_fitting_params.lb = lb;
single_diffusion_fitting_params.dx = dx;
single_diffusion_fitting_params.nmax = nmax;
single_diffusion_fitting_params.rlim = rlim;
single_diffusion_fitting_params.sigma = sigma;
cd(data_path)
save(date_of_experiment + "_single_diffusion_fitting_params.mat",'single_diffusion_fitting_params')

% Update lab notebook with results
f.fitMethod = 'diffusion_moving_beam.m';
cd(data_path);
save(f.dateString + "_single_diffusion_fitting","f")

obj = labarchivesCallObj('notebook','Matt Lab Notebook',...
    'folder','Experiments',...
    'page',page_title);
% microscope photo
figure(11)
obj = obj.updateFigureAttachment('caption','Microscope photo of the sample annotated with calculated values');
% uptake curve
figure(3)
obj = obj.updateFigureAttachment;
% single diffusion fitting result
figure(144)
caption = "Single diffusion coefficient fitting: ";
coeffs = coeffnames(f.diffusionFitResult.fobj);
units = ["um^2/s" "M" "s"];
if numel(units) ~= numel(coeffs)
    error("Cannot match all fitting parameters with a unit.")
end
for ii = 1:numel(coeffs)
   std_devs{ii} = (ci(2,ii) - ci(1,ii))/4;
   caption = caption + coeffs{ii} + " = " + f.diffusionFitResult.fobj.(coeffs{ii})...
       + " ± " + std_devs{ii} + " " + units(ii) + ", ";
end
obj = obj.updateFigureAttachment('caption',caption);
%% ---- END OF STEP 3 ----

%% Step 4: Fit for double diffusion coefficient

% Proceed to the double_diffusion_fitting.m script for this part. You will
% need to have run this entire script for it to work as it will draw from
% the .mat files you saved here.