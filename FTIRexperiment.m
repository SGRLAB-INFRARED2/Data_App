classdef FTIRexperiment
    properties
        data double
        freqAxis (:,1) double
        volume (:,1) double = 1 % must be in microliters
        pathLength (:,1) double = 12 % must be in micrometers
        radius (:,1) double
        timeInterval (:,1) double = 0 % must be in seconds
        timePts (:,1) = []
        sample string = ""
        finalSpectrum = [] % to get the concentration at saturation
        finalConc double = [] % concentration at saturation
        dateString string = ""
        name string = ""
        fittedSpectra struct = []
        fitMethod
        diffusionFitResult struct
    end
    methods
        % CONSTRUCTOR METHOD !!!
        function f = FTIRexperiment(data,freqAxis,volume,...
                pathLength,radius,timeInterval,sample,dateString,name)
            if nargin > 0
                f.data = data;
                f.freqAxis = freqAxis;
                f.volume = volume;
                f.pathLength = pathLength;
                f.radius = radius;
                f.sample = sample;
                f.timeInterval = timeInterval;
                f.dateString = dateString;
                f.name = name;
            end
        end
        function obj = timeAxis(obj,varargin)
            if nargin == 0
                filenames = uigetfile({'*.spa','Thermo Spectrum (*.spa)'}, ...
                    'MultiSelect','on','Select Spectra Files...');
            elseif nargin == 4
                pathname = varargin{1};
                fileroot = varargin{2};
                nums = varargin{3};
                
                cd(pathname);
                
                files = dir(fileroot + "*");
                filenames = {files(nums).name};
            else
                error("Error: Invalid set of arguments. ...If using arguments, enter (filepath,fileroot,nums)")
            end
            times = [];
            for ii = 1:numel(filenames)
                g = dir(filenames{ii});
                times = [times datetime(g.date)];
            end
            sec = [];
            for ii = 1:numel(times)
                sec = [sec seconds(times(ii) - times(1))];
            end
            obj.timePts = sec;
            %             axis = (0:(size(obj.data,2)-1)).*obj.timeInterval;
        end
        % % % % % ADD A METHOD TO FIND/PLOT DIFFERENCE BETWEEN CO2 PEAK AND REAL
        % SPECTRA
        % methods that require a fitted spectra set
        function concs = concOverTime(obj)
            if isempty(obj.fittedSpectra)
                error('You do not have any fitted spectra. Fit the gas lines out first.')
                return
            end
            n_spectra = numel(obj.fittedSpectra);
            OD = zeros(1,n_spectra);
            for ii = 1:n_spectra
                temp = obj.fittedSpectra(ii).fobj;
                fcn = co2GasLineFitFunction(obj.fittedSpectra(ii).x,...
                    temp.center,temp.w_g,temp.w_l,temp.a1,temp.a2,0,0,0);
                OD(ii) = max(fcn);
            end
            concs = OD./(1000*obj.pathLength*1e-4);
        end
        function obj = getFinalConc(obj)
            if isstruct(obj.finalSpectrum)
                temp = obj.finalSpectrum.fobj;
                fcn = co2GasLineFitFunction(obj.finalSpectrum.x,...
                    temp.center,temp.w_g,temp.w_l,temp.a1,temp.a2,0,0,0);
                OD = max(fcn);
                obj.finalConc = OD./(1000*obj.pathLength*1e-4);
            else
                temp = obj.concOverTime;
                obj.finalConc = temp(end);
            end
        end
        function obj = gasLineFit(obj,center,wg,wl,a1,a2,a3,c0,c1)
            
            n_spectra = size(obj.data,2); % number of columns
            
            %initial guess from inputs do this before calling function
            sp = [center,wg,wl,a1,a2,a3,c0,c1];
            %upper and lower bounds
            lb = [2300, 0.5, 0.5,   0, 0.0,   0, -10, -1];
            ub = [2400, 4,   4,   100, 0.2, inf,  10,  1];
            
            opts = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',lb,'Upper',ub,'StartPoint',sp,...
                'Display','Iter');
            ft = fittype(@(center,w_g,w_l,a1,a2,a3,c0,c1,w) co2GasLineFitFunction(w,center,w_g,w_l,a1,a2,a3,c0,c1),...
                'independent',{'w'},'dependent','absorbance',...
                'coefficients',{'center','w_g','w_l','a1','a2','a3','c0','c1'},...
                'options',opts);
            
            %clear out
            out(n_spectra) = struct('x',[],'ydata',[],'yfit',[],'res',[],...
                'fobj',[],'G',[],'O',[]);
            
            % start a timer
            tic
            
            % set the fit range
            range1 = [2290 2390];
            
            % fit each spectrum
            for ii = 1:n_spectra
                
                freq = flip(obj.freqAxis);
                s = flip(obj.data(:,ii));
                
                % update the fitting region (x and y)
                ind1 = find(freq>=range1(1) & freq<range1(2));
                x = freq(ind1);
                ydata = s(ind1);
                
                % do the fit
                [fobj, G, O] = fit(x,ydata,ft);
                
                
                % get fit result for plotting
                yfit = fobj(x);
                
                % pack up the data and results
                out(ii).x = x;
                out(ii).ydata = ydata;
                out(ii).yfit = yfit;
                out(ii).res = ydata - yfit;
                out(ii).fobj = fobj;
                out(ii).G = G;
                out(ii).O = O;
                
            end
            
            % stop the timer
            toc
            
            % check results
            for ii = 1:n_spectra
                if out(ii).O.exitflag < 1
                    warning('Spectrum %i did not converge!!! Results might not be trustworthy.',ii);
                end
            end
            
            obj.fittedSpectra = out;
            
            if ~isempty(obj.finalSpectrum) & ~isstruct(obj.finalSpectrum)
                n_spectra = 1; % number of columns
                
                %initial guess from inputs do this before calling function
                sp = [center,wg,wl,a1,a2,a3,c0,c1];
                %upper and lower bounds
                lb = [2300, 0.5, 0.5,   0, 0.0,   0, -10, -1];
                ub = [2400, 4,   4,   100, 0.2, inf,  10,  1];
                
                opts = fitoptions('Method','NonlinearLeastSquares',...
                    'Lower',lb,'Upper',ub,'StartPoint',sp,...
                    'Display','Iter');
                ft = fittype(@(center,w_g,w_l,a1,a2,a3,c0,c1,w) co2GasLineFitFunction(w,center,w_g,w_l,a1,a2,a3,c0,c1),...
                    'independent',{'w'},'dependent','absorbance',...
                    'coefficients',{'center','w_g','w_l','a1','a2','a3','c0','c1'},...
                    'options',opts);
                
                %clear out
                uh(n_spectra) = struct('data',[],'x',[],'ydata',[],'yfit',[],'res',[],...
                    'fobj',[],'G',[],'O',[]);
                
                % start a timer
                tic
                
                % set the fit range
                range1 = [2290 2390];
                
                % fit each spectrum
                
                freq = flip(obj.freqAxis);
                s = flip(obj.finalSpectrum);
                
                % update the fitting region (x and y)
                ind1 = find(freq>=range1(1) & freq<range1(2));
                x = freq(ind1);
                ydata = s(ind1);
                
                % do the fit
                [fobj, G, O] = fit(x,ydata,ft);
                
                
                % get fit result for plotting
                yfit = fobj(x);
                
                % pack up the data and results
                uh.data = obj.finalSpectrum;
                uh.x = x;
                uh.ydata = ydata;
                uh.yfit = yfit;
                uh.res = ydata - yfit;
                uh.fobj = fobj;
                uh.G = G;
                uh.O = O;
                
                
                
                % stop the timer
                toc
                
                % check results
                
                if uh.O.exitflag < 1
                    warning('Spectrum %i did not converge!!! Results might not be trustworthy.',ii);
                end
                obj.finalSpectrum = uh;
            end
        end
    end
end