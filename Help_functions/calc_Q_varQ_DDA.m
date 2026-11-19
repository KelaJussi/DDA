function [ff,success] = calc_Q_varQ_DDA(lambda, T, location,apprQ, plotQ,titlestr)

    %lambda, location should be in same units
    %lambda = vector of x values
    % locations: approximately the location of the spike
    % T: transmission data (y-values)
    % plotQ: false = no plot; true = plot in the figure tagged 'DDA_Q_fit'
    %        (one tile per fit, created when missing); or a tiledlayout or
    %        axes handle to plot into
    % intensitycutoff, peaks with power extinction lower than this will be
    % ignored
    % errorcutoff: the bigger this is, the worse fits get accepted
    if nargin < 6
        titlestr = "";
    end
    if nargin < 5
        plotQ = 1;
    end
    if nargin < 4
        apprQ = 1e3;
    end
    success = true;
    [~, minindex] = min(abs(lambda.' - location));
    slr_peak = lambda(minindex);             %(m)
    %slr_peak= 982.6e-9;
    apprWidth = slr_peak/apprQ/2;


    %% Calculate Q-factor etc. with Lorentz fit
    %reftrans = 1;               %reference transmittance for Q calculations (out of resonance transmittance)
    xmin = slr_peak-apprWidth*12; xmax= slr_peak+apprWidth*12;
    minindex = find(lambda>xmin,1);
    maxindex = find(lambda>xmax,1);
    if isempty(maxindex)
        maxindex = length(lambda);
    end
    if isempty(minindex)
        minindex = 1;
    end
    xdata = reshape(lambda(minindex:maxindex),[],1,1);
    ydata = reshape(T(minindex:maxindex),[],1,1);
    reftrans = min(2,max(ydata));
   % reftrans = 0;
    
    % kokonaan nouseva tai laskeva funktion ei voi sisältää piikkiä, ja fit
    % tarvitsee 5 datapistettä
    if any([isempty(ydata), min(diff(ydata))> 0 , max(diff(ydata))<0 , length(xdata)<5])
        success = false;
        ff = struct();

        return
    end

    if min(ydata)>0.998
        %% Small particles
        f = fano_fit(xdata,(reftrans-ydata)*1e3,[0, 0, 0.1, slr_peak, apprWidth]);
        ff = f;
            %% plot it
        if wants_plot(plotQ)
            ax = q_axes(plotQ);
            plot(ax, xdata,1e3*(-reftrans+ydata)+reftrans)
            hold(ax, 'on')
            plot(ax, xdata, reftrans-ff(xdata));
            %plot(xdata, ff(xdata));
            Q = abs(ff.x0/(ff.gamma*2));
           
            
            titleshow = ['Q: ', num2str(Q, '%.2f'),', a1:',num2str(ff.a1),' a2:',num2str(ff.a2),convertStringsToChars(titlestr)];
            title(ax, titleshow);
            legend(ax, {"extinction","fit"})
            
        end
        ff.I = 1e-3*(max(ff(xdata))-min(ff(xdata)));
    else
        %% Normal particles
          %f = lorentz_fit(xdata,reftrans-ydata,[1 1 slr_peak+ii]);
   % f = fano_fit(xdata,reftrans-ydata,[0, 0, 1, slr_peak, apprWidth]);
        ymax = max(abs(ydata));
       % ymax = 1;
     f = fano_fit2(xdata,reftrans-ydata); %./ymax
   % f  = fano_fit2(xdata,ydata, [1, slr_peak, apprWidth, 1e-5]);
    ff = f;
    %% plot it
        if wants_plot(plotQ)
            ax = q_axes(plotQ);
            plot(ax, xdata,ydata)
            hold(ax, 'on')
            plot(ax, xdata, reftrans-ff(xdata));
            %plot(xdata, ff(xdata));
            Q = abs(ff.x0/(ff.gamma));
             [~,~,promwidth,valleyProms] = findpeaks(ydata,'MinPeakProminence', 0.01);
             [~,index]=max(valleyProms);
             % width from findpeaks is in samples: peak position / (width * sample spacing)
             Qprom = slr_peak./(promwidth(index)*abs(xdata(2)-xdata(1)));
            titleshow = ['Q: ', num2str(Q, '%.2f'),', Qprom:',num2str(Qprom, '%.2f'),convertStringsToChars(titlestr)];%', sigma: ',num2str(ff.sigma),
            title(ax, titleshow);
            %legend({"extinction","fit"})
    
            
        end
         ff.I = max(ff(xdata))-min(ff(xdata));
    end

    if isempty(ff)
        success = false;
        return
    end
   warning('off','last');
function f = lorentz_fit(data_x,data_y,startpoint)

% startpoint is a vector with values [I gamma x0]
lorentz = @(I,gamma,x0,x) I.*(gamma.^2./((x-x0).^2+gamma.^2));
data_y = double(data_y);
f = fit(data_x,data_y,lorentz,'StartPoint', startpoint );

end

function f = fano_fit(data_x,data_y,startpoint)
%from: "High Q-factor multiple Fano resonances for high-sensitivity sensing
%in all-dielectric metamaterials", by Wang et al. (2019)
    fano = @(a1,a2,I,x0,gamma,x) abs(a1+a2*1i+I./(x-x0-1i*gamma)).^2;
    f = fit(data_x,double(data_y), fano, 'Startpoint', startpoint);
end

function f = fano_fit2(data_x,data_y)
%from: "High Q-factor multiple Fano resonances for high-sensitivity sensing
%in all-dielectric metamaterials", by Wang et al. (2019)
fano = @(sigma,x0, gamma, I, x) I.*4.*sin(sigma).^2.*(cot(sigma)+2*(x-x0)./gamma).^2./(1+4*(x-x0).^2./gamma.^2);
    options = fitoptions(fano);
   options.StartPoint = [pi/2, slr_peak,apprWidth, 1];
   options.TolFun = 1e-10;
   options.TolX = 1e-10;
   options.MaxIter = 2000;
   options.MaxFunEvals = 2e5;
   options.Display = 'off';
   wghts = ones(size(data_y));
   [~,mini] = min(data_y);
   [~,maxi] = max(data_y);
   %wghts([mini]) = 5;
   options.Weights = wghts;
   options.Lower = [-pi,slr_peak-apprWidth,apprWidth*0.001,0];
   options.Upper = [pi,slr_peak+apprWidth,apprWidth*10,inf];
    %options.Normalize = 'on';
    f = fit(data_x,double(data_y), fano, options);
end
end



function tf = wants_plot(plotQ)
% true, 1, or a graphics handle
if isa(plotQ, 'matlab.graphics.Graphics')
    tf = true;
else
    tf = ~isempty(plotQ) && (isnumeric(plotQ) || islogical(plotQ)) && all(plotQ(:) ~= 0);
end
end

function ax = q_axes(plotQ)
% Axes for the next fit plot. A tiledlayout or axes handle is used as given.
% Otherwise one figure (Tag 'DDA_Q_fit') with a 'flow' layout, which never
% runs out of tiles. (It used to be figure(10), which could be any figure,
% e.g. a full 2 x 3 dipole overview: "The layout does not have sufficient space".)
if isa(plotQ, 'matlab.graphics.layout.TiledChartLayout')
    ax = nexttile(plotQ);  return
elseif isa(plotQ, 'matlab.graphics.axis.Axes')
    ax = plotQ;  return
end
fig = findobj(groot, 'Type', 'figure', 'Tag', 'DDA_Q_fit');
if isempty(fig)
    if exist('dda_figure', 'file') == 2
        fig = dda_figure([100 100 700 500], 'Name', 'Q fit');
    else
        fig = figure('Name', 'Q fit');
    end
    fig.Tag = 'DDA_Q_fit';
end
fig = fig(1);
tl = findobj(fig, 'Type', 'tiledlayout');
if isempty(tl), tl = tiledlayout(fig, 'flow'); end
ax = nexttile(tl(1));
end

%fano_func = @(delta,x0,Gamma,x) 4*sin(delta).^2*(cot(delta)+2*(x-x0)/Gamma).^2./(1+(2*(x-x0)/Gamma).^2)


