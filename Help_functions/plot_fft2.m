function plot_fft2(X,Y,F,os)
    if nargin < 4
        os = 1;
    end
        Nx = os * ceil(sqrt(numel(X)));
        Ny = os * ceil(sqrt(numel(Y)));
        
        xmin = min(X);
        xmax = max(X);
        
        ymin = min(Y);
        ymax = max(Y);
        Lx = xmax - xmin;
        Ly = ymax - ymin;
        
        kx = 2*pi * (-Nx/2:Nx/2-1)/Lx;
        ky = 2*pi * (-Ny/2:Ny/2-1)/Ly;
        
        [KX,KY] = meshgrid(kx,ky);
        current_hFig = gcf;   % get current figure handle
        figNum = 24;

        hFig = findobj('Type','figure','Number',figNum);
        
        if isempty(hFig)
            hFig = figure(figNum);
            tiledlayout;
        else
            figure(hFig);   % make existing figure current
        end
        nexttile
        imagesc(kx,ky,abs(fftshift(F)).^2)
        axis xy
        xlabel('k_x')
        ylabel('k_y')
        colorbar
        figure(current_hFig); % make it current and bring it to front
end