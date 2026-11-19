function F = nufft_Gauss(X,Y,mode, os, pad, sigma, radius)
%Caluclates nufft by using gaussian oversampling
% X,Y Nx1 points for values mode   
    if nargin < 4
        % Oversampling factor
         os = 4;
         pad = 256; 
        sigma = 2;
        radius = 4;
    end
    
    
    % Choose grid size
    Nx = os * ceil(sqrt(numel(X)));
    Ny = os * ceil(sqrt(numel(Y)));
    
    % Coordinate limits
    xmin = min(X);
    xmax = max(X);
    
    ymin = min(Y);
    ymax = max(Y);
    
    % Uniform grid
    xg = linspace(xmin, xmax, Nx);
    yg = linspace(ymin, ymax, Ny);
    
    dx = xg(2)-xg(1);
    dy = yg(2)-yg(1);
    
    G = complex(zeros(Ny,Nx));
    
   
    for n = 1:numel(X)
    
        tx = (X(n)-xmin)/dx + 1;
        ty = (Y(n)-ymin)/dy + 1;
    
        ix0 = round(tx);
        iy0 = round(ty);
    
        for ix = max(1,ix0-radius):min(Nx,ix0+radius)
            for iy = max(1,iy0-radius):min(Ny,iy0+radius)
    
                w = exp(-((tx-ix)^2 + (ty-iy)^2)/(2*sigma^2));
    
                G(iy,ix) = G(iy,ix) + mode(n)*w;
            end
        end
    end
     % FFT
    F = fft2(G,pad,pad);
    if 0
        
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


end