function complex_freqs_plots_from_file(filename, plotting)
% COMPLEX_FREQS_PLOTS_FROM_FILE  Plots of a saved complex-frequency run.
%
%   complex_freqs_plots_from_file(file)            options saved in Input.plotting.complex
%   complex_freqs_plots_from_file(file, options)   override some of them, e.g.
%       complex_freqs_plots_from_file(f, struct('quiver', false, 'full_range', false, ...
%                                               're_lims', [1.12 1.18], 'im_lims', [25 30]))
%
% file: a name in dda_results_folder("complex") or a full path. options may be
% the complex options themselves or a whole plotting struct with .complex.
% Defaults (for files saved before these options existed) are in
% complex_plot_defaults at the end of this file, and in Make_predefined_Inputs.
%
% Mode labels (options.labels):
%   "sine" (default)  the eigenvector at each minimum of |lambda| is projected on
%                     the same basis as the real-frequency mode tracking
%                     (dda_mode_basis: sin((m+1)*pi*u) along open directions,
%                     Fourier along periodic ones, options.basis); the label is
%                     the basis function with the largest share of the mode's
%                     power, and that share is its confidence. Labels below
%                     options.min_confidence (0.25) get a '?' and are drawn grey.
%                     Complex (leaky) modes spread their power over neighbouring
%                     sines, so 0.3 - 0.6 is a typical confidence.
%   "nodes"           the number of sign changes of the mode (made as real as
%                     possible by a global phase) along x and y, through the row
%                     and column with the most power; Fourier/periodic directions
%                     keep the basis label. Confidence: how standing-wave-like
%                     the mode is (1 = real up to a phase, 0 = travelling).
%                     Use this for structures whose modes are not sine-like,
%                     e.g. half-shifted chains (array_shifts), where "sine"
%                     labels symmetric modes as nodes + 2.
%   "k"               the plane wave that couples most strongly to the mode:
%                     the |k| maximising |v.' * exp(i k r)|^2 (both signs of k,
%                     the unconjugated product of the reciprocal problem), mapped
%                     to a mode number with the standing waves sin((m+1) pi u)
%                     on the same positions as calibration. Confidence:
%                     1 - (second coupling peak)/(first); when it is low the
%                     tile title also shows the runner-up, e.g. (2|5).
%                     Modes of chains with a pi phase step (HG01: half-shifted
%                     arrays) couple to two |k| almost equally, k_m +- ~1.2 pi/L.
%                     options.k_step = true undoes the steps first (the sign of
%                     the dipoles flips at every half-period shift of the
%                     lattice), i.e. it labels by the step-matched (HG-shaped)
%                     field instead of a plane wave; each mode then has a single
%                     peak at its own carrier k_m.
%   "fft"             the peak of the padded spatial FFT (as before; plot_fft,
%                     nufft apply only here)
% options.label_map = true adds a figure with the label of the smallest-|lambda|
% eigenvector at every point of the complex plane, i.e. which mode is the
% nearest one where.
        %close all
        fpath = char(filename);
        if ~isfile(fpath), fpath = fullfile(dda_results_folder("complex"), fpath); end
        S = load(fpath, "C", "Input", "Array");
        [~, name, ext] = fileparts(fpath);  filename = [name ext];

        C = S.C;
        Input = S.Input;
        
        if isfield(S, "Array")
            Array = S.Array;
        else
            Array = [];   % or some default value
        end

        o = complex_plot_defaults();
        if isfield(Input, 'plotting') && isfield(Input.plotting, 'complex')
            o = merge_fields(o, Input.plotting.complex);
        end
        if nargin > 1 && ~isempty(plotting)
            if isfield(plotting, 'complex'), plotting = plotting.complex; end
            o = merge_fields(o, plotting);
        end
        docked = true;
        if isfield(Input, 'plotting') && isfield(Input.plotting, 'docked'), docked = Input.plotting.docked; end
        if nargin > 1 && isstruct(plotting) && isfield(plotting, 'docked'), docked = plotting.docked; end
        dda_figure('docking', docked);     % tabs in one window (Input.plotting.docked)
        is_nufft       = o.nufft;
        quiverplot     = o.quiver;
        imagescisZ     = o.imagesc_is_z;
        separate_pols  = o.separate_pols;
        surfplot1D     = o.surf_1D;
        tilesperfigure = o.tiles_per_figure;
        maxfftyplt     = o.max_fft_order;
        plot_fft       = o.plot_fft;
        use_sine       = any(strcmpi(string(o.labels), ["sine","nodes","k"]));
        use_nodes      = strcmpi(string(o.labels), "nodes");
        use_k          = strcmpi(string(o.labels), "k");
        last_alt       = [NaN NaN];             % runner-up label ("k" labels)
        LB             = struct('ok', false);   % sine basis (sine_setup)
        last_conf      = NaN;                   % confidence of the last sine label

        fullplot = o.full_range;
        % Real axis in energy (eV, default) or vacuum wavelength (nm). All plot
        % coordinates (heatmap, markers, limits) use the chosen unit; to_eV
        % converts back where the physics needs the energy.
        re_nm = strcmpi(string(o.real_axis), "nm");
        if re_nm, to_plot = @(e) 1239.84187./e; to_eV = @(x) 1239.84187./x;
        else,     to_plot = @(e) e;             to_eV = @(x) x;
        end
        ylims = o.im_lims;
        xlims = o.re_lims;
        if ~isempty(Array)
             ri = Array.r_i;
              X = ri(:,1); Y = ri(:,2);

        end
       

        % Axes in the same order as the data. The grid may run in either
        % direction (e.g. Im from -34 to -43 meV); unique() would sort the axes
        % and mirror the heatmap against the mode markers.
        if isfield(Input.range, 'real_parts') && isfield(Input.range, 'imag_parts')
            nre = numel(Input.range.real_parts);  nim = numel(Input.range.imag_parts);
        else
            nre = numel(unique(real(Input.range.eV)));  nim = numel(unique(imag(Input.range.eV)));
        end
        eV_grid = reshape(Input.range.eV, nim, nre);        % rows: Im, columns: Re
        real_freqs = real(eV_grid(1,:));
        imag_freqs = imag(eV_grid(:,1)).'*-1e3;
        Input.range.eV = to_plot(real(Input.range.eV))-imag(Input.range.eV)*1i*1e3;   % plot coordinates
        min_values = 1./abs(C.min_eig_values);
        shaped_min_values = reshape(min_values, length(imag_freqs), length(real_freqs));
        fig1=figure;
        ax1=axes;
        if re_nm
            heatmap_nm(ax1, real_freqs, imag_freqs, log10(shaped_min_values));   % uneven in nm: true cell edges
        else
            imagesc(ax1, real_freqs,imag_freqs,log10(shaped_min_values))
        end
        %colormap(flipud(pink))
        colormap(flipud(bone))
        set(gca,'YDir','normal')
        ax1.XColor = [0,0,0];ax1.YColor = [0,0,0];
        ax1.FontName = 'Arial';
        if re_nm, xlabel("Wavelength (nm)", 'Color' ,[0,0,0])
        else,     xlabel("Re[\omega] (eV)", 'Color' ,[0,0,0])
        end
        ylabel("-Im[\omega] (meV)", 'Color' ,[0,0,0])
        title(filename)
        
        if isfield(C,"fundindex") && 0
            fund_values = C.values(C.fundindex,:);
            fundlog = log10(1./abs(fund_values));
            shaped_fundlog = reshape(fundlog,length(imag_freqs), length(real_freqs));
            figure;
            
            imagesc(real_freqs,imag_freqs,shaped_fundlog)
            set(gca,'YDir','reverse')
            xlabel("Real frequency (eV)")
            ylabel("Imaginary frequency (eV)")
            title("Fundamental (0,0) mode")
            colorbar
        end
        set(0, 'currentfigure', fig1);
        
 
        
        if fullplot
            xlims = [min(real_freqs), max(real_freqs)];
            ylims = [min(imag_freqs), max(imag_freqs)];
        end
        xlims = sort(to_plot(xlims));          % re_lims and the grid are in eV
        ylim(ylims);
        xlim(xlims)
        %clim([-30.5,-28.3])
        a = colorbar;
        %ylabel(a,'log_{10}(1/\Lambda_{\mathrm{min}})','FontSize',12,'Rotation',270);
        if 0
            tic
            [r,pol] = aaa(min_values,Input.range.eV,'mmax',100);
            toc
              hold on
            plot(pol,'.r','markersize',10), 
            hold off
        else
            TF = islocalmax2(shaped_min_values);
            complex_arr = reshape(Input.range.eV,length(imag_freqs),length(real_freqs));
            TFasZ  = complex_arr(TF);

            % Interp2 for cubic search between data points
            %xvals = 1000; yvals = 1000;
            %xqvec = linspace(min(real(Input.range.eV),[],'all'),max(real(Input.range.eV),[],'all'), xvals);
            %yqvec = linspace(min(imag(Input.range.eV),[],'all'),max(imag(Input.range.eV),[],'all'), yvals);
            %[Xq, Yq] = meshgrid(xqvec,yqvec);
            %Zq = interp2(real(complex_arr), imag(complex_arr), shaped_min_values, Xq, Yq, 'makima');
            %interpTF = islocalmax2(Zq);
            %interpTFasZ = Xq(interpTF)+1i*Yq(interpTF);

            hold on
          %  plot(ax1, TFasZ,'.r','MarkerSize',10)

            %colorbar
            if any(Input.dimensionality == ["2D","cyl_TE"])
                mult = 2;
            elseif Input.dimensionality == "3D"
                mult = 3;
            else
                mult = 1;
            end
            if Input.array_type == "double_rectangular"
                mult = mult*2;
            end
            all_modes = reshape(C.min_eig_vectors,Input.arrayN(1)*Input.arrayN(2)*mult,length(imag_freqs),length(real_freqs));
            %wls = 1240./eVs;
            max_modes = all_modes(:,TF);
            max_freqs = complex_arr(TF);
            plotnmb =0;
            LB = sine_setup();
            if o.label_map && LB.ok
                plot_label_map();
            end
            if o.modes
                plot_all_modes_in_lims(xlims,ylims)
            end

        end



    function s = freq_label(z)
        % z in plot coordinates: real part eV or nm, imaginary part -Im in meV
        if re_nm, s = sprintf('%.2f nm, -Im %.2f meV', real(z), imag(z));
        else,     s = sprintf('%.4f eV - %.2fi meV', real(z), imag(z));
        end
    end

    function plot_all_modes()
        figure 
        tiledlayout
       for ii = 1:size(max_freqs,1)
           nexttile
           plot_mode(ii)      
       end

    end
    function plot_all_modes_in_lims(reallims, imaglims)
            fig2 = figure ;
            t= tiledlayout;
            tilesperfigure = tilesperfigure;
            is2d = 0;


           for ii = 1:size(max_freqs,1)
               if real(max_freqs(ii)) > reallims(1) && real(max_freqs(ii)) < reallims(2) ...
                       && imag(max_freqs(ii)) > imaglims(1) && imag(max_freqs(ii)) < imaglims(2)

                    
                     ax = findall(t, 'Type', 'Axes');
                    numAxes = numel(ax);
                    if numAxes == tilesperfigure
                  % if mod(plotnmb, tilesperfigure) == 0
                        fig2 = dda_figure([100,100,385,250]) ;
                        t = tiledlayout;
                    end
                   set(0, 'currentfigure', fig2);  %# for figures
                   ntext0 = numel(findobj(ax1, 'Type', 'text'));
                   last_conf = NaN;
                if any(Input.dimensionality == ["cyl_TM","cyl_TE"]) || Input.arrayN(1) == 1 || Input.arrayN(2) == 1
                        fftval = plot_mode_cyl_TM(ii);
                        set(0, 'currentfigure', fig1);  %# for figures
                        scatter(real(max_freqs(ii)),imag(max_freqs(ii)),30,'k','filled');%max(80-fftx*30,5)

                        %scatter(real(max_freqs(ii)),imag(max_freqs(ii)),5,'red','filled');
                        
                        text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d)",round(fftval,0)))
                elseif Input.dimensionality == "2D"
                      [dim, fftx,ffty] = plot_mode2D(ii);
                      set(0, 'currentfigure', fig1);  %# for figures
                      if isnan(fftx)
                        scatter(real(max_freqs(ii)),imag(max_freqs(ii)),'white','filled');
                      else
                            
                            if fftx+ffty < 30
                                if isfield(C,"degeneracy")
                                    
                                        %colors = [
                                           % 0.10 0.60 0.95;   % bright blue
                                          %  1.00 0.55 0.00;   % bright orange
                                         %   0.00 0.75 0.45;   % bright green
                                        %    0.90 0.35 0.70;   % bright magenta
                                       % 0.95, 0.20, 0.20       %bright red
                                       %0.20, 0.85, 0.90        %bright   cyan
                                       % ];

                                    degeneracy1 = reshape(C.degeneracy, length(imag_freqs), length(real_freqs));
                                    degeneracy = degeneracy1(TF);
                                    if degeneracy(ii) == 2
                                        scatter(real(max_freqs(ii)),imag(max_freqs(ii)),30,'red','filled');%, 'MarkerEdgeColor', 'k');
                                        text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d,%d)",floor(fftx),floor(ffty)),'Color','red')
                                    elseif dim == 1
                                        scatter(real(max_freqs(ii)),imag(max_freqs(ii)),30,[0.00 0.75 0.45],'filled');% 'MarkerEdgeColor', 'k');
                                        text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d,%d)",floor(fftx),floor(ffty)),'Color',[0.00 0.75 0.45])

                                    elseif dim == 2
                                        scatter(real(max_freqs(ii)),imag(max_freqs(ii)),30,[0.10 0.60 0.95],'filled');% 'MarkerEdgeColor', 'k');
                                        text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d,%d)",floor(fftx),floor(ffty)),'Color',[0.10 0.60 0.95])
                                    else
                                        scatter(real(max_freqs(ii)),imag(max_freqs(ii)),30,[0,0,0],'filled');% 'MarkerEdgeColor', 'k');
                                        text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d,%d)",floor(fftx),floor(ffty)),'Color',[0,0,0])
                                    end
                                    %text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d,%d)",floor(fftx),floor(ffty)))
                                else
                                    scatter(real(max_freqs(ii)),imag(max_freqs(ii)),max(30-(fftx+ffty)*5,5),'k','filled');
                                    text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d,%d)",floor(fftx),floor(ffty)))
                                end
                               
                            end
                      end
                   elseif Input.dimensionality == "3D"
                        
                        [dim, fftx,ffty] = plot_mode3D(ii);
                        if max(fftx,ffty) < 40
                            dotsize = 30;
                        else
                            dotsize = 8;
                        end
                        set(0, 'currentfigure', fig1);  %# for figures
                       
                          if isnan(fftx)
                            scatter(real(max_freqs(ii)),imag(max_freqs(ii)),'white','filled');
                          else
                              if isfield(C,"degeneracy")
                                    
                                            %colors = [
                                               % 0.10 0.60 0.95;   % bright blue
                                              %  1.00 0.55 0.00;   % bright orange
                                             %   0.00 0.75 0.45;   % bright green
                                            %    0.90 0.35 0.70;   % bright magenta
                                           % 0.95, 0.20, 0.20       %bright red
                                           %0.20, 0.85, 0.90        %bright   cyan
                                           % ];
    
                                        degeneracy1 = reshape(C.degeneracy, length(imag_freqs), length(real_freqs));
                                        degeneracy = degeneracy1(TF);
                                        if degeneracy(ii) == 3
                                            scatter(real(max_freqs(ii)),imag(max_freqs(ii)),30,'red','filled');%, 'MarkerEdgeColor', 'k');
                                            text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d,%d)",floor(fftx),floor(ffty)),'Color','red')
                                        elseif degeneracy(ii) == 2
                                            scatter(real(max_freqs(ii)),imag(max_freqs(ii)),30,[0.20, 0.85, 0.90],'filled');%, 'MarkerEdgeColor', 'k');
                                            text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d,%d)",floor(fftx),floor(ffty)),'Color',[0.20, 0.85, 0.90])
                                        elseif dim == 1
                                            scatter(real(max_freqs(ii)),imag(max_freqs(ii)),30,[0.00 0.75 0.45],'filled');% 'MarkerEdgeColor', 'k');
                                            text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d,%d)",floor(fftx),floor(ffty)),'Color',[0.00 0.75 0.45])
    
                                        elseif dim == 2
                                            scatter(real(max_freqs(ii)),imag(max_freqs(ii)),30,[0.10 0.60 0.95],'filled');% 'MarkerEdgeColor', 'k');
                                            text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d,%d)",floor(fftx),floor(ffty)),'Color',[0.10 0.60 0.95])
                                        elseif dim == 3
                                            scatter(real(max_freqs(ii)),imag(max_freqs(ii)),dotsize,[ 0.90 0.35 0.70],'filled');
                                            text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d,%d)",floor(fftx),floor(ffty)),'Color',[ 0.90 0.35 0.70])
                                        elseif dim == 4             %aka x==y
                                            scatter(real(max_freqs(ii)),imag(max_freqs(ii)),dotsize,'k','filled');
                                            text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d,%d)",floor(fftx),floor(ffty)))
                                        end
                                    else
                                        scatter(real(max_freqs(ii)),imag(max_freqs(ii)),30,[0,0,0],'filled');% 'MarkerEdgeColor', 'k');
                                        text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d,%d)",floor(fftx),floor(ffty)),'Color',[0,0,0])
                                end
                              
                              
                
                          end

                   
                   else

                        [fftx,ffty] = plot_mode(ii)  ;

                        set(0, 'currentfigure', fig1);  %# for figures
                         if isfield(C,"degeneracy")
                                    
                            %colors = [
                               % 0.10 0.60 0.95;   % bright blue
                              %  1.00 0.55 0.00;   % bright orange
                             %   0.00 0.75 0.45;   % bright green
                            %    0.90 0.35 0.70;   % bright magenta
                           % 0.95, 0.20, 0.20       %bright red
                           %0.20, 0.85, 0.90        %bright   cyan
                           % ];

                            degeneracy1 = reshape(C.degeneracy, length(imag_freqs), length(real_freqs));
                            degeneracy = degeneracy1(TF);
                            if degeneracy(ii) == 2
                                scatter(real(max_freqs(ii)),imag(max_freqs(ii)),30,'red','filled');%, 'MarkerEdgeColor', 'k');
                                text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d,%d)",floor(fftx),floor(ffty)),'Color','red')
                            else
                                scatter(real(max_freqs(ii)),imag(max_freqs(ii)),30,'k','filled');%, 'MarkerEdgeColor', 'k');
                                text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d,%d)",round(abs(fftx)),round(abs(ffty))),'Color','k')
                            end
                         else
        
                                hue = min(ffty/4,1);
                                
                                saturation = 1;%-min(1,ffty/3);% - fftx/Input.arrayN(1);
            
            
                                %if ffty<1                                         %3
                                %    scatter(real(max_freqs(ii)),imag(max_freqs(ii)),10,'r','filled');%max(80-fftx*30,5)
                                %elseif ffty<2                                        %3
                                %    scatter(real(max_freqs(ii)),imag(max_freqs(ii)),10,'g','filled');%max(80-fftx*30,5)
                                %elseif ffty<3                                        %3
                                %    scatter(real(max_freqs(ii)),imag(max_freqs(ii)),10,'blue','filled');%max(80-fftx*30,5)
                                %elseif ffty>3                                        %3
                                %    scatter(real(max_freqs(ii)),imag(max_freqs(ii)),10,'k','filled');%max(80-fftx*30,5)
                                %end
                                scatter(real(max_freqs(ii)),imag(max_freqs(ii)),10,'k','filled')
                                %scatter(real(max_freqs(ii)),imag(max_freqs(ii)),5,RGB,'filled');%max(80-fftx*30,5)
            
                                %scatter(real(max_freqs(ii)),imag(max_freqs(ii)),5,'red','filled');
                                if fftx+ffty < 100
                                    text(real(max_freqs(ii)),imag(max_freqs(ii)),sprintf("(%d,%d)",round(abs(fftx),0),round(abs(ffty),0)))
                                end
                         end
                   end
                   mark_uncertain(ntext0);
               end
           end
           %colorbar

    end

    function fftval = plot_mode_cyl_TM(modenmb)
        plotnmb = plotnmb+1;
        mode = max_modes(:,modenmb);
        mode_z = [];
        if Input.dimensionality == "cyl_TE"
            % line dipoles (p_x, p_z): label by the stronger component, plot both
            mode = permute(reshape(mode, 2,Input.arrayN(1),Input.arrayN(2)),[2,3,1]);
            modex = mode(:,:,1);
            mode_z = mode(:,:,2);
            if norm(mode_z(:)) > norm(modex(:)), fftval = label_1D(mode_z, modenmb, 2);
            else,                                fftval = label_1D(modex, modenmb, 1);
            end
            mode = modex;
        elseif Input.dimensionality == "2D"
            mode = permute(reshape(mode, 2,Input.arrayN(1),Input.arrayN(2)),[2,3,1]);
            modex = mode(:,:,1);
            modey = mode(:,:,2);
            fftval = label_1D(modex, modenmb, []);
            mode = modex;
        elseif Input.dimensionality == "3D"
            mode = permute(reshape(mode, 3,Input.arrayN(1),Input.arrayN(2)),[2,3,1]);
            modex = mode(:,:,1);
            modey = mode(:,:,2);
            modez = mode(:,:,3);
            fftval = label_1D(modex, modenmb, []);
            mode = modex;
        else
            mode = reshape(mode, Input.arrayN(1),Input.arrayN(2)).';
            fftval = label_1D(mode, modenmb, 1);
        end
        
        
        % Magnitude (grey dashed) and real part (solid) on the left axis, cumulative
        % phase on the right. The global phase is chosen to make the mode as real
        % as possible (1/2 arg(v.'v)), so Re p and the phase do not depend on the
        % arbitrary phase of the eigenvector; the mode is scaled to max |p| = 1.
        v  = mode(:);  w = mode_z(:);
        g  = exp(-0.5i*angle(sum([v; w].^2)));
        sc = max(abs([v; w]));  sc(sc == 0) = 1;
        v  = v*g/sc;  w = w*g/sc;
        if Input.arrayN(1) == 1 && Input.arrayN(min(2,end)) > 1, ldir = 'y'; else, ldir = 'x'; end
        ax = nexttile;
        yyaxis(ax, 'left');  hold(ax, 'on')
        plot(ax, abs(v),  '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        plot(ax, real(v), '-',  'Color', [0 0.447 0.741], 'LineWidth', 1.2);
        if ~isempty(w)                                   % cyl_TE: p_z dotted / dash-dotted
            plot(ax, abs(w),  ':',  'Color', [0.5 0.5 0.5], 'LineWidth', 1);
            plot(ax, real(w), '-.', 'Color', [0 0.447 0.741], 'LineWidth', 1.2);
        end
        ylim(ax, [-1.05 1.05]);  ax.YAxis(1).Color = [0 0.447 0.741];
        ylabel(ax, '|p| (dashed), Re p')
        yyaxis(ax, 'right')
        plot(ax, unwrap(angle(v))/pi, '-', 'Color', [0.85 0.325 0.098], 'LineWidth', 1.2);
        ax.YAxis(2).Color = [0.85 0.325 0.098];
        ylabel(ax, 'phase/\pi (cumulative)')
        yyaxis(ax, 'left')
        xlabel(ax, sprintf('%s (NPs)', ldir));  xlim(ax, [1 max(2, numel(v))]);
        box(ax, 'on')
        
        title(ax, sprintf("(%g%s)%s\n %s",round(fftval,1),alt_str(),conf_str(),freq_label(max_freqs(modenmb))))
   

    end



    function [fftx,ffty] = plot_mode(modenmb)
        if Input.array_type == "double_rectangular"
            mode = reshape(max_modes(1:size(max_modes,1)/2,modenmb), Input.arrayN(1), Input.arrayN(2)).';
        else
            mode = reshape(max_modes(:,modenmb), Input.arrayN(1),Input.arrayN(2)).';
        end
        [fftx, ffty] = label_xy(mode, modenmb, 1);
        %% plotting
        if fftx < maxfftyplt && ffty < maxfftyplt
            nexttile
           % figure('Position',[200,200,450,300])
        
            
            % Magnitude and phase
            mag = abs(mode);
            phase = angle(mode);
            
            % Normalize magnitude to [0,1]
            mag_norm = mag / max(mag(:));
            
            % Hue: map phase [-pi, pi] -> [0,1]
            H = (phase + pi) / (2*pi);
            H = mod(H+0.15,1);
            
            % Saturation: proportional to magnitude
            S = mag_norm;                       % use this for hsv type phase colormap
            %S = zeros(size(mag_norm));         %use this for grayscale
            % Value (brightness):
            % Make small values white (V=1, S=0), large more “colored”
            % You can tune gamma to control contrast
            gamma = 1;  
            
            Vmin = 0.6;  % minimum brightness
            V = 1 - (1 - Vmin) * (mag_norm).^gamma;

            
            % Combine HSV channels
            HSV = cat(3, H, S, V);
            
            % Convert to RGB
            RGB = hsv2rgb(HSV);
            
            % Plot
            imagesc(RGB,'Interpolation','bilinear');
            %axis image off;
            hold on 

           
            if surfplot1D
                surf(abs(mode)+0.1,angle(mode)/pi,'EdgeColor','none');
                surf(abs(mode)+0.1,2*H-1,'EdgeColor','none','FaceColor','interp');
                view(3)
               % zlabel('|p|')
                 zticks([0,0.05,0.1,0.15,0.2])
            end
             zticklabels(["","","0","1","2"])
             ylabel('y (NPs)')
             xlabel('x (NPs)')
             
             yticks([1, 10, 20]);xticks([1, 10, 20])
             %xticklabels([]);yticklabels([])
            colormap hsv
            %colorbar
            %view(2)
            clim([-1,1])
             a = colorbar;
             a.Ticks = [-1, -0.5, 0, 0.5, 1];
            ylabel(a,'Phase/\pi','FontSize',12,'Rotation',270);
            box on
            set(gca, 'LineWidth', 1.5)
   
            title(sprintf("Mode (%.0f,%.0f)%s, %s", ...
                fftx, ffty, conf_str(), freq_label(max_freqs(modenmb))));
            plotnmb = plotnmb + 1;
        end
       
    end

    function  [dim, fftx,ffty] = plot_mode2D(modenmb)
        nmbmodesx = Input.arrayN(1);
         nmbmodesy = Input.arrayN(2);
         [matrx,matry] = meshgrid([1:nmbmodesx],[1:nmbmodesy]);
        %  matrx = ([1:nmbmodesy].'*ones(nmbmodesx,1).');
        %matry = ([1:nmbmodesx].'*ones(nmbmodesy,1).').';
        mode = permute(reshape(max_modes(:,modenmb), 2,nmbmodesx,nmbmodesy),[1,3,2]);
        %charge = calculate_winding_number(mode);
        xpol = squeeze(mode(1,:,:));ypol = squeeze(mode(2,:,:));
        if separate_pols
             %[fftx, ffty] = get_fft_nmbrs(mode);
                nexttile
                surf(abs(xpol), angle(xpol), 'EdgeColor','none');
                %surf(modeabs,'EdgeColor','none')
                 ylabel('y')
                 title(freq_label(max_freqs(modenmb)))
                 clim([-pi,pi])
                 colormap hsv
                 nexttile
                surf(abs(ypol), angle(ypol), 'EdgeColor','none');
    
                colormap hsv
                %colorbar
                %view(2)
                clim([-pi,pi])
                
                title(freq_label(max_freqs(modenmb)))
                fftx = NaN; ffty = NaN;
        end
        if quiverplot
            
            nexttile
            %figure
            %Quiver plot
             %surf2Colormap = pink();
             %tc = colormapToTruecolor(surf2Colormap,-squeeze(rssq(mode(:,:,:),1)));
             
             %imagesc(permute(tc,[2,1,3]));%,'Interpolation','bilinear');

             imagesc(permute(squeeze(rssq(abs(mode),1)),[1,2,3]),'Interpolation','bilinear')
             %colormap(slanCM('freeze', 257))
             %colormap("pink")
             clim([0.*max(rssq(mode,1),[],"all"),max(rssq(mode,1),[],"all")*1])
             hold on
             %colorbar
            
            %circular polarization
            circularity = min(abs(angle(mode(1,:,:))-angle(mode(2,:,:))),abs(-angle(mode(1,:,:))+angle(mode(2,:,:))));
            circularity = squeeze(min(circularity,abs(pi-circularity)));
            S4 = -2*imag(xpol.*conj(ypol)); %stokes parameter, circularity [-1,1]
            U = squeeze(mode(1,:,:)); V = squeeze(mode(2,:,:));
           
            q = quiver(matrx,matry, U , V ,'k');
            %% check symmetry
             R_U = rot90(U);    R_V = rot90(V);
             U_sym_expected = -R_V; V_sym_expected =  R_U;
            err_sym = max(abs(U(:)- U_sym_expected(:))) + max(abs(V(:)-V_sym_expected(:)));
            U_anti_expected =  R_V; V_anti_expected = -R_U;
            
            err_anti = max(abs(U(:)-U_anti_expected(:))) + max(abs(V(:)-V_anti_expected(:)));
            ax = gca;
            if err_sym < 1e-4
                ax.Colormap = slanCM('freeze', 257);
            elseif err_anti < 1e-4
                ax.Colormap = slanCM('flamingo', 257);
            else 
                ax.Colormap = slanCM('hot', 257);
            end
           
            %// Get the current colormap
           % currentColormap = colormap("hot");
            
            %// Now determine the color to make each arrow using a colormap
            %[~, ~, ind] = histcounts(S4, size(currentColormap, 1));
            %[~, ~, ind] = histcounts(abs(S4), linspace(0,1,size(currentColormap,1)));

            %// Now map this to a colormap to get RGB
            %cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
            %cmap(:,:,4) = 255;
            %cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
            
            %// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
            %set(q.Head, ...
            %    'ColorBinding', 'interpolated', ...
            %    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'
            
            %// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
            %set(q.Tail, ...
            %    'ColorBinding', 'interpolated', ...
            %    'ColorData', reshape(cmap(1:2,:,:), [], 4).');
            xpolarisity = rssq(reshape(xpol,1,[]))/rssq(reshape(ypol,1,[]));  
            sums = sum(abs(mode(1:2,:,:)),[2,3]);
            
            polthresh = 2;
            divided = sums(1)./sums(2);
            if divided > polthresh
                dim = 1;
            elseif divided<1/polthresh
                dim =2 ;
            else 
                dim = 3;
            end


            if dim == 1
                polstr = "X-pol,";
                [fftx, ffty] = label_xy(xpol, modenmb, 1);
            elseif dim == 2
                polstr = "Y-pol,";
                [fftx, ffty] = label_xy(ypol, modenmb, 2);
            else
                [fftx, ffty] = label_xy(ypol, modenmb, 2);
                polstr = sprintf("Y(%d,%d),X",floor(fftx),floor(ffty));
                [fftx, ffty] = label_xy(xpol, modenmb, 1);
                %fftx= NaN;ffty= NaN;
            end

            title(sprintf("%s (%d,%d)%s\n %s",polstr,floor(fftx),floor(ffty), conf_str(), freq_label(max_freqs(modenmb))))
           %[fftx, ffty] = label_xy(xpol, modenmb, 1);
           box on 
           set(gca, 'LineWidth',2)


        end
       
    end
 
    function  [dim,fftx,ffty] = plot_mode3D(modenmb)
        
        nmbmodesx = Input.arrayN(1);
         nmbmodesy = Input.arrayN(2);
        %  matrx = ([1:nmbmodesy].*ones(nmbmodesx,1));
        %matry = ([1:nmbmodesx].*ones(nmbmodesy,1)).';
        %matrz = ones(nmbmodesy,1).*ones(nmbmodesx,1).';
        [matrx,matry,matrz] = meshgrid(1:nmbmodesx,1:nmbmodesy,1);
        mode = permute(reshape(max_modes(:,modenmb), 3,nmbmodesx,nmbmodesy),[1,3,2]);
        xpol = squeeze(mode(1,:,:));ypol = squeeze(mode(2,:,:)); zpol = squeeze(mode(3,:,:));
        zmax = max(abs(mode(3,:,:)),[],'all');
        xmax = max(abs(mode(1,:,:)),[],'all');
        ymax = max(abs(mode(2,:,:)),[],'all');
        
        %% Calculate dominating polarization component
        mag = sum(abs(reshape(mode(:,:,:),[],1)));
        zmag = sum(abs(reshape(mode(3,:,:),[],1)))/mag*100;
        xmag = sum(abs(reshape(mode(1,:,:),[],1)))/mag*100;
        ymag = sum(abs(reshape(mode(2,:,:),[],1)))/mag*100;

        
        [~,dim] = max(sum(abs(mode(1:2,:,:)),[2,3]));
        if ismembertol(xmag,ymag,0.1) % if x and y are alike
            dim =4;
        end
        if zmag > 1.3*max(xmag,ymag) %if z component cominates anyway
            dim = 3;
        end

        [xfftx, xffty, xc] = label_xy(xpol, modenmb, 1);
        [yfftx, yffty, yc] = label_xy(ypol, modenmb, 2);
        [zfftx, zffty, zc] = label_xy(zpol, modenmb, 3);
        
        if dim == 1 
            polstr = sprintf("\\textbf{X %.1f\\%% (%d,%d)},Y %.1f\\%% (%d,%d),Z %.1f\\%% (%d,%d)\n %s",xmag,floor(xfftx),floor(xffty),...
                ymag, floor(yfftx),floor(yffty),...
                zmag, floor(zfftx),floor(zffty), freq_label(max_freqs(modenmb)));
            fftx = xfftx;ffty = xffty;
        elseif dim == 2
            polstr = sprintf("X %.1f\\%% (%d,%d),\\textbf{Y %.1f\\%% (%d,%d)},Z %.1f\\%% (%d,%d)\n %s",xmag,floor(xfftx),floor(xffty),...
                ymag, floor(yfftx),floor(yffty),...
                zmag, floor(zfftx),floor(zffty), freq_label(max_freqs(modenmb)));
            fftx = yfftx;ffty = yffty;
        elseif dim == 3
            
            polstr = sprintf("X %.1f\\%% (%d,%d),Y %.1f\\%% (%d,%d),\\textbf{Z %.1f\\%% (%d,%d)}\n %s",xmag,floor(xfftx),floor(xffty),...
                ymag, floor(yfftx),floor(yffty),...
                zmag, floor(zfftx),floor(zffty), freq_label(max_freqs(modenmb)));
            fftx = zfftx;ffty = zffty;
            %[fftx, ffty] = get_fft_nmbrs(zpol);
            %fftx= NaN;ffty= NaN;
        else
            polstr = sprintf("\\textbf{X %.1f\\%% (%d,%d),Y %.1f\\%% (%d,%d)},Z %.1f\\%% (%d,%d)\n %s",xmag,floor(xfftx),floor(xffty),...
                ymag, floor(yfftx),floor(yffty),...
                zmag, floor(zfftx),floor(zffty), freq_label(max_freqs(modenmb)));
            fftx = xfftx;ffty = xffty;
        end

        cs = [xc yc zc];
        if dim == 4, last_conf = xc; else, last_conf = cs(dim); end
        if max(fftx,ffty) < 40
            boolplot = 1;
        else
            boolplot = 0;
        end
        if separate_pols
             %[fftx, ffty] = get_fft_nmbrs(mode);
                nexttile
                surf(abs(xpol), angle(xpol), 'EdgeColor','none');
                %surf(modeabs,'EdgeColor','none')
                 ylabel('y')
                 title(freq_label(max_freqs(modenmb)))
                 clim([-pi,pi])
                 colormap hsv
                 nexttile
                surf(abs(ypol), angle(ypol), 'EdgeColor','none');
    
                colormap hsv
                %colorbar
                %view(2)
                clim([-pi,pi])
                
                title(freq_label(max_freqs(modenmb)))
                fftx = NaN; ffty = NaN;
        end
        if quiverplot
            if boolplot
                plotnmb = plotnmb+1;
                nexttile
  
                M = max(abs(mode(3,:,:)),[],'all');   % symmetric range
                clim([-M M])                        % ensures 0 is exactly the center
    
         
                imagesc(squeeze(abs(mode(3,:,:)).*sign(real(mode(3,:,:)))),'Interpolation','bilinear')
                colormap(slanCM('bwr', 257))
                
                M = max(abs(mode(3,:,:)),[],'all');
                clim([-M M])
                
                %colorbar
    
                 hold on
                
                %circular polarization
                S4 = -2*imag(xpol.*conj(ypol)); %stokes parameter, circularity [-1,1]
                if imagescisZ
                    q = quiver(matrx,matry, squeeze(abs(mode(1,:,:)).*sign(real(mode(1,:,:)))),squeeze(abs(mode(2,:,:)).*sign(real(mode(2,:,:)))),'Color','k');
                else
                    q = quiver3(matrx,matry, matrz, squeeze(mode(1,:,:)),squeeze(mode(2,:,:)), squeeze(mode(3,:,:)));
                    imagesc()
                end
                xticks([]);yticks([])
            
            %axis equal
            %// Get the current colormap
            %{
            currentColormap = colormap("hot");
            
            %// Now determine the color to make each arrow using a colormap
            [~, ~, ind] = histcounts(S4, size(currentColormap, 1));
            [~, ~, ind] = histcounts(abs(S4), linspace(0,1,size(currentColormap,1)));

            %// Now map this to a colormap to get RGB
            cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
            cmap(:,:,4) = 255;
            cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
            
            %// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
            set(q.Head, ...
                'ColorBinding', 'interpolated', ...
                'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'
            
            %// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
            set(q.Tail, ...
                'ColorBinding', 'interpolated', ...
                'ColorData', reshape(cmap(1:2,:,:), [], 4).');
            xpolarisity = rssq(reshape(xpol,1,[]))/rssq(reshape(ypol,1,[]));  
            %}
           
            title(sprintf("%s",polstr),'Interpreter','latex')
            end


        end
       
    end





    %% ---------------------------------------------------------------
    %  Sine-basis mode labels (dda_mode_basis, as in Track_eigenmodes)
    function L = sine_setup()
        L = struct('ok', false);
        if ~use_sine, return; end
        if isempty(Array) || ~isfield(Array, 'r_i')
            fprintf('No particle positions in this file: FFT mode labels instead of sine labels.\n');
            return
        end
        rpos = double(Array.r_i);
        nx = Input.arrayN(1);  ny = Input.arrayN(min(2,end));
        L.M     = dda_mode_basis(rpos, nx, ny, o.basis, Input);
        L.ncomp = max(1, round(size(C.min_eig_vectors,1)/size(rpos,1)));
        L.ncell = nx*ny;
        % grid indices of the particles (for "nodes")
        tol = 1e-3*min(Input.p(1:2));
        [~, ~, L.ix] = unique(round(rpos(:,1)/tol));
        [~, ~, L.iy] = unique(round(rpos(:,2)/tol));
        L.nxv = max(L.ix);  L.nyv = max(L.iy);
        L.isgrid = L.nxv*L.nyv == size(rpos,1) && ...
                   numel(unique(L.iy + (L.ix-1)*L.nyv)) == size(rpos,1);
        L.per = startsWith([L.M.used], "fourier") * [1 0] + endsWith([L.M.used], "fourier") * [0 1];
        if use_nodes && ~L.isgrid
            fprintf('"nodes" labels need particles on a grid; sine labels instead.\n');
        end
        if use_k, L.K = k_setup(rpos, [nx ny]); end
        L.ok    = true;
    end

    function [m, conf] = sine_label(modenmb, comp)
        % m = [mx my component]; comp = [] picks the best component
        sc = dda_mode_basis(LB.M, double(max_modes(:,modenmb)), LB.ncomp);
        sc = reshape(sc, LB.ncell, LB.ncomp);
        if isempty(comp)
            [~, l] = max(sc(:));
            [~, comp] = ind2sub(size(sc), l);
        end
        comp = min(comp, LB.ncomp);
        [best, icell] = max(sc(:,comp));
        conf = best / max(sum(sc(:,comp)), eps);     % share of this component's power
        [ix, iy] = ind2sub([LB.M.nx LB.M.ny], icell);
        m = [LB.M.lab_x(ix), LB.M.lab_y(iy), comp];
        last_alt = [NaN NaN];
        if use_nodes && LB.isgrid
            [nxn, nyn, conf] = node_counts(double(max_modes(comp:LB.ncomp:end, modenmb)));
            if ~LB.per(1), m(1) = nxn; end
            if ~LB.per(2), m(2) = nyn; end
        elseif use_k
            [mk, conf, alt] = k_label(double(max_modes(comp:LB.ncomp:end, modenmb)));
            for dd = 1:2
                if ~LB.per(dd), m(dd) = mk(dd); last_alt(dd) = alt(dd); end
            end
        end
        last_conf = conf;
    end

    function K = k_setup(rpos, n)
        % |k| grids, calibration and step signs (shared with Track_eigenmodes)
        K = dda_mode_basis("k_setup", rpos, n);
    end

    function [mk, conf, alt] = k_label(v)
        % coupling |v.' e_k|^2 to the plane waves exp(i(+-kx x +- ky y)), summed over the signs
        v = v(:);
        if o.k_step, v = v.*LB.K.flip; end               % undo the pi steps of shifted lattices
        Ck = dda_mode_basis("k_spectrum", LB.K, v);      % nky x nkx
        [c1, l1] = max(Ck(:));  [iy, ix] = ind2sub(size(Ck), l1);
        kx = LB.K.k{1}(ix);  ky = LB.K.k{2}(iy);
        near = abs(LB.K.k{2} - ky) < LB.K.dk(2) & abs(LB.K.k{1}.' - kx) < LB.K.dk(1);
        Ck(near) = 0;
        [c2, l2] = max(Ck(:));  [iy2, ix2] = ind2sub(size(Ck), l2);
        conf = 1 - c2/max(c1, eps);
        mk  = [k_to_m(1, kx), k_to_m(2, ky)];
        alt = [k_to_m(1, LB.K.k{1}(ix2)), k_to_m(2, LB.K.k{2}(iy2))];
    end

    function m = k_to_m(dd, k)
        [~, i] = min(abs(LB.K.km{dd} - k));
        m = i - 1;
    end

    function [nxn, nyn, conf] = node_counts(v)
        % sign changes of the mode along x (row with the most power) and y (column)
        ph = 0.5*angle(sum(v.^2));                   % global phase making v as real as possible
        w  = v*exp(-1i*ph);
        conf = max(0, 2*sum(real(w).^2)/max(sum(abs(w).^2), eps) - 1);   % 1 standing, 0 travelling
        G = zeros(LB.nyv, LB.nxv);
        G(sub2ind(size(G), LB.iy, LB.ix)) = real(w);
        [~, rr] = max(sum(G.^2, 2));  [~, cc] = max(sum(G.^2, 1));
        nxn = sign_changes(G(rr,:));  nyn = sign_changes(G(:,cc));
    end

    function n = sign_changes(x)
        x = x(abs(x) > 1e-3*max(abs(x)));            % skip near-zero values at the nodes
        n = nnz(diff(sign(x)));
    end

    function [fx, fy, conf] = label_xy(compmap, modenmb, comp)
        % (mx, my) of one dipole component: sine basis, or the FFT peak
        conf = NaN;
        if LB.ok
            [m, conf] = sine_label(modenmb, comp);
            fx = m(1);  fy = m(2);
        else
            [fx, fy] = get_fft_nmbrs(compmap);
        end
    end

    function v = label_1D(vec, modenmb, comp)
        % mode number along a line array (or cylinders along x)
        if LB.ok
            m = sine_label(modenmb, comp);
            if LB.M.nx == 1 && LB.M.ny > 1, v = m(2); else, v = m(1); end
        else
            v = get_fft_nmbr_1D(vec);
        end
    end

    function s = alt_str()
        % runner-up of a "k" label along the line, when the label is unclear
        s = '';
        if ~(use_k && LB.ok) || isnan(last_conf) || last_conf >= o.min_confidence, return; end
        alts = last_alt;
        if LB.M.nx == 1 && LB.M.ny > 1, alts = last_alt(2); elseif LB.M.ny == 1, alts = last_alt(1); end
        alts = alts(~isnan(alts));
        if ~isempty(alts), s = sprintf('|%s', strjoin(string(alts), ',')); end
    end

    function s = conf_str()
        s = '';
        if LB.ok && ~isnan(last_conf), s = sprintf(' %.0f%%', 100*last_conf); end
    end

    function mark_uncertain(ntext0)
        % grey label with '?' when the sine label is unclear
        if ~LB.ok || isnan(last_conf) || last_conf >= o.min_confidence, return; end
        h = findobj(ax1, 'Type', 'text');
        if numel(h) > ntext0
            h(1).String = [char(h(1).String) '?'];
            h(1).Color  = [0.5 0.5 0.5];
        end
    end

    function plot_label_map()
        % Label of the smallest-|lambda| eigenvector at every grid point
        V  = C.min_eig_vectors;
        nk = size(V, 2);
        lab = zeros(1, nk);  conf = zeros(1, nk);
        nx = LB.M.nx;  ny = LB.M.ny;
        for c0 = 1:2000:nk                                   % in blocks: V can be large
            cols = c0:min(nk, c0+1999);
            sc = dda_mode_basis(LB.M, double(V(:,cols)), LB.ncomp);
            sc = reshape(sum(reshape(sc, LB.ncell, LB.ncomp, []), 2), LB.ncell, []);  % all components
            [best, icell] = max(sc, [], 1);
            lab(cols) = icell;  conf(cols) = best ./ max(sum(sc,1), eps);
        end
        [ix, iy] = ind2sub([nx ny], lab);
        maps = {};  names = {};
        if nx > 1, maps{end+1} = LB.M.lab_x(ix); names{end+1} = 'm_x'; end
        if ny > 1, maps{end+1} = LB.M.lab_y(iy); names{end+1} = 'm_y'; end
        if isempty(maps), return; end
        figm = dda_figure([100 100 700*numel(maps) 450], 'Name', 'mode label map');
        tlm  = tiledlayout(figm, 1, numel(maps), 'TileSpacing', 'compact');
        for q = 1:numel(maps)
            axm = nexttile(tlm);
            Mq  = reshape(double(maps{q}), length(imag_freqs), length(real_freqs));
            Aq  = reshape(conf, length(imag_freqs), length(real_freqs));
            if re_nm
                heatmap_nm(axm, real_freqs, imag_freqs, Mq);
            else
                hi = imagesc(axm, real_freqs, imag_freqs, Mq);
                hi.AlphaData = 0.35 + 0.65*min(1, Aq/max(o.min_confidence, eps));   % faded where unclear
            end
            set(axm, 'YDir', 'normal'); colormap(axm, turbo); hold(axm, 'on')
            cb = colorbar(axm); cb.Label.String = sprintf('mode number %s', names{q});
            plot(axm, real(max_freqs), imag(max_freqs), 'k.', 'MarkerSize', 12)   % the minima of |lambda|
            xlim(axm, xlims); ylim(axm, ylims)
            if re_nm, xlabel(axm, 'Wavelength (nm)'); else, xlabel(axm, 'Re[\omega] (eV)'); end
            ylabel(axm, '-Im[\omega] (meV)')
            title(axm, sprintf('%s of the smallest-|\\lambda| eigenvector (%s basis)', names{q}, LB.M.used))
        end
        set(0, 'currentfigure', fig1);
    end

    function [fftx,ffty] = get_fft_nmbrs(mode)
        pad = 2^nextpow2(sqrt(numel(mode))*8);
        %pad = 20;
        pad2 = ceil(pad/2);
        if is_nufft
            os = 1;
            ffted = nufft_Gauss(X,Y,reshape(mode,[],1),os,pad,2,4);
       
        else
            os = 1;
            ffted = fft2((mode),pad,pad).'; %unknown why permutation needed
            
        end
        if plot_fft
             plot_fft2(X,Y,ffted,os);
        end
        fftedshifted = fftshift(ffted);
       % figure
        %imagesc(abs(ffted))
        %ffted = ffted();'
        dx = Input.p(1);                 % was 700e-9
        L = dx*size(mode,1);
        dy = Input.p(2);
        ks = 2*pi/dx; %Sampling frequency as wavenumber
       %
        %imagesc(ks*(0:pad-1)/pad, ks*(0:pad-1)/pad, abs(fftedshifted))
       % figure(20)
       % imagesc([-pad/2:pad/2-1]/pad*2*20, [-pad/2:pad/2-1]/pad*2*20, abs(fftedshifted))
       % xlim([0,5]);ylim([0,5]);clim([0,max(abs(ffted),[],"all")])
        %colormap parula
        %colorbar
        [maxffted,maxrow] = max(abs(fftedshifted(1:pad2,1:pad2)),[],1);
        [maxval,maxcol] = max(maxffted);
        maxrow = maxrow(maxcol);
        fftytemp = pad2-maxcol;fftxtemp = pad2-maxrow;
        ffty = fftytemp*size(mode,1)/pad2*os; %change of ffty to fftx happens here because mode should be permuted
        fftx = fftxtemp*size(mode,2)/pad2*os;
        kx = fftx*ks/size(mode,1);
        ky = ffty*ks/size(mode,2);

        F = fftshift(fft2(mode,pad,pad));

        [~,idx] = max(abs(F(:)));
        [row,col] = ind2sub(size(F),idx);
        
        kx_axis = 2*pi*(-pad/2:pad/2-1)/(pad*dx);
        ky_axis = 2*pi*(-pad/2:pad/2-1)/(pad*dy);
        
        kx_peak = kx_axis(col);
        ky_peak = ky_axis(row);

        m = max(abs(kx_peak)*L/pi - 1,0);
        n = max(abs(ky_peak)*L/pi - 1,0);
        A = abs(F).^2;

        % extract line through peak
        profile = A(maxrow,:);
        
        halfmax = max(profile)/2;
        
        idx = find(profile > halfmax);
        
        FWHM_pixels = idx(end)-idx(1)+1;
        dkx = 2*pi/(pad*dx);
        sigma_kx = FWHM_pixels * dkx /(2*sqrt(2*log(2)));

        return


    end
      function [fftval] = get_fft_nmbr_1D(mode)
        pad = 2^nextpow2(length(mode)*8);
        pad2 = ceil(pad/2);
        ffted = squeeze(fftshift(fft((mode.'),pad)));
        %figure
       %plot(abs(ffted))
        %ff
        % ted = ffted();
        %imagesc(abs(ffted))
        [maxffted,maxindex] = max(ffted(1:pad2),[],"all");
        %[~,maxcol] = max(maxffted);
        ffttemp = pad2-maxindex;
        fftval = ffttemp*length(mode)/pad2; %normalize
        return


      end

    function [charge] = calculate_winding_number(mode)
        %% Extract dipoles

            px = squeeze(mode(1,:,:));
            py = squeeze(mode(2,:,:));
            
            %% Fourier transform to momentum space
            nk = 64;
            Px = fftshift(fft2(real(px),nk,nk));
            Py = fftshift(fft2(real(py),nk,nk));
            S1 = abs(Px).^2 - abs(Py).^2;
            S2 = 2*real(Px.*conj(Py));
            %% Polarization angle
            
            %phi = atan2(real(Py), real(Px));
            
            % For complex dipoles, a better choice is:
            %phi = 0.5 * angle((Px + 1i*Py).^2);
            % polarisation angle from the Stokes parameters S1, S2
            phi = 0.5 * angle(S1+1i*S2);

            %% Gamma point
            
            N = size(phi,1);
            cx = floor(N/2)+1;
            cy = floor(N/2)+1;
            
            %% Square contour size
            % half-way from Gamma to X/Y point
            
            R = floor(N/8);
            
            %% Build contour (counterclockwise)
            
            x = [];
            y = [];
            
            % bottom edge
            x = [x, cx-R:cx+R-1];
            y = [y, (cy-R)*ones(1,2*R)];
            
            % right edge
            x = [x, (cx+R)*ones(1,2*R)];
            y = [y, cy-R:cy+R-1];
            
            % top edge
            x = [x, cx+R: -1 : cx-R+1];
            y = [y, (cy+R)*ones(1,2*R)];
            
            % left edge
            x = [x, (cx-R)*ones(1,2*R)];
            y = [y, cy+R:-1:cy-R+1];
            
            %% Extract phase on contour
            
            phaseContour = phi(sub2ind(size(phi),y,x));
            
            %% Winding number
            
            wrap = @(a) mod(a+pi,2*pi)-pi;
            
            dphi = wrap(diff([phaseContour phaseContour(1)]));
            
            winding = sum(dphi)/(2*pi);
            
            fprintf('Winding number around Gamma = %.6f\n',winding);

            
            %% Function for wrapping angles to [-pi, pi]
            
            wrap = @(x) mod(x + pi, 2*pi) - pi;
            
            %% Compute winding number on each plaquette
            
            [Nx, Ny] = size(phi);
            
            winding = zeros(Nx-1,Ny-1);
            
            for ix = 1:Nx-1
                for iy = 1:Ny-1
            
                    p1 = phi(ix,iy);
                    p2 = phi(ix+1,iy);
                    p3 = phi(ix+1,iy+1);
                    p4 = phi(ix,iy+1);
            
                    circulation = ...
                        wrap(p2-p1) + ...
                        wrap(p3-p2) + ...
                        wrap(p4-p3) + ...
                        wrap(p1-p4);
            
                    winding(ix,iy) = circulation/(2*pi);
            
                end
            end
            
            %% Round to nearest integer charge
            
            charge = round(winding);
            
            %% Find nonzero charges
            
            [row,col] = find(abs(charge) > 0);
            
            fprintf('Found %d singularities\n',length(row));
            
            for n = 1:length(row)
                fprintf('(%d,%d) charge = %+d\n',...
                    row(n), col(n), charge(row(n),col(n)));
            end
            
            %% Plot
            
            figure;
            imagesc(charge);
            axis image;
            colorbar;
            title('Polarization winding number');
            xlabel('k_x');
            ylabel('k_y');
            figure;

            quiver(real(Px),real(Py),'k');
            
            axis equal;
            title('Momentum-space polarization');

    end
     



end

function heatmap_nm(ax, re_eV, im_meV, C)
% Heatmap on a wavelength axis. The grid is uniform in energy, so in nm the
% cells have uneven widths: draw each cell between the midpoints of its
% neighbours (as imagesc does for even spacing), converted to nm.
edges = @(v) [v(1) - (v(2)-v(1))/2, (v(1:end-1) + v(2:end))/2, v(end) + (v(end)-v(end-1))/2];
if numel(re_eV) > 1, xe = 1239.84187./edges(re_eV(:).'); else, xe = 1239.84187./re_eV*[0.999 1.001]; end
if numel(im_meV) > 1, ye = edges(im_meV(:).');          else, ye = im_meV + [-0.5 0.5]; end
[X, Y] = meshgrid(xe, ye);
Cp = [C, C(:,end); C(end,:), C(end,end)];              % 'flat' faces use the first corner
surface(ax, X, Y, zeros(size(X)), Cp, 'EdgeColor', 'none', 'FaceColor', 'flat');
view(ax, 2); axis(ax, 'tight');
end

function o = complex_plot_defaults()
% Values used when the file has no Input.plotting.complex (same as Make_predefined_Inputs)
o = struct('full_range', true, 're_lims', [1.12 1.18], 'im_lims', [25 30], 'modes', true, ...
           'tiles_per_figure', 20, 'max_fft_order', 1000, 'quiver', true, 'imagesc_is_z', true, ...
           'separate_pols', false, 'surf_1D', true, 'plot_fft', true, 'nufft', false, ...
           'real_axis', "eV", 'labels', "sine", 'basis', "auto", 'min_confidence', 0.25, ...
           'label_map', false, 'k_step', false);
end

function a = merge_fields(a, b)
f = fieldnames(b);
for n = 1:numel(f), a.(f{n}) = b.(f{n}); end
end
