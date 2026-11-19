function [edit, ok] = DDA_array_editor(Input, before)
% DDA_ARRAY_EDITOR  Show the array the inputs construct and edit its particles.
%
%   [edit, ok] = DDA_array_editor(Input)
%   [edit, ok] = DDA_array_editor(Input, before)   before: Input before GUI changes
%                                                   (derived inputs are updated first)
%
% Opened from the Array tab of DDA_parameter_GUI. Builds the array with
% Make_generic_Array, exactly as the run will, and shows every particle at its
% position with its radius, coloured by z.
%
% Mouse
%   click                 select a particle (and only it)
%   Shift- or Ctrl-click  add / remove a particle from the selection
%   drag on empty space   select the particles inside the box (Shift/Ctrl: add)
%   drag a selected one   move the whole selection
% Panel
%   x, y, z, radius (nm)  set for all selected particles (empty = mixed values)
%   move by dx, dy, dz    shift the selection
%   undo, reset selected, reset all, select all, clear selection
%
% Returns edit = [] when nothing differs from the constructed array, otherwise
% the struct stored in Input.array_edit (applied by dda_apply_array_edit):
%   enable, r_i (N x 3, m), radius (N x 1, m, NaN = Input.particle.r),
%   base_hash, n_moved, n_resized, hash, base_inputs
% ok = false if the window was cancelled or closed.

ok = false;
edit = [];
if isfield(Input, 'array_edit'), edit = Input.array_edit; end
if nargin > 1 && ~isempty(before), Input = dda_update_derived_inputs(Input, before); end

%% The array the inputs construct
Ib = rmfield_if(Input, 'array_edit');
Ib.N = prod(Ib.arrayN);
tic; Ib.toc1 = toc;
[A0, ~] = Make_generic_Array(Ib);
base_r = double(A0.r_i);
if size(base_r, 2) < 3, base_r(:,3) = 0; end
N = size(base_r, 1);
base_hash = dda_checksum(A0.r_i);
r_default = Input.particle.r;
lorentz = any(strcmpi(string(Input.sample.type), ["Lorentzian","lorentzian Au","nonidentical"]));
walls   = any(string(Input.array_type) == ["Outerwall","HG_wall","Gaussian","Inverted_Gaussian"]);
radius_editable = ~lorentz && ~walls;
spacing = min(Input.p(1:2));

pos = base_r;  rad = nan(N, 1);  notes = strings(0);
if ~isempty(edit) && isfield(edit, 'base_hash') && string(edit.base_hash) == base_hash && size(edit.r_i,1) == N
    pos = double(edit.r_i);  rad = edit.radius(:);
elseif ~isempty(edit) && isfield(edit, 'enable') && edit.enable
    notes(end+1) = "Earlier layout edits were made for a different array and are not shown; Done replaces them.";
end
if ~radius_editable
    if lorentz, notes(end+1) = "Radius editing: not for Lorentzian particles (they have no radius).";
    else,       notes(end+1) = "Radius editing: not for wall/Gaussian lattices (they set their own radii)."; end
end
if contains(string(Input.inversion_tool), "fft")
    notes(end+1) = "Note: FFT-based solvers need the regular lattice; an edited layout needs eig or another solver.";
end
sel = false(N, 1);
undo_stack = {};
drag_mode = "";  p0 = [0 0];  pos0 = pos;  moved_in_drag = false;

%% Window
fig = uifigure('Name', sprintf('Array layout: %d particles', N), 'Position', [80 60 980 640], ...
    'WindowStyle', 'modal');
fig.CloseRequestFcn = @(~,~) finish(false);
g = uigridlayout(fig, [1 2]);  g.ColumnWidth = {'1x', 270};
ax = uiaxes(g);
hold(ax, 'on'); axis(ax, 'equal'); box(ax, 'on');
xlabel(ax, 'x (\mum)'); ylabel(ax, 'y (\mum)');
try, disableDefaultInteractivity(ax); catch, end      % clicks select instead of panning
colormap(ax, parula);
cb = colorbar(ax); cb.Label.String = 'z (nm)';

nv = 24;  th = linspace(0, 2*pi, nv+1);  th(end) = [];
Fall = reshape(1:N*nv, nv, N).';
hp = patch(ax, 'Faces', Fall, 'Vertices', zeros(N*nv, 2), 'FaceVertexCData', zeros(N,1), ...
    'FaceColor', 'flat', 'EdgeColor', [0.25 0.25 0.25], 'PickableParts', 'none');
hs = patch(ax, 'Faces', Fall(1,:), 'Vertices', zeros(N*nv, 2), 'FaceColor', 'none', ...
    'EdgeColor', [0.85 0.1 0.1], 'LineWidth', 2, 'PickableParts', 'none', 'Visible', 'off');
hbox = plot(ax, NaN, NaN, '--', 'Color', [0.1 0.1 0.1], 'PickableParts', 'none');

side = uigridlayout(g, [17 3]);
side.ColumnWidth = {80, '1x', '1x'};
side.RowHeight = [repmat({24}, 1, 15), {'1x', 30}];
side.RowSpacing = 5;
info_lbl = uilabel(side, 'Text', '', 'WordWrap', 'on');  place(info_lbl, 1, [1 3]);  side.RowHeight{1} = 40;
place(uilabel(side, 'Text', 'Selected particles (nm)', 'FontWeight', 'bold'), 2, [1 3]);
names = ["x","y","z","radius"];
fld = gobjects(1, 4);
for q = 1:4
    place(uilabel(side, 'Text', names(q)), 2+q, 1);
    fld(q) = uieditfield(side, 'text', 'HorizontalAlignment', 'left', ...
        'ValueChangedFcn', @(src,~) set_value(q, src));
    place(fld(q), 2+q, [2 3]);
end
place(uilabel(side, 'Text', 'Move selected by (nm)', 'FontWeight', 'bold'), 7, [1 3]);
dfld = gobjects(1, 3);
dl = ["dx","dy","dz"];
for q = 1:3
    place(uilabel(side, 'Text', dl(q)), 7+q, 1);
    dfld(q) = uieditfield(side, 'numeric', 'Value', 0, 'HorizontalAlignment', 'left');
    place(dfld(q), 7+q, [2 3]);
end
place(uibutton(side, 'Text', 'Move', 'ButtonPushedFcn', @(~,~) move_by()), 11, [2 3]);
place(uibutton(side, 'Text', 'Select all', 'ButtonPushedFcn', @(~,~) select_set(true(N,1))), 12, [1 2]);
place(uibutton(side, 'Text', 'Clear selection', 'ButtonPushedFcn', @(~,~) select_set(false(N,1))), 12, 3);
place(uibutton(side, 'Text', 'Undo', 'ButtonPushedFcn', @(~,~) undo()), 13, [1 2]);
place(uibutton(side, 'Text', 'Reset selected', 'ButtonPushedFcn', @(~,~) reset_sel()), 13, 3);
place(uibutton(side, 'Text', 'Reset all particles', 'ButtonPushedFcn', @(~,~) reset_all()), 14, [1 3]);
place(uilabel(side, 'Text', ['Click: select. Shift/Ctrl-click: add or remove. Drag on empty ' ...
    'space: box select. Drag a selected particle: move the selection.'], ...
    'WordWrap', 'on', 'FontColor', [0.35 0.35 0.35]), 15, [1 3]);
side.RowHeight{15} = 48;
nl = uilabel(side, 'Text', strjoin(notes, newline), 'WordWrap', 'on', 'FontColor', [0.6 0.25 0]);
place(nl, 16, [1 3]);
place(uibutton(side, 'Text', 'Cancel', 'ButtonPushedFcn', @(~,~) finish(false)), 17, 1);
place(uibutton(side, 'Text', 'Done', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~,~) finish(true)), 17, [2 3]);
if ~radius_editable, fld(4).Enable = 'off'; end

fig.WindowButtonDownFcn   = @(~,~) mouse_down();
fig.WindowButtonMotionFcn = @(~,~) mouse_move();
fig.WindowButtonUpFcn     = @(~,~) mouse_up();

draw_all();
pad = max(spacing, 2*r_default)*1e6;
xlim(ax, [min(pos(:,1))*1e6 - pad, max(pos(:,1))*1e6 + pad + eps]);
ylim(ax, [min(pos(:,2))*1e6 - pad, max(pos(:,2))*1e6 + pad + eps]);
show_selection();
uiwait(fig);

%% ========================================================================
    function place(h, row, col), h.Layout.Row = row; h.Layout.Column = col; end

    function R = radius_now()
        R = rad;  R(isnan(R)) = r_default;
    end

    function V = vertices()
        Rm = radius_now()*1e6;
        X = pos(:,1)*1e6 + Rm.*cos(th);
        Y = pos(:,2)*1e6 + Rm.*sin(th);
        V = [reshape(X.', [], 1), reshape(Y.', [], 1)];
    end

    function draw_all()
        V = vertices();
        z = pos(:,3)*1e9;
        hp.Vertices = V;  hp.FaceVertexCData = z;
        if max(z) - min(z) < 1e-9, clim(ax, min(z) + [-1 1]); else, clim(ax, [min(z) max(z)]); end
        draw_selection(V);
    end

    function draw_selection(V)
        if nargin < 1, V = hp.Vertices; end
        hs.Vertices = V;
        if any(sel)
            hs.Faces = Fall(sel, :);  hs.Visible = 'on';
        else
            hs.Visible = 'off';
        end
    end

    function show_selection()
        ns = nnz(sel);
        moved = any(abs(pos - base_r) > 1e-15, 2);
        resized = ~isnan(rad) & abs(rad - r_default) > 1e-15;
        txt = sprintf('%d of %d selected. Edited: %d moved, %d resized.', ns, N, nnz(moved), nnz(resized));
        ov = overlaps();
        if ov > 0, txt = sprintf('%s\nWarning: %d overlapping pair(s).', txt, ov); end
        info_lbl.Text = txt;
        vals = {pos(:,1)*1e9, pos(:,2)*1e9, pos(:,3)*1e9, radius_now()*1e9};
        for qq = 1:4
            if ns == 0
                fld(qq).Value = '';  fld(qq).Placeholder = '';
                fld(qq).Enable = 'off';
                continue
            end
            v = vals{qq}(sel);
            if max(v) - min(v) < 1e-6
                fld(qq).Value = num2str(v(1), 8);  fld(qq).Placeholder = '';
            else
                fld(qq).Value = '';  fld(qq).Placeholder = 'mixed';
            end
            fld(qq).Enable = (qq < 4 || radius_editable);
        end
    end

    function n = overlaps()
        % pairs of particles closer than the sum of their radii
        R = radius_now();
        if N > 4000, n = 0; return; end
        D = hypot(pos(:,1) - pos(:,1).', pos(:,2) - pos(:,2).');
        D = hypot(D, pos(:,3) - pos(:,3).');
        S = R + R.';
        M = D < S*(1 - 1e-9);
        M(1:N+1:end) = false;
        n = nnz(M)/2;
    end

    function push_undo()
        undo_stack{end+1} = {pos, rad};
        if numel(undo_stack) > 50, undo_stack(1) = []; end
    end

    function undo()
        if isempty(undo_stack), return; end
        s = undo_stack{end};  undo_stack(end) = [];
        pos = s{1};  rad = s{2};
        draw_all();  show_selection();
    end

    function set_value(q, src)
        if ~any(sel), return; end
        txt = strtrim(src.Value);
        if isempty(txt), show_selection(); return; end
        v = str2double(txt);
        if isnan(v), uialert(fig, sprintf('"%s" is not a number.', txt), 'Invalid value'); show_selection(); return; end
        if q == 4 && v <= 0, uialert(fig, 'The radius must be positive.', 'Invalid value'); show_selection(); return; end
        push_undo();
        if q <= 3, pos(sel, q) = v*1e-9; else, rad(sel) = v*1e-9; end
        draw_all();  show_selection();
    end

    function move_by()
        if ~any(sel), return; end
        d = [dfld(1).Value, dfld(2).Value, dfld(3).Value]*1e-9;
        if all(d == 0), return; end
        push_undo();
        pos(sel, :) = pos(sel, :) + d;
        draw_all();  show_selection();
    end

    function select_set(s)
        sel = s;  draw_selection();  show_selection();
    end

    function reset_sel()
        if ~any(sel), return; end
        push_undo();
        pos(sel, :) = base_r(sel, :);  rad(sel) = NaN;
        draw_all();  show_selection();
    end

    function reset_all()
        push_undo();
        pos = base_r;  rad(:) = NaN;
        draw_all();  show_selection();
    end

    %% Mouse --------------------------------------------------------------
    function p = cursor()
        cp = ax.CurrentPoint;  p = cp(1, 1:2);
    end

    function tf = inside(p)
        xl = xlim(ax); yl = ylim(ax);
        tf = p(1) >= xl(1) && p(1) <= xl(2) && p(2) >= yl(1) && p(2) <= yl(2);
    end

    function j = hit(p)
        % the particle under the cursor (radius, at least 1% of the view), or []
        d = hypot(pos(:,1)*1e6 - p(1), pos(:,2)*1e6 - p(2));
        tol = max(radius_now()*1e6, 0.01*diff(xlim(ax)));
        [m, j] = min(d - tol);
        if m > 0, j = []; end
    end

    function mouse_down()
        p = cursor();
        if ~inside(p), return; end
        typ = string(fig.SelectionType);
        j = hit(p);
        if any(typ == ["extend","alt"])
            if ~isempty(j)
                sel(j) = ~sel(j);  draw_selection();  show_selection();
            else
                drag_mode = "box_add";  p0 = p;
            end
            return
        end
        if typ ~= "normal", return; end
        if ~isempty(j)
            if ~sel(j), sel(:) = false; sel(j) = true; draw_selection(); show_selection(); end
            drag_mode = "drag";  p0 = p;  pos0 = pos;  moved_in_drag = false;
        else
            drag_mode = "box";  p0 = p;
        end
    end

    function mouse_move()
        if drag_mode == "", return; end
        p = cursor();
        switch drag_mode
            case "drag"
                d = (p - p0)*1e-6;
                if any(d ~= 0)
                    moved_in_drag = true;
                    pos(sel, 1:2) = pos0(sel, 1:2) + d;
                    V = vertices();  hp.Vertices = V;  draw_selection(V);
                end
            case {"box","box_add"}
                hbox.XData = [p0(1) p(1) p(1) p0(1) p0(1)];
                hbox.YData = [p0(2) p0(2) p(2) p(2) p0(2)];
        end
    end

    function mouse_up()
        if drag_mode == "", return; end
        p = cursor();
        switch drag_mode
            case "drag"
                if moved_in_drag
                    undo_stack{end+1} = {pos0, rad};
                    if numel(undo_stack) > 50, undo_stack(1) = []; end
                end
                show_selection();
            case {"box","box_add"}
                x = sort([p0(1) p(1)]);  y = sort([p0(2) p(2)]);
                in = pos(:,1)*1e6 >= x(1) & pos(:,1)*1e6 <= x(2) & ...
                     pos(:,2)*1e6 >= y(1) & pos(:,2)*1e6 <= y(2);
                if drag_mode == "box", sel = in; else, sel = sel | in; end
                hbox.XData = NaN;  hbox.YData = NaN;
                draw_selection();  show_selection();
        end
        drag_mode = "";
    end

    function finish(done)
        if done
            moved = any(abs(pos - base_r) > 1e-15, 2);
            resized = ~isnan(rad) & abs(rad - r_default) > 1e-15;
            rad_out = rad;  rad_out(~resized) = NaN;
            if any(moved) || any(resized)
                edit = struct('enable', true, 'r_i', pos, 'radius', rad_out, ...
                    'base_hash', base_hash, 'n_moved', nnz(moved), 'n_resized', nnz(resized), ...
                    'hash', dda_checksum(pos, rad_out), 'base_inputs', array_inputs(Input));
            else
                edit = [];
            end
            ok = true;
        end
        delete(fig);
    end
end

%% ------------------------------------------------------------------------
function s = array_inputs(Input)
% The inputs that decide the constructed array (for the GUI's "edits outdated" note)
names = {'p','arrayN','array_type','array_shifts_x','array_shifts_y','second_array_offset', ...
         'PC_mirror_thickness','PC_mirror_p','PC_mirror_offset','Outerwallxy','wall_thickness', ...
         'wall_increase_in_radius','Gaussian_sigma','diag_shift','mirror_A','periodic_boundary'};
s = struct();
for n = 1:numel(names)
    if isfield(Input, names{n}), s.(names{n}) = Input.(names{n}); end
end
end

function S = rmfield_if(S, f)
if isfield(S, f), S = rmfield(S, f); end
end
