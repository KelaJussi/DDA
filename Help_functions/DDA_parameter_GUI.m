function [Input, ok] = DDA_parameter_GUI(Input, code_defaults)
% DDA_PARAMETER_GUI  Choose the run options in a tabbed window before a run.
%
%   [Input, ok] = DDA_parameter_GUI(Input)
%   [Input, ok] = DDA_parameter_GUI(Input, code_defaults)
%
% code_defaults: the values written in Make_predefined_Inputs, before any
% saved defaults (DDA_saved_defaults.m) were applied. Buttons:
%   Reset          back to the values the window opened with
%   Code values    back to the values written in Make_predefined_Inputs
%   Save defaults  write the settings that differ from the code values to
%                  DDA_saved_defaults.m (next to DDA_Master.m); Make_predefined_Inputs
%                  applies that file when Input.use_saved_defaults is true
%
% Opened by DDA_Master when Input.open_parameter_GUI is true. The window
% starts from the defaults of Make_predefined_Inputs; "Run" returns the edited
% Input (ok = true), "Cancel" or closing the window returns it unchanged
% (ok = false). "Reset" goes back to the defaults. On "Run" the options that
% differ from the defaults are printed, so the command window records what
% was run.
%
% Options that only matter for a chosen mode (Lorentzian particles, gmres,
% threads/background/cluster, a Hermite or Gaussian field, enabled plots ...)
% are shown only when that mode is chosen. inversion_tool appears on both the
% Physics and the Solver tab; the two stay in sync.
%
% Vector fields accept any MATLAB expression, e.g. [700 700 100] or
% linspace(1.1,1.5,40). Units are shown in the labels.

if nargin < 2 || isempty(code_defaults), code_defaults = Input; end
I0 = Input;                 % starting values, for Reset and for the change report
W  = Input;                 % working copy
ok = false;
opts = option_list();
tab_names = ["Array","Physics","Solver","Tracking","Incident field","Plotting"];

%% Window
fig = uifigure('Name', 'DDA run options', 'Position', [100 80 572 470]);   % width 1.1 x 520
fig.CloseRequestFcn = @(~,~) finish(false);
main = uigridlayout(fig, [2 1]);
main.RowHeight = {'1x', 34};
tg = uitabgroup(main);

ctrl  = cell(1, numel(opts));       % control of each option
lab   = cell(1, numel(opts));       % its label
grids = gobjects(1, numel(tab_names));
rowof = zeros(1, numel(opts));
for t = 1:numel(tab_names)
    tab = uitab(tg, 'Title', tab_names(t));
    k_in = find([opts.tab] == tab_names(t));
    g = uigridlayout(tab, [max(numel(k_in),1) 2]);
    g.ColumnWidth = {'1x', '1x'};
    g.RowHeight = repmat({22}, 1, max(numel(k_in),1));
    g.Scrollable = 'on';
    g.RowSpacing = 6;
    grids(t) = g;
    for r = 1:numel(k_in)
        k = k_in(r);
        rowof(k) = r;
        make_row(k, g, r);
    end
end

%% Buttons
bar = uigridlayout(main, [1 6]);
bar.ColumnWidth = {60, 85, 95, '1x', 65, 65};
bar.Padding = [0 0 0 0];
bar.ColumnSpacing = 6;
uibutton(bar, 'Text', 'Reset', 'ButtonPushedFcn', @(~,~) reset_all(I0), ...
    'Tooltip', 'back to the values the window opened with');
uibutton(bar, 'Text', 'Code values', 'ButtonPushedFcn', @(~,~) reset_all(code_defaults), ...
    'Tooltip', 'back to the values written in Make_predefined_Inputs');
uibutton(bar, 'Text', 'Save defaults', 'ButtonPushedFcn', @(~,~) save_defaults(), ...
    'Tooltip', 'save the settings that differ from the code values to DDA_saved_defaults.m');
uilabel(bar, 'Text', '');
uibutton(bar, 'Text', 'Cancel', 'ButtonPushedFcn', @(~,~) finish(false));
uibutton(bar, 'Text', 'Run', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~,~) finish(true));

%% Favorites tab: save the current settings under a name, load saved ones
ftab = uitab(tg, 'Title', 'Favorites');
fg = uigridlayout(ftab, [8 2]);
fg.ColumnWidth = {110, '1x'};
fg.RowHeight = {48, 22, 22, 26, 22, '1x', 52, 26};
fl = uilabel(fg, 'Text', ['A favorite stores every option of this window, and the edited array ' ...
    'layout. Loading one replaces the values here; Run then uses them.'], 'WordWrap', 'on', ...
    'FontColor', [0.35 0.35 0.35]);
fl.Layout.Row = 1; fl.Layout.Column = [1 2];
fav_lbl = uilabel(fg, 'Text', 'Name');  fav_lbl.Layout.Row = 2; fav_lbl.Layout.Column = 1;
fav_name = uieditfield(fg, 'text', 'Placeholder', 'e.g. 3D 20x20 Roihu');
fav_name.Layout.Row = 2; fav_name.Layout.Column = 2;
fav_lbl = uilabel(fg, 'Text', 'Note');  fav_lbl.Layout.Row = 3; fav_lbl.Layout.Column = 1;
fav_note = uieditfield(fg, 'text', 'Placeholder', 'optional');
fav_note.Layout.Row = 3; fav_note.Layout.Column = 2;
fav_btn = uibutton(fg, 'Text', 'Save current settings as favorite', 'ButtonPushedFcn', @(~,~) fav_save());
fav_btn.Layout.Row = 4; fav_btn.Layout.Column = [1 2];
fav_lbl = uilabel(fg, 'Text', 'Saved favorites', 'FontWeight', 'bold');  fav_lbl.Layout.Row = 5; fav_lbl.Layout.Column = [1 2];
fav_list = uilistbox(fg, 'Items', {}, 'ValueChangedFcn', @(~,~) fav_show());
fav_list.Layout.Row = 6; fav_list.Layout.Column = [1 2];
fav_info = uilabel(fg, 'Text', '', 'WordWrap', 'on');
fav_info.Layout.Row = 7; fav_info.Layout.Column = [1 2];
bl = uibutton(fg, 'Text', 'Load', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~,~) fav_load());
bl.Layout.Row = 8; bl.Layout.Column = 1;
bd = uibutton(fg, 'Text', 'Delete', 'ButtonPushedFcn', @(~,~) fav_delete());
bd.Layout.Row = 8; bd.Layout.Column = 2;
favs = [];
fav_fill();

refresh();
uiwait(fig);
if isvalid(fig), delete(fig); end
if ok
    Input = W;
    report_changes(I0, W, opts);
end

%% ========================================================================
%  Nested helpers (share W, ctrl, lab, ...)
%% ========================================================================
    function make_row(k, g, r)
        o = opts(k);
        switch o.kind
            case {'heading','note','info'}
                l = uilabel(g, 'Text', o.label, 'WordWrap', 'on');
                if o.kind == "heading", l.FontWeight = 'bold'; end
                if o.kind == "note",    l.FontColor = [0.35 0.35 0.35]; end
                l.Layout.Row = r; l.Layout.Column = [1 2];
                lab{k} = l; ctrl{k} = [];
                return
        end
        l = uilabel(g, 'Text', o.label);
        l.Layout.Row = r; l.Layout.Column = 1;
        switch o.kind
            case 'num'
                c = uieditfield(g, 'numeric', 'Limits', o.limits, 'RoundFractionalValues', o.integer, ...
                    'HorizontalAlignment', 'left');
            case {'vec','text'}
                c = uieditfield(g, 'text', 'Placeholder', o.placeholder);
            case 'drop'
                c = uidropdown(g, 'Items', cellstr(o.items));
            case 'check'
                c = uicheckbox(g, 'Text', '');
            case 'list'
                c = uilistbox(g, 'Items', cellstr(o.items), 'Multiselect', 'on');
            case 'button'
                l.Text = '';
                c = uibutton(g, 'Text', o.label, 'ButtonPushedFcn', @(~,~) on_button(k));
        end
        c.Layout.Row = r; c.Layout.Column = 2;
        if o.kind ~= "button", c.ValueChangedFcn = @(src,~) on_change(k, src); end
        if strlength(o.tip) > 0, l.Tooltip = o.tip; c.Tooltip = o.tip; end
        lab{k} = l; ctrl{k} = c;
        show_value(k);
    end

    function show_value(k)
        % Put the current value of option k (from W) into its control
        o = opts(k); c = ctrl{k};
        if isempty(c) || o.kind == "button", return; end
        v = o.get(W);
        switch o.kind
            case 'num'
                v = double(v);
                if isempty(v) || isnan(v(1)), v = max(0, o.limits(1)); end   % numeric fields cannot be empty
                c.Value = min(max(v(1), o.limits(1)), o.limits(2));
            case 'vec',   c.Value = vec2str(v);
            case 'text',  c.Value = char(string(v));
            case 'check', c.Value = logical(v);
            case 'drop'
                v = char(string(v));
                if ~ismember(v, c.Items), c.Items = [c.Items, {v}]; end   % keep unusual values
                c.Value = v;
            case 'list'
                v = cellstr(string(v));
                c.Value = v(ismember(v, c.Items));
        end
    end

    function on_change(k, src)
        o = opts(k);
        try
            switch o.kind
                case 'vec'
                    v = parse_vec(src.Value, o);
                case 'text'
                    v = string(src.Value);
                case 'list'
                    v = string(src.Value);
                otherwise
                    v = src.Value;
            end
            W = o.set(W, v);
        catch err
            uialert(fig, err.message, 'Invalid value');
            show_value(k);                              % back to the last valid value
            return
        end
        for j = find([opts.key] == o.key & (1:numel(opts)) ~= k)
            show_value(j);                              % e.g. inversion_tool on two tabs
        end
        refresh();
    end

    function on_button(k)
        % A button changes several options at once: show them all again.
        % Actions take (I) or (I, starting values), e.g. the array editor.
        a = opts(k).action;
        if nargin(a) >= 2, W = a(W, I0); else, W = a(W); end
        for j = 1:numel(opts), show_value(j); end
        refresh();
    end

    function refresh()
        % Show the rows whose condition holds, collapse the others
        for t = 1:numel(tab_names)
            k_in = find([opts.tab] == tab_names(t));
            h = grids(t).RowHeight;
            for k = k_in
                o = opts(k);
                vis = o.visible(W);
                if vis
                    h{rowof(k)} = o.height;
                    if o.kind == "info", lab{k}.Text = o.text(W); end
                    if any(o.kind == ["note","info"])
                        % about 65 characters per line at this window width
                        h{rowof(k)} = max(22, 16*ceil(strlength(string(lab{k}.Text))/65) + 6);
                    end
                else
                    h{rowof(k)} = 0;
                end
                lab{k}.Visible = vis;
                if ~isempty(ctrl{k}), ctrl{k}.Visible = vis; end
            end
            grids(t).RowHeight = h;
        end
    end

    %% Favorites --------------------------------------------------------
    function fav_fill(select)
        favs = dda_favorites('list');
        if isempty(favs)
            fav_list.Items = {};  fav_info.Text = 'No favorites saved yet.';
            return
        end
        fav_list.Items = cellstr([favs.name]);
        if nargin > 0 && any([favs.name] == select), fav_list.Value = char(select); end
        fav_show();
    end

    function fav_show()
        if isempty(favs) || isempty(fav_list.Value), fav_info.Text = ''; return; end
        f = favs([favs.name] == string(fav_list.Value));
        txt = sprintf('Saved %s, %d settings%s.', char(f.saved, 'yyyy-MM-dd HH:mm'), f.n, ...
            ternary_s(f.layout, ', with an edited array layout', ''));
        if strlength(f.note) > 0, txt = sprintf('%s\n%s', txt, f.note); end
        fav_info.Text = txt;
    end

    function fav_save()
        name = strtrim(string(fav_name.Value));
        if strlength(name) == 0
            uialert(fig, 'Give the favorite a name first.', 'Save favorite'); return
        end
        if ~isempty(favs) && any(lower([favs.name]) == lower(name))
            a = uiconfirm(fig, sprintf('Replace the favorite "%s"?', name), 'Save favorite', ...
                'Options', {'Replace', 'Cancel'}, 'DefaultOption', 2, 'CancelOption', 2);
            if a ~= "Replace", return; end
        end
        paths = option_paths(opts);
        keep = false(size(paths));  vals = cell(size(paths));
        for q = 1:numel(paths)
            [vals{q}, keep(q)] = getp_has(W, paths(q));
        end
        ae = [];
        if isfield(W, 'array_edit') && ~isempty(W.array_edit) && W.array_edit.enable, ae = W.array_edit; end
        try
            dda_favorites('save', name, fav_note.Value, paths(keep), vals(keep), ae);
            fav_fill(name);
            uialert(fig, sprintf('Saved "%s" (%d settings%s).', name, nnz(keep), ...
                ternary_s(isempty(ae), '', ' and the array layout')), 'Save favorite', 'Icon', 'success');
        catch err
            uialert(fig, err.message, 'Could not save the favorite');
        end
    end

    function fav_load()
        if isempty(favs) || isempty(fav_list.Value), return; end
        try
            W = dda_favorites('apply', W, string(fav_list.Value));
        catch err
            uialert(fig, err.message, 'Could not load the favorite'); return
        end
        for k = 1:numel(opts), show_value(k); end
        refresh();
        fav_name.Value = fav_list.Value;
        uialert(fig, sprintf('Loaded "%s". Check the tabs, then press Run.', fav_list.Value), ...
            'Favorite loaded', 'Icon', 'success');
    end

    function fav_delete()
        if isempty(favs) || isempty(fav_list.Value), return; end
        a = uiconfirm(fig, sprintf('Delete the favorite "%s"?', fav_list.Value), 'Delete favorite', ...
            'Options', {'Delete', 'Cancel'}, 'DefaultOption', 2, 'CancelOption', 2);
        if a ~= "Delete", return; end
        dda_favorites('delete', string(fav_list.Value));
        fav_fill();
    end

    function reset_all(values)
        W = values;
        for k = 1:numel(opts), show_value(k); end
        refresh();
    end

    function save_defaults()
        file = saved_defaults_file();
        layout_file = fullfile(fileparts(file), 'DDA_saved_layout.mat');
        [paths, scales] = changed_paths(code_defaults, W, opts);
        has_layout = isfield(W, 'array_edit') && ~isempty(W.array_edit) && W.array_edit.enable;
        if isempty(paths) && ~has_layout
            if isfile(file)
                answer = uiconfirm(fig, ['All settings equal the code values. Delete ' file ...
                    ' so that the code values are used?'], 'Save defaults', ...
                    'Options', {'Delete', 'Cancel'}, 'DefaultOption', 2, 'CancelOption', 2);
                if answer == "Delete"
                    delete(file); clear('DDA_saved_defaults');
                    if isfile(layout_file), delete(layout_file); end
                end
            else
                uialert(fig, 'All settings equal the code values; nothing to save.', 'Save defaults', 'Icon', 'info');
            end
            return
        end
        listed = paths;
        if has_layout
            listed(end+1) = sprintf('array layout (%d moved, %d resized)', ...
                W.array_edit.n_moved, W.array_edit.n_resized);
        end
        answer = uiconfirm(fig, sprintf('Save these %d settings as defaults?\n\n%s\n\nto %s', ...
            numel(listed), strjoin(listed, ', '), file), 'Save defaults', ...
            'Options', {'Save', 'Cancel'}, 'DefaultOption', 1, 'CancelOption', 2);
        if answer ~= "Save", return; end
        try
            if has_layout
                array_edit = W.array_edit; %#ok<NASGU>
                save(layout_file, 'array_edit');
            elseif isfile(layout_file)
                delete(layout_file);                     % saved defaults without a layout
            end
            write_saved_defaults(file, paths, scales, W, has_layout);
            clear('DDA_saved_defaults');                % use the new file next time
            uialert(fig, sprintf('Saved %d settings to %s', numel(listed), file), 'Save defaults', 'Icon', 'success');
        catch err
            uialert(fig, err.message, 'Could not save defaults');
        end
    end

    function finish(run)
        ok = run;
        uiresume(fig);
        if ~run, delete(fig); end
    end
end

%% ========================================================================
%  The options
%% ========================================================================
function O = option_list()
O = opt('heading', 'Array', 'Lattice');
O(end+1) = opt('vec', 'Array', 'Period p = [x y z] (nm)', 'path', "p", 'scale', 1e-9, 'numel', 3);
O(end+1) = opt('vec', 'Array', 'Particles arrayN = [nx ny nz]', 'path', "arrayN", 'numel', 3, 'integer', true);
O(end+1) = opt('drop', 'Array', 'Lattice type (array_type)', 'path', "array_type", ...
    'items', ["1Dsquare_halfshifts","1Dsquare_halfshifts2","square","double_rectangular","PC_mirrors", ...
              "PC_mirror_reflection_test","Outerwall","Gaussian","Inverted_Gaussian","shifted_diagonals","hexagonal"], ...
    'tip', 'PC_mirror_reflection_test also sets the incident field to that test');
O(end+1) = opt('num', 'Array', 'Shift x (array_shifts_x)', 'path', "array_shifts_x", ...
    'visible', @(I) atype(I, ["1Dsquare_halfshifts","1Dsquare_halfshifts2","double_rectangular", ...
        "PC_mirrors","PC_mirror_reflection_test","PC_middle","Outerwall","HG_wall","Gaussian","Inverted_Gaussian"]));
O(end+1) = opt('num', 'Array', 'Shift y (array_shifts_y)', 'path', "array_shifts_y", ...
    'visible', @(I) atype(I, ["1Dsquare_halfshifts","1Dsquare_halfshifts2","double_rectangular", ...
        "PC_mirrors","PC_mirror_reflection_test","PC_middle","Outerwall","HG_wall","Gaussian","Inverted_Gaussian"]));
O(end+1) = opt('vec', 'Array', 'Second lattice offset [x y z] (nm)', 'path', "second_array_offset", ...
    'scale', 1e-9, 'numel', 3, 'visible', @(I) atype(I, "double_rectangular"));
O(end+1) = opt('num', 'Array', 'PC mirror thickness (particles)', 'path', "PC_mirror_thickness", ...
    'integer', true, 'limits', [0 Inf], 'visible', @is_pc);
O(end+1) = opt('vec', 'Array', 'PC mirror period [x y] (nm)', 'path', "PC_mirror_p", 'scale', 1e-9, ...
    'numel', 2, 'visible', @is_pc);
O(end+1) = opt('vec', 'Array', 'PC mirror offset [x y] (nm)', 'path', "PC_mirror_offset", 'scale', 1e-9, ...
    'numel', 2, 'visible', @is_pc);
O(end+1) = opt('num', 'Array', 'Wall on (Outerwallxy)', 'path', "Outerwallxy", 'integer', true, ...
    'limits', [0 2], 'tip', '0: x and y, 1: only x, 2: only y', 'visible', @(I) atype(I, "Outerwall"));
O(end+1) = opt('num', 'Array', 'Wall thickness (particles)', 'path', "wall_thickness", 'integer', true, ...
    'limits', [0 Inf], 'visible', @(I) atype(I, "Outerwall"));
O(end+1) = opt('num', 'Array', 'Radius increase (nm)', 'path', "wall_increase_in_radius", 'scale', 1e-9, ...
    'tip', 'Outerwall: added at the walls; Gaussian types: at the centre/edges', ...
    'visible', @(I) atype(I, ["Outerwall","Gaussian","Inverted_Gaussian"]));
O(end+1) = opt('num', 'Array', 'Gaussian width sigma (um)', 'path', "Gaussian_sigma", 'scale', 1e-6, ...
    'limits', [0 Inf], 'visible', @(I) atype(I, ["Gaussian","Inverted_Gaussian"]));
O(end+1) = opt('num', 'Array', 'Diagonal shift (nm)', 'path', "diag_shift", 'scale', 1e-9, ...
    'visible', @(I) atype(I, "shifted_diagonals"));
O(end+1) = opt('heading', 'Array', 'Layout');
O(end+1) = opt('button', 'Array', 'Edit array layout...', 'action', @open_array_editor, ...
    'tip', 'show the array the run will construct; move particles, change radii and z');
O(end+1) = opt('info', 'Array', '', 'text', @layout_info, 'visible', @has_layout_edit);
O(end+1) = opt('button', 'Array', 'Clear layout edits', 'action', @clear_layout_edit, ...
    'visible', @has_layout_edit);
O(end+1) = opt('heading', 'Array', 'Particles');
O(end+1) = opt('drop', 'Array', 'Substrate / host', 'path', "substrate", 'items', ["BK7","const"], ...
    'tip', 'BK7: dispersive glass; const: constant index n');
O(end+1) = opt('num', 'Array', 'Host refractive index n', 'path', "n", 'limits', [1 Inf], ...
    'visible', @(I) string(I.substrate) == "const");
O(end+1) = opt('drop', 'Array', 'Particle model (sample.type)', 'path', "sample.type", ...
    'items', ["Au_Drude_Mie","Au_Mie","glass_Mie","Lorentzian"], ...
    'tip', 'Mie theory from the radius, or MJH Lorentzian dipoles');
O(end+1) = opt('num', 'Array', 'Radius particle.r (nm)', 'path', "particle.r", 'scale', 1e-9, ...
    'limits', [0 Inf], 'visible', @(I) ~is_lorentz(I));
O(end+1) = opt('num', 'Array', 'Glass index (Mie_glass_n)', 'path', "Mie_glass_n", ...
    'limits', [1 Inf], 'visible', @(I) string(I.sample.type) == "glass_Mie");
O(end+1) = opt('heading', 'Array', 'MJH Lorentzian particles', 'visible', @is_lorentz);
O(end+1) = opt('note', 'Array', ['alpha(w) = A0/(w0 - w - i*gamma) per axis. Each value: a scalar, ' ...
    '[x y z], or N x 3 (one row per particle).'], 'visible', @is_lorentz, 'height', 34);
O(end+1) = opt('vec', 'Array', 'Strength A0 (scatter_strength)', 'path', "particle.scatter_strength", ...
    'visible', @is_lorentz);
O(end+1) = opt('vec', 'Array', 'Lifetime tau (fs)', 'path', "particle.tau", 'scale', 1e-15, ...
    'visible', @is_lorentz);
O(end+1) = opt('vec', 'Array', 'Resonance wavelength (nm)', 'path', "particle.lambda_peak_vec", ...
    'scale', 1e-9, 'visible', @is_lorentz);
O(end+1) = opt('drop', 'Array', 'Units of A0', 'path', "particle.lorentz_units", ...
    'items', ["gaussian","SI"], 'default', "gaussian", 'visible', @is_lorentz, ...
    'tip', 'gaussian: values as in the old MJH code (converted to SI); SI: A0 already in C m^2/V');
O(end+1) = opt('drop', 'Array', 'Correction', 'path', "particle.lorentz_correction", ...
    'items', ["none","radiative","MLWA"], 'default', "none", 'visible', @is_lorentz, ...
    'tip', 'none: plain Lorentzian (as before); radiative: energy conserving; MLWA: + depolarization');
O(end+1) = opt('num', 'Array', 'MLWA particle size a (nm)', 'path', "particle.lorentz_mlwa_a", ...
    'scale', 1e-9, 'default', 100e-9, 'limits', [0 Inf], ...
    'visible', @(I) is_lorentz(I) && strcmpi(string(getp(I,"particle.lorentz_correction","none")), "MLWA"));

% ---------------------------------------------------------------- Physics
O(end+1) = opt('heading', 'Physics', 'Dipoles and boundaries');
O(end+1) = opt('drop', 'Physics', 'Dimensionality', 'path', "dimensionality", ...
    'items', ["3D","2D","1Dx","1Dy","1Dz","cyl_TM","cyl_TE"], ...
    'tip', ['cyl_TM / cyl_TE: infinite cylinders along y, array along x (arrayN = [Nx 1 1]); ' ...
            'TM: field along the cylinders, TE: field in the xz plane (2 unknowns per cylinder)']);
O(end+1) = opt('heading', 'Physics', 'Layered background (cylinders)', 'visible', @is_cyl);
O(end+1) = opt('check', 'Physics', 'Use layers (layers.enable)', 'path', "layers.enable", ...
    'default', false, 'visible', @is_cyl, ...
    'tip', 'the cylinders sit in a host layer (the host index of the Array tab) between stacks of layers');
O(end+1) = opt('note', 'Physics', ['Layers are stacked along z. Indices are listed from the host layer ' ...
    'outwards; the last one is the half-space, so give one thickness less than indices.'], ...
    'visible', @layers_on, 'height', 34);
O(end+1) = opt('num', 'Physics', 'Host layer thickness (nm)', 'path', "layers.host_thickness", ...
    'scale', 1e-9, 'default', 1e-6, 'limits', [0 Inf], 'visible', @layers_on);
O(end+1) = opt('num', 'Physics', 'Cylinder height z_cyl (nm)', 'path', "layers.z_cyl", ...
    'scale', 1e-9, 'default', 100e-9, 'limits', [0 Inf], 'visible', @layers_on, ...
    'tip', 'height of the cylinder axes above the lower interface of the host layer');
O(end+1) = opt('vec', 'Physics', 'Indices above (nearest first)', 'path', "layers.above_n", ...
    'default', 1, 'visible', @layers_on, 'tip', 'e.g. 1 for air, or [1.5 1] for a 1.5 layer and then air');
O(end+1) = opt('vec', 'Physics', 'Thicknesses above (nm)', 'path', "layers.above_d", 'scale', 1e-9, ...
    'default', [], 'allow_empty', true, 'placeholder', 'one less than the indices', 'visible', @layers_on);
O(end+1) = opt('vec', 'Physics', 'Indices below (nearest first)', 'path', "layers.below_n", ...
    'default', 1, 'visible', @layers_on, 'tip', 'e.g. [1.33 1.52]: a 1.33 layer on a 1.52 substrate');
O(end+1) = opt('vec', 'Physics', 'Thicknesses below (nm)', 'path', "layers.below_d", 'scale', 1e-9, ...
    'default', [], 'allow_empty', true, 'placeholder', 'one less than the indices', 'visible', @layers_on);
O(end+1) = opt('check', 'Physics', 'Benchmark the layer cost', 'path', "layers.benchmark", ...
    'default', true, 'visible', @layers_on, ...
    'tip', ['time the layered Green''s function for this array on this computer for the time estimate ' ...
            '(cached; up to a few seconds after a change of the array, layers or frequency range)']);
O(end+1) = opt('check', 'Physics', 'Old cyl_TM polarizability', 'path', "particle.cyl_alpha_legacy", ...
    'default', false, 'visible', @(I) string(I.dimensionality) == "cyl_TM", ...
    'tip', 'the expression used before the fix (too small); only to reproduce old results');
O(end+1) = opt('vec', 'Physics', 'Mirrors mirror_A = [x0 xL y0 yL]', 'path', "mirror_A", 'numel', 4, ...
    'tip', ['0 none; 2D/3D: +1 field mirror-symmetric (PMC-like), -1 antisymmetric (PEC-like); ' ...
            '1D: p(-x) = +-p(x). Not together with a periodic boundary in the same direction']);
O(end+1) = opt('check', 'Physics', 'Periodic in x', 'key', "periodic_x", 'save_path', "periodic_boundary", ...
    'get', @(I) logical(getp(I,"periodic_boundary",[0 0])*[1;0]), ...
    'set', @(I,v) set_periodic(I, 1, v));
O(end+1) = opt('check', 'Physics', 'Periodic in y', 'key', "periodic_y", 'save_path', "periodic_boundary", ...
    'get', @(I) logical(getp(I,"periodic_boundary",[0 0])*[0;1]), ...
    'set', @(I,v) set_periodic(I, 2, v));
O(end+1) = opt('vec', 'Physics', 'Periodic image cells [x y]', 'path', "periodic_Ncells", ...
    'default', [10 10], 'numel', 2, 'integer', true, 'visible', @any_periodic, ...
    'tip', 'Image cells on EACH side of the central one: 0 = none, 1 = 3 cells in all, 10 = 21 cells');
O(end+1) = opt('vec', 'Physics', 'Periodic cell size [x y] (nm)', 'path', "periodic_cell", ...
    'scale', 1e-9, 'default', [], 'allow_empty', true, 'placeholder', 'empty = from the positions', ...
    'visible', @any_periodic);
O(end+1) = opt('heading', 'Physics', 'Frequency range');
O(end+1) = range_opt('Real part from (eV)',  "real_parts", 1, 1);
O(end+1) = range_opt('Real part to (eV)',    "real_parts", 2, 1);
O(end+1) = range_opt('Real part points',     "real_parts", 3, 1);
O(end+1) = range_opt('Imag. part from (meV)', "imag_parts", 1, 1e-3);
O(end+1) = range_opt('Imag. part to (meV)',   "imag_parts", 2, 1e-3);
O(end+1) = range_opt('Imag. part points',     "imag_parts", 3, 1e-3);
O(end+1) = opt('button', 'Physics', 'Make the run real (Im = 0, 1 point)', ...
    'action', @make_real, 'visible', @is_complex, ...
    'tip', 'sets range.imag_parts = 0: a real-frequency run');
O(end+1) = opt('info', 'Physics', '', 'text', @range_info);
O(end+1) = opt('heading', 'Physics', 'Solver');
O(end+1) = inversion_opt('Physics');
O(end+1) = opt('num', 'Physics', 'Gain (direct solver only)', 'path', "gain", 'visible', @(I) ~is_eig(I), ...
    'tip', 'Old gain parameter of the direct solver; prefer complex frequencies');

% ---------------------------------------------------------------- Solver
O(end+1) = opt('heading', 'Solver', 'Solver');
O(end+1) = inversion_opt('Solver');
O(end+1) = opt('num', 'Solver', 'gmres restart', 'path', "gmres.restart", 'integer', true, ...
    'limits', [1 Inf], 'visible', @is_gmres);
O(end+1) = opt('num', 'Solver', 'gmres max. iterations', 'path', "gmres.maxit", 'integer', true, ...
    'limits', [1 Inf], 'visible', @is_gmres);
O(end+1) = opt('num', 'Solver', 'gmres tolerance', 'path', "gmres.tol", 'limits', [0 Inf], ...
    'visible', @is_gmres);
O(end+1) = opt('check', 'Solver', 'Nearest neighbours only', 'path', "nearest_neighbour", ...
    'default', false, 'visible', @(I) contains(string(I.inversion_tool), "fft_gmres"), ...
    'tip', 'fft_gmres: only nearest-neighbour interactions');
O(end+1) = opt('num', 'Solver', 'Degeneracy tolerance', 'path', "degeneracy_tol", ...
    'default', 1e29, 'limits', [0 Inf], 'visible', @is_eig, 'tip', '|lambda - lambda_min| below this counts as degenerate');
O(end+1) = opt('drop', 'Solver', 'Driven response P', 'path', "P_method", ...
    'items', ["direct","eigenmodes"], 'default', "direct", 'visible', @(I) is_eig(I) && ~is_complex(I), ...
    'tip', 'direct: P = A\E (exact); eigenmodes: P = V*(c./lambda), the sum over eigenmodes');
O(end+1) = opt('check', 'Solver', 'Legacy eigen solvers (validation)', 'path', "use_legacy_eig_solvers", ...
    'default', false, 'visible', @is_eig, 'tip', 'old per-dimension files; only with parallel none/threads');
O(end+1) = opt('check', 'Solver', 'Recompute existing results', 'path', "recompute", ...
    'tip', 'solve again even if a result file with these inputs exists');
O(end+1) = opt('heading', 'Solver', 'Where the solve runs');
O(end+1) = opt('drop', 'Solver', 'Parallel', 'path', "parallel", ...
    'items', ["none","threads","background","cluster"], ...
    'tip', ['none/threads: in this MATLAB (busy until done); background: local pool, MATLAB free; ' ...
            'cluster: batch job (or "Processes" = this computer), MATLAB free']);
O(end+1) = opt('info', 'Solver', '', 'text', @time_estimate, 'height', 52);
O(end+1) = opt('info', 'Solver', '', 'text', @memory_info, 'visible', @is_eig);
O(end+1) = opt('vec', 'Solver', 'Thread pool workers', 'path', "threads.num_workers", 'default', [], ...
    'allow_empty', true, 'integer', true, 'placeholder', 'empty = all cores', 'visible', @(I) par(I,"threads"));
O(end+1) = opt('check', 'Solver', 'Keep the pool open between runs', 'path', "threads.keep_pool", ...
    'default', true, 'visible', @(I) par(I,"threads"));
O(end+1) = opt('drop', 'Solver', 'Background pool', 'path', "background.pool", ...
    'items', ["Processes","Threads"], 'default', "Processes", 'visible', @(I) par(I,"background"));
O(end+1) = opt('vec', 'Solver', 'Background pool workers', 'path', "background.num_workers", ...
    'default', [], 'allow_empty', true, 'integer', true, 'placeholder', 'empty = default size', ...
    'visible', @(I) par(I,"background"));
O(end+1) = opt('text', 'Solver', 'Cluster profile (cluster.name)', 'path', "cluster.name", ...
    'visible', @(I) par(I,"cluster"), 'tip', '"Processes" runs the job on this computer');
O(end+1) = opt('check', 'Solver', 'Wait for the job (block MATLAB)', 'path', "cluster.wait", ...
    'default', false, 'visible', @(I) par(I,"cluster"));
O(end+1) = opt('text', 'Solver', 'Wall time (h:mm:ss)', 'path', "cluster.wall_time", ...
    'visible', @(I) par(I,"cluster") && ~is_local_profile(I));
O(end+1) = opt('info', 'Solver', '', 'text', @(I) "(1 thread, per freq/worker) 1Dx20x20 0.09s, 2Dx20x20 0.59s, 3Dx20x20 2.8s, 2Dx30x30 10.5s, 2Dx40x40 57s; + 0.7s per worker, + 6s start-up", ...
    'visible', @(I) par(I,"cluster") && ~is_local_profile(I));
O(end+1) = opt('info', 'Solver', '', 'text', @walltime_check, ...
    'visible', @(I) par(I,"cluster") && ~is_local_profile(I));
O(end+1) = opt('num', 'Solver', 'Workers', 'path', "cluster.num_workers", 'integer', true, ...
    'limits', [1 Inf], 'visible', @(I) par(I,"cluster"));
O(end+1) = opt('num', 'Solver', 'Threads per worker', 'path', "cluster.num_threads", 'integer', true, ...
    'limits', [1 Inf], 'visible', @(I) par(I,"cluster"));
O(end+1) = opt('text', 'Solver', 'CSC project (cluster.project)', 'path', "cluster.project", ...
    'placeholder', 'project_2001234', 'visible', @(I) par(I,"cluster") && ~is_local_profile(I));
O(end+1) = opt('text', 'Solver', 'Memory per CPU (e.g. 4g)', 'path', "cluster.mem_per_cpu", ...
    'placeholder', 'empty = profile default', 'visible', @(I) par(I,"cluster") && ~is_local_profile(I));
O(end+1) = opt('info', 'Solver', '', 'text', @cores_info, 'visible', @(I) par(I,"cluster"));
O(end+1) = opt('button', 'Solver', 'Use a full node', 'action', @use_full_node, ...
    'visible', @(I) par(I,"cluster") && ~is_local_profile(I), ...
    'tip', 'partition medium (whole node), workers to fill 384 cores with the chosen threads, 1900 MB per core (fits a node)');
O(end+1) = opt('text', 'Solver', 'Partition', 'path', "cluster.partition", 'placeholder', 'empty = profile default', ...
    'visible', @(I) par(I,"cluster") && ~is_local_profile(I));
O(end+1) = opt('text', 'Solver', 'E-mail for Slurm notifications', 'path', "cluster.email", ...
    'visible', @(I) par(I,"cluster") && ~is_local_profile(I));
O(end+1) = opt('text', 'Solver', 'Extra sbatch arguments', 'path', "cluster.submit_args", ...
    'visible', @(I) par(I,"cluster") && ~is_local_profile(I));
O(end+1) = opt('check', 'Solver', 'Delete the job after collecting', 'path', "cluster.delete_job", ...
    'default', true, 'visible', @(I) par(I,"cluster"));
O(end+1) = opt('check', 'Solver', 'Send the code with the job', 'path', "cluster.auto_attach_files", ...
    'default', true, 'visible', @(I) par(I,"cluster"));
O(end+1) = opt('check', 'Solver', 'Save these settings in the profile', 'path', "cluster.save_profile", ...
    'default', false, 'visible', @(I) par(I,"cluster"));
O(end+1) = opt('check', 'Solver', 'Plot when collected', ...
    'path', "jobs.plot_on_collect", 'default', false, ...
    'visible', @(I) par(I,"cluster") || par(I,"background"), ...
    'tip', 'plot submitted runs when collect_dda_jobs saves them');
O(end+1) = opt('heading', 'Solver', 'Saving');
O(end+1) = opt('text', 'Solver', 'File name note (save_message)', 'path', "save_message", ...
    'default', "", 'placeholder', 'optional');
O(end+1) = opt('check', 'Solver', 'Save mode dipoles (P_eigs)', 'path', "save.P_eigs", ...
    'default', true, 'visible', @(I) is_eig(I) && ~is_complex(I));
O(end+1) = opt('check', 'Solver', 'Keep all eigenvalues (complex)', 'path', "save.complex_all", ...
    'default', false, 'visible', @(I) is_eig(I) && is_complex(I), ...
    'tip', 'off: keep only what the complex plots use (about 8x smaller results)');
O(end+1) = opt('num', 'Solver', 'Drop P_eigs above file size (GB)', 'path', "save.max_file_GB", ...
    'default', 4, 'limits', [0 Inf], 'visible', @(I) is_eig(I) && ~is_complex(I));

% ---------------------------------------------------------------- Tracking
O(end+1) = opt('note', 'Tracking', ['Tracking concerns real-frequency runs: the labelled modes feed the ' ...
    'real-frequency plots. Complex-frequency runs plot the eigenvalues and need no tracking.'], 'height', 34);
O(end+1) = opt('note', 'Tracking', 'Mode tracking is used only by the eigen solver (inversion_tool = eig).', ...
    'visible', @(I) ~is_eig(I));
O(end+1) = opt('drop', 'Tracking', 'Tracking (track.enable)', 'key', "track.enable", 'save_path', "track.enable", ...
    'items', ["auto","on","off"], 'visible', @is_eig, ...
    'get', @track_enable_get, 'set', @track_enable_set, ...
    'tip', 'auto: real-frequency runs only');
O(end+1) = opt('drop', 'Tracking', 'Method', 'path', "track.method", ...
    'items', ["basis","overlap","hybrid"], 'default', "hybrid", 'visible', @tracking_on, ...
    'tip', 'basis: label by spatial pattern at each frequency (band diagrams); overlap/hybrid: follow eigenvectors');
O(end+1) = opt('num', 'Tracking', 'Overlap weight', 'path', "track.overlap_weight", 'default', 0.5, ...
    'limits', [0 1], 'visible', @(I) tracking_on(I) && string(getp(I,"track.method","hybrid")) == "hybrid");
O(end+1) = opt('drop', 'Tracking', 'Basis', 'path', "track.basis", 'items', ["auto","sine","fourier","k"], ...
    'default', "auto", 'visible', @tracking_on, ...
    'tip', 'auto: Fourier along periodic directions, sine otherwise; k: the plane wave that couples best (as the complex "k" labels)');
O(end+1) = opt('check', 'Tracking', 'k labels: undo the pi steps of shifted lattices', 'path', "track.k_step", ...
    'default', false, 'visible', @(I) tracking_on(I) && string(getp(I,"track.basis","auto")) == "k", ...
    'tip', 'label by the step-matched (HG-shaped) field instead of a plane wave');
O(end+1) = opt('check', 'Tracking', 'Keep mode dipoles (P_eigs)', 'path', "track.store_P_eigs", ...
    'default', true, 'visible', @tracking_on);
O(end+1) = opt('num', 'Tracking', 'Minimum confidence', 'path', "track.min_confidence", 'default', 0.5, ...
    'limits', [0 1], 'visible', @tracking_on);
O(end+1) = opt('num', 'Tracking', 'Memory limit for eigenvectors (GB)', 'path', "track.max_memory_GB", ...
    'default', 8, 'limits', [0 Inf], 'visible', @tracking_on);
O(end+1) = opt('check', 'Tracking', 'Print a summary', 'path', "track.verbose", 'default', true, ...
    'visible', @tracking_on);

% ---------------------------------------------------------------- Incident field
O(end+1) = opt('note', 'Incident field', ['The incident field concerns real-frequency runs (driven response, ' ...
    'extinction, plots). Complex-frequency eigenmodes do not depend on it.'], 'height', 34);
O(end+1) = opt('heading', 'Incident field', 'Plane wave');
O(end+1) = opt('drop', 'Incident field', 'Polarization', 'path', "polarization", 'items', ["linear","RC","LC"], ...
    'tip', 'RC/LC (circular) only for normal incidence');
O(end+1) = opt('num', 'Incident field', 'theta (deg)', 'path', "theta", 'tip', 'incident angle, 0 = normal incidence');
O(end+1) = opt('num', 'Incident field', 'phi (deg)', 'path', "phi", 'tip', 'angle in the xy-plane, 0 = along x');
O(end+1) = opt('num', 'Incident field', 'psi (deg)', 'path', "psi", 'tip', 'polarization rotation; (0,0,90) = x-polarized along z');
O(end+1) = opt('heading', 'Incident field', 'Field profile');
O(end+1) = opt('drop', 'Incident field', 'Field profile', 'path', "preset_field", ...
    'items', ["none","Gaussian"], 'default', "none", ...
    'tip', 'none: plane wave; Gaussian: Hermite-Gaussian beam centred on the array');
O(end+1) = opt('num', 'Incident field', 'Hermite order x', 'path', "Hermite_order_x", 'integer', true, ...
    'limits', [0 Inf], 'visible', @(I) preset(I, "Gaussian"));
O(end+1) = opt('num', 'Incident field', 'Hermite order y', 'path', "Hermite_order_y", 'integer', true, ...
    'limits', [0 Inf], 'visible', @(I) preset(I, "Gaussian") && ~is_cyl(I), ...
    'tip', 'e.g. x 0, y 1: HG01, a node along y through the array centre');
O(end+1) = opt('num', 'Incident field', 'Beam waist w0 (um)', 'path', "Gaussian_w0", 'scale', 1e-6, ...
    'limits', [0 Inf], 'visible', @(I) preset(I, "Gaussian"), ...
    'tip', 'waist radius (1/e^2 of the intensity); intensity FWHM = w0*sqrt(2 ln 2)');
O(end+1) = opt('num', 'Incident field', 'Waist distance z (um)', 'path', "Gaussian_z", 'scale', 1e-6, ...
    'visible', @(I) preset(I, "Gaussian"), ...
    'tip', 'distance from the array to the waist; > 0: waist beyond the array (converging beam), 0: waist on the array');

% ---------------------------------------------------------------- Plotting
O(end+1) = opt('note', 'Plotting', ['Real-frequency runs: real_freqs_plots_from_file. Complex-frequency ' ...
    'runs: complex_freqs_plots_from_file. The options are saved with the results; a file can be ' ...
    're-plotted with other options later.'], 'height', 48);
O(end+1) = opt('check', 'Plotting', 'Plot after the run', 'path', "plotting.enable");
O(end+1) = opt('check', 'Plotting', 'Figures as tabs in one window', 'path', "plotting.docked", ...
    'default', true, 'visible', @plots_on, 'tip', 'off: each figure in its own window');
O(end+1) = opt('text', 'Plotting', 'Export folder for PNGs', 'path', "plotting.export_folder", ...
    'default', "", 'placeholder', 'empty = no export', 'visible', @plots_real);
O(end+1) = opt('heading', 'Plotting', 'Dipole moments', 'visible', @plots_real);
O(end+1) = opt('check', 'Plotting', 'Show', 'path', "plotting.dipoles.enable", 'visible', @plots_real);
O(end+1) = opt('vec', 'Plotting', 'Eigenmodes [mx my; ...] (0-based)', 'path', "plotting.dipoles.modes", ...
    'allow_empty', true, 'visible', @(I) sub_on(I,"dipoles"));
O(end+1) = opt('num', 'Plotting', '+ most excited modes', 'path', "plotting.dipoles.n_strongest", ...
    'integer', true, 'limits', [0 Inf], 'visible', @(I) sub_on(I,"dipoles"));
O(end+1) = opt('vec', 'Plotting', 'Profile wavelengths (nm)', 'path', "plotting.dipoles.profile_wls", ...
    'allow_empty', true, 'visible', @(I) sub_on(I,"dipoles"));
O(end+1) = opt('vec', 'Plotting', 'Map wavelength (nm)', 'path', "plotting.dipoles.map_wl", ...
    'allow_empty', true, 'placeholder', 'empty = first wavelength', 'visible', @(I) sub_on(I,"dipoles"));
O(end+1) = opt('vec', 'Plotting', 'Surface wavelengths (nm)', 'path', "plotting.dipoles.surf_wls", ...
    'allow_empty', true, 'placeholder', 'empty = none', 'visible', @(I) sub_on(I,"dipoles"));
O(end+1) = opt('heading', 'Plotting', 'Q factor', 'visible', @plots_real);
O(end+1) = opt('check', 'Plotting', 'Show', 'path', "plotting.Q.enable", 'visible', @plots_real);
O(end+1) = opt('vec', 'Plotting', 'Wavelength (nm)', 'path', "plotting.Q.wl", 'allow_empty', true, ...
    'placeholder', 'empty = extinction peak', 'visible', @(I) sub_on(I,"Q"));
O(end+1) = opt('num', 'Plotting', 'Initial guess', 'path', "plotting.Q.guess", 'visible', @(I) sub_on(I,"Q"));
O(end+1) = opt('heading', 'Plotting', 'Incident field', 'visible', @plots_real);
O(end+1) = opt('check', 'Plotting', 'Show', 'path', "plotting.incidence.enable", 'visible', @plots_real);
O(end+1) = opt('num', 'Plotting', 'Wavelength (nm)', 'path', "plotting.incidence.wl", ...
    'visible', @(I) sub_on(I,"incidence"));
O(end+1) = opt('heading', 'Plotting', 'Apertures', 'visible', @plots_real);
O(end+1) = opt('check', 'Plotting', 'Show', 'path', "plotting.apertures.enable", 'visible', @plots_real);
O(end+1) = opt('vec', 'Plotting', 'Sizes (particles per side)', 'path', "plotting.apertures.sizes", ...
    'integer', true, 'visible', @(I) sub_on(I,"apertures"));
O(end+1) = opt('heading', 'Plotting', 'Emitted / transmitted field', 'visible', @plots_real);
O(end+1) = opt('check', 'Plotting', 'Show', 'path', "plotting.emitted.enable", 'visible', @plots_real);
O(end+1) = opt('list', 'Plotting', 'Types', 'path', "plotting.emitted.types", ...
    'items', ["fraunhofer","transmitted","far_field","quiver_ball"], 'height', 76, ...
    'visible', @(I) sub_on(I,"emitted"));
O(end+1) = opt('vec', 'Plotting', 'Wavelengths (nm)', 'path', "plotting.emitted.wls", 'allow_empty', true, ...
    'placeholder', 'empty = extinction peak', 'visible', @(I) sub_on(I,"emitted"));
O(end+1) = opt('num', 'Plotting', 'Angle range (deg)', 'path', "plotting.emitted.angle_lim", ...
    'limits', [0 90], 'visible', @(I) sub_on(I,"emitted"));
O(end+1) = opt('drop', 'Plotting', 'Fraunhofer quantity', 'path', "plotting.emitted.fraunhofer_quantity", ...
    'items', ["P/alpha","P"], 'visible', @(I) sub_on(I,"emitted"));
O(end+1) = opt('vec', 'Plotting', 'Far-field view angles x (rad)', 'path', "plotting.emitted.anglex", ...
    'visible', @(I) sub_on(I,"emitted") && any(string(getp(I,"plotting.emitted.types","")) == "far_field"));
O(end+1) = opt('vec', 'Plotting', 'Far-field view angles y (rad)', 'path', "plotting.emitted.angley", ...
    'allow_empty', true, 'placeholder', 'empty = same as x', ...
    'visible', @(I) sub_on(I,"emitted") && any(string(getp(I,"plotting.emitted.types","")) == "far_field"));
O(end+1) = opt('heading', 'Plotting', 'Eigenmodes', 'visible', @plots_real);
O(end+1) = opt('check', 'Plotting', 'Show coefficients', 'path', "plotting.coefficients.enable", ...
    'visible', @plots_real);
O(end+1) = opt('num', 'Plotting', 'Number of modes', 'path', "plotting.coefficients.n_modes", ...
    'integer', true, 'limits', [1 Inf], 'visible', @(I) sub_on(I,"coefficients"));
O(end+1) = opt('check', 'Plotting', 'Show band diagram', 'path', "plotting.eig_band.enable", ...
    'visible', @plots_real);
% ---- complex-frequency plots
O(end+1) = opt('heading', 'Plotting', 'Complex-frequency plots', 'visible', @plots_complex);
O(end+1) = opt('drop', 'Plotting', 'Real axis', 'path', "plotting.complex.real_axis", ...
    'items', ["eV","nm"], 'default', "eV", 'visible', @plots_complex, ...
    'tip', 'eV: real part of the frequency; nm: vacuum wavelength (Re range still entered in eV)');
O(end+1) = opt('check', 'Plotting', 'Whole computed range', 'path', "plotting.complex.full_range", ...
    'default', true, 'visible', @plots_complex);
O(end+1) = opt('vec', 'Plotting', 'Re range [from to] (eV)', 'path', "plotting.complex.re_lims", ...
    'default', [1.12 1.18], 'numel', 2, ...
    'visible', @(I) plots_complex(I) && ~logical(getp(I,"plotting.complex.full_range",true)));
O(end+1) = opt('vec', 'Plotting', '-Im range [from to] (meV)', 'path', "plotting.complex.im_lims", ...
    'default', [25 30], 'numel', 2, ...
    'visible', @(I) plots_complex(I) && ~logical(getp(I,"plotting.complex.full_range",true)));
O(end+1) = opt('check', 'Plotting', 'Plot the modes at the minima', 'path', "plotting.complex.modes", ...
    'default', true, 'visible', @plots_complex);
O(end+1) = opt('num', 'Plotting', 'Tiles per figure', 'path', "plotting.complex.tiles_per_figure", ...
    'default', 20, 'integer', true, 'limits', [1 Inf], 'visible', @cmodes_on);
O(end+1) = opt('num', 'Plotting', 'Highest spatial order shown', 'path', "plotting.complex.max_fft_order", ...
    'default', 1000, 'limits', [0 Inf], 'visible', @cmodes_on);
O(end+1) = opt('check', 'Plotting', 'Arrows (quiver) for 2D/3D modes', 'path', "plotting.complex.quiver", ...
    'default', true, 'visible', @(I) cmodes_on(I) && any(string(I.dimensionality) == ["2D","3D"]), ...
    'tip', 'needs slanCM (File Exchange) on the path');
O(end+1) = opt('check', 'Plotting', 'Colour map shows p_z (3D)', 'path', "plotting.complex.imagesc_is_z", ...
    'default', true, 'visible', @(I) cmodes_on(I) && string(I.dimensionality) == "3D");
O(end+1) = opt('check', 'Plotting', 'One tile per polarization', 'path', "plotting.complex.separate_pols", ...
    'default', false, 'visible', @(I) cmodes_on(I) && any(string(I.dimensionality) == ["2D","3D"]));
O(end+1) = opt('check', 'Plotting', 'Surfaces for scalar modes', 'path', "plotting.complex.surf_1D", ...
    'default', true, 'visible', @(I) cmodes_on(I) && ~any(string(I.dimensionality) == ["2D","3D"]));
O(end+1) = opt('drop', 'Plotting', 'Mode labels', 'path', "plotting.complex.labels", ...
    'items', ["sine","k","nodes","fft"], 'default', "sine", 'visible', @cmodes_on, ...
    'tip', ['sine: projection on the basis of the real-frequency mode tracking; k: the plane wave ' ...
            'that couples most strongly; nodes: sign changes of the mode; fft: peak of the spatial FFT (as before)']);
O(end+1) = opt('check', 'Plotting', 'k labels: undo the pi steps of shifted lattices', 'path', "plotting.complex.k_step", ...
    'default', false, 'visible', @(I) cmodes_on(I) && clabels(I, "k"), ...
    'tip', 'label by the step-matched (HG-shaped) field instead of a plane wave: one coupling peak per mode');
O(end+1) = opt('drop', 'Plotting', 'Label basis', 'path', "plotting.complex.basis", ...
    'items', ["auto","sine","fourier"], 'default', "auto", ...
    'visible', @(I) cmodes_on(I) && ~clabels(I, "fft"), ...
    'tip', 'auto: Fourier along periodic directions, sine otherwise (as the mode tracking)');
O(end+1) = opt('num', 'Plotting', 'Label confidence for "?"', 'path', "plotting.complex.min_confidence", ...
    'default', 0.25, 'limits', [0 1], 'visible', @(I) cmodes_on(I) && ~clabels(I, "fft"), ...
    'tip', 'labels below this are drawn grey with a ?');
O(end+1) = opt('check', 'Plotting', 'Label map of the complex plane', 'path', "plotting.complex.label_map", ...
    'default', false, 'visible', @(I) cmodes_on(I) && ~clabels(I, "fft"), ...
    'tip', 'mode number of the smallest-|lambda| eigenvector at every grid point');
O(end+1) = opt('check', 'Plotting', 'Show the spatial FFT of each mode', 'path', "plotting.complex.plot_fft", ...
    'default', true, 'visible', @(I) cmodes_on(I) && clabels(I, "fft"));
O(end+1) = opt('check', 'Plotting', 'Non-uniform FFT for mode labels', 'path', "plotting.complex.nufft", ...
    'default', false, 'visible', @(I) cmodes_on(I) && clabels(I, "fft"));
end

%% ------------------------------------------------------------------------
function o = opt(kind, tab, label, varargin)
% One option row. Name-value pairs:
%   path      field of Input, e.g. "particle.r"     key   sync key (default: path)
%   scale     SI value = shown value * scale        default  value if the field is missing
%   get/set   custom accessors (shown value)        visible  @(I) condition for showing
%   items     dropdown/list items                   numel, integer, allow_empty, limits
%   tip, placeholder, height, text (@(I) for 'info' rows), action (@(I) -> I for 'button' rows)
o = struct('kind', string(kind), 'tab', string(tab), 'label', string(label), 'path', "", ...
    'key', "", 'scale', 1, 'default', [], 'get', [], 'set', [], 'visible', @(I) true, ...
    'items', strings(0), 'numel', [], 'integer', false, 'allow_empty', false, 'limits', [-Inf Inf], ...
    'tip', "", 'placeholder', "", 'height', 22, 'text', @(I) "", 'action', [], 'save_path', "");
for n = 1:2:numel(varargin)
    o.(varargin{n}) = varargin{n+1};
end
if o.key == "", o.key = o.path; end
if o.key == "", o.key = o.kind + ":" + o.tab + ":" + o.label; end      % headings, notes
if isempty(o.get) && o.path ~= ""
    p = o.path; s = o.scale; d = o.default;
    o.get = @(I) scaled(getp(I, p, d), 1/s);
    o.set = @(I, v) setp(I, p, convert_like(scaled(v, s), getp(I, p, d)));
end
end

function o = range_opt(label, field, which, unit)
% Frequency range as from / to / points (range.real_parts, range.imag_parts)
o = opt('num', 'Physics', label, 'key', "range." + field + which, 'save_path', "range." + field, ...
    'get', @(I) range_get(I, field, which, unit), 'set', @(I, v) range_set(I, field, which, unit, v), ...
    'integer', which == 3);
if which == 3, o.limits = [1 Inf]; end
end

function v = range_get(I, field, which, unit)
x = I.range.(field);
switch which
    case 1, v = x(1)/unit;
    case 2, v = x(end)/unit;
    case 3, v = numel(x);
end
end

function I = range_set(I, field, which, unit, v)
x = I.range.(field);
a = x(1); b = x(end); n = numel(x);
switch which
    case 1, a = v*unit;
    case 2, b = v*unit;
    case 3, n = round(v);
end
I.range.(field) = linspace(a, b, n);
end

function o = inversion_opt(tab)
o = opt('drop', tab, 'Inversion tool', 'path', "inversion_tool", ...
    'items', ["eig","gmres","fft_gmres","mldivide","linsolve","inv","gmres_gpu","fft_gmres_gpu"], ...
    'tip', 'eig: eigenmode solver (needed for complex frequencies); others: direct solver');
end

%% ------------------------------------------------------------------------
%  Conditions and info texts
function tf = is_lorentz(I),  tf = any(strcmpi(string(I.sample.type), ["Lorentzian","lorentzian Au"])); end
function tf = any_periodic(I), tf = any(getp(I, "periodic_boundary", [0 0])); end
function tf = is_cyl(I),      tf = any(string(I.dimensionality) == ["cyl_TM","cyl_TE"]); end
function tf = layers_on(I),   tf = is_cyl(I) && logical(getp(I, "layers.enable", false)); end
function tf = is_eig(I),      tf = string(I.inversion_tool) == "eig"; end
function tf = is_gmres(I),    tf = contains(string(I.inversion_tool), "gmres"); end
function tf = par(I, m),      tf = string(getp(I, "parallel", "none")) == m; end
function tf = is_local_profile(I), tf = strcmpi(string(getp(I, "cluster.name", "")), "Processes"); end
function tf = preset(I, names), tf = any(strcmpi(string(I.preset_field), names)); end
function tf = plots_on(I),    tf = logical(getp(I, "plotting.enable", true)); end
function tf = is_complex(I),  tf = any(I.range.imag_parts ~= 0); end
function tf = plots_real(I),  tf = plots_on(I) && ~is_complex(I); end
function tf = plots_complex(I), tf = plots_on(I) && is_complex(I); end
function tf = cmodes_on(I),   tf = plots_complex(I) && logical(getp(I, "plotting.complex.modes", true)); end
function tf = clabels(I, m),  tf = strcmpi(string(getp(I, "plotting.complex.labels", "sine")), m); end
function tf = sub_on(I, name), tf = plots_real(I) && logical(getp(I, "plotting." + name + ".enable", false)); end
function tf = atype(I, names), tf = any(string(I.array_type) == names); end
function tf = is_pc(I), tf = atype(I, ["PC_mirrors","PC_mirror_reflection_test","PC_middle"]); end
function tf = tracking_on(I), tf = is_eig(I) && track_enable_get(I) ~= "off"; end

function v = track_enable_get(I)
e = getp(I, "track.enable", "auto");
if isstring(e) || ischar(e), v = "auto";
elseif logical(e),           v = "on";
else,                        v = "off";
end
end

function I = track_enable_set(I, v)
switch string(v)
    case "auto", I.track.enable = "auto";
    case "on",   I.track.enable = true;
    otherwise,   I.track.enable = false;
end
end

function I = open_array_editor(I, I0)
[ed, ok] = DDA_array_editor(I, I0);
if ~ok, return; end
if isempty(ed)
    if isfield(I, 'array_edit'), I = rmfield(I, 'array_edit'); end
else
    I.array_edit = ed;
end
end

function I = clear_layout_edit(I)
if isfield(I, 'array_edit'), I = rmfield(I, 'array_edit'); end
end

function tf = has_layout_edit(I)
tf = isfield(I, 'array_edit') && ~isempty(I.array_edit) && I.array_edit.enable;
end

function s = layout_info(I)
E = I.array_edit;
s = string(sprintf('Edited layout: %d particle(s) moved, %d resized.', E.n_moved, E.n_resized));
if isfield(E, 'base_inputs')
    f = fieldnames(E.base_inputs);
    changed = f(cellfun(@(n) ~isfield(I, n) || ~isequal(I.(n), E.base_inputs.(n)), f));
    if ~isempty(changed)
        s = s + " OUTDATED: " + strjoin(changed, ', ') + " changed since the edit; the run will " + ...
            "stop. Open the editor again or clear the edits.";
    end
end
end

function I = make_real(I)
I.range.imag_parts = 0;
end

function I = set_periodic(I, d, v)
pb = getp(I, "periodic_boundary", [0 0]);
pb(d) = double(v);
I.periodic_boundary = pb;
end

function s = range_info(I)
nr = numel(I.range.real_parts); ni = numel(I.range.imag_parts);
wl = 1239.84187./I.range.real_parts([1 end]);
s = string(sprintf('%d x %d = %d frequencies; real part %.1f - %.1f nm', nr, ni, nr*ni, min(wl), max(wl)));
if all(I.range.imag_parts == 0), s = s + " (real frequencies)"; end
end

function s = time_estimate(I)
% Solve-time estimate from the shared model (dda_time_estimate)
E = dda_time_estimate(I);
s = string(sprintf('Estimated (1 thread/worker, Ndof = %d): %.3g s per frequency x %d = %s serial', ...
    E.ndof, E.per_freq, E.nk, fmt_dur(E.per_freq*E.nk)));
if isfinite(E.workers) && E.workers > 1
    s = s + sprintf(', about %s on %d workers', fmt_dur(E.solve), E.workers);
end
s = s + ". Excludes tracking and queueing.";
if strlength(E.layer_note) > 0
    s = s + " Includes " + sprintf('%.3g s per frequency', E.layers.per_freq_used) + " for the layers (" + E.layer_note + ").";
end
end

function s = walltime_check(I)
E = dda_time_estimate(I);
req = parse_walltime(getp(I, "cluster.wall_time", ""));
est = E.run;                                   % solve + pool start-up
note = "";
w = string(getp(I, "cluster.wall_time", ""));
if ~isnan(req) && ~isempty(regexp(w, ':\d+:\d+$|^\d+:\d+$', 'once'))
    secs = str2double(regexp(w, '(\d+)$', 'tokens', 'once'));
    if secs > 0, note = sprintf(' Slurm counts whole minutes: %s becomes %s.', w, fmt_dur(req)); end
end
if isnan(req)
    s = sprintf('Estimated %s on %d workers (incl. ~%d s start-up).', fmt_dur(est), E.workers, E.job_startup);
elseif est > req
    s = sprintf('Estimated %s on %d workers EXCEEDS the wall time %s.', fmt_dur(est), E.workers, fmt_dur(req));
else
    s = sprintf('Estimated %s on %d workers (incl. ~%d s start-up), %.0f%% of the wall time %s.', ...
        fmt_dur(est), E.workers, E.job_startup, 100*est/req, fmt_dur(req));
end
s = string(s) + note;
end

function sec = parse_walltime(w)
% as Slurm applies it: whole minutes, rounded up (dda_walltime_seconds)
sec = dda_walltime_seconds(w);
end

function s = fmt_dur(sec)
if sec < 60,       s = sprintf('%.0f s', sec);
elseif sec < 3600, s = sprintf('%.1f min', sec/60);
else,              s = sprintf('%.1f h', sec/3600);
end
end

function s = cores_info(I)
w = getp(I, "cluster.num_workers", 1); t = getp(I, "cluster.num_threads", 1);
c = (w+1)*t;
s = string(sprintf('Cores reserved: (%d workers + 1) x %d threads = %d', w, t, c));   % string: + appends
if ~is_local_profile(I)
    node = 384;                                     % cores of a Roihu CPU node
    if c > node
        s = s + sprintf(' - MORE than one node (%d cores); the job will not start with -N 1', node);
    else
        s = s + sprintf(' of %d on a Roihu CPU node', node);
        part = lower(string(getp(I, "cluster.partition", "")));
        if any(part == ["medium","large"]) && c < node
            s = s + sprintf(' (the whole node is reserved: %d cores unused)', node - c);
        end
        % total memory of the request against a node: 768 GiB, of which about
        % 744 GiB is usable for jobs (CSC); Slurm memory units are binary
        mem = parse_mem(getp(I, "cluster.mem_per_cpu", ""));
        usable = 744*2^30;
        if isfinite(mem)
            total = c*mem;
            s = s + sprintf('; memory %.0f GiB in total', total/2^30);
            if total > usable
                s = s + sprintf(' - MORE than a node has (~744 GiB): lower mem_per_cpu to at most %.0fmb', ...
                    floor(usable/c/2^20));
            end
        elseif c > 0
            s = s + "; memory per core: the profile's default (MemPerCPU) - check it fits";
        end
    end
end
end

function I = use_full_node(I)
% Fill one Roihu CPU node (384 cores): whole-node partition, workers from the
% chosen threads per worker ((workers + 1) x threads <= 384), node memory
node = 384;
t = max(1, getp(I, "cluster.num_threads", 1));
I.cluster.partition   = "medium";
I.cluster.num_workers = floor(node/t) - 1;
% An empty mem_per_cpu would mean the profile's MemPerCPU default (e.g. 4gb),
% and 384 x 4 GiB = 1536 GiB is more than a node has: Slurm then refuses the
% job ("Requested node configuration is not available"). About 744 GiB of a
% 768 GiB node is usable: 384 x 1900 MiB = 712.5 GiB fits.
I.cluster.mem_per_cpu = "1900mb";
end

function s = memory_info(I)
% Memory estimate (same model as dda_size_check, which checks it before submitting)
J = I;
J.N = prod(J.arrayN);
if isfield(J, 'array_type') && string(J.array_type) == "double_rectangular", J.N = 2*J.N; end
J.range.k_N = numel(J.range.real_parts)*numel(J.range.imag_parts);
R = dda_result_size(J);
arr = 16*R.nk*(9 + 2) + 24*J.N;                    % compact polarizability, k and n, positions
inp = 4*16*R.nk;                                   % range vectors in Input
worker = inp + arr + R.eig_workspace;
client = inp + arr + R.results + R.vectors;
s = string(sprintf('Memory: results ~%s%s; each worker ~%s; job client ~%s', human_bytes(R.results), ...
    ternary_s(R.lean, ' (lean)', ''), human_bytes(worker), human_bytes(client)));
if R.track, s = s + sprintf(' (incl. %s of eigenvectors for tracking)', human_bytes(R.vectors)); end
if par(I, "cluster") && ~is_local_profile(I)
    lim = parse_mem(getp(I, "cluster.mem_per_cpu", "")) * max(1, getp(I, "cluster.num_threads", 1));
    if isfinite(lim)
        s = s + sprintf('; %s per worker requested', human_bytes(lim));
        if client > lim || worker > lim, s = s + " - DOES NOT FIT: raise mem_per_cpu or threads, or split the grid"; end
    end
end
if R.results > 1.9e9
    s = s + " - results over the 2 GB .mat limit: split the frequency grid";
end
s = s + ".";
end

function b = parse_mem(m)
% Slurm memory such as "4g", "4gb", "4000m" in bytes (Slurm units are binary:
% 1 mb = 1 MiB, 1 gb = 1 GiB; a plain number is MiB); NaN if unset
b = NaN;
m = strtrim(lower(string(m)));
if strlength(m) == 0, return; end
t = regexp(m, '^([\d.]+)\s*([kmgt]?)b?$', 'tokens', 'once');
if isempty(t), return; end
f = struct('k', 2^10, 'm', 2^20, 'g', 2^30, 't', 2^40);
u = t{2}; if u == "", u = "m"; end
b = str2double(t{1}) * f.(char(u));
end

function s = human_bytes(b)
if b >= 1e9, s = sprintf('%.2f GB', b/1e9); elseif b >= 1e6, s = sprintf('%.0f MB', b/1e6); else, s = sprintf('%.0f kB', b/1e3); end
end

%% ------------------------------------------------------------------------
%  Field access and parsing
function v = getp(I, path, default)
parts = split(string(path), ".");
v = I;
for p = parts.'
    if isstruct(v) && isfield(v, p), v = v.(p); else, v = default; return; end
end
end

function I = setp(I, path, v)
parts = cellstr(split(string(path), "."));
I = setfield(I, parts{:}, v);
end

function v = scaled(v, factor)
% Unit conversion for numbers only (text and logical values pass through)
if isnumeric(v) && factor ~= 1, v = v .* factor; end
end

function v = convert_like(v, old)
% Keep the type of the default: char stays char, string stays string
if ischar(old) && (isstring(v) || ischar(v)),    v = char(v);
elseif isstring(old) && ischar(v),               v = string(v);
elseif islogical(old) && isnumeric(v) && isscalar(v), v = logical(v);
elseif isnumeric(old) && islogical(v),           v = double(v);
end
end

function v = parse_vec(txt, o)
txt = strtrim(char(txt));
if isempty(txt) || strcmp(txt, '[]')
    if ~o.allow_empty, error('%s cannot be empty.', o.label); end
    v = [];
    return
end
[v, okp] = str2num(txt); %#ok<ST2NM>  evaluates expressions such as linspace(1,2,5)
if ~okp || ~isnumeric(v)
    error('"%s" is not a number, vector or MATLAB expression.', txt);
end
if ~isempty(o.numel) && numel(v) ~= o.numel
    error('%s needs %d values, got %d.', o.label, o.numel, numel(v));
end
if o.integer && any(v(:) ~= round(v(:)))
    error('%s must be whole numbers.', o.label);
end
end

function s = vec2str(v)
if isempty(v), s = ''; elseif isscalar(v), s = num2str(v, 10); else, s = mat2str(v, 8); end
end

function paths = option_paths(opts)
% every Input field the window edits (option paths, and the fields behind the
% range, periodic-boundary and tracking rows)
paths = strings(0);
for k = 1:numel(opts)
    p = opts(k).path;
    if p == "", p = opts(k).save_path; end
    if p ~= "" && ~any(paths == p), paths(end+1) = p; end %#ok<AGROW>
end
end

function s = ternary_s(c, a, b)
if c, s = a; else, s = b; end
end

function file = saved_defaults_file()
% DDA_saved_defaults.m lives next to DDA_Master.m
root = fileparts(which('DDA_Master'));
if isempty(root), root = pwd; end
file = fullfile(root, 'DDA_saved_defaults.m');
end

function [paths, scales] = changed_paths(I0, W, opts)
% Input fields edited in the GUI whose value differs from I0
paths = strings(0); scales = [];
for k = 1:numel(opts)
    o = opts(k);
    p = o.path; sc = o.scale;
    if p == "", p = o.save_path; sc = 1; end
    if p == "" || any(paths == p), continue; end
    [a, ha] = getp_has(I0, p); [b, hb] = getp_has(W, p);
    if hb && (~ha || ~isequal(a, b))
        paths(end+1) = p; scales(end+1) = sc; %#ok<AGROW>
    end
end
end

function [v, has] = getp_has(I, path)
v = []; has = false;
parts = split(string(path), ".");
s = I;
for q = parts.'
    if isstruct(s) && isfield(s, q), s = s.(q); else, return; end
end
v = s; has = true;
end

function write_saved_defaults(file, paths, scales, W, has_layout)
if nargin < 5, has_layout = false; end
lines = strings(0);
for n = 1:numel(paths)
    v = getp(W, paths(n), []);
    lines(end+1) = sprintf('Input.%s = %s;', paths(n), value_expr(v, scales(n))); %#ok<AGROW>
end
if has_layout
    % positions and radii are too many for text: they are in DDA_saved_layout.mat
    lines(end+1) = "";
    lines(end+1) = "% Array layout edited in the array editor (in DDA_saved_layout.mat, same folder)";
    lines(end+1) = "L = load(fullfile(fileparts(mfilename('fullpath')), 'DDA_saved_layout.mat'), 'array_edit');";
    lines(end+1) = "Input.array_edit = L.array_edit;";
    paths(end+1) = "array layout";
end
txt = [ ...
"function Input = DDA_saved_defaults(Input)"
"% DDA_SAVED_DEFAULTS  Defaults saved from the parameter GUI (""Save defaults"")."
"% Saved " + string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')) + "."
"%"
"% Applied by Make_predefined_Inputs (DDA_Master.m) right after its own values,"
"% when Input.use_saved_defaults is true, so these settings win. Edit by hand,"
"% or delete this file to use the values written in DDA_Master.m."
""
lines(:)
""
"fprintf('Saved defaults applied (DDA_saved_defaults.m): %s\n', '" + strjoin(paths, ', ') + "');"
"end"];
fid = fopen(file, 'w');
if fid < 0, error('Cannot write %s', file); end
fprintf(fid, '%s\n', txt);
fclose(fid);
end

function s = value_expr(v, scale)
% A MATLAB expression for v, in the GUI units (e.g. [700 700 100]*1e-9)
if isstring(v)
    q = """" + replace(v, """", """""") + """";
    if isscalar(q), s = q; else, s = "[" + strjoin(q, ", ") + "]"; end
elseif ischar(v)
    s = "'" + replace(string(v), "'", "''") + "'";
elseif islogical(v)
    if isscalar(v), s = string(mat2str(v)); else, s = string(mat2str(v)); end
elseif isnumeric(v)
    if isempty(v), s = "[]"; return; end
    x = double(v)/scale;
    d = diff(x(:));
    if isvector(x) && numel(x) > 3 && all(abs(d - d(1)) <= 1e-9*max(abs(x(:))))
        s = sprintf('linspace(%s, %s, %d)', num2str(x(1), 12), num2str(x(end), 12), numel(x));
        if size(x,1) > 1, s = s + "'"; end
    else
        s = string(mat2str(x, 10));
    end
    if scale ~= 1
        s = s + "*" + regexprep(sprintf('%.0e', scale), 'e([+-])0*(\d)', 'e$1$2');
    end
else
    error('Cannot save a value of class %s.', class(v));
end
end

function report_changes(I0, W, opts)
% Print the options that differ from the defaults
lines = strings(0);
done = strings(0);
for k = 1:numel(opts)
    o = opts(k);
    if isempty(o.get) || any(done == o.key), continue; end
    done(end+1) = o.key; %#ok<AGROW>
    a = o.get(I0); b = o.get(W);
    if ~isequal(a, b)
        lines(end+1) = sprintf('  %-40s %s  ->  %s', o.tab + ": " + o.label, show(a), show(b)); %#ok<AGROW>
    end
end
if isfield(W, 'array_edit') && ~isempty(W.array_edit) && ...
        ~(isfield(I0, 'array_edit') && isequal(I0.array_edit, W.array_edit))
    lines(end+1) = sprintf('  %-40s %d moved, %d resized (#%s)', 'Array: edited layout', ...
        W.array_edit.n_moved, W.array_edit.n_resized, W.array_edit.hash);
elseif isfield(I0, 'array_edit') && ~isfield(W, 'array_edit')
    lines(end+1) = "  Array: layout edits cleared";
end
if isempty(lines)
    fprintf('Run options: defaults of Make_predefined_Inputs.\n');
else
    fprintf('Run options changed in the GUI:\n%s\n', strjoin(lines, newline));
end
end

function s = show(v)
if isstring(v) || ischar(v), s = sprintf('"%s"', strjoin(string(v), '", "'));
elseif islogical(v) && isscalar(v), s = mat2str(v);
elseif isempty(v), s = '[]';
else, s = vec2str(v);
end
end
