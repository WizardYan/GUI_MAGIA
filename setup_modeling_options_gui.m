function setup_modeling_options_gui()
    f = figure('Name','Modeling Options Setup','Position',[300 300 460 360]);

    fontSize = 9;
    ctrlHeight = 22;
    xInput = 220;
    inputWidth = 200;
    y = 300;
    spacing = 40;

    % Helper to create auto-sized label
    function makeLabel(text, ypos)
        temp = uicontrol(f, 'Style', 'text', ...
            'String', text, ...
            'FontSize', fontSize, ...
            'Units', 'pixels', ...
            'Position', [0 0 10 10], ...
            'Visible', 'off');
        extent = get(temp, 'Extent');
        delete(temp);
        uicontrol(f, 'Style', 'text', ...
            'String', text, ...
            'Position', [20 ypos extent(3) ctrlHeight], ...
            'FontSize', fontSize, ...
            'HorizontalAlignment', 'left');
    end

    % Model
    makeLabel('Model:', y);
    modelPopup = uicontrol(f, 'Style', 'popupmenu', ...
        'String', {'srtm','logan'}, ...
        'Position', [xInput y inputWidth ctrlHeight], ...
        'FontSize', fontSize, 'Value', 1);

    % Lower bounds
    y = y - spacing;
    makeLabel('Lower Bounds [R1, k2, BP]:', y);
    lbEdit = uicontrol(f, 'Style', 'edit', ...
        'Position', [xInput y inputWidth ctrlHeight], ...
        'String', '[0,0,0]', 'FontSize', fontSize);

    % Upper bounds
    y = y - spacing;
    makeLabel('Upper Bounds [R1, k2, BP]:', y);
    ubEdit = uicontrol(f, 'Style', 'edit', ...
        'Position', [xInput y inputWidth ctrlHeight], ...
        'String', '[ , ,]', 'FontSize', fontSize);

    % theta3 lower
    y = y - spacing;
    makeLabel('Lower Bound [theta3]:', y);
    theta3lbEdit = uicontrol(f, 'Style', 'edit', ...
        'Position', [xInput y inputWidth ctrlHeight], ...
        'String', '', 'FontSize', fontSize);

    % theta3 upper
    y = y - spacing;
    makeLabel('Upper Bound [theta3]:', y);
    theta3ubEdit = uicontrol(f, 'Style', 'edit', ...
        'Position', [xInput y inputWidth ctrlHeight], ...
        'String', '', 'FontSize', fontSize);

    % Number of bases
    y = y - spacing;
    makeLabel('Number of Bases:', y);
    nbasesEdit = uicontrol(f, 'Style', 'edit', ...
        'Position', [xInput y inputWidth ctrlHeight], ...
        'String', '300', 'FontSize', fontSize);

    % Save button
    uicontrol(f, 'Style', 'pushbutton', ...
        'String', 'Save Options', ...
        'Position', [150 20 160 30], ...
        'FontSize', fontSize + 1, ...
        'Callback', @(~,~)save_callback());

    %% Save logic (TXT only)
    function save_callback()
        modeling_options.model = modelPopup.String{modelPopup.Value};
        modeling_options.lb = str2num(get(lbEdit, 'String')); %#ok<ST2NM>
        modeling_options.ub = str2num(get(ubEdit, 'String')); %#ok<ST2NM>
        modeling_options.theta3_lb = str2double(get(theta3lbEdit, 'String'));
        modeling_options.theta3_ub = str2double(get(theta3ubEdit, 'String'));
        modeling_options.nbases = str2double(get(nbasesEdit, 'String'));

        [file, path] = uiputfile('modeling_options.txt', 'Save Options As');
        if isequal(file, 0)
            msgbox('Save canceled. No file written.','Notice');
            return;
        end

        txtfile = fullfile(path, file);
        fid = fopen(txtfile, 'w');
        fprintf(fid, 'model: %s\n', modeling_options.model);
        fprintf(fid, 'lb: %s\n', mat2str(modeling_options.lb));
        fprintf(fid, 'ub: %s\n', mat2str(modeling_options.ub));
        fprintf(fid, 'theta3_lb: %.6f\n', modeling_options.theta3_lb);
        fprintf(fid, 'theta3_ub: %.6f\n', modeling_options.theta3_ub);
        fprintf(fid, 'nbases: %d\n', modeling_options.nbases);
        fclose(fid);

        msgbox('Modeling options saved to .txt file.','Success');
    end
end
