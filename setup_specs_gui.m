function setup_specs_gui()
    %% Field definitions (template removed from here)
    fields = {
        'tracer', '';
        'scanner', '';
        'model', '';
        'dc', '1';
        'rc', '0';
        'cut_time', '90';
        'roi_set', '';
        'input_type', '';
        'roi_type', '';
        'fwhm_pre', '5';
        'fwhm_post', '5';
        'fwhm_roi', '0';
        'cpi', '1';
        'norm_method', 'pet';        
        'mc_fwhm', '7';
        'mc_excluded_frames', '1';
        'mc_ref_frame', '0';
        'mc_rtm', '1';
        'gu', '';
    };

    comments = {
        'Radiotracer used, e.g., [11c]nnc';
        'PET scanner type, e.g., hrrt';
        'Kinetic model to use, e.g., srtm';
        'Already Decay Corrected: 1=yes, 0=no';
        'Run ROI correction: 1=yes, 0=no';
        'Analysis cut time in minutes';
        'ROI set name, e.g., tracer_default';
        'Input type: ref or plasma';
        'ROI type: use atlas or freesurfer-based ROI definition';
        'Pre-smoothing FWHM (mm)';
        'Post-smoothing FWHM (mm)';
        'ROI smoothing FWHM (mm)';
        'Run Voxel model: 1=yes;0=no';
        'Normalization method, e.g., pet or mri';
        'Motion correction: FWHM (mm)';
        'Motion correction: Frames to skip motion correction';
        'Motion correction: Reference frame';
        'Motion correction: Realign to mean, 1=yes;0=no';
        'Glucose';
    };

    %% Layout settings
    numFields = size(fields, 1);
    vSpacing = 24;
    topMargin = 40;
    ctrlHeight = 18;
    baseHeight = (numFields + 3) * vSpacing + 140;
    winW = 820;
    winH = baseHeight;

    % Center on screen, lower slightly
    screenSize = get(0, 'ScreenSize');
    winX = round((screenSize(3) - winW) / 2);
    winY = max(round((screenSize(4) - winH) / 2), 50);

    %% Create figure
    f = figure('Name', 'Set Specs', ...
               'Position', [winX winY winW winH], ...
               'Resize', 'off');

    controls = cell(numFields, 2);

    %% Create input fields
    for i = 1:numFields
        y_pos = winH - topMargin - vSpacing * (i - 1);

        uicontrol(f, 'Style', 'text', 'String', fields{i,1}, ...
            'Position', [20 y_pos 100 ctrlHeight], ...
            'FontSize', 8, 'HorizontalAlignment', 'left');

        if strcmp(fields{i,1},'tracer')
            controls{i,2} = uicontrol(f, 'Style', 'popupmenu', 'String', {'Click&Select','[11c]nnc', '[11c]raclopride'},'value',1,...
                'Position',[120 y_pos 280 ctrlHeight],'FontSize', 8);

        elseif strcmp(fields{i,1},'scanner')
            controls{i,2} = uicontrol(f, 'Style', 'popupmenu', 'String', {'Unknown','HRRT', 'PET/CT'},'value',1,...
                'Position',[120 y_pos 280 ctrlHeight],'FontSize', 8);

        elseif strcmp(fields{i,1},'model')
            controls{i,2} = uicontrol(f, 'Style', 'popupmenu', 'String', {'Click&elect','srtm','logan'},'value',1,...
                'Position',[120 y_pos 280 ctrlHeight],'FontSize', 8);

        elseif strcmp(fields{i,1},'roi_type')
            controls{i,2} = uicontrol(f, 'Style', 'popupmenu', 'String', {'atlas', 'freesurfer'},'value',1,...
                'Position',[120 y_pos 280 ctrlHeight],'FontSize', 8);

        elseif strcmp(fields{i,1},'input_type')
            controls{i,2} = uicontrol(f, 'Style', 'popupmenu', 'String', {'ref', 'plasma'},'value',1,...
                'Position',[120 y_pos 280 ctrlHeight],'FontSize', 8);

        elseif strcmp(fields{i,1},'norm_method')
            controls{i,2} = uicontrol(f, 'Style', 'popupmenu', 'String', {'pet', 'mri'},'value',1,...
                'Position',[120 y_pos 280 ctrlHeight],'FontSize', 8);

        else
        controls{i,2} = uicontrol(f, 'Style', 'edit', 'String', fields{i,2}, ...
            'Position', [120 y_pos 280 ctrlHeight], ...
            'FontSize', 8);
        end
        uicontrol(f, 'Style', 'text', 'String', comments{i}, ...
            'Position', [410 y_pos 390 ctrlHeight], ...
            'FontSize', 7, 'ForegroundColor', [0.4 0.4 0.4], ...
            'HorizontalAlignment', 'left');

        controls{i,1} = fields{i,1};
    end





    % %% ROI Type dropdown
    % base_y = y_pos - vSpacing;
    % uicontrol(f, 'Style', 'text', 'String', 'roi_type', ...
    %     'Position', [20 base_y 100 ctrlHeight], ...
    %     'FontSize', 8, 'HorizontalAlignment', 'left');
    % roiTypeMenu = uicontrol(f, 'Style', 'popupmenu', 'String', {'atlas', 'freesurfer'}, ...
    %     'Position', [120 base_y 200 ctrlHeight], ...
    %     'FontSize', 8, 'Value', 1);
    % uicontrol(f, 'Style', 'text', 'String', 'ROI type selection method', ...
    %     'Position', [330 base_y 300 ctrlHeight], ...
    %     'FontSize', 7, 'ForegroundColor', [0.4 0.4 0.4], ...
    %     'HorizontalAlignment', 'left');

    %% mni_roi_mask_dir input with browse
    base_y = y_pos - vSpacing;
    uicontrol(f, 'Style', 'text', 'String', 'mni_roi_mask_dir', ...
        'Position', [20 base_y 100 ctrlHeight], ...
        'FontSize', 8, 'HorizontalAlignment', 'left');
    roiPathEdit = uicontrol(f, 'Style', 'edit', ...
        'String', 'atlas_rois', ...
        'Position', [120 base_y 280 ctrlHeight], ...
        'FontSize', 8);
    uicontrol(f, 'Style', 'pushbutton', 'String', 'Browse...', ...
        'Position', [410 base_y 70 ctrlHeight], ...
        'FontSize', 8, ...
        'Callback', @(~,~)browse_dir());
    uicontrol(f, 'Style', 'text', 'String', 'Path to ROI mask NIfTI directory', ...
        'Position', [490 base_y 300 ctrlHeight], ...
        'FontSize', 7, 'ForegroundColor', [0.4 0.4 0.4], ...
        'HorizontalAlignment', 'left');

    %% template input with browse
    base_y = base_y - vSpacing;
    uicontrol(f, 'Style', 'text', 'String', 'template', ...
        'Position', [20 base_y 100 ctrlHeight], ...
        'FontSize', 8, 'HorizontalAlignment', 'left');
    templatePathEdit = uicontrol(f, 'Style', 'edit', ...
        'String', 'magia-master/templates/[11c]nnc.nii', ...
        'Position', [120 base_y 280 ctrlHeight], ...
        'FontSize', 8);
    uicontrol(f, 'Style', 'pushbutton', 'String', 'Browse...', ...
        'Position', [410 base_y 70 ctrlHeight], ...
        'FontSize', 8, ...
        'Callback', @(~,~)browse_template());
    uicontrol(f, 'Style', 'text', 'String', 'Path to the template NIfTI', ...
        'Position', [490 base_y 300 ctrlHeight], ...
        'FontSize', 7, 'ForegroundColor', [0.4 0.4 0.4], ...
        'HorizontalAlignment', 'left');

    %% Save Button
    uicontrol(f, 'Style', 'pushbutton', 'String', 'Save to .txt', ...
        'Position', [winW/2 - 50 20 100 26], ...
        'FontSize', 9, ...
        'Callback', @(~,~)save_callback());

    %% Nested functions
    function browse_dir()
        folder = uigetdir;
        if folder ~= 0
            set(roiPathEdit, 'String', folder);
        end
    end

    function browse_template()
        [file, path] = uigetfile('*.nii', 'Select Template NIfTI File');
        if isequal(file, 0)
            return;
        end
        set(templatePathEdit, 'String', fullfile(path, file));
    end

    function save_callback()
        [file, path] = uiputfile('*.txt', 'Save Specs As');
        if isequal(file, 0)
            msgbox('Save canceled.', 'Notice');
            return;
        end

        fid = fopen(fullfile(path, file), 'w');
        for i = 1:numFields
            key = controls{i,1};
            % val = get(controls{i,2}, 'String');

            if strcmp(controls{i,1}, 'tracer')
                items = get(controls{i,2}, 'String');
                val = items{get(controls{i,2}, 'Value')};

            elseif strcmp(controls{i,1}, 'scanner')
                items = get(controls{i,2}, 'String');
                val = items{get(controls{i,2}, 'Value')};


            elseif strcmp(controls{i,1},'roi_type')
                items = get(controls{i,2}, 'String');
                val = items{get(controls{i,2}, 'Value')};

            elseif strcmp(controls{i,1},'input_type')
                items = get(controls{i,2}, 'String');
                val = items{get(controls{i,2}, 'Value')};

            elseif strcmp(controls{i,1},'model')
                items = get(controls{i,2}, 'String');
                val = items{get(controls{i,2}, 'Value')};

            elseif strcmp(controls{i,1},'norm_method')
                items = get(controls{i,2}, 'String');
                val = items{get(controls{i,2}, 'Value')};
                
            else
                val = get(controls{i,2}, 'String');
            end


            fprintf(fid, '%s: %s\n', key, val);
        end
        % fprintf(fid, 'roi_type: %s\n', roiTypeMenu.String{roiTypeMenu.Value});
        fprintf(fid, 'mni_roi_mask_dir: %s\n', get(roiPathEdit, 'String'));
        fprintf(fid, 'template: %s\n', get(templatePathEdit, 'String'));
        fclose(fid);

        msgbox('Parameters saved to file.', 'Success');
    end
end