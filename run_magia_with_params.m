function run_magia_with_params(data_parent_directory, sbj, pet_tr, workfolder, sessions, modeling_options_file, specs_file)

    modeling_options = load_modeling_options_txt(modeling_options_file);
    disp(sbj)
    ses={'ses-1'  'ses-2'};

    cd(workfolder);
    addpath('./spm12')
    addpath(workfolder);
    addpath(genpath(strcat(workfolder,'/magia-master/')));


    function modeling_options = load_modeling_options_txt(modeling_options_file)
        if ~isfile(modeling_options_file)
            error('File does not exist: %s', modeling_options_file);
        end
    
        % Initialize struct
        modeling_options = struct();
    
        % Open and read file
        fid = fopen(modeling_options_file, 'r');
        if fid == -1
            error('Failed to open file: %s', modeling_options_file);
        end
    
        % Read and parse each line
        while ~feof(fid)
            line = strtrim(fgetl(fid));
            if isempty(line)
                continue;
            end
            tokens = strsplit(line, ':');
            if numel(tokens) < 2
                continue;
            end
            key = strtrim(tokens{1});
            valStr = strtrim(strjoin(tokens(2:end), ':'));
    
            % Parse and assign
            switch lower(key)
                case 'model'
                    modeling_options.model = valStr;
                case 'lb'
                    modeling_options.lb = str2num(valStr); %#ok<ST2NM>
                case 'ub'
                    modeling_options.ub = str2num(valStr); %#ok<ST2NM>
                case 'theta3_lb'
                    modeling_options.theta3_lb = str2double(valStr);
                case 'theta3_ub'
                    modeling_options.theta3_ub = str2double(valStr);
                case 'nbases'
                    modeling_options.nbases = str2double(valStr);
                otherwise
                    warning('Unknown key: %s', key);
            end
        end
        fclose(fid);
    end


    for kl = sessions %for nnc only sess1; for rac sess1 and sess2
    
        sess = ses{kl}

        session_dir = fullfile(data_parent_directory, sbj, sess, pet_tr);
    
        % Check if the session directory exists
        if ~isfolder(session_dir)
            fprintf('Session folder does not exist: %s. Skipping.\n', session_dir);
            continue;
        end
    
        fprintf('Processing subject %s, session %s\n', sbj, sess);

        setenv('DATA_DIR',strcat(workfolder,'/',pet_tr,'/',sess));
        setenv('MAGIA_ARCHIVE',strcat(workfolder,'/MAGIA_ARCHIEVE/',pet_tr,'/',sess))
        setenv('SPM_DIR',strcat(workfolder, '/spm12'))
        setenv('MAGIA_PATH',strcat(workfolder,'/magia-master/'))
        setenv('FS_LICENSE','/data/yanw4/license.txt');
        setenv('MRI_DIR',strcat(workfolder,'/MRI'))
        setenv('FREESURFER_DIR','/usr/local/apps/freesurfer/7.4.1/bin/freesurfer') %needed if freesurfer is run from Magia
        setenv('FREESURFER_HOME','/vf/users/LNI/Project_114/freesurfer_out') %needed if freesurfer is run from Magia
        
        setenv('FS_FINAL_DIR','/vf/users/LNI/Project_114/freesurfer_out')
        setenv('FS_TEMP_DIR',strcat(workfolder,'/FREESURFER_TEMP'))
        magia_check_envs


        % if ~exist(strcat(workfolder,'/MRI/',sbj,'/T1/mri_',sbj,'.nii'),'file')
        %     if exist(strcat('/data/LNI/Project_114/freesurfer_out/',sbj,'/mri/T1.mgz'), 'file')
        %         display("Copying T1w from FREESURFER FOLDER")            
        %         mkdir(strcat(workfolder,'/MRI/',sbj,'/T1'))
        %         copyfile(['/vf/users/LNI/Project_114/freesurfer_out/',sbj,'/mri/T1.mgz'],[workfolder,'/MRI/',sbj,'/T1/'],'f')
        %         z_cmd=['mri_convert ' , workfolder,'/MRI/',sbj,'/T1/T1.mgz ',workfolder,'/MRI/',sbj,'/T1/T1.nii'  ]
        %         system(z_cmd)
        %         z_cmd=['mv ', workfolder,'/MRI/',sbj,'/T1/T1.nii ', workfolder,'/MRI/',sbj,'/T1/mri_',sbj,'.nii']
        %         system(z_cmd)
        % 
        %         %Import basic specifications
        %         % specs_file=strcat(workfolder, '/modeling_options_',pet_tr,'_srtm/',pet_tr,'_specs_srtm.txt');
        %         specs_file=strcat(workfolder, '/modeling_options_',pet_tr,'_srtm_yan/',pet_tr,'_specs_srtm.txt'); % no smoothing, higher low bound
        % 
        %         specs=magia_read_specs(specs_file);
        % 
        %     else
        %         % Import basic specifications
        %         % specs_file=strcat(workfolder, '/modeling_options_',pet_tr,'_srtm_pet_normalization/',pet_tr,'_specs_srtm.txt');
        %         specs_file=strcat(workfolder, '/modeling_options_',pet_tr,'_srtm_pet_normalization_yan/',pet_tr,'_specs_srtm.txt');
        %         specs=magia_read_specs(specs_file);
        %     end
        % 
        %     specs=magia_read_specs(specs_file);
        % 
        % else
        % 
        %     specs=magia_read_specs(specs_file);
        % 
        % end
    
    
        specs=magia_read_specs(specs_file);


        %% DOUBLE CHECK MODEL AND RADOTRACER'S CONSISTENCY;
        % Check model consistency
        if ~strcmp(specs.magia.model, modeling_options.model)
            error('Mismatch in model types: specs.magia.model = %s, modeling_options.model = %s', ...
                specs.magia.model, modeling_options.model);
        end
        
        % Extract template tracer name
        template_name = specs.magia.template;
        tracer_from_template = regexp(template_name, '\[.*?\][a-zA-Z0-9_-]+', 'match', 'once');
        
        % Check tracer consistency
        if ~strcmp(tracer_from_template, specs.study.tracer)
            error('Mismatch in tracer: template = %s, specs.study.tracer = %s', ...
                template_name, specs.study.tracer);
        end


        %% DOUBLE CHECK SPECIFICATIONS IF HAVE JSON FILES
        %Edit subject
        specs.study.mri_code=sbj;
        
        fname = dir(strcat(data_parent_directory, '/', sbj,'/',sess,'/',pet_tr,'/*json'));
    
        if ((size(fname, 1) == 1))
    
        %Read json frame durations
        fname = [data_parent_directory, '/', sbj,'/',sess,'/',pet_tr,'/',sbj,'_',sess,'_pet.json'];
        fid = fopen(fname);
        raw = fread(fid,inf);
        str = char(raw');
        fclose(fid);
        val = jsondecode(str);
    
        frames=val.FrameDuration;
        numFrames=length(frames)
        new_frames=[0 0];
            for fr=1:(numFrames)
                tmp=[new_frames(fr,2) new_frames(fr,2)+frames(fr)];
                new_frames=[new_frames;tmp];
            end
    
        new_frames=new_frames(2:end,:);
    
        %Replace new dose and weight values and fix frame times turn them into minutes
        specs.study.weight='NaN'; %weight;
        specs.study.dose=val.RadionuclideTotalDose; %'NaN'; %dose;
    
        specs.study.frames = new_frames/60 %specs.study.frames/60 %Enter frame durations
        specs.study.scanner = val.ManufacturersModelName ;%enter HRRT/PETCT info
        % spec.study.mc_excluded_frames= 1 ; %[1 numFrames] %Exclude 1st first frame from the analysis
        end    
    
        %  MAKE A COPY OF THE PET IMAGE
        fname = dir(strcat(data_parent_directory, '/',sbj,'/',sess,'/',pet_tr,'/*nii'));
    
        if ((size(fname, 1) == 1))
    
            mkdir([workfolder,'/',pet_tr,'/', sess,'/',sbj,'/PET/nii']);
            disp('Copying PET IMAGE');
            cp_cmd = ['cp -rf ',data_parent_directory, '/', sbj,'/',sess,'/',pet_tr,'/' ,fname.name,...
                ' ' workfolder,'/',pet_tr,'/', sess,'/',sbj,'/PET/nii/pet_',sbj,'.nii'];
            disp(cp_cmd)
            system(cp_cmd); disp("RUNNING MAGIA, it takes around 10 minutes....")
    
            run_magia(sbj,specs,modeling_options)
    
            %  z_cmd=['gunzip -f ' workfolder,'/',pet_tr,'/', sess,'/',sbj,'/PET/nii/pet_',sbj,'.nii.gz']
            %   disp(z_cmd)
            % system(z_cmd);

        else
            disp(['There are more or less than 1 NII PET image...Check subject folder: ',...
            data_parent_directory, '/', sbj,'/',sess,'/',pet_tr,'/']);
    
        end
    
        %md_op=magia_write_modeling_options(model_options_file)
        %run_magia(sbj,specs,modeling_options)
    
    
    end

end