%--------------------------------------------------------------------------
% Generating multi-angle rfdata and its coorisponding target 
% Transducer L7-4 (modified), Number of elements increase from 128 to 164*2
% Trasnducer's element width, height, and kert are decreased by 2
%--------------------------------------------------------------------------
%-- Selecting Phantom type 
Ph_type_name = 'SmallS_phantom_';
% Ph_type_name = 'MediumS_phantom_';
% Ph_type_name = 'LargeS_phantom_';
%--------------------------------------------------------------------------

S_ind = 451;                  %-- Starting index
E_ind = 600;                  %-- Ending index
%--------------------------------------------------------------------------
    % Attenuation configuration 
%     set_field('Freq_att', 0.5*100/1e6);
%     set_field('use_att',1)


for ph = S_ind:E_ind 
    
    % ********************************************************************
    %-- Generating Rfdata ------------------------
    % ********************************************************************
    
    %-- Load Phantom 
%     Phantom_name = [Ph_type_name,num2str(ph),'.mat']
%     load(Phantom_name)
    [phantom_positions, phantom_amplitudes] = cyst_move(round((rand*100) + 15));
    
    phantom_amplitudes(-0.006<phantom_positions(:,1) & phantom_positions(:,1)<0.006) = [];
    phantom_positions(-0.006<phantom_positions(:,1) & phantom_positions(:,1)<0.006,:) = [];
%     phantom_amplitudes(phantom_positions(:,3)<(20/1000)) =  phantom_amplitudes(phantom_positions(:,3)<(20/1000)) / 5;
    
    %-- Init. simulation parameters 
    f0=5.208e6;             %-- Transducer center frequency [Hz]
    fs=f0*4;                %-- Sampling frequency [Hz]
    c=1540;                 %-- Speed of sound [m/s]
    lambda=c/f0;            %-- Wavelength [m]
    width=(25/1)/100000;    %-- Width of element
    element_height=7/1000;  %-- Height of element [m]
    kerf=(0.048/1)/1000;     %-- Kerf [m]
    focus=[0 0 3e4]/1000;   %-- Fixed focal point [m]
    Rfocus = 25/1000;
    N_elements=128;       %-- Number of pc active elements
    transducer_width=N_elements*width+((N_elements-1)*kerf);    %-- Aperture size
    angle_width = 9;  %-- in degrees 
    N_angles = 100;          %-- Number of Angels



    %-- Set the sampling frequency

    set_sampling(fs);

    %-- Generate aperture for emission
    xmit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1,focus);
%     xmit_aperture = xdc_focused_array (N_elements, width, element_height, kerf, Rfocus, 1, 10,focus);


    %-- Set the impulse response and excitation of the xmit aperture
    impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
    impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
    xdc_impulse (xmit_aperture, impulse_response);
    excitation=sin(2*pi*f0*(0:1/fs:2/f0));
    xdc_excitation (xmit_aperture, excitation);

    %-- Generate aperture for reception
    receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1,focus);
%     receive_aperture = xdc_focused_array (N_elements, width, element_height, kerf, Rfocus, 1, 10,focus);

    %-- Set the impulse response for the receive aperture
    xdc_impulse (receive_aperture, impulse_response);

    %-- plane wave steering 
    lateral_center = transducer_width/2-width/2;
    ele_pos = linspace(-lateral_center, lateral_center, N_elements); % element positions
    steer_ang = linspace(-angle_width, angle_width, N_angles);  % steering angle
    CPWC_RF = zeros(N_angles,2706,N_elements);
    
    for a = 1:N_angles

        fprintf('Phantom:%2d ,Simulating angle = %2d/%2d \n ', ph, a, N_angles);
        %-- delays with angle
        delays = (ele_pos*sind(steer_ang(a)) )/c;               %-- delay in time
        xdc_focus_times (xmit_aperture, 0, delays);             %-- set transmit delays for steering

        %-- Calculate the received response
        [rf_data, tstart]=calc_scat_multi(xmit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);

        %-- Add data
        no_added_sampling=tstart*fs;
        added_data=zeros(round(no_added_sampling),N_elements);
        final_data=[added_data' rf_data'];

        %-- Stack with added data
        CPWC_rf(1:max(size(final_data)),1:N_elements)=final_data';
        CPWC_RF(a,1:length(CPWC_rf),:) = CPWC_rf;
        

    end
    CPWC_RF = single(CPWC_RF);
    
    %-- RF data is saved only if there are more than 20000 scatteres 
%     if Np > 2e4
%         save(['RFcpwc_',Ph_type_name , num2str(ph),'.mat'],'CPWC_RF')
%     end
    
%     if strcmp(Ph_type_name, 'SmallS_phantom_') == 1
%         fprintf('Skip');
%     else 
%         save(['RFcpwc_',Ph_type_name , num2str(ph),'.mat'],'CPWC_RF');
%     end 
%     CPWC_RF = zeros(size(CPWC_RF));
    
    % ********************************************************************
    %-- Beamforming Session ------------------------
    % ********************************************************************

    clc
    fprintf('Generating beamformed target ... \n');

    scan = linear_scan();
    scan.read_file('resolution_distorsion_expe_scan.hdf5');
    dataset = us_dataset();
    dataset.read_file('resolution_distorsion_expe_dataset_rf.hdf5');

    st = single(steer_ang' * (pi/180));
    dataset.angles = st;

    xpos = [(-transducer_width/2)+(width/2):(kerf+width):(transducer_width/2)+(width/2)]';
    a = zeros(N_elements,3);
    a(:,1) = xpos;
    dataset.probe_geometry = single(a);
    %-- Replacing the scan part
    rf_depth = size(dataset.data,1);
    z_delta = dataset.c0/dataset.sampling_frequency/2;
    z_axis_temp = [0:z_delta:rf_depth*z_delta];
    scan.z_axis = z_axis_temp';
    scan.x_axis = dataset.probe_geometry(:,1);

    %-- Beamforming loop and dataloader  
    
    %-- Rearrange and assign rf-data
    CPWC_RF = permute(CPWC_RF,[2,3,1]);
    RFDATA = zeros(3328,N_elements,N_angles);
    n_s = size(CPWC_RF,1);
    RFDATA(1:n_s,:,:) = CPWC_RF;
    dataset.data = single(RFDATA);
    CPWC_RF = zeros(size(CPWC_RF));

    %-- Check whether data is RF or IQ
    if(isempty(dataset.modulation_frequency)||dataset.modulation_frequency==0)
        dataset.data=reshape(hilbert(dataset.data(:,:)),size(dataset.data));
    end
    %
    % ---------------------------------------------
    %-- Define new recon positions
    % ---------------------------------------------
    x_cp = scan.x;
    z_cp = scan.z;

    starting_index = find(x_cp<=-0.0094);
    x_cp(1:starting_index(end)) = [];
    z_cp(1:starting_index(end)) = [];
    ending_index = find(x_cp>=0.0094);
    x_cp(ending_index(1):end) = [];
    z_cp(ending_index(1):end) = [];


%     starting_index = find(z_cp>=0.05);
%     z_cp(1:starting_index(end)) = [];
%     x_cp(1:starting_index(end)) = [];


%     %-- Define value of each x element
%     ele_val = x_cp(1:size(scan.x_matrix,1):end);
%     del_ele_val = ele_val(1:2:end);
%     ind = 0;
%     Index = [];
%     for kk = 1:length(del_ele_val)
%         del_ind = find(x_cp==del_ele_val(kk));
%         Index = [Index;del_ind];
%     end 
%     x_cp(Index) = [];
%     z_cp(Index) = [];
    pixel = size(x_cp,1);
    % ---------------------------------------------

    %-- receive apodization:
    %-- dynamically expanding receive aperture with Tukey 25% apodization
    rx_f_number = 1.6;
    rx_aperture = z_cp/rx_f_number;
    rx_aperture_distance = abs(x_cp*ones(1,dataset.channels)-ones(pixel,1)*dataset.probe_geometry(:,1).');
    receive_apodization = tools.apodization(rx_aperture_distance,rx_aperture*ones(1,dataset.channels),'tukey25');
    receive_apodization = bsxfun(@rdivide,receive_apodization,sum(receive_apodization,2));

    angular_apodization = ones(pixel,dataset.firings);
    angular_apodization  = bsxfun(@rdivide,angular_apodization,sum(angular_apodization,2));

    %-- beamforming loop
    beamformed_data = zeros(pixel,1);
    time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
    w0 = 2*pi*dataset.modulation_frequency;
    tic
    %
    for pw= 1:N_angles
        fprintf('Beamforming..., Angle = %2d/%2d \n ', pw, N_angles);
        %-- transmit delay
        transmit_delay = z_cp*cos(dataset.angles(pw))+x_cp*sin(dataset.angles(pw));
        for nrx=1:dataset.channels
            %-- progress bar
            step=(nrx + (pw-1)*dataset.channels)/length(dataset.angles)/dataset.channels;
            %-- receive delay
            receive_delay = sqrt((dataset.probe_geometry(nrx,1)-x_cp).^2+(dataset.probe_geometry(nrx,3)-z_cp).^2);
            %-- total delay
            delay =  ((transmit_delay+receive_delay)/dataset.c0);
            %-- phase shift
            phase_shift = exp(1i.*w0*(delay-2*z_cp/dataset.c0));

            beamformed_data = beamformed_data + phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
        end
    end
    t = toc
    envelope_beamformed_data = abs(reshape(beamformed_data,[numel(scan.z_axis) 64]));
    fprintf('Saving target... \n ------------------------------------------');
    save(['Target_',Ph_type_name ,num2str(ph),'.mat'],'envelope_beamformed_data','phantom_amplitudes', 'phantom_positions','t')
    
    
    % ---------------------------------------------
    %-- Plot result 
    f = figure;
    f.Position = [50 120 200 500];
    setdB = 60;     
    log_env = 20*log10(abs(envelope_beamformed_data)/max(max(abs(envelope_beamformed_data))));
    log_env(log_env<-setdB) = -setdB;
    imagesc(log_env(1:1300,:)), colormap gray
    title('Target')
end

%%
 for i = 1:100
     load(['Target_SmallS_phantom_',num2str(i)]);
     log_com_plot(envelope_beamformed_data,1,1200);
     pause 
 end 
%%

[phantom_positions, phantom_amplitudes] = cyst_move(round((rand*100) + 100));
phantom_amplitudes(-0.005<phantom_positions(:,1) & phantom_positions(:,1)<0.005) = [];
phantom_positions(-0.005<phantom_positions(:,1) & phantom_positions(:,1)<0.005,:) = [];
scatter(phantom_positions(:,1),phantom_positions(:,3))

