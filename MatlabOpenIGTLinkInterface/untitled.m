%--------------------------------------------------------------------------
% Generating Planewave rfdata  
% Transducer L7-4 (modified), Number of elements 128
%--------------------------------------------------------------------------
%-- Selecting Phantom type 
Ph_type_name = 'Target_SmallS_phantom_';
% Ph_type_name = 'MediumS_phantom_';
% Ph_type_name = 'LargeS_phantom_';
%--------------------------------------------------------------------------

S_ind = 301;                  %-- Starting index
E_ind = 600;                  %-- Ending index
%--------------------------------------------------------------------------

for num = S_ind:E_ind
    f0=5.208e6;             %  Transducer center frequency [Hz]
    fs=f0*4;                %  Sampling frequency [Hz]
    c=1540;                 %  Speed of sound [m/s]
    lambda=c/f0;            %  Wavelength [m]
    width=25/100000;        %  Width of element
    element_height=7/1000;  %  Height of element [m]
    kerf=0.048/1000;        %  Kerf [m]
    focus=[0 0 3e4]/1000;   %  Fixed focal point [m]
    N_elements=64;         %  Number of pc active elements 
    transducer_width=N_elements*width+((N_elements-1)*kerf);

    %--  Set the sampling frequency
    set_sampling(fs);
    %--  Generate aperture for emission
    xmit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1,focus);
    %--  Set the impulse response and excitation of the xmit aperture
    impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
    impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
    xdc_impulse (xmit_aperture, impulse_response);
    excitation=sin(2*pi*f0*(0:1/fs:2/f0));
    xdc_excitation (xmit_aperture, excitation);

    %--  Generate aperture for reception
    receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1,focus);
    %--  Set the impulse response for the receive aperture

    xdc_impulse (receive_aperture, impulse_response);
    no_lines=1;


    %-- Load Phantom 
    Phantom_name = [Ph_type_name,num2str(num),'.mat']
    load(Phantom_name)

    for k=[1:no_lines]

    %--   Calculate the received response
    [rf_data, tstart]=calc_scat_multi(xmit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);

    %--   Add data
    no_added_sampling=tstart*fs;
    added_data=zeros(round(no_added_sampling),N_elements);
    final_data=[added_data' rf_data'];


    save_folder = ['PlanewaveRF'];
    if ~exist(save_folder, 'dir')
        mkdir(save_folder)
    end

    %--   Save data
    pw_rf(1:max(size(final_data)),1:N_elements)=final_data';
    name_RF = append(save_folder,'/PWRF_',num2str(num),'.mat');

    save(name_RF,'pw_rf');
    pw_rf = zeros(size(pw_rf));
    final_data = zeros(size(final_data));
    rf_data = zeros(size(rf_data));
        
end
end
% n3_GAIN
%%
for i = 1:300
    disp(i)
    load(['PWRF_',num2str(i+300)])
    RF(i,1:1100,:) = pw_rf(1:1100,:);
end

%%
for i = 1:300
    disp(i)
    load(['Target_SmallS_phantom_',num2str(i+300)])
    Target(i,1:1100,:) = envelope_beamformed_data(1:1100,:);
end
%%
for i = 1:300
    disp(i)
    rf = squeeze(RF(i,:,:));
    tar = squeeze(Target(i,:,:));
    
    subplot 121
    log_com_plot(rf,1,1080);
    subplot 122
    log_com_plot(tar,1,1080);
    title(num2str(i))
    pause
end



