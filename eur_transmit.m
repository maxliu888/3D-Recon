clear all
close all
s = serialport("COM5", 115200,"ByteOrder","big-endian");
% figure

igtlConnection = igtlConnect('127.0.0.1',18944);
transform.name = 'IMU';

startTime = igtlTimestampNow();

%     for t = 1:10000
while(1)
%   t=igtlTimestampNow()-startTime;

    data = read(s,5,"int16");
    hex_data = dec2hex(data);
    ST = hex_data(1,:);
    roll_hex = hex_data(2,:);
    pitch_hex = hex_data(3,:);
    yaw_hex = hex_data(4,:);
    CKS = hex_data(5,:);
    
    roll = hex2dec(strcat('0x',roll_hex,'s16'))/100;
    pitch = hex2dec(strcat('0x',pitch_hex,'s16'))/100;
    yaw = hex2dec(strcat('0x',yaw_hex,'s16'))/100;
    rot_ang = eul2rotm(deg2rad([roll pitch yaw]));
    t_matrix = zeros(4,4);
    t_matrix(4,4) = 1;
    t_matrix(1:3,1:3) = rot_ang;
    transform.matrix = t_matrix;
    a = strcat(string(roll/1),' / ', string(pitch/1),' / ', string(yaw/1)) ;
    disp(a)





%   transform.matrix = [ 1 0 0 12+30*sin(t*0.5); 0 1 0 -5; 0 0 1 20; 0 0 0 1 ];




  transform.timestamp = igtlTimestampNow();
  transform
  igtlSendTransform(igtlConnection, transform);
%   pause(0.1)
end

igtlDisconnect(igtlConnection);



% 
% while(1)
%     n = n + 1;
%     data = read(s,8,"int16");
%     hex_data = dec2hex(data);
%     ST = hex_data(1,:);
%     roll_hex = hex_data(2,:);
%     pitch_hex = hex_data(3,:);
%     yaw_hex = hex_data(4,:);
%     dx_hex = hex_data(5,:);
%     dy_hex = hex_data(6,:);
%     dz_hex = hex_data(7,:);
%     CKS = hex_data(8,:);
% 
%     roll = hex2dec(strcat('0x',roll_hex,'s16'));
%     pitch = hex2dec(strcat('0x',pitch_hex,'s16'));
%     yaw = hex2dec(strcat('0x',yaw_hex,'s16'));
% %     dx = hex2dec(strcat('0x',dx_hex,'s16'));
% %     dy = hex2dec(strcat('0x',dy_hex,'s16'));
% %     dz = hex2dec(strcat('0x',dz_hex,'s16'));
% 
%     ROLL(n) = roll/100;
%     PITCH(n) = pitch/100;
%     YAW(n) = yaw/100;
% %     DX(n) = dx;
% %     DY(n) = dy;
% %     DZ(n) = dz;
%     N(n) = n;
% 
% %     scatter3(DX,DY,DZ)
% %     scatter(N,DX)
%     
% 
% % %%% RPY--------------------
% %     subplot 311
% %     scatter(N(n), ROLL(n)/100,'.')
% %     title('ROLL')
% %     hold on
% % 
% %     subplot 312
% %     scatter(N(n), PITCH(n)/100,'.')
% %     title('PITCH')
% %     hold on
% %     
% %     subplot 313
% %     scatter(N(n), YAW(n)/100,'.')
% %     title('YAW')
% %     hold on
% % %%%--------------------
% %%%--------------------
% 
% % %%% XYZ--------------------
% %     
% %     subplot 311
% %     scatter(n, dx/1000,'.')
% %     title('X')
% %     hold on
% % 
% %     subplot 312
% %     scatter(n, dy/1000,'.')
% %     title('Y')
% %     hold on
% %     
% %     subplot 313
% %     scatter(n, dz/1000,'.')
% %     title('Z')
% %     hold on
% %%%--------------------
% %%%--------------------
% 
% 
% %%% RPY--------------------
%     subplot 311
%     scatter(N(n), ROLL(n)/100)
%     title('ROLL')
%     hold on
% 
%     subplot 312
%     scatter(N(n), PITCH(n)/100)
%     title('PITCH')
%     hold on
%     
%     subplot 313
%     scatter(N(n), YAW(n)/100)
%     title('YAW')
%     hold on
% %%%--------------------
% % %%%--------------------
% % 
% % %%% XYZ--------------------
% %     
% %     subplot 322
% %     scatter(N(n), DX(n))
% %     title('X')
% %     hold on
% % 
% %     subplot 324
% %     scatter(N(n), DY(n))
% %     title('Y')
% %     hold on
% %     
% %     subplot 326
% %     scatter(N(n), DZ(n))
% %     title('Z')
% %     hold on
% % %%%--------------------
% % %%%--------------------
% 
% %     a = strcat(string(roll/100),' / ', string(pitch/100),' / ', string(yaw/100),' / ', string(dx/100),' / ', string(dy/100),' / ', string(dz/100),' / ') ;
% %     disp(a)
%  a = strcat(string(roll/100),' / ', string(pitch/100),' / ', string(yaw/100)) ;
%     disp(a)
% 
% %         a = strcat(string(roll_hex),' / ', string(pitch_hex),' / ', string(yaw_hex), ' / ', string(dx_hex),' / ', string(dy_hex),' / ', string(dz_hex),' / ') ;
% %     disp(a)
% 
% end
% %%


