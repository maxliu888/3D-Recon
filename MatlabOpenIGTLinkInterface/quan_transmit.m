clear all
close all
s = serialport("COM5", 115200,"ByteOrder","big-endian");
% figure

igtlConnection = igtlConnect('127.0.0.1',18944);
transform.name = 'IMU';clear all
close all
s = serialport("COM5", 115200,"ByteOrder","big-endian");
% figure

igtlConnection = igtlConnect('127.0.0.1',18944);
transform.name = 'IMU';

startTime = igtlTimestampNow();

%     for t = 1:10000
while(1)
%   t=igtlTimestampNow()-startTime;

    data = read(s,6,"int16");
    hex_data = dec2hex(data);
    ST = hex_data(1,:);
    a = hex_data(2,:);
    b = hex_data(3,:);
    c = hex_data(4,:);
    d = hex_data(5,:);

    CKS = hex_data(6,:);
    
    A = hex2dec(strcat('0x',a,'s16'));
    B = hex2dec(strcat('0x',b,'s16'));
    C = hex2dec(strcat('0x',c,'s16'));
    D = hex2dec(strcat('0x',d,'s16'));

    rot_ang = quat2rotm([A B C D]);
    t_matrix = zeros(4,4);
    t_matrix(4,4) = 1;
    t_matrix(1:3,1:3) = rot_ang;
    transform.matrix = t_matrix;
    a = strcat(string(A/1),' / ', string(B/1),' / ', string(C/1),' / ', string(D/1)) ;
    disp(a)





%   transform.matrix = [ 1 0 0 12+30*sin(t*0.5); 0 1 0 -5; 0 0 1 20; 0 0 0 1 ];




  transform.timestamp = igtlTimestampNow();
  transform
  igtlSendTransform(igtlConnection, transform);
%   pause(0.1)
end

igtlDisconnect(igtlConnection);







startTime = igtlTimestampNow();

%     for t = 1:10000
while(1)
%   t=igtlTimestampNow()-startTime;

    data = read(s,6,"int16");
    hex_data = dec2hex(data);
    ST = hex_data(1,:);
    a = hex_data(2,:);
    b = hex_data(3,:);
    c = hex_data(4,:);
    d = hex_data(5,:);

    CKS = hex_data(6,:);
    
    A = hex2dec(strcat('0x',a,'s16'));
    B = hex2dec(strcat('0x',b,'s16'));
    C = hex2dec(strcat('0x',c,'s16'));
    D = hex2dec(strcat('0x',d,'s16'));

    rot_ang = quat2rotm([A B C D]);
    t_matrix = zeros(4,4);
    t_matrix(4,4) = 1;
    t_matrix(1:3,1:3) = rot_ang;
    transform.matrix = t_matrix;
    a = strcat(string(A/1),' / ', string(B/1),' / ', string(C/1),' / ', string(D/1)) ;
    disp(a)





%   transform.matrix = [ 1 0 0 12+30*sin(t*0.5); 0 1 0 -5; 0 0 1 20; 0 0 0 1 ];




  transform.timestamp = igtlTimestampNow();
  transform
  igtlSendTransform(igtlConnection, transform);
%   pause(0.1)
end

igtlDisconnect(igtlConnection);





