load('Testing_64e_experiment.mat')
Target_data = Target;
load('Target_pw.mat')

load('Result_EXP_sc_396')
Result_cut = results;
load('old/Result_EXP_aug_410_old.mat')
Result_aug_rand = results;

load('old/Result_ft_196_old.mat')
Result_pretrain_simu = results;

load('Result_EXPPretrained_500.mat')
Result_pretrain_exp = results;
%%
z = 1/27.4;
st_depth = 1;
depth = 1080;
v = VideoWriter('video_simu1.avi','Uncompressed AVI')
v.FrameRate = 1;
open(v)
% FrameRate = 3;

for i = 1:76
    
%     v.FrameRate = 3;
    target = squeeze(Target_data(i,1:depth,:));
    target_pw = squeeze(Target_pw(i,1:depth,:));
    result_cut = squeeze(Result_cut(i,1:depth,:));
    result_new = squeeze(Result_aug_rand(i,1:depth,:));
    result_pt_simu = squeeze(Result_pretrain_simu(i,1:depth,:));
    result_pt_exp = squeeze(Result_pretrain_exp(i,1:depth,:));
    rf_data = squeeze(input(i,1:depth,:));

    subplot 171
    log_com_plot(target_pw,st_depth,depth);
    title('Target 64t 64r pw')
    
    subplot 172
    log_com_plot(target,st_depth,depth);
    title('Target 64t 64r')

    subplot 173
    log_com_plot(result_cut, st_depth,depth);
    title('exp more multi')

    subplot 174
    log_com_plot(result_new, st_depth,depth);
    title(['exp rand aug, #:',num2str(i)]) 

    subplot 175
    log_com_plot(result_pt_simu, st_depth,depth);
    title('with simu_ft')

    subplot 176
    log_com_plot(result_pt_exp, st_depth,depth);
    title('with exp pretrained')

    subplot 177
    log_com_plot(rf_data, st_depth,depth);
    title('RF data (log scale)')





    
    pause
    frame = getframe(gcf);
    writeVideo(v,frame)
end 
close(v)

%%

for i = 1:10: 500
    a = squeeze(tar(i,:,:));
    log_com_plot(a,1,1300);
    title(num2str(i))
    pause 
end
%%
log_com_plot(squeeze(Target(1,:,:)),1,1000);

