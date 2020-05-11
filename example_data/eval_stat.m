function statistics = eval_stat(data,out)
% data : Original simulated data
% out : output structure obtained from running VI_BVAR
% statistics: vector of FPR,FNR,Accuracy,F_1
% True Positive
g_1_TP =sum(data.gamma_sim.one(out.q_gamma(1,:)>0.9)==1);
% False Positive
g_1_FP =sum(data.gamma_sim.one(out.q_gamma(1,:)>0.9)==0);
% True Negative
g_1_TN = sum(data.gamma_sim.one(out.q_gamma(1,:)<0.9)==0);
% False Negative
g_1_FN =sum(data.gamma_sim.one(out.q_gamma(1,:)<0.9)==1);
% True Positive
g_2_TP = sum(data.gamma_sim.two(out.q_gamma(2,:)>0.9)==1);
% False Positive
g_2_FP= sum(data.gamma_sim.two(out.q_gamma(2,:)>0.9)==0);
% True Negative
g_2_TN=sum(data.gamma_sim.two(out.q_gamma(2,:)<0.9)==0);
% False Negative
g_2_FN= sum(data.gamma_sim.two(out.q_gamma(2,:)<0.9)==1);
% Now compute statistics
g_1_FPR = g_1_FP/(g_1_FP +g_1_TN);
g_1_FNR = g_1_FN/(g_1_FN +g_1_TP);
g_1_Accuracy = (g_1_TP + g_1_TN)/(g_1_TP + g_1_TN +g_1_FP + g_1_FN);
g_1_F_1 = 2*(g_1_TP/(g_1_TP+g_1_FP)*g_1_TP/(g_1_TP+g_1_FN))/...
    (g_1_TP/(g_1_TP+g_1_FP)+g_1_TP/(g_1_TP+g_1_FN));
g_2_FPR = g_2_FP/(g_2_FP +g_2_TN);
g_2_FNR = g_2_FN/(g_2_FN +g_2_TP);
g_2_Accuracy = (g_2_TP + g_2_TN)/(g_2_TP + g_2_TN +g_2_FP + g_2_FN);
g_2_F_1 = 2*(g_2_TP/(g_2_TP+g_2_FP)*g_2_TP/(g_2_TP+g_2_FN))/...
    (g_2_TP/(g_2_TP+g_2_FP)+g_2_TP/(g_2_TP+g_2_FN));

statistics =[g_1_FPR,g_1_FNR,g_1_Accuracy,g_1_F_1,g_2_FPR,g_2_FNR,g_2_Accuracy,g_2_F_1];
end