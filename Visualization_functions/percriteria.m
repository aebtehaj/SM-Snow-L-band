% Emission Model
% Divya Kumawat, 09/2022
%% Calculates the error metrics

function [R2,NSE,RMSE,Bias,ubRMSE] = percriteria(sim,obs)
 R2=sum((obs(:)-mean(obs(:))).*((sim(:)-mean(sim(:)))))^2/(sum((obs(:)-mean(obs(:))).^2)*sum((sim(:)-mean(sim(:))).^2));
 NSE=1 - sum((sim(:)-obs(:)).^2)/sum((obs(:)-mean(obs(:))).^2);
 RMSE=(1/length(obs)*sum((sim(:)-obs(:)).^2))^0.5;
 Bias=sum(sim-obs)/length(obs);
 ubRMSE = sqrt(((1/length(obs))*sum((sim(:)-obs(:)).^2))-Bias.^2);
end
