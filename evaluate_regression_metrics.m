function metrics = evaluate_regression_metrics(y_pred,y_true)
    [y_true,y_pred] = prepareCurveData(y_true(:),y_pred(:));   
    error = y_true - y_pred;

    

    fitres=fitlm(y_true,y_pred);    
     R2 = fitres.Rsquared.Ordinary;
    RMSE = sqrt(mean(error.^2));
    RRMSE = RMSE / mean(y_true);
    MAE = mean(abs(error));
    MAPE = mean(abs(error ./ y_true)) * 100;
    Bias = mean(error);
    
    MSE = mean(error.^2);
    NRMSE = RMSE / (max(y_true) - min(y_true));
    SMAPE = mean(2 * abs(error) ./ (abs(y_true) + abs(y_pred))) * 100;
    MSLE = mean((log1p(y_pred) - log1p(y_true)).^2);
    MedAE = median(abs(error));
    explainedVar = 1 - var(error) / var(y_true);    
    metrics = struct();
    metrics.R2 = R2;
    metrics.RMSE = RMSE;
    metrics.RRMSE = RRMSE;
    metrics.MAE = MAE;
    metrics.MAPE = MAPE;
    metrics.Bias = Bias;
    metrics.MSE = MSE;
    metrics.NRMSE = NRMSE;
    metrics.SMAPE = SMAPE;
    metrics.MSLE = MSLE;
    metrics.MedAE = MedAE;
    metrics.ExplainedVariance = explainedVar;
    metrics.fitres=fitres;

  
end
