
####################################
### Data loading/manipulating funcs
####################################
load_data <- function(data_name, file=NULL){
    if (data_name == "bh"){
        data(BostonHousing)
        bh <- BostonHousing
        # print("DATASET SUMMARY:")
        # summary(bh)

        # set dependent variable
        y = bh$medv
        # set predictor variables
        X = bh[, !names(bh) %in% c("medv", "chas")]
        # Scale input
        y = scale(y)
        X = scale(X)
        data = NULL
        cat(sprintf("X dims: [%i, %i], y dims: [%i, %i]",
                    dim(X)[1], dim(X)[2], dim(y)[1], dim(y)[2]))
    } else if (data_name == "bb_syndrome"){
        data(eyedata)
        X = x
        # Scale input
        y = scale(y)
        X = scale(X)
        data = NULL

    } else if (data_name == "synthetic"){
        stopifnot(!is.null(file))
        load(file)
        data = get(ls()[1])
        cat(sprintf("Num datasets: %d\n\t X dims: [%d, %d], y dims: [%d, %d]", length(data), dim(data[[1]][["X"]])[1], 
                    dim(data[[1]][["X"]])[2], dim(data[[1]][["y"]])[1], dim(data[[1]][["y"]])[2]), "\n")
        X = NULL
        y = NULL
    } 
    if (!is.null(data)){
        cat(sprintf("Num datasets: %d\n\t", length(data)))
    }
    cat(sprintf("X dims: [%i, %i], y dims: [%i, %i]",
                dim(X)[1], dim(X)[2], dim(y)[1], dim(y)[2]), "\n")
    return (list(X=X, y=y, data=data, dataset_name=data_name))
}


split_data <- function(data, split_pct, i){
    if (data$dataset_name == "synthetic"){
        # use the i^th dataset
        X = data[["data"]][[i]][["X"]]
        y = data[["data"]][[i]][["y"]]
        y_train = scale(y)
        X_train = scale(X)
        y_test = NULL
        X_test = NULL
        theta = data[["data"]][[i]][["theta"]] 
    }else{
        X = data$X
        y = data$y
        theta = NULL
        n = dim(X)[1]
        test_size = round(split_pct*n, digits=0)
        test_idxs = sample(1:n, test_size, replace=FALSE)

        X_test = X[test_idxs,] 
        X_train = X[!(1:n %in% test_idxs),]
        y_test = as.matrix(y[test_idxs])
        y_train = as.matrix(y[!(1:n %in% test_idxs)])
    }

    return (list(X_train=X_train, X_test=X_test, y_train=y_train, y_test=y_test, theta=theta))
}

add_implicit_bias <- function(X){
    X = cbind( as.matrix(rep(1, dim(X)[1])), X)
    return (X)
}


##############################
### metric computing funcs
##############################
mse_func <- function(actual, prediction, normalizer=NULL){
    stopifnot(is.matrix(actual) && dim(actual)[2] == 1)
    stopifnot(is.matrix(prediction) && dim(prediction)[2] == 1)
    if (is.null(normalizer)){
        mse = colMeans((actual-prediction)**2)
        mse_sd = sqrt(colVars(actual-prediction))
        mse_se = mse_sd/sqrt(dim(actual)[1])
        return (list(mse=mse, mse_sd=mse_sd, mse_se=mse_se))
    }else{
        mse = colSums((actual-prediction)**2)/normalizer
        return (mse)
    }
}

cos_sim_func <- function(actual, prediction){
    stopifnot(is.matrix(actual) && dim(actual)[2] == 1)
    stopifnot(is.matrix(prediction) && dim(prediction)[2] == 1)
    cos_sim = t(actual) %*% prediction 
    denom = (norm(actual, type="2")*norm(prediction, type="2"))
    if (denom ==0.0){
        return (0.0)
    }else{
        return (cos_sim/denom)
    }
}

mcc_func <- function(tp, tn, fp, fn){
    denom = (tp+fp)*(tp + fn)*(tn + fp)*(tn + fn)
    if (denom == 0 || is.na(denom)){
        return (0)
    }else{
        denom = sqrt(denom)
    }
    return ((tp * tn - fp *fn)/denom)
}

fp_func <- function(tn, fp){
    denom = fp+tn
    if (denom == 0) return (0)
    return (fp/denom)

}


get_tfpn_info <- function(actual, prediction){
    stopifnot(is.matrix(actual) && dim(actual)[2] == 1)
    stopifnot(is.matrix(prediction) && dim(prediction)[2] == 1)
    actual_nz = abs(actual) > .00000001
    prediction_nz = abs(prediction) > .00000001
    tp = sum(actual_nz & prediction_nz)
    tn = sum((!actual_nz) & (!prediction_nz))
    fp = sum((!actual_nz) & (prediction_nz))
    fn = sum((actual_nz) & (!prediction_nz))
    return (list(tp=tp, fp=fp, tn=tn, fn=fn))
}


ess_func <- function(beta_hat, mcmc_samps){
    "
        mcmc_samps -- type: matrix, dims: [num_samples, num_coefficients]
    "
    beta_vars = apply(mcmc_samps, 2, var)
    most_sig_idxs = order(beta_vars, decreasing=TRUE)[1:10]
    most_sig_coefs = mcmc_samps[,most_sig_idxs,drop=FALSE]
    n = dim(mcmc_samps)[1]
    accum = 0
    for (i in most_sig_idxs){
        ret = acf(mcmc_samps[,i], plot=FALSE)#, lag.max=(dim(most_sig_coefs)[1]-2))
        denom = 1+ 2*sum(ret$acf[,1,1])
        accum = accum + n/denom
    }
    ess = accum/10
    return (ess)
}

##############################
### eval computing funcs
##############################
non_synthetic_eval <- function(X, beta_hat, y, stats_record){
    # add implicit bias for evaluation
    if (stats_record[["dataset_name"]] != "synthetic"){
        X = add_implicit_bias(X)
    }
    y_hat = X %*% beta_hat
    mspe_info = mse_func(y, y_hat)
    cos_sim = cos_sim_func(y, y_hat)
    model_size = length(beta_hat[abs(beta_hat) > .00000001])
    if (is.null(stats_record[["mspe"]])){
        stats_record = list(mspe=list(as.numeric(mspe_info$mse)), mspe_sd=list(as.numeric(mspe_info$mse_sd)),
                            mspe_se=list(as.numeric(mspe_info$mse_se)),cos_sim=list(as.numeric(cos_sim)), 
                            model_size=list(as.numeric(model_size)), dataset_name=stats_record[["dataset_name"]], 
                            true_theta=NULL)
    }else{
        stats_record[["mspe"]] = append(stats_record[["mspe"]], as.numeric(mspe_info$mse))
        stats_record[["mspe_sd"]] = append(stats_record[["mspe_sd"]], as.numeric(mspe_info$mse_sd))
        stats_record[["mspe_se"]] = append(stats_record[["mspe_se"]], as.numeric(mspe_info$mse_se))
        stats_record[["cos_sim"]] = append(stats_record[["cos_sim"]], as.numeric(cos_sim))
        stats_record[["model_size"]] = append(stats_record[["model_size"]], as.numeric(model_size))
    }
    return (stats_record)
}

synthetic_eval <- function(X, beta_hat, y, stats_record, true_beta, mcmc_samps=NULL){
    # "
    #     mcmc_samps -- type: matrix, dims: [num_coefficients, num_samples]
    # "
    # remove intercept index. drop=FALSE keeps it as single column matrix
    beta_hat = beta_hat[2:dim(beta_hat)[1],,drop=FALSE]
    tfpn_info = get_tfpn_info(true_beta, beta_hat)
    mcc = mcc_func(tfpn_info$tp, tfpn_info$tn, tfpn_info$fp, tfpn_info$fn)
    fp = tfpn_info$fp
    true_p = dim(true_beta)[1]/10
    mse = mse_func(true_beta, beta_hat, normalizer=1)
    cos_sim = cos_sim_func(true_beta, beta_hat)
    model_size = length(beta_hat[abs(beta_hat) > .00000001])
    if (!is.null(mcmc_samps)){
        ess = ess_func(beta_hat, mcmc_samps)
    }else{
        ess = NULL   
    }
    if (is.null(stats_record[["mcc"]])){
        stats_record = list(dataset_name=stats_record[["dataset_name"]], mcc=mcc, mse=mse,
                            cos_sim=cos_sim, fp=fp, ess=ess, model_size=model_size)
    }else{
        stats_record[["mcc"]]     = append(stats_record[["mcc"]], as.numeric(mcc))
        stats_record[["mse"]]     = append(stats_record[["mse"]], as.numeric(mse))
        stats_record[["cos_sim"]] = append(stats_record[["cos_sim"]], as.numeric(cos_sim))
        stats_record[["fp"]]      = append(stats_record[["fp"]], as.numeric(fp))
        stats_record[["ess"]]     = append(stats_record[["ess"]], as.numeric(ess))
        stats_record[["model_size"]] = append(stats_record[["model_size"]], as.numeric(model_size))
    }
    return (stats_record)
}

my_eval <- function(X, beta_hat, y, stats_record, true_beta=NULLk, mcmc_samps=NULL){
    if (stats_record[["dataset_name"]] == "synthetic"){
        return (synthetic_eval(X, beta_hat, y, stats_record, true_beta, mcmc_samps))
    } else{
        return (non_synthetic_eval(X, beta_hat, y, stats_record))
    }
}

####################################
### Organizing predictions
####################################
get_lasso_beta_hat <- function(lasso_fitted_results, lambda_cv){
    # includes bias param by default
    lasso_beta_hat = as.matrix(coef(lasso_fitted_results, s=lambda_cv))#[2:13]
    # lasso_beta_hat dims: [p+1, 1]  bc of intercept
    return (lasso_beta_hat)
}

get_blasso_beta_hat <- function(fitted_results, burnin, dataset_name, X_train, y_train){
    betas = fitted_results$beta
    icepts = fitted_results$mu

    # Concat the intercept to betas
    bl_icept = as.matrix(mean(icepts[burnin:length(icepts)[1]]))
    betas = betas[burnin:dim(betas)[1], ]
    bl_beta_hat = rbind(bl_icept, as.matrix(colMeans(betas)))

    # ###  Set threshold for betas 
    X_train_bias = add_implicit_bias(X_train)
    # y_train used to calculate error variance
    beta_hat = set_zeros(X_train_bias, bl_beta_hat, y_train, fitted_results$beta)
    return (beta_hat)
}

get_threshold <- function(X_train, beta_hat, y_train, mcmc_samps){
    train_y_hat = X_train %*% beta_hat
    train_error_var = colVars(train_y_hat) 
    # train_error_var = mean(colVars(mcmc_samps))
    thresh = .1* sqrt(train_error_var)
    return (thresh)
}

set_zeros <- function(X_train, beta_hat, y_train, mcmc_samps){
    thresh = get_threshold(X_train, beta_hat, y_train, mcmc_samps)
    # Set coefficients below threshold to 0.0
    beta_mags = abs(beta_hat)
    beta_hat[beta_mags < thresh] = 0.0
    return (beta_hat)
}



####################################
### Display funcs
####################################
display_stats <- function(stats){
    if (stats[["dataset_name"]] == "synthetic"){
        synth_display_stats(stats)
    }else{
        non_synth_display_stats(stats)
    }
}

display_stat <- function(stat_vec, name, stat_sd=NULL, stat_se=NULL, latex=FALSE){
    stat_mean = mean(stat_vec)
    if (is.null(stat_sd)){
        stat_sd = sd(stat_vec)
    }
    if (is.null(stat_se)){
        stat_se = stat_sd/sqrt(length(stat_vec))
    }
    if (!latex){
        cat(sprintf("%-15s %10.4f  (sd: %7.4f, se: %7.4f)\n", paste(name, "mean:"),  stat_mean, stat_sd, stat_se))
    }else{
        cat(sprintf("&%10.3f(%.3f)", stat_mean, stat_se))
    }
}
non_synth_display_stats <- function(stats){
    display_stat(unlist(stats[["mspe"]]), "mspe")#, stat_sd=mean(unlist(stats[["mspe_sd"]])), stat_se=mean(unlist(stats[["mspe_se"]])))
    # display_stat(unlist(stats[["mspe_sd"]]), "mspe_sd")
    # display_stat(unlist(stats[["mspe_se"]]), "mspe_se")
    display_stat(unlist(stats[["cos_sim"]]), "cos_sim")
    display_stat(unlist(stats[["model_size"]]), "ms")
    display_stat(unlist(stats[["mspe"]]), "mspe", latex=TRUE)
    display_stat(unlist(stats[["cos_sim"]]), "cos_sim", latex=TRUE)
    display_stat(unlist(stats[["model_size"]]), "ms", latex=TRUE)
    cat("\\\\\n")
}

synth_display_stats <- function(stats){
    display_stat(unlist(stats[["mse"]]), "mse")
    display_stat(unlist(stats[["cos_sim"]]), "cos_sim")
    display_stat(unlist(stats[["mcc"]]), "mcc")
    display_stat(unlist(stats[["fp"]]), "fp")
    display_stat(unlist(stats[["ess"]]), "ess")
    display_stat(unlist(stats[["model_size"]]), "ms")
    display_stat(unlist(stats[["mse"]]), "mse", latex=TRUE)
    display_stat(unlist(stats[["cos_sim"]]), "cos_sim", latex=TRUE)
    display_stat(unlist(stats[["mcc"]]), "mcc", latex=TRUE)
    display_stat(unlist(stats[["fp"]]), "fp", latex=TRUE)
    display_stat(unlist(stats[["ess"]]), "ess", latex=TRUE)
    display_stat(unlist(stats[["model_size"]]), "ms", latex=TRUE)
    cat("\\\\\n")
}


