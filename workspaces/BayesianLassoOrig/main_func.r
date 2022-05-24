main <- function(dataset_name, dataset_fname=NULL, reps=20){
    # Experiment settings
    nsamps = 10000
    burnin = 2000 + 1
    thinning = 1
    split_pct = .1
    if (dataset_name == "synthetic"){
        fit_intercept = FALSE
    }else{
        fit_intercept = TRUE
        # fit_intercept = FALSE
    }

    # Load data
    source("../my_funcs.R")
    cat("Dataset:", dataset_name,"\n")
    if (dataset_name == "synthetic"){
        cat("\tFilename:", dataset_fname,"\n")
    }
    data = load_data(dataset_name, file=dataset_fname)


    # Init stat records
    bl_stats = list(dataset_name=dataset_name)
    l_stats = list(dataset_name=dataset_name)

    # Repeat experiment and collect stats 
    for (i in seq(1, reps)){
        cat(" Iteration: ", i, "\r")

        # split data
        data_parts = split_data(data, split_pct, i)
        X_train = data_parts$X_train
        X_test = data_parts$X_test
        y_train = data_parts$y_train 
        y_test = data_parts$y_test
        true_beta = data_parts$theta


        # Lasso cv estimate and lambda_cv param
        lasso_fitted_results = cv.glmnet(X_train, y_train, alpha=1, intercept=fit_intercept, lambda=seq(0,.5, length.out=5000))
        lambda_cv = as.numeric(lasso_fitted_results["lambda.min"])

        # Get param estimate 
        lasso_beta_hat =  get_lasso_beta_hat(lasso_fitted_results, lambda_cv)

        # evaluate param estimate (Sythetic datasets don't use X_test/y_test)
        l_stats = my_eval(X_test, lasso_beta_hat, y_test, l_stats, true_beta=true_beta)

        # Bayesian Lasso estimate 
        # "rd"=FALSE causes the lambda2 parameter to fixed at it's initial value
        blasso_fitted_results <- blasso(X_train, y_train, T=nsamps,  thin=thinning, icept=fit_intercept, 
                                        lambda2=(lambda_cv), 
                                        verb=0, 
                                        rd=FALSE,
                                        # RJ=FALSE,
                                        rao.s2=FALSE
                                        ) 
        bl_beta_hat = get_blasso_beta_hat(blasso_fitted_results, burnin, dataset_name, X_train, y_train)
        # bl_beta_hat dims: [p+1, 1]  bc of intercept

        # evaluate param estimate (Sythetic datasets don't use X_test/y_test)
        bl_stats = my_eval(X_test, bl_beta_hat, y_test, bl_stats, mcmc_samps=blasso_fitted_results$beta[burnin:nsamps,,drop=FALSE], true_beta=true_beta)
    }

    # Display results
    cat("\n==========\nLasso\n----------\n")
    display_stats(l_stats)

    cat("\n==========\nBayesian Lasso\n----------\n")
    display_stats(bl_stats)
}
