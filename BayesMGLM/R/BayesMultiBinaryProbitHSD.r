#library(MGLM)
#library(reshape)
#library(Matrix)
#rm(list=ls(all=TRUE))

#load("/Users/kjlee/Research/KJLEE_Papers/MultivariateLongitudinalModelSelection/MultivariateProbitModel/Data/KoGES.rda")

#source("/Users/kjlee/Research/KJLEE_Papers/MultivariateLongitudinalModelSelection/MultivariateProbitModel/Code/Rpackage/MGLM/R/HSModelDesignMatrix.r")


BayesMultiBinaryProbitHSD = function(fixed, data, random, HS.model.within, HS.model.cross, hyper.params, num.of.iter)
{

# process data: reponse, fixed and random effects matrices. 

    cl = match.call()
    mf = match.call(expand.dots = FALSE)  

    #aa = call("BayesMultiBinaryProbitHSD", fixed = hyp+hdlc~drink+smoke+age, data = KoGES, random = ~1) #

    #match.call(aa)

    #cat("1 mf", length(mf), "\n")

    #print(mf)
    #print(attributes(mf))
    #for(i in 1:length(mf))
    #print(mf[i])



    m = match(c("fixed", "data", "subset", "na.action"), names(mf), 0L)

    #cat("m = ", m, "\n")

    mf = mf[c(1L, m)]

    #cat("2 mf\n")
    #print(mf)



    mf$drop.unused.levels = TRUE
    mf[[1L]] = quote(model.frame)


    names(mf)[2] = "formula"

    #cat("3 mf\n")
    #for(i in 1:length(mf))
    #print(mf[i])


    mf2 = eval(mf, parent.frame())

    yy = model.response(mf2, "numeric")

    fixed.eff = all.vars.character(fixed[-2])$m[[2]]

    #cat("fixed.eff = ", fixed.eff, "\n")
    #fixed.eff.intercept.included = !any(grepl("-1", fixed.eff))
    
    random.eff = all.vars.character(random)$m[[2]]

    #cat("random.eff = ", random.eff, "\n")

    #HS.model.cov = all.vars.character(HS.model)$m[[2]]
    #cat("HS.model.cov = ", HS.model.cov, "\n")



    #print(head(yy))

    Terms = attr(mf2, "terms")
    fixed.eff = colnames(model.matrix(Terms, mf2))


    #fixed.eff = fixed.eff[-1]


    mf[[2L]] = update(fixed, as.formula(paste("~.+", paste(random.eff, collapse="+") )))
    
    #cat("4 mf\n")
    #print(mf[[2L]])

    mf[[2L]] = update(mf[[2L]], ~.+id)

    #cat("5 mf\n")
    #print(mf[[2L]])

    mf  = eval(mf, parent.frame())

    m.design.mat = attr(mf, "terms")

    #cat("mfixed.design.mat = \n")
    #print(mfixed.design.mat)

    yy = model.response(mf, "numeric") #model.response(mf, "numeric")
    xx = model.matrix(m.design.mat, mf)

    
    #yy = model.response(mf, "numeric") #model.response(mf, "numeric")

    #cat("yy = \n")
    #print(head(yy))
    #print(head(xx))

    #print(rownames(yy))

    #xx <- model.matrix(m.design.mat, mf)
    random.eff[random.eff=="1"] = "(Intercept)"

    #cat("random.eff = ", random.eff, "\n")

    x.fixed = xx[, colnames(xx)%in%fixed.eff, drop=FALSE]


    z.random = xx[, colnames(xx)%in%random.eff, drop=FALSE]

#==========================================================
    #HS.model.cov = all.vars.character(HS.model)$m[[2]]
    #cat("HS.model.cov = ", HS.model.cov, "\n")

    #TimeOrder = sort(gsub("IndTime", "", HS.model.cov[HS.model.cov %in% paste0("IndTime", 1:10)]))
    #cat("TimeOrder = ", TimeOrder, "\n")

    #DiffTime = sort(gsub("DiffTime", "", HS.model.cov[HS.model.cov %in% paste0("DiffTime", 1:10)]))
    #cat("DiffTime = ", DiffTime, "\n")
    
    #cat("HS.model = \n")
    #print(as.formula(HS.model))

    #interaction.terms = attr(terms.formula(as.formula(HS.model)), "term.labels")
   
#==========================================================



    id = xx[, colnames(xx)%in%"id"]

    P = dim(x.fixed)[2]
    Q = dim(z.random)[2]
    N =length(table(id))
    T = range(table(id))[2]
    K = ncol(yy)

    #cat("T = ", T, "\n")

    TimePointsAvailable = as.vector(table(id))

    sub.id = as.numeric(names(table(id)))

 
        #id = as.numeric(names(table(KoGES$id)))
        #id = as.numeric(names(table(KoGES$id)))
    Y.Binary = matrix(NA, T*K, N) #NULL #y((subject,time), attributes))
    for(i in 1:N)
        Y.Binary[1:(K*TimePointsAvailable[i]), i]  = matrix(t(yy[id==sub.id[i], ]))

    X.all = array(NA, c(N, T*K, K*P))
    Z.all = array(NA, c(N, T*K, K*Q))


    for(i in 1:N){
        X.all.tmp = NULL
        x.tmp = x.fixed[id==sub.id[i], ,drop=F]
        for(t in 1:(dim(x.tmp)[1]))
            X.all[i, ((t-1)*K+1):(K*t), ] = diag(K)%x%x.tmp[t, , drop=F]
        #X.all[i,  , ] = X.all.tmp
    }

    #print(X.all[2000, , ])

    for(i in 1:N){
        Z.all.tmp = NULL
        z.tmp = z.random[id==sub.id[i], ,drop=F]
        for(t in 1:(dim(z.tmp)[1]))
            Z.all[i, ((t-1)*K+1):(K*t), ] = rbind(Z.all.tmp, diag(K)%x%z.tmp[t, , drop=F])
        #X.all[i,  , ] = X.all.tmp
    }
    #X.all.out <<-X.all
#===========================================================================#
    
    HS.model.cov.within = all.vars.character(HS.model.within)$m[[2]]
    HS.model.cov.cross = all.vars.character(HS.model.cross)$m[[2]]


    if(all(HS.model.cov.within ==0)){
        UpdateAlpha = FALSE
        HS.model.Within = NULL
        alpha.ini = NULL
        tuning.alpha = NULL
        sigma2.alpha = NULL
        
        
    }
    else{
        HS.model.Within = HSModelDesignMatrix(HS.model= HS.model.within, Num.Of.Subjects = N, Time.Points = T)
        UpdateAlpha = TRUE
        S = dim(HS.model.Within)[1]/N
        #cat("dim(HS.model.Within)\n")
        #print(dim(HS.model.Within))
        alpha.ini = rep(0, S)
        tuning.alpha = 0.01
        sigma2.alpha = 1
        
    }

    if(all(HS.model.cov.cross ==0)){
        UpdateDelta  = FALSE
        HS.model.Cross = NULL
        delta.ini = NULL
        tuning.delta = NULL
        sigma2.delta = NULL
        
    }
    else if(HS.model.cov.cross==2){
        UpdateDelta = TRUE

        u = array(0, c(K, K, N, K))
        for(i in 1:N)
            for(k in 1:K){
                u[k, , i, k] = 1
                u[, k, i, k] = 1
                u[k, k, i, k] = 1
            }

        HS.model.Cross = aperm(u, c(4, 3, 1, 2))
        #print(dim(u))
        dim(HS.model.Cross) = c(N*K, K, K)

        A = K
        #cat("dim(HS.model.Cross)\n")
        #print(dim(HS.model.Cross))
        delta.ini = rep(0, A)

        tuning.delta = 0.1
        sigma2.delta = 1

    }
    else{       
        HS.model.Cross = HSModelDesignMatrix(HS.model = HS.model.cross, Num.Of.Subjects = N, Time.Points = K)
        UpdateDelta = TRUE

        A = dim(HS.model.Cross)[1]/N
        #cat("dim(HS.model.Cross)\n")
        #print(dim(HS.model.Cross))
        delta.ini = rep(0, A)

        tuning.delta = 0.1
        sigma2.delta = 1
        
    }




    #InitialValues = list(beta.ini = beta.ini, b.ini = bi.ini, Sigma.ini = Sigma.ini, alpha.ini = alpha.ini, delta.ini = delta.ini)
    #HyperPara = list(sigma2.beta = sigma2.beta, Vb = Vb, sigma2.alpha = sigma2.alpha, sigma2.delta = sigma2.delta, Lambda = Lambda)
    #UpdatePara = list(UpdateYstar = UpdateYstar, UpdateBeta=UpdateBeta, Updateb =Updateb, UpdateSigma= UpdateSigma, UpdateDelta=UpdateDelta, UpdateAlpha=UpdateAlpha)
    #TuningPara = list(tuning.alpha= tuning.alpha, tuning.delta = tuning.delta)
        


 #==============================================================#
    UpdateYstar = TRUE
    UpdateBeta = TRUE
    Updateb = TRUE
    UpdateSigma = TRUE
    #UpdateAlpha = TRUE
    #UpdateDelta = TRUE

    
    Data = list(Y = Y.Binary, X=X.all, Z=Z.all, TimePointsAvailable = TimePointsAvailable, W=HS.model.Within, V=HS.model.Cross)


    beta.ini = matrix(0, K*P)
    bi.ini = matrix(0, N, K*Q)


    Sigma = rWishart(K,Q,diag(Q)) 

    Sigma.ini =  0.01*as.matrix(do.call("bdiag", lapply(seq(dim(Sigma)[3]), function(x) Sigma[ , , x])))


    #alpha.ini = rep(0, S)
    #delta.ini = rep(0, A)

    #tuning.alpha = 0.1
    #tuning.delta = 0.1

    sigma2.beta = 1
    Lambda = diag(K*Q)
    Vb = 10
    #sigma2.alpha = 1
    #sigma2.delta = 1
    Lambda = diag(K*Q)

    InitialValues = list(beta.ini = beta.ini, b.ini = bi.ini, Sigma.ini = Sigma.ini, alpha.ini = alpha.ini, delta.ini = delta.ini)
    HyperPara = list(sigma2.beta = sigma2.beta, Vb = Vb, sigma2.alpha = sigma2.alpha, sigma2.delta = sigma2.delta, Lambda = Lambda)
    UpdatePara = list(UpdateYstar = UpdateYstar, UpdateBeta=UpdateBeta, Updateb =Updateb, UpdateSigma= UpdateSigma, UpdateDelta=UpdateDelta, UpdateAlpha=UpdateAlpha)
    TuningPara = list(tuning.alpha= tuning.alpha, tuning.delta = tuning.delta)


    start.time <- Sys.time()

    PosteriorSamples = MultiLinearModelMCMC(num.of.iter, Data, InitialValues, HyperPara, UpdatePara, TuningPara)

    end.time <- Sys.time()




    #cat("\nCall:\n", printCall(cl), "\n\n", sep = "")
    cat("\nData Descriptives:\n")
    cat("Longitudinal Data Information:")
    cat("\nNumber of Observations: ", sum(TimePointsAvailable), "\tNumber of Covariates: ", P-1)
    cat("\nNumber of subjects:", N, "\n\n")

    #cat("HSModelDesignMatrix = \n")
    #print(HS.model.Cross[1, , ])
    #print(HS.model.Cross[2, , ])
    #print(HS.model.Cross[3, , ])
    #print(HS.model.Cross[4, , ])
    #print(HS.model.Cross[5, , ])

    #cat("HS.model.Within = \n")
    #print(HS.model.Within[1, , ])
    #print(HS.model.Within[2, , ])

    out <- list(Posterior.Samples = PosteriorSamples, Fixed.Effects.Names = fixed.eff, 
                Random.Effects.Names = random.eff, 
                Response = Y.Binary, Fixed.Effects.Mat = X.all, Random.Effects.Mat = Z.all, 
                HS.model.Within.Mat = HS.model.Within,  HS.model.Cross.Mat = HS.model.Cross, 
                call = cl, Num.of.Iter = num.of.iter)
    
    #class(out)
    out 



}

#BayesMultiBinaryProbitHSD(fixed = cbind(hyp, hdlc)~drink+smoke+age+smet+meta+exf+sex, 
#    data = KoGES, random = ~1,  HS.model.within = ~IndTime1+IndTime2, 
#    HS.model.cross = ~1, num.of.iter = 10)
#"hyp", "hdlc", "tri", "hfg", "abd"

#drink", "smoke", "age", "smet", "meta", "exf", "sex")
