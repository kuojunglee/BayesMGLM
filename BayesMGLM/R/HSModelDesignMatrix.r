#create the design matrix for HS model
#library(abind)
#rm(list=ls(all=TRUE))



HSModelDesignMatrix = function(HS.model, Num.Of.Subjects, Time.Points)
{
	N = Num.Of.Subjects
	T = Time.Points
	#==========================================================
    HS.model.cov = all.vars.character(HS.model)$m[[2]]
    #cat("HS.model.cov = ", HS.model.cov, "\n")

    TimeOrder = sort(gsub("IndTime", "", HS.model.cov[HS.model.cov %in% paste0("IndTime", 1:10)]))
    #cat("TimeOrder = ", TimeOrder, "\n")

    DiffTime = sort(gsub("DiffTime", "", HS.model.cov[HS.model.cov %in% paste0("DiffTime", 1:10)]))
    #cat("DiffTime = ", DiffTime, "\n")
    
    #cat("HS.model = \n")
    #print(as.formula(HS.model))

    interaction.terms = attr(terms.formula(as.formula(HS.model)), "term.labels")
   
   	#cat("interaction.terms = ", interaction.terms, "\n")
   	#====================================================================#
    a = length(interaction.terms) 

    #cat("a = ", a, "\n")
    u = NULL

    delta.num = 0
    # for HSD model 

    TimeOrder = unique(TimeOrder)
    DiffTime = unique(DiffTime)

 
    if(a>0){
        u = array(0, c(T, T, N, a))
        if(length(TimeOrder)>0){
            for(t in 1:length(TimeOrder))
                u[,,1:N, t][as.matrix(dist(1:T, method="euclidean", diag = TRUE, upper = TRUE))==t] = 1
            delta.num = delta.num + length(TimeOrder)
            cat("length(TimeOrder)>0 delta.num = ", delta.num, "\n")
        }
        #cat("length(DiffTime) = ", length(DiffTime)>0, "\n")
        if(length(DiffTime)>0){
            for(t in 1:length(DiffTime))
                u[,,1:N, delta.num+t] = (as.matrix(dist(1:T, method="euclidean", diag = TRUE, upper = TRUE)))^t
            delta.num = delta.num + length(DiffTime)
        }

        main.terms = interaction.terms[-grep("Time",interaction.terms)]
        int.terms = interaction.terms[grep(":",interaction.terms)]

        #cat("main = ", main.terms, "\n")
        #cat("int.terms = ", int.terms, "\n")
        
        
        if(length(main.terms)>0){
            #cat("=============== 1 ============\n")
            for(t in 1:length(main.terms)){
                uu = data[, names(data)%in%c("id", main.terms[t]), drop=FALSE]
                HSD.cov = unique(uu)[, , drop=FALSE]
                HSD.cov = HSD.cov[HSD.cov$id%in%id, ]
                HSD.cov = as.matrix(HSD.cov[, -1])
                for(sub in 1:N)
                    u[,,sub, delta.num+t] = matrix(HSD.cov[sub], T, T)
            }
            delta.num = delta.num + length(main.terms)
        }
        if(length(int.terms)>0){
            #cat("=============== 2 ============\n")
            for(t in 1:length(int.terms)){
                #cat("int.terms[t]= ", int.terms[t], "\n")
                int.terms.tmp = strsplit(int.terms[t], ":")[[1]]
                #cat("int.terms.tmp = ", int.terms.tmp, "\n")
                int.terms.tmp.IndTime.TF = int.terms.tmp %in% paste0("IndTime", 1:10)
                int.terms.tmp.DiffTime.TF = int.terms.tmp %in% paste0("DiffTime", 1:10)
                int.terms.tmp.IndTime = int.terms.tmp[int.terms.tmp %in% paste0("IndTime", 1:10)]
                int.terms.tmp.DiffTime = int.terms.tmp[int.terms.tmp %in% paste0("DiffTime", 1:10)]
                #cat("any(int.terms.tmp.IndTime) = ", length(int.terms.tmp.IndTime), "\n")
                #cat("any(int.terms.tmp.DiffTime) = ", length(int.terms.tmp.DiffTime), "\n")

                if(length(int.terms.tmp.IndTime)>0){
                    #cat("=============== 3 ============\n")
                    #cat("int.terms.tmp.IndTime = ", int.terms.tmp.IndTime, "\n")
                    IndTime.tmp = as.numeric(gsub("\\D", "", int.terms.tmp.IndTime))
                    #as.numeric(int.terms.tmp.IndTime) #as.numeric(gsub("IndTime", "", int.terms.tmp.IndTime))
                    #cat("IndTime.tmp = ", IndTime.tmp, "\n")
                    HSD.cov.tmp = int.terms.tmp[!int.terms.tmp.IndTime.TF]
                    #cat("HSD.cov.tmp = ", HSD.cov.tmp, "\n")
                    
                    uu.tmp = matrix(0, T, T)
                    uu.tmp[as.matrix(dist(1:T, method="euclidean", diag = TRUE, upper = TRUE))==IndTime.tmp] = 1

                    uu = data[, names(data)%in%c("id", HSD.cov.tmp), drop=FALSE]

                    uu = uu[complete.cases(uu), ]
                    #cat("uu = \n")
                    #print(uu)

                    HSD.cov = unique(uu)[, , drop=FALSE]

                    #cat("HSD.cov 1= \n")
                    #print(HSD.cov)

                    HSD.cov = HSD.cov[HSD.cov$id%in%id, ]
                    HSD.cov = as.matrix(HSD.cov[, -1])

                    #cat("HSD.cov 2= \n")
                    #print(HSD.cov)

                    for(sub in 1:N){
                        u[,, sub, delta.num+1] = uu.tmp*matrix(HSD.cov[sub], T, T)
                    }

                    delta.num = delta.num + 1

                }
                else if(length(int.terms.tmp.DiffTime)>0){
                    #cat("=============== 4 ============\n")
                    DiffTime.tmp = as.numeric(gsub("DiffTime", "", int.terms.tmp.DiffTime))
                    HSD.cov.tmp = int.terms.tmp[!int.terms.tmp.DiffTime.TF]
                    uu.tmp = matrix(0, T, T)
                    uu.tmp[as.matrix(dist(1:T, method="euclidean", diag = TRUE, upper = TRUE))==DiffTime.tmp] = 1

                    uu = data[, names(data)%in%c("id", HSD.cov.tmp), drop=FALSE]
                    HSD.cov = unique(uu)[, , drop=FALSE]
                    HSD.cov = HSD.cov[HSD.cov$id%in%id, ]
                    HSD.cov = as.matrix(HSD.cov[, -1])

                    for(sub in 1:N)
                        u[,, sub, delta.num + 1] = uu.tmp*matrix(HSD.cov[sub], T, T)

                    delta.num = delta.num + 1
                }
                else{
                    #cat("=============== 5 ============\n")
                    #cat("int.terms.tmp = ", int.terms.tmp, "\n")
                    uu = data[, names(data)%in%c("id", int.terms.tmp), drop=FALSE]
                    HSD.cov1 = unique(uu)[, 1, drop=FALSE]
                    #print(head(HSD.cov1))
                    HSD.cov2 = unique(uu)[, 2, drop=FALSE]
                    #print(head(HSD.cov2))
                    HSD.cov = cbind(HSD.cov1, HSD.cov2)
                    HSD.cov = HSD.cov[HSD.cov$id%in%id, ]
                    #print(head(HSD.cov))
                    HSD.cov = as.matrix(HSD.cov[, -1])
                    for(sub in 1:N)
                        u[,, sub ,delta.num + 1] = matrix(prod(HSD.cov[sub, ]), T, T)


                }


            }

        }

    }

    #print(dim(u))

    #cat("delta.num = ", delta.num, "\n")
    #cat("=============== 6 ============\n")
    if(a != delta.num)
        stop("Something wrong to specify the design matrix in HSD model.\n")
    if(any(HS.model.cov==1)){
        a = a+1
        u = abind(array(1, c(T, T, N)), u, along=4)
    }

    #print(dim(u))
    #print(u[, , , 1])
    #print("================")
    #print(u[, , , 2])
    uu = u

    #cat("dim(u) = ", dim(u), "\n")
    #cat("u = \n")

    u = aperm(u, c(4, 3, 1, 2))
    #print(dim(u))
    if(a>0)
        dim(u) = c(N*a, T, T)


    return(u)


}

#AA = HSModelDesignMatrix(HS.model = ~IndTime1+IndTime2, Num.Of.Subjects = 10, Time.Points = 4)
