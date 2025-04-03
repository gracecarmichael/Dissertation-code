################################################################################
# R code used to implement headcount-based optimisation
# author: Grace Carmichael
################################################################################

# load packages 
library(gurobi) #(need a licence to use)
library(tmap)

# create vector of different values of D (max travel time) to try
D_opt <- c(5,10,20,30,40,50,60,90)

# create empty list to store results
no_restrict_sens_res <- list()


#-------------------------------------------------------------------------------
# implement unrestricted optimisation (no limit on pup number)
#-------------------------------------------------------------------------------

# loop through values of D
for(i in 1:8){
  # set D for the iteration
  D <- D_opt[i]
  # read in travel times matrix
  tt_mat <- read.table("tt_mat.txt",header=TRUE,row.names=1, check.names=FALSE)
  # if travel time between i and j <= D then TRUE (1), else FALSE (0)
  N <- tt_mat<=D
  
  # create model as a list
  model <- list()
  # name model
  model$modelname <- 'pup_allocation_unlim_pup'
  # tell gurobi that we are maximising objective
  model$modelsense <- 'max'
  
  #create LHS constraint coefficient matrix (A)
  Y_const <- diag(x=-1,nrow = nrow(N),ncol = ncol(N)) # Y portion of constraints
  # full constraint LHS matrix
  A <- cbind(N,Y_const) # first constraint we multiply N by X and Y_const by Y 
 
  # set variable types ("B" for binary)
  model$vtype <- rep("B", ncol(A))
  # set objective function coefficients
  model$obj <- c(rep(0,ncol(N)),h)
  # name the decision variables
  model$varnames <- c(rep("allocated",ncol(N)), rep("served",ncol(Y_const)))
  # set LHS constrain coefficients to previously created matrix A
  model$A <- as.matrix(A)
  # set RHS values of constraints
  model$rhs <- c(rep(0,ncol(N)))
  # set which constraints should be <= or >=
  model$sense <- c(rep(">",nrow(N)))
  # name the constraints
  model$constrnames <- c(rep("coverage",nrow(N)))
  
  # set up model
  gurobi_write(model,'pup_allocation_base.lp')
  
  # set parameters
  params <- list()
  params$method <- 1
  
  # optimise
  res <- gurobi(model, params)
  # store results
  no_restrict_sens_res[[i]]<-res
    
  # extract ward ids of wards assigned a pup in optimisation
  wards_with_pup <- sort(wards)[which(res$x[1:ncol(N)] == 1)]
  # extract rows of wards dataset that were assigned pups
  wards_assigned <- wards[which(wards$WardID %in% wards_with_pup),]
  # find centers of wards that were assigned pup to plot on map
  centers <- st_centroid(st_transform(wards_assigned,2263))
    
  # create map - plot clinic locations(lighter grey), current pup locations 
  # (darker grey) and suggested new pup locations (red)
  tmap_mode("plot")
  tmap_options(check.and.fix = TRUE)
  m <- tm_shape(wards) + tm_fill("ProvinceName", n = 8, palette = pal[-1], legend.show = FALSE) + tm_borders() + tm_shape(pup_sf) + tm_dots(col = "grey21") + tm_shape(centers) + tm_dots(col = "red2", size=0.08)
  # save map with new pup locations
  tmap_save(m, filename = paste("base_","D",D_opt[i],"_","punlim",".pdf", sep = ""), height = 5.45, width = 6.5, units = "in")
}

# save results
saveRDS(no_restrict_sens_res, file="no_restrict_sens_res.RData")
no_restrict_sens_res <- readRDS("no_restrict_sens_res.RData")


#-------------------------------------------------------------------------------
# implement optimisation with different max pup numbers and max travel time
#-------------------------------------------------------------------------------

# create vector of different values of D (max travel time) to try
D_opt <- c(5,10,20,30,40,50,60,90)
# create vector of different values of p (max pup number) to try
pup_num <- c(1500,1000,500,200,100,50,20,10)

# create a list to store results
base_sens_res <- list()

# first loop through values of D
for(i in 1:6){
  # create vector to store results for different values of p
  p_sens_res <- list()
  # set value of D for this iteration
  D <- D_opt[i]
  # read in travel times matrix
  tt_mat <- read.table("tt_mat.txt",header=TRUE,row.names=1, check.names=FALSE)
  # if travel time between i and j <= D then TRUE (1), else FALSE (0)
  N <- tt_mat<=D
  for(j in 1:8){
    # set maximum number of pup to assign
    p <- pup_num[j]
    
    # create model as a list
    model <- list()
    # name model
    model$modelname <- 'pup_allocation'
    # tell gurobi that we are maximising objective
    model$modelsense <- 'max'
    
    #create LHS constrain coefficient matrix (A)
    Y_const <- diag(x=-1,nrow = nrow(N),ncol = ncol(N)) # Y portion of constraints
    
    C1 <- cbind(N,Y_const)# first constraint multiply N by X and Y_const by Y 
    C2 <- c(rep(1,ncol(N)), rep(0,ncol(Y_const)))# second constraint sums all X
    # bind 2 constraints into 1 constraint matrix
    A <- rbind(C1,C2) 
    
    # set variable types ("B" for binary)
    model$vtype <- rep("B", ncol(A))
    # set objective function coefficients
    model$obj <- c(rep(0,ncol(N)),h)
    # name the decision variables
    model$varnames <- c(rep("allocated",ncol(N)), rep("served",ncol(Y_const)))
    # set LHS constrain coefficients to previously created matrix A
    model$A <- as.matrix(A)
    # set RHS values of constraints
    model$rhs <- c(rep(0,ncol(N)),p)
    # set which constraints should be <= or >=
    model$sense <- c(rep(">",nrow(N)),"<")
    # name the constraints
    model$constrnames <- c(rep("coverage",nrow(N)), "facility_num")
    
    # set up model
    gurobi_write(model,'pup_allocation_base.lp')
    
    # set parameters
    params <- list()
    params$method <- 1
    
    # optimise
    res <- gurobi(model, params)
    # store results
    p_sens_res[[j]]<-res
    
    # extract ward ids of wards assigned a pup in optimisation
    wards_with_pup <- sort(wards)[which(res$x[1:ncol(N)] == 1)]
    # extract rows of wards dataset that were assigned pups
    wards_assigned <- wards[which(wards$WardID %in% wards_with_pup),]
    # find centers of wards that were assigned pup to plot on map
    centers <- st_centroid(st_transform(wards_assigned,2263))
    

    # create map - plot current pup locations (dark grey) 
    # and suggested new pup locations (red)
    library(tmap)
    tmap_mode("plot")
    tmap_options(check.and.fix = TRUE)
    m <- tm_shape(wards) + tm_fill("ProvinceName", n = 8, palette = pal[-1], legend.show = FALSE) + tm_borders() + tm_shape(pup_sf) + tm_dots(col = "grey21") + tm_shape(centers) + tm_dots(col = "red2", size=0.08)
    # save map with new pup locations
    tmap_save(m, filename = paste("base_","D",D_opt[i],"_","p",pup_num[j],".pdf", sep = ""), height = 5.45, width = 6.5, units = "in")
  }
  # store results for all p for D[i]
  base_sens_res[[i]] <- p_sens_res
}


# save results
saveRDS(base_sens_res, file="base_sens_res.RData")
base_sens_res <- readRDS("base_sens_res.RData")

