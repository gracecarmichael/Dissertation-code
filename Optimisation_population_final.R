################################################################################
# R code used to implement population-based optimisation
# author: Grace Carmichael
################################################################################

# load packages 
library(gurobi) #(need a licence to use)
library(tmap)

### Removing wards with NAs for population from all data ###
pop <- pop_by_ward$Population

# get index of missing populations
NA_ind <- which(is.na(pop)==TRUE)

# remove from N matrix
N <- N[-NA_ind,-NA_ind]
# remove from travel times matrix
tt_mat <- tt_mat[-NA_ind,-NA_ind]
# remove from populations
pop <- na.omit(pop)
# remove from wards
wards_pop <- wards[order(wards$WardID),]
wards_pop <- wards[-which(wards$WardID %in% NA_ind),]


# create vector of different values of D (max travel time) to try
D_opt <- c(5,10,20,30,40,50,60,90)

# create empty list to store results
pop_no_restrict_sens_res <- list()


#-------------------------------------------------------------------------------
# implement unrestricted optimisation (no limit on pup number)
#-------------------------------------------------------------------------------

# loop through values of D
for(i in 1:8){
  # set D for the iteration
  D <- D_opt[i]
  # if travel time between i and j <= D then TRUE (1), else FALSE (0)
  N <- tt_mat<=D
  
  # create model as a list
  model <- list()
  # name model
  model$modelname <- 'pop_allocation_unlim_pup'
  # tell gurobi that we are maximising objective
  model$modelsense <- 'max'
  
  # create LHS constraint coefficient matrix (A)
  Y_const <- diag(x=-1,nrow = nrow(N),ncol = ncol(N)) # Y portion of constraints
  # full constraint LHS matrix
  A <- cbind(N,Y_const) # multiply N by X and Y_const by Y 
  
  # set variable types ("B" for binary)
  model$vtype <- rep("B", ncol(A))
  # set objective function coefficients
  model$obj <- c(rep(0,ncol(N)),pop)
  # name the decision variables
  model$varnames <- c(rep("allocated",ncol(N)), rep("served",ncol(Y_const)))
  # set LHS constraint coefficients to previously created matrix A
  model$A <- as.matrix(A)
  # set RHS values of constraints
  model$rhs <- c(rep(0,ncol(N)))
  # set which constraints should be <= or >=
  model$sense <- c(rep(">",nrow(N)))
  # name the constraints
  model$constrnames <- c(rep("coverage",nrow(N)))
  
  # set up model
  gurobi_write(model,'pop_allocation_pup.lp')
  
  # set parameters
  params <- list()
  # set method
  params$method <- 1
  
  # optimise
  res <- gurobi(model, params)
  # store results
  pop_no_restrict_sens_res[[i]]<-res
  
  # extract ward ids of wards assigned a pup in optimisation
  wards_with_pup <- sort(wards_pop$WardID)[which(res$x[1:ncol(N)] == 1)]
  # extract rows of wards dataset that were assigned pups
  wards_assigned <- wards[which(wards_pop$WardID %in% wards_with_pup),]
  # find centers of wards that were assigned pup to plot on map
  centers <- st_centroid(st_transform(wards_assigned,2263))
  
  # create map - plot current pup locations (dark grey) 
  # and suggested new pup locations (red)
  tmap_mode("plot")
  tmap_options(check.and.fix = TRUE)
  m <- tm_shape(wards) + tm_fill("ProvinceName", n = 8, palette = pal[-1], legend.show = FALSE) + tm_borders() + tm_shape(pup_sf) + tm_dots(col = "grey21") + tm_shape(centers) + tm_dots(col = "red2", size=0.08)
  # save map with new pup locations
  tmap_save(m, filename = paste("pop_","D",D_opt[i],"_","punlim",".pdf", sep = ""), height = 5.45, width = 6.5, units = "in")
}

# save results
saveRDS(pop_no_restrict_sens_res, file="pop_no_restrict_sens_res.RData")
pop_no_restrict_sens_res <- readRDS("pop_no_restrict_sens_res.RData")


#-------------------------------------------------------------------------------
# implement optimisation with different max pup numbers
#-------------------------------------------------------------------------------

# create vector of different values of D (max travel time) to try
D_opt <- c(5,10,20,30,40,50,60,90)
# create vector of different values of p (max pup number) to try
pup_num <- c(1500,1000,500,200,100,50,20,10)

# create a list to store results
pop_sens_res <- list()

# first loop through values of D
for(i in 1:8){
  # create vector to store results for different values of p
  p_sens_res <- list()
  # set D for the iteration
  D <- D_opt[i]
  # if travel time between i and j <= D then TRUE (1), else FALSE (0)
  N <- tt_mat<=D
  
  # loop through values of p
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
    
    C1 <- cbind(N,Y_const) # first constraint multiply N by X and Y_const by Y 
    C2 <- c(rep(1,ncol(N)), rep(0,ncol(Y_const))) # second constraint we sum all 
    # X on LHS to get number of pups added
    # bind 2 constraints into 1 constraint matrix
    A <- rbind(C1,C2) 
    
    # set variable types ("B" for binary)
    model$vtype <- rep("B", ncol(A))
    # set objective function coefficients
    model$obj <- c(rep(0,ncol(N)),pop)
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
    gurobi_write(model,'pop_allocation_pup2.lp')
    
    # set parameters
    params <- list()
    params$method <- 1
    
    # optimise
    res <- gurobi(model, params)
    p_sens_res[[j]]<-res
    
    # extract ward ids of wards assigned a pup in optimisation
    wards_with_pup <- sort(wards_pop$WardID)[which(res$x[1:ncol(N)] == 1)]
    # extract rows of wards dataset that were assigned pups
    wards_assigned <- wards[which(wards_pop$WardID %in% wards_with_pup),]
    # find centers of wards that were assigned pup to plot on map
    centers <- st_centroid(st_transform(wards_assigned,2263))
    
    
    # create map - plot current pup locations (dark grey) 
    # and suggested new pup locations (red)
    library(tmap)
    tmap_mode("plot")
    tmap_options(check.and.fix = TRUE)
    m <- tm_shape(wards) + tm_fill("ProvinceName", n = 8, palette = pal[-1], legend.show = FALSE) + tm_borders() + tm_shape(pup_sf) + tm_dots(col = "grey21") + tm_shape(centers) + tm_dots(col = "red2", size=0.08)
    # save map with new pup locations
    tmap_save(m, filename = paste("pop_","D",D_opt[i],"_","p",pup_num[j],".pdf", sep = ""), height = 5.45, width = 6.5, units = "in")
  }
  # store results for all p for D[i]
  pop_sens_res[[i]] <- p_sens_res
}

# save results
saveRDS(pop_sens_res, file="pop_sens_res.RData")
pop_sens_res <- readRDS("pop_sens_res.RData")
