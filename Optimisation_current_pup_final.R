################################################################################
# R code used to implement headcount-based optimisation with current pup 
# locations
# author: Grace Carmichael
################################################################################

# load packages 
library(gurobi) #(need a licence to use)
library(tmap)

# create dataset with current pick up points and which ward they fall into

pups <- st_transform(pup_sf, 2263)
# vector of zeros to populate with wards
pups$Ward <- rep(0, nrow(pups))
# get logical matrix with pups on rows and wards on columnsindicating whether a
# pup is inside a ward (true) or not (false)
int <- st_intersects(pups, st_transform(wards,2263), sparse = FALSE)

# populate Ward column with the ward that pup i is located in
for(i in 1:nrow(int)){
  j <- which(int[i,] == TRUE)
  pups$Ward[i] <- wards[j,]$WardID
}

# save as file
write.csv(pups, file = "pups.csv")
pups <- read.csv("pups.csv",header=TRUE,row.names=1, check.names=FALSE)


# create a vector that has 1 entry for each ward indicating how many pups are 
# in the ward
current_pups <- rep(0, length(wards_no_WC))

for (i in 1:length(wards_no_WC)) {
  current_pups[i] <- length(which(pups$Ward %in% wards[i]))
}

# create a vector indicating whether a ward has at least 1 pup or not
current_pups_ind <- current_pups>=1


#-------------------------------------------------------------------------------
# implement optimisation with different max pup numbers and different c
#-------------------------------------------------------------------------------

# create vector of different values of p (max pup number) to try
pup_num <- c(1500,1000,500,200,100,50,20,10)
# create vector of different values of c to try
c_opt <- c(0.2,0.5,0.8)

# create a list to store results
c_sens_res <- list()

# create N matrix with set D (D=5)
N <- tt_mat<=5

# first loop through values of c
for(i in 1:3){
  # create vector to store results for different values of p
  p_sens_res <- list()
  # loop through values of p
  for(j in 1:8){
    # set maximum number of pup to assign
    p <- pup_num[j]
    # set max no of pups that can be in a ward with a current pup
    c <- c_opt[i]*p
    # create model as a list
    model <- list()
    # name model
    model$modelname <- 'pup_allocation'
    # tell gurobi that we are maximising objective
    model$modelsense <- 'max'
    
    #create LHS constraint coefficient matrix (A)
    Y_const <- diag(x=-1,nrow = nrow(N),ncol = ncol(N)) # Y portion of constraints
    
    C1 <- cbind(N,Y_const) # first constraint multiply N by X and Y_const by Y
    C2 <- c(rep(1,ncol(N)), rep(0,ncol(Y_const))) # second constraint sum all X 
    C3 <- c(current_pups_ind, rep(0,ncol(Y_const))) # current pup constrain
    
    # bind 3 constraints into 1 constraint matrix
    A <- rbind(C1,C2,C3) 
    
    # set variable types ("B" for binary)
    model$vtype <- rep("B", ncol(A))
    # set objective function coefficients
    model$obj <- c(rep(0,ncol(N)),h)
    # name the decision variables
    model$varnames <- c(rep("allocated",ncol(N)), rep("served",ncol(Y_const)))
    # set LHS constrain coefficients to previously created matrix A
    model$A <- as.matrix(A)
    # set RHS values of constraints
    model$rhs <- c(rep(0,ncol(N)),p,c)
    # set which constraints should be <= or >=
    model$sense <- c(rep(">",nrow(N)),"<","<")
    # name the constraints
    model$constrnames <- c(rep("coverage",nrow(N)), "facility_num", "current_pup")
    
    # set up model
    gurobi_write(model,'pup_allocation_base.lp')
    
    # set parameters
    params <- list()
    #set method
    params$method <- 1
    
    # optimise
    res <- gurobi(model, params)
    # store results
    p_sens_res[[j]]<-res
    
    # extract ward ids of wards assigned a pup in optimisation
    wards_with_pup <- sort(wards_no_WC)[which(res$x[1:ncol(N)] == 1)]
    # extract rows of wards dataset that were assigned pups
    wards_assigned <- wards[which(wards$WardID %in% wards_with_pup),]
    # find centers of wards that were assigned pup to plot on map
    centers <- st_centroid(st_transform(wards_assigned,2263))
    
    # create map - plot current pup locations (dark grey) 
    # and suggested new pup locations (red)
    tmap_mode("plot")
    tmap_options(check.and.fix = TRUE)
    m <- tm_shape(wards) + tm_fill("ProvinceName", n = 8, palette = pal[-1], legend.show = FALSE) + tm_borders() + tm_shape(pup_sf) + tm_dots(col = "grey21") + tm_shape(centers) + tm_dots(col = "red2", size=0.08)
    # save map with new pup locations
    tmap_save(m, filename = paste("base_","c",c_opt[i],"_","p",pup_num[j],".pdf", sep = ""), height = 5.45, width = 6.5, units = "in")
  }
  # store results for all p for c[i]
  c_sens_res[[i]] <- p_sens_res
}

# save results
saveRDS(c_sens_res, file="c_sens_res.RData")
c_sens_res <- readRDS("c_sens_res.RData")

