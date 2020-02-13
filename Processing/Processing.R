setwd("C:/Users/beacr/Politecnico di Milano/Simone Panzeri - Bayesian Project/data_extraction/dati_iniziali")
library(readxl)
library(coda)
#### IMPORT DATASET ####

#Intensity values for every pixel at every frame
dataset <- read_excel("Dataset preliminari Tr1.xlsx", sheet = "SM+LM")
head(dataset)
n_frame <- dim(dataset)[1]    # rows <-> values at each frame
n_pixel <- dim(dataset)[2]    # columns <-> pixels


# Coordinates
Coord_data <- read_excel("Coord.xlsx")
head(Coord_data)      #(y,x)
Coord_data = rev(Coord_data)
#head(Coord_data)     #(x,y)
#x11()
#plot(Coord_data)

#180° rotation
Coord_data[,2] <- max(Coord_data[,2])-Coord_data[,2]
#x11()
#plot(Coord_data)

# write.table(Coord_data, file="Coordinates.xlsx")

Coord <- as.matrix(Coord_data)
v <- as.matrix(dataset)

#### INITIALISATION OF SOME OBJECTS ####
#SET the time when the laser is switched off (fixed technical parameter)
laser_off <- 185

#initialisation of some matrices -- will define them later --
v_nospark <- v # dataset senza scintille in forma matriciale
c1 <- matrix(0, n_frame, n_pixel) # distanza tra un passaggio                                  # del laser e il successivo > 25 frame
c2 <- c1 # c1 raffinata con picco maggiore della soglia

N_pass <- matrix(NA,n_frame, n_pixel)
N_on <- matrix(NA,n_frame, n_pixel)
mean_intensity <- matrix(NA,n_frame,n_pixel)
T_raff <- matrix(NA,n_frame, n_pixel) # Cooling times

N_intervals <- rep(NA,n_pixel)
mean_T <- rep(NA,n_pixel) # mean cooling time
last_laser <- rep(0,n_pixel) # last_laser[p] = last frame for which the laser is on (on pixel p)
rif <- rep(0, n_pixel) # mean intensity when there isn't the laser
mean_val <- rif  # mean intensity (with zeros)
count <- rep(0,n_pixel)
max_int <- 0

#### PROCESSING ####


# Simulation of the missing values (sparkles)
trunc_norm <- function(n, mu, sigma_sq, a1, a2){
  results <- c()
  for(i in 1:n){
    F1 <- pnorm (a1,mu,sqrt(sigma_sq))
    F2 <- pnorm (a2,mu,sqrt(sigma_sq))
    u <- runif(1,min=F1,max=F2)   # U(F1,F2)
    results[i] <- qnorm(u,mu,sqrt(sigma_sq))
  }
  return(results)
}

gibbs <- function (niter,data,mu0,sigma0,tau) # n0,nu0)
{
  mis <- which(is.na(data))
  n_mis <- length(mis)
  
  PHI <- matrix(nrow=niter,ncol=1+n_mis)
  PHI[1,1] <- mu0
  #PHI[1,2] <- sigma0
  for (i in 1:n_mis)
    PHI[1,i+1] <- data[mis[i]-1]
  
  n <- length(data)
  
  Mu <- mu0
  Sigma <- sigma0
  
  for (iter in 2:niter)
  {
    for (i in 1:n_mis)
    {
      # Sample each missing value
      PHI[iter,i+1] <- data[mis[i]] <- trunc_norm(1,Mu,Sigma,50,200)
    }
    # update the parameters
    mu1 <- (n/sigma0*mean(data[which(data!=0)])+mu0/tau)/(n/sigma0+1/tau)
    taun <- n*tau/(n*tau+sigma0)+sigma0*mu0/(n*tau+sigma0)
    Mu <- rnorm(1,mu1,sqrt(taun))
    PHI[iter,1] <- Mu
  }
  
  return (PHI)
}

for (p in 1:n_pixel)
{
  mean_val[p]<-mean(v_nospark[which(v_nospark[,p]!=0),p])
}
contatore<-rep(0,n_pixel)
for (p in 1:n_pixel)
{
  V <- var(v_nospark[which(v_nospark[,p]!=0),p])
  Mu <- mean(v_nospark[which(v_nospark[,p]!=0),p])
  for (f in 2:(n_frame-1))
    if (v_nospark[f,p]==0 & v_nospark[f+1,p]!=0 & v_nospark[f-1,p]!=0)
    {
      contatore[p] <-contatore[p] +1
      v_nospark[f,p] <- NA
      contatore[p] <- contatore[p]+1
    }
  if (contatore[p] > 0)
  {
    miss <- which(is.na(v_nospark[,p]))
    simul <- gibbs(1000,v_nospark[,p],Mu,V,var(mean_val)) # length(which(v_nospark[,p]!=0))-length(miss),1)
    for (i in 1:length(miss))
      v_nospark[miss[i],p] <- simul[1000,i+1]
    burnin <- 200
    simulazione <- mcmc(simul, start=burnin+1, end=1000)
    # plot(simulazione)
    # x11()
    # plot( simul[burnin:1000,1:2],type="l",
    #       lty=1,col="gray",xlab=expression(mu),ylab=expression(y_miss))
    # summary(simulazione)
    # acfplot(simulazione)
    # effectiveSize(simulazione)
  }
}

mean_val <- colMeans(v_nospark)
#alternative (laser included) :
#rif[p] <- ifelse( (mean(v_ss[[p]][(n_frame-100):n_frame])) < 85 , mean(v_ss[[p]][(n_frame-100):n_frame]) , 85)

# Creation of the final dataset
# write.table(v_nospark, file="Dataset_definitivo_gibbs.xlsx")

# Cut the dataset (we will only consider the top-right angle)
cut <- matrix()
cut <- v_nospark[,which(Coord[,1]>=60)]
# write.table(cut, file="Dataset_cut_gibbs.xlsx")
Coord_cut <- matrix()
Coord_cut <- Coord[which(Coord[,1]>=60),]
# write.table(Coord_cut, file="Coord_cut.xlsx")

#### directly import the final data: ####
# v_nospark <- read_table("Dataset_cut.xlsx")
# v_nospark <- read_table("Dataset_definitivo.xlsx")
# rif <- colMeans(v_nospark[which(v_nospark!=0)])
# mean_val <- colMeans(v_nospark)
