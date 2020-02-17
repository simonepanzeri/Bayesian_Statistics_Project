setwd("C:/Users/beacr/Politecnico di Milano/Simone Panzeri - Bayesian Project/data_extraction")
library(readxl)
library(mvtnorm)
library(coda)
library("RColorBrewer", lib.loc="~/R/win-library/3.5")

#### IMPORT DATASET ####

# We only work with the top right angle (420 pixels, 329 frames)
dataset <- read.table("Dataset_cut_gibbs.xlsx")

# if you want to work on the full image (996 pixels, 329 frames)
# dataset <- read.table("Dataset_definitivo_gibbs.xlsx")

head(dataset)
n_frame <- dim(dataset)[1]    # rows <-> values at each frame
n_pixel <- dim(dataset)[2]    # columns <-> pixels

Coord_data <- read.table("Coord_cut.xlsx")

# if you want to work on the full image
# Coord_data <- read.table("Coordinates.xlsx")

Coord <- as.matrix(Coord_data)
v_nospark <- as.matrix(dataset)


#### DECLARATION OF SOME OBJECTS ####

#SET the time when the laser is switched off (fixed technical parameter)
laser_off <- 185

c1 <- matrix(0, n_frame, n_pixel) # distanza tra un passaggio                                  # del laser e il successivo > 25 frame
c2 <- c1 # c1 raffinata con picco maggiore della soglia

N_pass <- matrix(NA,n_frame, n_pixel)
N_on <- matrix(NA,n_frame, n_pixel)
mean_intensity <- matrix(NA,n_frame,n_pixel)
T_raff <- matrix(NA,n_frame, n_pixel) # Cooling times

N_intervals <- rep(NA,n_pixel)
mean_T <- rep(NA,n_pixel) # mean cooling time
last_laser <- rep(0,n_pixel) # last_laser[p] = last frame for which the laser is on (on pixel p)
max_int <- 0

rif <- rep(0, n_pixel) # mean intensity when there isn't the laser
mean_val <- rif  # mean intensity (with zeros)



#### VARIABLES DEFINITION ####
mean_val <- colMeans(v_nospark)

for (p in 1:n_pixel)
{
  rif[p] <- mean(v_nospark[which(v_nospark[,p]!=0),])
  
  int <- 0
  flag <- 0  # indicator: 1 if pos1 has been updated
  # extrema of the interval between two following passages of the laser
  pos1 <- min(which(v_nospark[,p]==0))  
  pos2 <- pos1
  # definition of the variables while the laser is on
  for (f in pos1:laser_off) {
    if(!is.na(v_nospark[f,p])){
      # if I find a laser, I update the interval
      if (v_nospark[f,p]==0)
      {
        # if pos1 hasn't been updated yet, I fix pos1 at f, otherwise I fix pos2 at f
        if (flag==0) {
          pos1 <- f
          # flag <- 1
        } 
        else
        {
          pos2<-f
          # if the distance between two following passaes of the laser
          # is greater than 25, register it in c1 and update the number of
          # suspect intevals for p
          # then refine better by controlling the peak value and register in c2
          if (pos2 - pos1 > 2)
          {  
            flag1 <- 0
            int <- int+1
            c1[(pos1+1):(pos2-1),p] <- int
            picco <- pos1+min(which(v_nospark[(pos1+1):min(pos2-1,pos1+10),p]==max(na.omit(v_nospark[(pos1+1):min(pos2-1,pos1+10),p]))))
            k <- ifelse (sum(na.omit(v_nospark[(picco+1):(pos2-1),p])<rif[p])>0, picco+min(which(na.omit(v_nospark[(picco+1):(pos2-1),p])<rif[p])), pos2)
            flag1 <- 1
            T_raff[int,p] <- k - picco
            c2[picco:k,p] <- int
            
            N_on[int,p] <- sum(v_nospark[1:pos2,p]==0)
            
            mean_intensity[int+1,p] <- mean(v_nospark[pos1:pos2,p])
            
            N_pass[int,p] <- int
          }
          pos1 <- pos2
          flag <- 0
        }
      }
      else
      {
        if (pos1 >= pos2)
          flag <- 1
      }
    }
  }
  
  # definition of the variables when the laser is switched off
  last_laser[p] <- pos1
  int <- int+1
  c1[pos1:n_frame,p] <- int
  
  flag2 <- 0
  picco <- last_laser[p]+min(which(v_nospark[(last_laser[p]+1):(last_laser[p]+10),p]==max(na.omit(v_nospark[(last_laser[p]+1):(last_laser[p]+10),p]))))
  k <- ifelse (sum(na.omit(v_nospark[(picco+1):n_frame,p])<rif[p])>0, picco+min(which(na.omit(v_nospark[(picco+1):n_frame,p])<rif[p])), n_frame)
  flag1 <- 1
  T_raff[int,p] <- k - picco
  c2[picco:k,p] <- int
  
  mean_intensity[int+1,p] <- mean(v_nospark[pos1:n_frame,p])
  N_pass[int,p] <- int
  N_on[int,p] <- sum(v_nospark[1:pos1,p]==0)
  
  N_intervals[p] <- int
}
# maximum number of intervals of the pixels
max_int <- max(N_intervals)

N_on <- N_on[1:max_int,]
N_pass <- N_pass[1:max_int,]
mean_intensity <- mean_intensity[1:(max_int+1),]
T_raff <- T_raff[1:max_int,]

for (p in 1:n_pixel)
  mean_intensity[1,p]=mean(v_nospark[1:min(which(v_nospark[,p]==0)),p])

T_raff2 <- T_raff # without zeros
T_raff <- rbind(0,T_raff)

# save(N_on, file='N_on.dat')
# save(N_pass, file='N_pass.dat')
# save(mean_intensity, file='mean_intensity.dat')
# save(T_raff, file='T_raff.dat')

# mean cooling time of the pixel
mean_T2 <- rep(0,n_pixel)
for (p in 1:n_pixel) 
  mean_T2[p] <- mean(na.omit(T_raff2[,p])) # without zeros
mean_T <- colMeans(na.omit(T_raff))

# save(N_intervals, file='N_intervals.dat')

#### SUMMARIES ####
p <- 256 #set pixel (corresponding to p <- 519 in the full triangle)
summary(T_raff[,p])
summary(N_on[,p])
summary(N_pass[,p])
summary(mean_intensity[,p])

#intervals <- which(!is.na(N_pass[,p]))

variables <- cbind(mean_intensity[,p],T_raff[,p],N_on[,p],N_pass[,p])
colnames(variables) <- c('mean intensity','T raff','N on', 'N pass')
x11()
plot(as.data.frame(variables))

x11()
par(mfrow=c(2,2))
plot(T_raff[1:N_intervals[p],p],xlab='intervals',ylab='cooling times')
plot(N_on[1:N_intervals[p],p],type='l',xlab='intervals',ylab='number of frames with laser on')
plot(N_pass[1:N_intervals[p],p],xlab='intervals',ylab='number of passages of the laser')
plot(mean_intensity[1:N_intervals[p],p],xlab='intervals',ylab='mean intensities')

#### SOME GRAPHICAL REPRESENTATIONS OF THE DATA####
# x11()
# par(mfrow=c(2,1))
# plot(1:n_pixel,mean_val,xlab='pixel',ylab='intenisty',main='mean intensity over the frames without laser',type='l',ylim=c(70,170))
# plot(1:n_pixel,rif,xlab='pixel',ylab='intensity',main='mean intensity over all the frames',type='l',ylim=c(70,170))

p <- 77 # choose a pixel to visualize
x11()
plot(1:n_frame,v_nospark[,p],xlab='Frames',ylab='Intensity', type='l',xlim=c(0,n_frame),ylim=c(0,255),main=paste("Intensity (without spatters) of pixel ", p))
abline(h=rif[p],col='red')

x11()
par(mfrow=c(2,1))
plot(1:n_frame,c1[,p], type='l',col='black',xlab='Frames',ylab='c1',main=paste("Detecting Cooling Problem 1 for the pixel ",p))
plot(1:n_frame,c2[,p],type='l',col='red',xlab='Frames',ylab='c2',main=paste("Detecting Cooling Problem 2 for the pixel",p))

x11()
par(mfrow=c(1,3))
plot(1:n_pixel,mean_T[1:n_pixel],col='black',xlab='Pixels',ylab='mean_cool_time[p]',xlim=c(0,n_pixel),main="Mean cooling times")
hist(mean_T)

X <- seq(50,110)
Y <- seq(15,50)
W <- matrix(NA,110,50)
for(i in 1:n_pixel){
  W[Coord[i,1],Coord[i,2]] <- mean_T2[i]  
}
my_colors = brewer.pal( 9, 'Blues')
image(X,Y,W[50:110,15:50], col=my_colors, xlab='pixels',main=paste('Mean cooling times'))
legend(title  = 'Mean cooling times','bottomright',rev(c('high values (9.37-10.43)','low values')) , 
       col = cbind(my_colors[1],my_colors[9]), pch = 15)


# High mean intensity
#x11()
#plot(Coord, col=ifelse(mean(mean_intensity[1:n_pixel])<10,'black','red'),pch=16)
# we observe that the highest intensity pixels are in the top right angle
# we know (from the video) that the top right angle will burn

W <- matrix(NA,110,50)
for(i in 1:n_pixel){
  W[Coord[i,1],Coord[i,2]] <- N_intervals[i]  
}
par(mfrow=c(1,2))
my_colors = brewer.pal( 9, 'Blues')
image(X,Y,W[50:110,15:50], col=my_colors, xlab='pixels',main=paste('Number of intervals'))
hist(N_intervals)

x11()
par(mfrow=c(1,2))
W <- matrix(NA,110,50)
for(i in 1:n_pixel){
  W[Coord[i,1],Coord[i,2]] <- mean(mean_intensity[1:N_intervals[i],i])
}
my_colors = brewer.pal( 9, 'Blues')
image(X,Y,W[50:110,15:50], col=my_colors, xlab='pixels',main=paste('Mean of the mean intensities'))
hist(colMeans(na.omit(mean_intensity)), main=paste('Histogram of the mean of mean intensities'))

#### GEOMETRICAL COVARIATES ####

# WORK ON THE FULL IMAGE !!
Coord_full <- read.table("Coordinates.xlsx")
Coord_full <- as.matrix(Coord_full)

# Perimeter extraction
per <- matrix(NA,n_pixel,2)
i <- 1
h <- 0
for (cont in 45:0)
{
  h <- h+1
  row_pixels <- which(Coord_full[,2]==cont)
  #only pixels for which the previous and following x are not in Coord are in the perimeter
  for (j in row_pixels)
  {
    if (sum(Coord_full[row_pixels,1]==Coord_full[j,1]-1)==0 || sum(Coord_full[row_pixels,1]==Coord_full[j,1]+1)==0
        || sum(Coord_full[which(Coord_full[,1]==Coord_full[j,1]),2]==cont-1)==0 || sum(Coord_full[which(Coord_full[,1]==Coord_full[j,1]),2]==cont+1)==0)
    {
      per[i,1] <- Coord_full[j,1]
      per[i,2] <- Coord_full[j,2]
      i <- i+1
    }
  }
}
per <- per[which(!is.na(per[,1])),]
per_cut <- per[which(per[,1]>=60),]

x11()
plot(Coord,pch=0,xlab='x',ylab='y')
points(per_cut,pch=15)

# Angle amplitude
# to each pixel we assign the fraction of
# his circular neighbourhood in the triangle

geom_full <- rep(NA,996)  # full triangle
geom <- rep(NA, n_pixel)  # top rigth angle

rsq <- mean(c((4.5/2)^2, (5/2)^2+(4/2)^2, (10/2)^2+(4/2)^2))

for(p in 1:996)
{
  int_p <- 0
  for (pix in 1:996)
  {
    if ((Coord_full[pix,1]-Coord_full[p,1])^2+(Coord_full[pix,2]-Coord_full[p,2])^2 <= rsq)
      int_p <- int_p+1
  }
  geom_full[p] <- int_p/(pi*rsq)
}

geom <- geom_full[which(Coord_full[,1]>=60)]

x11()
X <- seq(50,110)
Y <- seq(15,50)
W <- matrix(NA,110,50)
for(i in 1:n_pixel){
  W[Coord[i,1],Coord[i,2]] <- geom[i]  
}
image(X,Y,W[50:110,15:50], xlab='pixels',main=paste('Geometrical property'))

summary(geom)

x11()
X_full <- seq(0,110)
Y_full <- seq(0,50)
W_full <- matrix(NA,110,50)
for(i in 1:996){
  W_full[Coord_full[i,1],Coord_full[i,2]] <- geom_full[i]  
}
image(X_full,Y_full,W_full, xlab='pixels',main=paste('Geometrical property'))

summary(geom_full)

# Second geometrical covariate: distance from the vertices

lato_up <- as.numeric(which(Coord_full[,2]>=40))
lato_dx <- as.numeric(which(Coord_full[,1]>60 & Coord_full[,2]<40))
lato_sx <- as.numeric(which(Coord_full[,1]<=60 & Coord_full[,2]<40))

g_full <- rep(NA, 996)

v1_up <- Coord_full[lato_up[which.min(Coord_full[lato_up,1])],]
v2_up <- Coord_full[lato_up[which.max(Coord_full[lato_up,1])],]
v1_sx <- v1_up
v2_sx <- Coord_full[lato_sx[which.min(Coord_full[lato_sx,2])],]
v1_dx <- v2_up
v2_dx <- v2_sx

l_up <- dist(rbind(v1_up,v2_up))
l_sx <- dist(rbind(v1_sx,v2_sx))
l_dx <- dist(rbind(v1_dx,v2_dx))

for (p in lato_up)
  g_full[p] <- sqrt(max((Coord_full[p,1]-v1_up[1])^2+(Coord_full[p,2]-v1_up[2])^2,
                        (Coord_full[p,1]-v2_up[1])^2+(Coord_full[p,2]-v2_up[2])^2))/l_up
for (p in lato_dx)
  g_full[p] <- sqrt(max((Coord_full[p,1]-v1_dx[1])^2+(Coord_full[p,2]-v1_dx[2])^2,
                   (Coord_full[p,1]-v2_dx[1])^2+(Coord_full[p,2]-v2_dx[2])^2))/l_dx
for (p in lato_sx)
  g_full[p] <- sqrt(max((Coord_full[p,1]-v1_sx[1])^2+(Coord_full[p,2]-v1_sx[2])^2,
                        (Coord_full[p,1]-v2_sx[1])^2+(Coord_full[p,2]-v2_sx[2])^2))/l_sx

X_full <- seq(0,110)
Y_full <- seq(0,50)
W_full <- matrix(NA,110,50)
for(i in 1:996){
  W_full[Coord_full[i,1],Coord_full[i,2]] <- g_full[i]  
}
image(X_full,Y_full,W_full, xlab='pixels',main=paste('Geometrical property'))

g <- g_full[which(Coord_full[,1]>=60)]

X <- seq(50,110)
Y <- seq(15,50)
W <- matrix(NA,110,50)
for(i in 1:n_pixel){
  W[Coord[i,1],Coord[i,2]] <- g[i]  
}
image(X,Y,W[50:110,15:50], xlab='pixels',main=paste('Geometrical property'))

# save(g, file='g.dat')
# write.table(geom, 'geom_cut.xlsx')
# write.table(geom_full, 'geom_full.xlsx')

#### "PREVIOUS" TIME INTERVALS & NEIGHBOURS ####

#fix an interval for pixel p i in 2:N_intervals[p]
i <- 3

# Neighbours of each pixel
B <- matrix(NA,5,n_pixel)

# Previous intervals corresponding to each neighbour
t <- matrix(NA,5,n_pixel)

for (p in 1:n_pixel)
{
  # Definition of the neighbourhood of pixel p at interval iC
  temp <- which((Coord[,1]==Coord[p,1] & Coord[,2]==Coord[p,2]+1) |
                  (Coord[,1]==Coord[p,1] & Coord[,2]==Coord[p,2]-1) | 
                  (Coord[,2]==Coord[p,2] & Coord[,1]==Coord[p,1]+1) |
                  (Coord[,2]==Coord[p,2] & Coord[,1]==Coord[p,1]-1))
  B[1:length(temp),p] <- temp
  B[length(temp)+1,p] <- p
  # Definition of the "previous intervals" for each neighbour of pixel p a interval i
  pos1 <- min(which(c1[,p]==i))-1
  pos2 <- max(which(c1[,p]==i))+1
  index <- 0
  neighs <- na.omit(B[,p])
  for (q in neighs)
  {
    if(q!=p){
      f <- pos1
      if(v_nospark[f,q]==0 & c1[f,q]>1)
      {
        index <- index+1
        t[index,p] <- c1[f,q]-1
      }
      else
      {
        while (v_nospark[f,q]==0)
          f = f-1
        if (f!=pos1 & c1[f,q]>0)
        {
          index <- index+1
          t[index,p] <- c1[f,q]
        }
        else
        {
          B[(index+1):4,p] <- B[(index+2):5,p]
        }
      }
    }
  }
  t[index+1,p] <- i-1
}

#### DELTA DATA for the LOGIT model ####

# frequence of diagonal covariances
delta <- rep (0, n_pixel)

# How to determine if a matrix is almost diagonal
# (and the diagonal values can be considered almost equal)
is_diag <- function (Q)
{
  nr <- dim(Q)[1]
  nc <- dim(Q)[2]
  
  d <- median(diag(Q))
  
  sum <- 0
  
  for (i in 1:nr)
    for (j in 1:nc)
    {
      if (j==i & (abs(Q[i,j]/d) < 1 | abs(Q[i,j]/d) > 10)) # diagonal values of the same order
        sum <- sum + 1
      else if (j!=i & abs(Q[i,j]/d) > 1) # non-diagonal values of inferior order than the diagonal 
        sum <- sum + 1
    }
  
  return (sum/(nr*nc) < 0.1)
}

min_length<-rep(n_frame, n_pixel)

for (p in 1:n_pixel)
{
  temp <- NA
  for (n in 1:N_intervals[p])
  {
    pos1 <- min(which(c1[,p] == n))
    pos2 <- max(which(c1[,p] == n))
    int <- pos1:pos2
    if (pos2 - pos1 >= 5)
    {
      if (pos2 - pos1 < min_length[p])
        min_length[p] <- pos2 - pos1
      if(is.na(temp))
        temp <- v_nospark[int,p]
      else
        temp <- cbind(temp, v_nospark[int[1:min_length[p]],p])
        # WARNING: la seganalazione è qui
        #          ci dice che gli intervalli più lunghi del primo vengono tagliati
        #          --> usiamo sempre lo stesso numero di dati !!
    }
  }
  temp <- temp[1:min_length[p],]
  if (is_diag(var(temp)))
    delta[p] <- 1
}

x11()
plot(Coord)
points(Coord[which(delta==1),],pch=19)


# save(delta, file='delta.dat')

delta_cut <- delta[which(Coord[,1]>=60)]
# save(delta_cut, file='delta_cut.dat')

#### HYPERPARAMETERS FOR THE MODELS ####

# Model 2

# logistic regression for the parameter delta with geometrical covariates

fit2 <- glm(delta ~ geom + g, family = binomial (link="logit"))

mub <- fit2$coefficients
B <- vcov(fit2)

# linear regression for the parameter mean_intensity depending on its previous values

S <- rep(NA, n_pixel)
for (p in 1:n_pixel)
{
  fit2ar <- lm(mean_intensity[2:13,p] ~ mean_intensity[1:12,p])
  S[p] <- sum(residuals(fit2ar)^2)/fit2ar$df
}

x11()
plot(S, type='l')
x11()
plot(density(S))
x11()
par(mfrow=c(2,2))
plot(fit2ar)

# save(S, file='sigma2.dat')
# save(mub, file='mub.dat')
# save(B, file='B_matrix.dat')

# Model 1

mean_intensity <- mean_intensity[2:13,]

mi <- lev <- non <- np <- t <- NA
for (p in 1:n_pixel)
  for (f in 1:12)
    if(f <= N_intervals[p])
    {
      mi <- c(mi,mean_intensity[f,p])
      non <- c(non, N_on[f,p])
      np <- c(np, N_pass[f,p])
      t <- c(t, T_raff2[f,p])
      lev <- c(lev,p)
    }

tab <- as.data.frame(cbind(mi,lev, non, np, t))
attach(tab)
# save(tab, file='table.dat')

fitaov <- aov(mi ~ lev)
s2_err <- sum(residuals(fitaov)^2)/(fitaov$df.residual) # variance within group
# save(s2_err, file='sigma0.dat')

fit <- lm (mi ~ non + np + t)
mub <- fit$coefficients
B <- vcov(fit)
# save(mub, file='mu.dat')
# save(B, file='matrice_B.dat')
