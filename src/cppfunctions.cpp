
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector linpredcompute(NumericMatrix X, const int nsites, const int p, 
                          NumericVector beta, NumericVector offset)
{
//Create new objects
// Compute the linear predictor
NumericVector linpred(nsites);
double temp; 


//  Compute the linear predictor via a double for loop
     for(int j = 0; j < nsites; j++)
     {
     temp = 0;
      
          for(int l = 0; l < p; l++) temp = temp + X(j,l) * beta[l];     
          
     linpred[j] = temp + offset[j];  
     }


// Return the result
return linpred;
}







// [[Rcpp::export]]
double binomialbetaupdate(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
               NumericVector proposal, NumericVector offset, NumericVector y, 
               NumericVector failures, NumericVector prior_meanbeta,
               NumericVector prior_varbeta)
{
// Compute the acceptance probability for beta
//Create new objects
double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit=0;
NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites);


// Create the log likelihood acceptance probability component
lp_current = linpredcompute(X, nsites, p, beta, offset);
lp_proposal = linpredcompute(X, nsites, p, proposal, offset);     
     for(int j = 0; j < nsites; j++)     
     {
     p_current[j] = exp(lp_current[j]) / (1 + exp(lp_current[j]));
     p_proposal[j] = exp(lp_proposal[j]) / (1 + exp(lp_proposal[j]));
     oldlikebit = oldlikebit + y[j] * log(p_current[j]) + failures[j] * log((1-p_current[j]));
     newlikebit = newlikebit + y[j] * log(p_proposal[j]) + failures[j] * log((1-p_proposal[j]));
     }
likebit = newlikebit - oldlikebit;


// Create the prior acceptance component
     for(int j = 0; j < p; j++)     
     {
     priorbit = priorbit + 0.5 * pow((beta[j]-prior_meanbeta[j]),2) / prior_varbeta[j] - 0.5 * pow((proposal[j]-prior_meanbeta[j]),2) / prior_varbeta[j];
     }


// Compute the acceptance probability and return the value
acceptance = exp(likebit + priorbit);
return acceptance;
}





// [[Rcpp::export]]
double quadform(IntegerVector W_duplet1, IntegerVector W_duplet2, const int n_duplet, const int nsites, 
                    NumericVector phi, NumericVector nneighbours, double diagonal, double offdiagonal)
{
// Compute a quadratic form for the random effects
// Create new objects 
double tau2_posteriorscale;
double tau2_quadform = 0, tau2_phisq = 0;
int row, col;

    
// Compute the off diagonal elements of the quadratic form
     for(int l = 0; l < n_duplet; l++)
     {
     row = W_duplet1[l] - 1;
     col = W_duplet2[l] - 1;
     tau2_quadform = tau2_quadform + phi[row] * phi[col]; 
     }
 
 
// Compute the diagonal elements of the quadratic form          
     for(int l = 0; l < nsites; l++)
     {
     tau2_phisq =  tau2_phisq + pow(phi[l],2) * (diagonal * nneighbours[l] + 1 - diagonal);
     }
           
     
// Compute the quadratic form
tau2_posteriorscale = 0.5 * (tau2_phisq - offdiagonal * tau2_quadform);

 
// Return the simulated value
return tau2_posteriorscale;
}






// [[Rcpp::export]]
List binomialcarupdate(List W_list, const int nsites, NumericVector phi, NumericVector nneighbours, 
                         double tau2, const NumericVector y, const NumericVector failures, 
                         const double phi_tune, double rho_num, double rho_den, 
                         NumericVector offset)
{
// Update the spatially correlated random effects 
//Create new objects
int accept=0;
double acceptance, sumphi;
double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
double priorvardenom, priormean, priorvar;
double propphi, pold, pnew;
NumericVector phinew(nsites);
 
   
//  Update each random effect in turn
phinew = phi;
     for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     IntegerVector neighbourvec = W_list[j];
     int m = neighbourvec.size();
     sumphi = 0;
          for(int l = 0; l < m; l++) sumphi += phinew[(neighbourvec[l]-1)];
      priorvardenom = (nneighbours[j] * rho_den + (1-rho_den));
      priorvar = tau2 / priorvardenom;
      priormean = rho_num * sumphi / priorvardenom; 
      
      // propose a value  
      propphi = rnorm(1, phinew[j], sqrt(priorvar*phi_tune))[0];
      
      // Accept or reject it
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
      pold = exp(offset[j] + phinew[j]) / (1 + exp(offset[j] + phinew[j]));
      pnew = exp(offset[j] + propphi) / (1 + exp(offset[j] + propphi));
      oldlikebit = y[j] * log(pold) + failures[j] * log((1-pold));
      newlikebit = y[j] * log(pnew) + failures[j] * log((1-pnew));
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
          if(runif(1)[0] <= acceptance) 
          {
          phinew[j] = propphi;
          accept = accept + 1;
          }
          else
          { 
          }
    }


List out(2);
out[0] = phinew;
out[1] = accept;
return out;
}



// [[Rcpp::export]]
List binomialindepupdate(const int nsites, NumericVector theta, double sigma2, const NumericVector y, 
               const NumericVector failures, const double theta_tune,  NumericVector offset)
{
// Update the spatially correlated random effects 
//Create new objects
int accept=0;
double acceptance;
double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
double proptheta, pold, pnew;
NumericVector thetanew(nsites);
 
   
//  Update each random effect in turn
thetanew = theta;
     for(int j = 0; j < nsites; j++)
     {
      // propose a value  
      proptheta = rnorm(1, thetanew[j], theta_tune)[0];
      
      // Accept or reject it
      newpriorbit = (0.5/sigma2) * pow(proptheta, 2); 
      oldpriorbit = (0.5/sigma2) * pow(thetanew[j], 2);
      pold = exp(offset[j] + thetanew[j]) / (1 + exp(offset[j] + thetanew[j]));
      pnew = exp(offset[j] + proptheta) / (1 + exp(offset[j] + proptheta));
      oldlikebit = y[j] * log(pold) + failures[j] * log((1-pold));
      newlikebit = y[j] * log(pnew) + failures[j] * log((1-pnew));
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
          if(runif(1)[0] <= acceptance) 
          {
          thetanew[j] = proptheta;
          accept = accept + 1;
          }
          else
          { 
          }
    }


List out(2);
out[0] = thetanew;
out[1] = accept;
return out;
}



// [[Rcpp::export]]
List poissonindepupdate(const int nsites, NumericVector theta, double sigma2, const NumericVector y, 
               const double theta_tune,  NumericVector offset)
{
// Update the spatially correlated random effects 
//Create new objects
int accept=0;
double acceptance;
double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
double proptheta, fittedold, fittednew;
NumericVector thetanew(nsites);
 
   
//  Update each random effect in turn
thetanew = theta;
     for(int j = 0; j < nsites; j++)
     {
      // propose a value  
      proptheta = rnorm(1, thetanew[j], theta_tune)[0];
      
      // Accept or reject it
      newpriorbit = (0.5/sigma2) * pow(proptheta, 2); 
      oldpriorbit = (0.5/sigma2) * pow(thetanew[j], 2);
      fittedold = exp(offset[j] + thetanew[j]);
      fittednew = exp(offset[j] + proptheta);
      oldlikebit = y[j] * log(fittedold) - fittedold;
      newlikebit = y[j] * log(fittednew) - fittednew;
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
          if(runif(1)[0] <= acceptance) 
          {
          thetanew[j] = proptheta;
          accept = accept + 1;
          }
          else
          { 
          }
    }


List out(2);
out[0] = thetanew;
out[1] = accept;
return out;
}


// [[Rcpp::export]]
List poissoncarupdate(List W_list, const int nsites, NumericVector phi, NumericVector nneighbours, 
                         double tau2, const NumericVector y, const double phi_tune, 
                         double rho_num, double rho_den, NumericVector offset)
{
// Update the spatially correlated random effects 
//Create new objects
int accept=0;
double acceptance, sumphi;
double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
double priorvardenom, priormean, priorvar;
double propphi, fittedold, fittednew;
NumericVector phinew(nsites);
 
   
//  Update each random effect in turn
phinew = phi;
     for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     IntegerVector neighbourvec = W_list[j];
     int m = neighbourvec.size();
     sumphi = 0;
          for(int l = 0; l < m; l++) sumphi += phinew[(neighbourvec[l]-1)];
      priorvardenom = (nneighbours[j] * rho_den + (1-rho_den));
      priorvar = tau2 / priorvardenom;
      priormean = rho_num * sumphi / priorvardenom; 
      
      // propose a value  
      propphi = rnorm(1, phinew[j], sqrt(priorvar*phi_tune))[0];
      
      // Accept or reject it
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
      fittedold = exp(offset[j] + phinew[j]);
      fittednew = exp(offset[j] + propphi);
      oldlikebit = y[j] * log(fittedold) - fittedold;
      newlikebit = y[j] * log(fittednew) - fittednew;
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
          if(runif(1)[0] <= acceptance) 
          {
          phinew[j] = propphi;
          accept = accept + 1;
          }
          else
          { 
          }
    }


List out(2);
out[0] = phinew;
out[1] = accept;
return out;
}




// [[Rcpp::export]]
double poissonbetaupdate(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
               NumericVector proposal, NumericVector offset, NumericVector y, 
               NumericVector prior_meanbeta, NumericVector prior_varbeta)
{
// Compute the acceptance probability for beta
//Create new objects
double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit=0;
NumericVector lp_current(nsites), lp_proposal(nsites), fitted_current(nsites), fitted_proposal(nsites);


// Create the log likelihood acceptance probability component
lp_current = linpredcompute(X, nsites, p, beta, offset);
lp_proposal = linpredcompute(X, nsites, p, proposal, offset);     
     for(int j = 0; j < nsites; j++)     
     {
     fitted_current[j] = exp(lp_current[j]);
     fitted_proposal[j] = exp(lp_proposal[j]);
     oldlikebit = oldlikebit + y[j] * log(fitted_current[j]) - fitted_current[j];
     newlikebit = newlikebit + y[j] * log(fitted_proposal[j]) - fitted_proposal[j];
     }
likebit = newlikebit - oldlikebit;


// Create the prior acceptance component
     for(int j = 0; j < p; j++)     
     {
     priorbit = priorbit + 0.5 * pow((beta[j]-prior_meanbeta[j]),2) / prior_varbeta[j] - 0.5 * pow((proposal[j]-prior_meanbeta[j]),2) / prior_varbeta[j];
     }


// Compute the acceptance probability and return the value
acceptance = exp(likebit + priorbit);
return acceptance;
}



// [[Rcpp::export]]
NumericVector gaussiancarupdate(List W_list, const int nsites, NumericVector phi, NumericVector nneighbours, 
                         double tau2, double rho_num, double rho_den, double nu2, NumericVector offset)
{
// Update the spatially correlated random effects 
//Create new objects
double sumphi;
double fcprecision, fcsd, fcmean;
double priorvardenom, priormean, priorvar;
NumericVector phinew(nsites);
 
   
//  Update each random effect in turn
phinew = phi;
     for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     IntegerVector neighbourvec = W_list[j];
     int m = neighbourvec.size();
     sumphi = 0;
          for(int l = 0; l < m; l++) sumphi += phinew[(neighbourvec[l]-1)];
      priorvardenom = (nneighbours[j] * rho_den + (1-rho_den));
      priorvar = tau2 / priorvardenom;
      priormean = rho_num * sumphi / priorvardenom; 
      
      // propose a value  
      fcprecision = (1/nu2) + (1/priorvar);
      fcsd = pow((1/fcprecision),0.5);
      fcmean = (priormean / priorvar + offset[j]) / fcprecision;
      phinew[j] = rnorm(1, fcmean, fcsd)[0];      
      }


return phinew;
}







// [[Rcpp::export]]
List poissondissimilaritycarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin,
                         const int nsites, NumericVector phi, double tau2, 
                         const NumericVector y, const double phi_tune, 
                         double rho, NumericVector offset)
{
// Update the spatially correlated random effects 
//Create new objects
int accept=0;
double acceptance, sumphi;
double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
double priorvardenom, priormean, priorvar;
double propphi, fittedold, fittednew;
NumericVector phinew(nsites);
int rowstart, rowend, rowtotal;
   
//  Update each random effect in turn
phinew = phi;


     for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     rowtotal = 0;
     sumphi = 0;
          for(int l = rowstart; l < rowend; l++) 
          {
          int r = Wtriplet(l,1) - 1;
          sumphi += Wtriplet(l, 2) * phinew[r];
          rowtotal += Wtriplet(l,2);
          }
     priorvardenom = (rowtotal * rho + (1-rho));
     priorvar = tau2 / priorvardenom;
     priormean = rho * sumphi / priorvardenom; 
     
      // propose a value  
      propphi = rnorm(1, phinew[j], sqrt(priorvar*phi_tune))[0];
      
      // Accept or reject it
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
      fittedold = exp(offset[j] + phinew[j]);
      fittednew = exp(offset[j] + propphi);
      oldlikebit = y[j] * log(fittedold) - fittedold;
      newlikebit = y[j] * log(fittednew) - fittednew;
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
          if(runif(1)[0] <= acceptance) 
          {
          phinew[j] = propphi;
          accept = accept + 1;
          }
          else
          { 
          }
    }




     
List out(2);
out[0] = phinew;
out[1] = accept;
return out;
}



// [[Rcpp::export]]
List binomialdissimilaritycarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin,
                         const int nsites, NumericVector phi, double tau2, 
                         const NumericVector y, const NumericVector failures,
                         const double phi_tune, double rho, NumericVector offset)
{
// Update the spatially correlated random effects 
//Create new objects
int accept=0;
double acceptance, sumphi;
double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
double priorvardenom, priormean, priorvar;
double propphi,  pold, pnew;
NumericVector phinew(nsites);
int rowstart, rowend, rowtotal;
   
//  Update each random effect in turn
phinew = phi;


     for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     rowtotal = 0;
     sumphi = 0;
          for(int l = rowstart; l < rowend; l++) 
          {
          int r = Wtriplet(l,1) - 1;
          sumphi += Wtriplet(l, 2) * phinew[r];
          rowtotal += Wtriplet(l,2);
          }
     priorvardenom = (rowtotal * rho + (1-rho));
     priorvar = tau2 / priorvardenom;
     priormean = rho * sumphi / priorvardenom; 
     
      // propose a value  
      propphi = rnorm(1, phinew[j], sqrt(priorvar*phi_tune))[0];
      
      // Accept or reject it
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
      pold = exp(offset[j] + phinew[j]) / (1 + exp(offset[j] + phinew[j]));
      pnew = exp(offset[j] + propphi) / (1 + exp(offset[j] + propphi));
      oldlikebit = y[j] * log(pold) + failures[j] * log((1-pold));
      newlikebit = y[j] * log(pnew) + failures[j] * log((1-pnew));
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
          if(runif(1)[0] <= acceptance) 
          {
          phinew[j] = propphi;
          accept = accept + 1;
          }
          else
          { 
          }
    }




     
List out(2);
out[0] = phinew;
out[1] = accept;
return out;
}




// [[Rcpp::export]]
NumericVector gaussiandissimilaritycarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin,
                         const int nsites, NumericVector phi, double tau2,  
                         double rho, double nu2, NumericVector offset)
{
// Update the spatially correlated random effects 
//Create new objects
double sumphi;
double fcprecision, fcsd, fcmean;
double priorvardenom, priormean, priorvar;
NumericVector phinew(nsites);
int rowstart, rowend, rowtotal;

//  Update each random effect in turn
phinew = phi;


     for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     rowtotal = 0;
     sumphi = 0;
          for(int l = rowstart; l < rowend; l++) 
          {
          int r = Wtriplet(l,1) - 1;
          sumphi += Wtriplet(l, 2) * phinew[r];
          rowtotal += Wtriplet(l,2);
          }
      priorvardenom = (rowtotal * rho + (1-rho));
      priorvar = tau2 / priorvardenom;
      priormean = rho * sumphi / priorvardenom; 
     
      // propose a value  
      fcprecision = (1/nu2) + (1/priorvar);
      fcsd = pow((1/fcprecision),0.5);
      fcmean = (priormean / priorvar + offset[j]) / fcprecision;
      phinew[j] = rnorm(1, fcmean, fcsd)[0];      
     }

return phinew;
}





// [[Rcpp::export]]
double quadformW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, const int n_triplet, const int nsites, 
                    NumericVector phi, double rho)
{
// Compute a quadratic form for the random effects
// Create new objects 
double tau2_posteriorscale;
double tau2_quadform = 0, tau2_phisq = 0;
int row, col;
int rowstart, rowend, rowtotal;
   
   
// Compute the off diagonal elements of the quadratic form
     for(int l = 0; l < n_triplet; l++)
     {
     row = Wtriplet(l,0) - 1;
     col = Wtriplet(l,1) - 1;
     tau2_quadform = tau2_quadform + phi[row] * phi[col] * Wtriplet(l,2); 
     }
 
 
 // Compute the diagonal elements of the quadratic form          
     for(int l = 0; l < nsites; l++)
     {
     rowstart = Wbegfin(l,0) - 1;
     rowend = Wbegfin(l,1);
     rowtotal = 0;
     for(int j = rowstart; j < rowend; j++) rowtotal += Wtriplet(j,2);
     tau2_phisq =  tau2_phisq + pow(phi[l],2) * (rho * rowtotal + 1 - rho);
     }
           
     
// Compute the quadratic form
tau2_posteriorscale = 0.5 * (tau2_phisq - rho * tau2_quadform);

 
// Return the simulated value
return tau2_posteriorscale;
}




// [[Rcpp::export]]
NumericVector binomialxupdate(const NumericVector y, const NumericVector failures, NumericVector offset,
               NumericVector x, const int nsites, NumericVector beta, NumericVector proposal){
  
  //Create new objects
  NumericVector xnew(nsites);
  int current, prop;
  double pold, pnew;
  double oldlikebit, newlikebit;  
  double acceptance;
  

  //  outer loop for each random effect
  xnew = x;
    for(int j = 0; j < nsites; j++){
       // Accept or reject the proposed value
     current = x[j] - 1;
     prop = proposal[j] - 1;
     pold = exp(offset[j] + beta[current]) / (1 + exp(offset[j] + beta[current]));
     pnew = exp(offset[j] + beta[prop]) / (1 + exp(offset[j] + beta[prop]));
     oldlikebit = y[j] * log(pold) + failures[j] * log((1-pold));
     newlikebit = y[j] * log(pnew) + failures[j] * log((1-pnew));
     acceptance = exp(newlikebit - oldlikebit);
      if(acceptance >= 1){
        xnew[j] = proposal[j];
      }else 
        if(runif(1)[0] <= acceptance) {
          xnew[j] = proposal[j];
        }
      else{ 
      }
    }
  
  return xnew;
}






// [[Rcpp::export]]
NumericVector poissonxupdate(const NumericVector y, NumericVector offset, NumericVector x, 
     const int nsites, NumericVector beta, NumericVector proposal)
     {
  
  //Create new objects
  NumericVector xnew(nsites);
  int current, prop;
  double muold, munew;
  double oldlikebit, newlikebit;  
  double acceptance;
  

  //  outer loop for each random effect
  xnew = x;
  
    for(int j = 0; j < nsites; j++)
    {
    // Accept or reject the proposed value
     current = x[j] - 1;
     prop = proposal[j] - 1;
     muold = offset[j] * exp(beta[current]);
     munew = offset[j] * exp(beta[prop]);
     oldlikebit = y[j] * log(muold) - muold;
     newlikebit = y[j] * log(munew) - munew;
     acceptance = exp(newlikebit - oldlikebit);
      if(acceptance >= 1){
        xnew[j] = proposal[j];
      }else 
        if(runif(1)[0] <= acceptance) {
          xnew[j] = proposal[j];
        }
      else{ 
      }
    }

  return xnew;
}




// [[Rcpp::export]]
NumericVector gaussianxupdate(const NumericVector y, NumericVector offset, double nu2,
     NumericVector x, const int nsites, NumericVector beta, NumericVector proposal)
     {
  
  //Create new objects
  NumericVector xnew(nsites);
  int current, prop;
  double oldlikebit, newlikebit;  
  double acceptance;
  

  //  outer loop for each random effect
  xnew = x;
    for(int j = 0; j < nsites; j++){
       // Accept or reject the proposed value
     current = x[j] - 1;
     prop = proposal[j] - 1;
     oldlikebit = pow((y[j] - offset[j] - beta[current]),2) / (2*nu2);
     newlikebit = pow((y[j] - offset[j] - beta[prop]),2) / (2*nu2);
     acceptance = exp(oldlikebit - newlikebit);
      if(acceptance >= 1){
        xnew[j] = proposal[j];
      }else 
        if(runif(1)[0] <= acceptance) {
          xnew[j] = proposal[j];
        }
      else{ 
      }
    }

  return xnew;
}
