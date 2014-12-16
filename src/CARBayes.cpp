#include <Rcpp.h>
using namespace Rcpp;



// This file contains the following functions:

// linpredcompute - computing the linear predictor for covariates.
// quadform - computing quadratic forms phi %*% Q %*% theta.
// binomialbetaupdate - update regression parameters in the binomial model
// binomialcarupdate - update random effects in the binomial model
// binomialindepupdate - update the independent effects in the binomial model
// poissonbetaupdate - update regression parameters in the poisson model
// poissoncarupdate - update random effects in the poisson model
// poissonindepupdate - update the independent effects in the poisson model
// gaussiancarupdate - update random effects in the Gaussian model


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
double quadform(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, const int nsites, 
                    NumericVector phi, NumericVector theta, double rho)
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
     tau2_quadform = tau2_quadform + phi[(Wtriplet(l,0) - 1)] * theta[(Wtriplet(l,1) - 1)] * Wtriplet(l,2); 
     }
 
 
 // Compute the diagonal elements of the quadratic form          
     for(int l = 0; l < nsites; l++)
     {
     tau2_phisq = tau2_phisq + phi[l] * theta[l] * (rho * Wtripletsum[l] + 1 - rho);    
     }
           
     
// Compute the quadratic form
tau2_posteriorscale = 0.5 * (tau2_phisq - rho * tau2_quadform);

 
// Return the simulated value
return tau2_posteriorscale;
}




// [[Rcpp::export]]
List binomialcarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum,const int nsites, NumericVector phi, double tau2, 
     const NumericVector y, const NumericVector failures, const double phi_tune, 
     double rho, NumericVector offset)
{
// Update the spatially correlated random effects 
//Create new objects
int accept=0, rowstart=0, rowend=0;
double acceptance, sumphi;
double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
double priorvardenom, priormean, priorvar;
double propphi, pold, pnew;
NumericVector phinew(nsites);


//  Update each random effect in turn
phinew = phi;
     for(int j = 0; j < nsites; j++)
     {
     // Calculate prior variance
     priorvardenom = rho * Wtripletsum[j] + 1 - rho;
     priorvar = tau2 / priorvardenom;
     
     // Calculate the prior mean
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     sumphi = 0;
          for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
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
List poissoncarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, NumericVector phi, 
     double tau2, const NumericVector y, const double phi_tune, 
     double rho, NumericVector offset)
{
// Update the spatially correlated random effects 
//Create new objects
int accept=0,rowstart=0, rowend=0;
double acceptance, sumphi;
double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
double priorvardenom, priormean, priorvar;
double propphi, lpold, lpnew;
NumericVector phinew(nsites);
   
   
//  Update each random effect in turn
phinew = phi;

     for(int j = 0; j < nsites; j++)
     {
     // Calculate prior variance
     priorvardenom = rho * Wtripletsum[j] + 1 - rho;
     priorvar = tau2 / priorvardenom;
     
     // Calculate the prior mean
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     sumphi = 0;
          for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
     priormean = rho * sumphi / priorvardenom; 
     
      // propose a value  
      propphi = rnorm(1, phinew[j], sqrt(priorvar*phi_tune))[0];
      
      // Accept or reject it
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
      lpold = offset[j] + phinew[j];
      lpnew = offset[j] + propphi;
      oldlikebit = y[j] * lpold - exp(lpold);
      newlikebit = y[j] * lpnew - exp(lpnew);
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
NumericVector gaussiancarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, NumericVector phi, double tau2, 
     double rho, double nu2, NumericVector offset)
{
// Update the spatially correlated random effects 
//Create new objects
int rowstart=0, rowend=0;
double sumphi;
double fcprecision, fcsd, fcmean;
double priorvardenom, priormean, priorvar;
NumericVector phinew(nsites);


//  Update each random effect in turn
phinew = phi;
     for(int j = 0; j < nsites; j++)
     {
     // Calculate prior variance
     priorvardenom = rho * Wtripletsum[j] + 1 - rho;
     priorvar = tau2 / priorvardenom;
     
     // Calculate the prior mean
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     sumphi = 0;
          for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
     priormean = rho * sumphi / priorvardenom; 
      
      // propose a value  
      fcprecision = (1/nu2) + (1/priorvar);
      fcsd = pow((1/fcprecision),0.5);
      fcmean = (priormean / priorvar + offset[j]) / fcprecision;
      phinew[j] = rnorm(1, fcmean, fcsd)[0];      
      }


return phinew;
}
