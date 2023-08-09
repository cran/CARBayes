#include <Rcpp.h>
using namespace Rcpp;

// This file contains the following functions:

// linpredcompute - computing the linear predictor for covariates.
// quadform - computing quadratic forms phi %*% Q %*% theta.
// binomialbetaupdateMALA - update regression parameters in the binomial model using MALA
// binomialbetaupdateRW - update regression parameters in the binomial model using RW
// binomialcarupdateRW - update random effects in the binomial model using RW
// binomialindepupdateRW - update the independent effects in the binomial model using RW
// poissonbetaupdateMALA - update regression parameters in the poisson model using MALA
// poissonbetaupdateRW - update regression parameters in the poisson model using RW
// poissoncarupdateRW - update random effects in the poisson model using RW
// poissonindepupdateRW - update the independent effects in the poisson model using RW
// zipcarupdateRW - update the random effects in the zip models using RW
// zipindepupdateRW - update the independent random effects in the zip models using RW
// gaussiancarupdate - update random effects in the Gaussian model
// binomialmcarupdateRW - update random effects in the binomial MCAR model using RW
// poissonmcarupdateRW - update random effects in the poisson MCAR model using RW
// gaussianmcarupdateRW - update random effecs in the Gaussian MCAR model usng RW
// multinomialbetaupdateRW - update beta in the multinomial model using RW
// poissoncarmultilevelupdate - Poisson spatial random effects updates
// binomialcarmultilevelupdate - binomial spatial random effects updates
// gaussiancarmultilevelupdate - gaussian spatial random effects updates
// gaussiancarmultilevelupdateindiv - Gaussian indep random effect updates
// poissoncarmultilevelupdateindiv - Poisson indep random effect updates
// binomialcarmultilevelupdateindiv - binomial indep random effect updates



///////////////////////////////////////////////////////////////////////////////
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



///////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double quadform(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, const int nsites, 
                    NumericVector phi, NumericVector theta, double rho)
{
// Compute a quadratic form for the random effects
// Create new objects 
double tau2_posteriorscale;
double tau2_quadform = 0, tau2_phisq = 0;

// Compute the off diagonal elements of the quadratic form
     for(int l = 0; l < n_triplet; l++)
     {
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



//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List binomialindepupdateRW(const int nsites, NumericVector theta, double sigma2, const NumericVector y, 
                         const NumericVector failures, const double theta_tune,  NumericVector offset)
{
    // Update the independent random effects 
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
    // Full conditional ratio
    newpriorbit = (0.5/sigma2) * pow(proptheta, 2); 
    oldpriorbit = (0.5/sigma2) * pow(thetanew[j], 2);
            
    pold = exp(offset[j] + thetanew[j]) / (1 + exp(offset[j] + thetanew[j]));
    pnew = exp(offset[j] + proptheta) / (1 + exp(offset[j] + proptheta));
    oldlikebit = y[j] * log(pold) + failures[j] * log((1-pold));
    newlikebit = y[j] * log(pnew) + failures[j] * log((1-pnew));
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            

    // Acceptace or reject the proposal
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



///////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List binomialcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
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
    double propphi, pold, pnew, proposal_var;
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
        proposal_var = priorvar * phi_tune;
        propphi = rnorm(1, phinew[j], sqrt(proposal_var))[0];
            
        // Accept or reject it
        // Full conditional ratio
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
        pold = exp(offset[j] + phinew[j]) / (1 + exp(offset[j] + phinew[j]));
        pnew = exp(offset[j] + propphi) / (1 + exp(offset[j] + propphi));
        oldlikebit = y[j] * log(pold) + failures[j] * log((1-pold));
        newlikebit = y[j] * log(pnew) + failures[j] * log((1-pnew));
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            

            // Acceptace or reject the proposal
            if(runif(1)[0] <= acceptance) 
            {
                phinew[j] = propphi;
                accept = accept + 1;
            }
            else
            { 
            }    
        }
    
    // Return the results
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}



//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List binomialbetaupdateRW(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                          NumericVector offset, NumericVector y,  NumericVector failures,
                          NumericVector prior_meanbeta, NumericVector prior_varbeta, 
                          const int nblock,double beta_tune, List block_list)
{
  // Compute the acceptance probability for beta
  //Create new objects
  int accept=0;
  double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
  double acceptance;
  NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites);
  
  // Create two beta vectors
  NumericVector beta_old(p);
  NumericVector beta_new(p);
  for(int g=0; g<p; g++)
  {
    beta_old[g] = beta[g];
    beta_new[g] = beta[g];
  }
  
  // Update each block in turn
  for(int r=0; r<nblock; r++)
  {
    // Determine the block to update
    IntegerVector idx = block_list[r];
    int len = block_list[(nblock+r)];
    
    // Propose a value
    for(int g=0; g<len; g++)
    {
      beta_new[idx[g]] = rnorm(1, beta_old[idx[g]], beta_tune)[0];
    }
    
    
    // Compute the acceptance ratio - full conditionals  
    oldlikebit = 0;
    newlikebit=0;
    lp_current = linpredcompute(X, nsites, p, beta_old, offset);
    lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
    for(int j = 0; j < nsites; j++)     
    {
      p_current[j] = exp(lp_current[j]) / (1 + exp(lp_current[j]));
      p_proposal[j] = exp(lp_proposal[j]) / (1 + exp(lp_proposal[j]));
      oldlikebit = oldlikebit +  y[j] * log(p_current[j]) + failures[j] * log((1-p_current[j]));
      newlikebit = newlikebit +  y[j] * log(p_proposal[j]) + failures[j] * log((1-p_proposal[j]));
    }
    likebit = newlikebit - oldlikebit;
    
    priorbit = 0;
    for(int g = 0; g < len; g++)     
    {
      priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
    }
    
    
    // Accept or reject hte proposal      
    acceptance = exp(likebit + priorbit);
    if(runif(1)[0] <= acceptance) 
    {
      for(int g=0; g<len; g++)
      {
        beta_old[idx[g]] = beta_new[idx[g]];  
      }
      accept = accept + 1;
    }
    else
    { 
      for(int g=0; g<len; g++)
      {
        beta_new[idx[g]] = beta_old[idx[g]];  
      }   
    }
  }
  
  
  // Compute the acceptance probability and return the value
  List out(2);
  out[0] = beta_new;
  out[1] = accept;
  return out;    
}



//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List binomialbetaupdateMALA(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                            NumericVector offset, NumericVector y,  NumericVector failures,
                            NumericVector trials, NumericVector prior_meanbeta, 
                            NumericVector prior_varbeta, const int nblock,
                            double beta_tune, List block_list)
{
    // Compute the acceptance probability for beta
    //Create new objects
    int accept=0;
    double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
    double acceptance;
    NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites), mala_temp1(nsites);
    
    // Create two beta vectors
    NumericVector beta_old(p);
    NumericVector beta_new(p);
    for(int g=0; g<p; g++)
    {
        beta_old[g] = beta[g];
        beta_new[g] = beta[g];
    }
    
    // Update each block in turn
    for(int r=0; r<nblock; r++)
    {
        // Determine the block to update
        IntegerVector idx = block_list[r];
        int len = block_list[(nblock+r)];
        
        // Propose a value
        lp_current = linpredcompute(X, nsites, p, beta_old, offset);
        mala_temp1 = y -  trials * exp(lp_current) / (1 + exp(lp_current));
        NumericVector mala_temp2(len), mala_old(len);
        for(int g=0; g<len; g++)
        {
            mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
            mala_old[g] = beta_old[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_old[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
            beta_new[idx[g]] = rnorm(1, mala_old[g], beta_tune)[0];
        }
        
        // Compute the acceptance ratio - full conditionals  
        oldlikebit = 0;
        newlikebit=0;
        lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
        for(int j = 0; j < nsites; j++)     
        {
            p_current[j] = exp(lp_current[j]) / (1 + exp(lp_current[j]));
            p_proposal[j] = exp(lp_proposal[j]) / (1 + exp(lp_proposal[j]));
            oldlikebit = oldlikebit + y[j] * log(p_current[j]) + failures[j] * log((1-p_current[j]));
            newlikebit = newlikebit + y[j] * log(p_proposal[j]) + failures[j] * log((1-p_proposal[j]));
        }
        likebit = newlikebit - oldlikebit;
        
        
        for(int g = 0; g < len; g++)     
        {
            priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
        }
        
        // Compute the acceptance ratio - proposal distributions
        mala_temp1 = y -  trials * exp(lp_proposal) / (1 + exp(lp_proposal));
        NumericVector mala_new(len);
        double prop_accept=0;
        for(int g=0; g<len; g++)
        {
            mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
            mala_new[g] = beta_new[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_new[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
            prop_accept = prop_accept +   pow((beta_new[idx[g]] - mala_old[g]), 2) -  pow((beta_old[idx[g]] - mala_new[g]), 2); 
        }
        
        // Accept or reject hte proposal      
        acceptance = exp(0.5 * prop_accept / pow(beta_tune,2) + likebit + priorbit);
        if(runif(1)[0] <= acceptance) 
        {
            for(int g=0; g<len; g++)
            {
                beta_old[idx[g]] = beta_new[idx[g]];  
            }
            accept = accept + 1;
        }
        else
        { 
            for(int g=0; g<len; g++)
            {
                beta_new[idx[g]] = beta_old[idx[g]];  
            }   
        }
    }
    
    
    // Compute the acceptance probability and return the value
    //acceptance = exp(likebit + priorbit);
    List out(2);
    out[0] = beta_new;
    out[1] = accept;
    return out;    
}



//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List poissonindepupdateRW(const int nsites, NumericVector theta, double sigma2, const NumericVector y, 
                            const double theta_tune,  NumericVector offset)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0;
    double acceptance;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double proptheta, lpold, lpnew;
    NumericVector thetanew(nsites);
    
    
    //  Update each random effect in turn
    thetanew = theta;
    for(int j = 0; j < nsites; j++)
    {
    // propose a value
    proptheta = rnorm(1, thetanew[j], theta_tune)[0];
    
    // Accept or reject it
    // Full conditional ratio
    newpriorbit = (0.5/sigma2) * pow(proptheta, 2); 
    oldpriorbit = (0.5/sigma2) * pow(thetanew[j], 2);
    lpold = offset[j] + thetanew[j];
    lpnew = offset[j] + proptheta;
    oldlikebit = y[j] * lpold - exp(lpold);
    newlikebit = y[j] * lpnew - exp(lpnew);
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);

    // Acceptace or reject the proposal
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



//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List poissoncarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                          NumericVector Wtripletsum, const int nsites, NumericVector phi, 
                          double tau2, const NumericVector y, const double phi_tune, 
                          double rho, NumericVector offset)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0,rowstart=0, rowend=0;
    double acceptance, sumphi, proposal_var;
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
        proposal_var = priorvar * phi_tune;
        propphi = rnorm(1, phinew[j], sqrt(proposal_var))[0];
            
        // Accept or reject it
        // Full conditional ratio
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
        lpold = offset[j] + phinew[j];
        lpnew = offset[j] + propphi;
        oldlikebit = y[j] * lpold - exp(lpold);
        newlikebit = y[j] * lpnew - exp(lpnew);
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            
        // Acceptace or reject the proposal
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



//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List poissonbetaupdateMALA(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                       NumericVector offset, NumericVector y, NumericVector prior_meanbeta, 
                       NumericVector prior_varbeta, const int nblock,
                       double beta_tune, List block_list)
{
    // Compute the acceptance probability for beta
    //Create new objects
    int accept=0;
    double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
    double acceptance;
    NumericVector lp_current(nsites), lp_proposal(nsites), mala_temp1(nsites);
    
    // Create two beta vectors
    NumericVector beta_old(p);
    NumericVector beta_new(p);
    for(int g=0; g<p; g++)
    {
        beta_old[g] = beta[g];
        beta_new[g] = beta[g];
    }
    
    // Update each block in turn
    for(int r=0; r<nblock; r++)
    {
        // Determine the block to update
        IntegerVector idx = block_list[r];
        int len = block_list[(nblock+r)];
        
        // Propose a value
        lp_current = linpredcompute(X, nsites, p, beta_old, offset);
        mala_temp1 = y - exp(lp_current);
        NumericVector mala_temp2(len), mala_old(len);
        for(int g=0; g<len; g++)
        {
            mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
            mala_old[g] = beta_old[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_old[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
            beta_new[idx[g]] = rnorm(1, mala_old[g], beta_tune)[0];
        }
        
        // Compute the acceptance ratio - full conditionals 
        oldlikebit = 0;
        newlikebit=0;
        lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
        for(int j = 0; j < nsites; j++)     
        {
            oldlikebit = oldlikebit + y[j] * lp_current[j] - exp(lp_current[j]);
            newlikebit = newlikebit + y[j] * lp_proposal[j] - exp(lp_proposal[j]);
        }
        likebit = newlikebit - oldlikebit;
        
        for(int g = 0; g < len; g++)     
        {
            priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
        }
        
        // Compute the acceptance ratio - proposal distributions
        mala_temp1 = y - exp(lp_proposal);
        NumericVector mala_new(len);
        double prop_accept=0;
        for(int g=0; g<len; g++)
        {
            mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
            mala_new[g] = beta_new[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_new[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
            prop_accept = prop_accept +   pow((beta_new[idx[g]] - mala_old[g]), 2) -  pow((beta_old[idx[g]] - mala_new[g]), 2); 
        }
        
        // Accept or reject hte proposal      
        acceptance = exp(0.5 * prop_accept / pow(beta_tune,2) + likebit + priorbit);
        if(runif(1)[0] <= acceptance) 
        {
            for(int g=0; g<len; g++)
            {
                beta_old[idx[g]] = beta_new[idx[g]];  
            }
            accept = accept + 1;
        }
        else
        { 
            for(int g=0; g<len; g++)
            {
                beta_new[idx[g]] = beta_old[idx[g]];  
            }   
        }
    }
    
    
    
    // Compute the acceptance probability and return the value
    //acceptance = exp(likebit + priorbit);
    List out(2);
    out[0] = beta_new;
    out[1] = accept;
    return out;    
}



//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List poissonbetaupdateRW(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                         NumericVector offset, NumericVector y, NumericVector prior_meanbeta, 
                         NumericVector prior_varbeta, const int nblock, double beta_tune, 
                         List block_list)
{
  // Compute the acceptance probability for beta
  //Create new objects
  int accept=0;
  double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
  double acceptance;
  NumericVector lp_current(nsites), lp_proposal(nsites);
  
  // Create two beta vectors
  NumericVector beta_old(p);
  NumericVector beta_new(p);
  for(int g=0; g<p; g++)
  {
    beta_old[g] = beta[g];
    beta_new[g] = beta[g];
  }

  
// Update each block in turn
  for(int r=0; r<nblock; r++)
  {
    // Determine the block to update
    IntegerVector idx = block_list[r];
    int len = block_list[(nblock+r)];
    
    // Propose a value
    for(int g=0; g<len; g++)
    {
      beta_new[idx[g]] = rnorm(1, beta_old[idx[g]], beta_tune)[0];
    }

    
    // Compute the acceptance ratio - likelihood part
    lp_current = linpredcompute(X, nsites, p, beta_old, offset);
    lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
    oldlikebit = 0;
    newlikebit=0;
    for(int j = 0; j < nsites; j++)     
    {
      oldlikebit = oldlikebit + y[j] * lp_current[j] - exp(lp_current[j]);
      newlikebit = newlikebit + y[j] * lp_proposal[j] - exp(lp_proposal[j]);
    }
    likebit = newlikebit - oldlikebit;
    
    
    // Compute the acceptance ratio - prior part
    priorbit = 0;
    for(int g = 0; g < len; g++)     
    {
      priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
    }

    
    // Accept or reject the proposal    
    acceptance = exp(likebit + priorbit);
    if(runif(1)[0] <= acceptance) 
    {
      for(int g=0; g<len; g++)
      {
        beta_old[idx[g]] = beta_new[idx[g]];  
      }
      accept = accept + 1;
    }
    else
    { 
      for(int g=0; g<len; g++)
      {
        beta_new[idx[g]] = beta_old[idx[g]];  
      }   
    }
  }


  // Return the value
  List out(2);
  out[0] = beta_new;
  out[1] = accept;
  return out;    
}



//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List zipcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                    NumericVector Wtripletsum, const int nsites, NumericVector phi, 
                    double tau2, const NumericVector y, const double phi_tune, 
                    double rho, NumericVector offset, NumericVector poiind)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0,rowstart=0, rowend=0;
    double acceptance, sumphi, proposal_var;
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
        
        // Different updates depending on whether the y[j] is missing or not.
        if(poiind[j]==1)
        {
            // propose a value
            proposal_var = priorvar * phi_tune;
            propphi = rnorm(1, phinew[j], sqrt(proposal_var))[0];
            
            // Accept or reject it
            // Full conditional ratio
            newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
            oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
            lpold = offset[j] + phinew[j];
            lpnew = offset[j] + propphi;
            oldlikebit = y[j] * lpold - exp(lpold);
            newlikebit = y[j] * lpnew - exp(lpnew);
            acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            
            // Acceptace or reject the proposal
            if(runif(1)[0] <= acceptance) 
            {
                phinew[j] = propphi;
                accept = accept + 1;
            }
            else
            { 
            }    
        }else
        {
            phinew[j] = rnorm(1, priormean, sqrt(priorvar))[0];    
        }
    }
    
    
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}



//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List zipindepupdateRW(const int nsites, NumericVector theta, double sigma2, const NumericVector y, 
                      const double theta_tune, NumericVector offset, NumericVector poiind)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0;
    double acceptance;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double proptheta, lpold, lpnew;
    NumericVector thetanew(nsites);
    
    
    //  Update each random effect in turn
    thetanew = theta;
    
    for(int j = 0; j < nsites; j++)
    {
        // Different updates depending on whether the y[j] is missing or not.
        if(poiind[j]==1)
        {
            // propose a value
            proptheta = rnorm(1, thetanew[j], theta_tune)[0];
            
            // Accept or reject it
            // Full conditional ratio
            newpriorbit = (0.5/sigma2) * pow(proptheta, 2); 
            oldpriorbit = (0.5/sigma2) * pow(thetanew[j], 2);
            lpold = offset[j] + thetanew[j];
            lpnew = offset[j] + proptheta;
            oldlikebit = y[j] * lpold - exp(lpold);
            newlikebit = y[j] * lpnew - exp(lpnew);
            acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            
            // Acceptace or reject the proposal
            if(runif(1)[0] <= acceptance) 
            {
                thetanew[j] = proptheta;
                accept = accept + 1;
            }
            else
            { 
            }    
        }else
        {
            thetanew[j] = rnorm(1, 0, sqrt(sigma2))[0];    
        }
    }
    
    
    List out(2);
    out[0] = thetanew;
    out[1] = accept;
    return out;
}



//////////////////////////////////////////////////////////////////////////////
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







//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List binomialmcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                        const int nsites,  const int nvar, NumericMatrix phi, 
                        NumericMatrix Y, NumericMatrix failures,
                        NumericMatrix phioffset, NumericVector denoffset,  
                        NumericMatrix Sigmainv, double rho, double phi_tune, 
                        NumericMatrix innovations)
{
    // Update the spatially correlated random effects 
    //Create new objects
    NumericMatrix fcprec(nvar, nvar);
    int rowstart=0, rowend=0, accept=0;
    NumericVector sumphi(nvar), fcmean(nvar), propphi(nvar);
    NumericVector diffcurrent(nvar), diffprop(nvar);        
    NumericVector quadcurrent(nvar), quadprop(nvar);  
    NumericVector lpold(nvar), lpnew(nvar), pold(nvar), pnew(nvar);
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;
    NumericMatrix phinew(nsites, nvar);
    phinew = clone(phi);
    
  
    //  Update each random effect in turn
    for(int j = 0; j < nsites; j++)
    {      
        // Calculate the prior precision
        for(int r=0; r<nvar; r++)
        {
        fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
        }

        // Calculate the prior mean
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = rep(0,nvar);
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew((Wtriplet(l,1) - 1),_);
        fcmean = rho * sumphi / denoffset[j]; 
 
        // Generate the proposal distribution mean and propose a value
        for(int r=0; r<nvar; r++)
        {
        propphi[r] = phinew(j,r) + innovations(j, r);
        }

        // Compute the prior ratio
        diffcurrent = phinew(j,_) - fcmean;
        diffprop = propphi - fcmean;
            for(int r=0; r<nvar; r++)
            {
            quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
            quadprop[r] = sum(diffprop * fcprec(_,r));  
            }
        oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
        newpriorbit = 0.5 * sum(quadprop * diffprop);      
        
        // Likelihood ratio
        lpold = phioffset(j,_) + phinew(j,_);
        lpnew = phioffset(j,_) + propphi;
        pold = exp(lpold) / (1 + exp(lpold));
        pnew = exp(lpnew) / (1 + exp(lpnew));
        oldlikebit = sum(Y(j,_) * log(pold) + failures(j,_) * log(1 - pold));
        newlikebit = sum(Y(j,_) * log(pnew) + failures(j,_) * log(1 - pnew));
        

        // Accept or reject the value
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            phinew(j,_) = propphi;
            accept = accept + 1;
        }
        else
        { 
        }
    }     
    
    
    // Return the results
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}






//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List poissonmcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                           const int nsites,  const int nvar, NumericMatrix phi, 
                           NumericMatrix Y, NumericMatrix phioffset, 
                           NumericVector denoffset, NumericMatrix Sigmainv, double rho, 
                           double phi_tune, NumericMatrix innovations)
{
    // Update the spatially correlated random effects 
    //Create new objects
    NumericMatrix fcprec(nvar, nvar);
    int rowstart=0, rowend=0, accept=0;
    NumericVector sumphi(nvar), fcmean(nvar), propphi(nvar);
    NumericVector diffcurrent(nvar), diffprop(nvar);        
    NumericVector quadcurrent(nvar), quadprop(nvar);  
    NumericVector lpold(nvar), lpnew(nvar);
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;
    NumericMatrix phinew(nsites, nvar);
    phinew = clone(phi);
    
    //  Update each random effect in turn
    for(int j = 0; j < nsites; j++)
    {     
        // Calculate the prior precision and mean
        for(int r=0; r<nvar; r++)
        {
        fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
        }
        
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = rep(0,nvar);
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew((Wtriplet(l,1) - 1),_);
        fcmean = rho * sumphi / denoffset[j]; 
        
        // Generate the proposal distribution mean and propose a value
            for(int r=0; r<nvar; r++)
            {
            propphi[r] = phinew(j,r) + innovations(j, r);
            }
                
        // Compute the prior ratio
        diffcurrent = phinew(j,_) - fcmean;
        diffprop = propphi - fcmean;
            for(int r=0; r<nvar; r++)
            {
            quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
            quadprop[r] = sum(diffprop * fcprec(_,r));  
            }
        oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
        newpriorbit = 0.5 * sum(quadprop * diffprop);      
        
        // Likelihood ratio
        lpold = phioffset(j,_) + phinew(j,_);
        lpnew = phioffset(j,_) + propphi;
        oldlikebit = sum(Y(j,_) * lpold - exp(lpold));
        newlikebit = sum(Y(j,_) * lpnew - exp(lpnew));
        

        // Accept or reject the value
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            phinew(j,_) = propphi;
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



//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List gaussianmcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                         const int nsites,  const int nvar, NumericMatrix phi, 
                         NumericMatrix phioffset, NumericVector denoffset, 
                         NumericMatrix Sigmainv, double rho, NumericVector nu2,
                         double phi_tune, NumericMatrix innovations)
{
    // Update the spatially correlated random effects 
    //Create new objects
    NumericMatrix fcprec(nvar, nvar);
    int rowstart=0, rowend=0, accept=0;
    NumericVector sumphi(nvar), fcmean(nvar), propphi(nvar);
    NumericVector diffcurrent(nvar), diffprop(nvar);        
    NumericVector quadcurrent(nvar), quadprop(nvar);  
    NumericVector lpold(nvar), lpnew(nvar);
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;
    NumericMatrix phinew(nsites, nvar);
    phinew = clone(phi);
    
    //  Update each random effect in turn
    for(int j = 0; j < nsites; j++)
    {     
        // Calculate the prior precision and mean
        for(int r=0; r<nvar; r++)
        {
            fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
        }
        
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = rep(0,nvar);
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew((Wtriplet(l,1) - 1),_);
        fcmean = rho * sumphi / denoffset[j]; 
        
        // Generate a proposal value
        for(int r=0; r<nvar; r++)
        {
        propphi[r] = phinew(j,r) + innovations(j, r);
        }
                    
                    
        // Compute the prior ratio
        diffcurrent = phinew(j,_) - fcmean;
        diffprop = propphi - fcmean;
            for(int r=0; r<nvar; r++)
            {
            quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
            quadprop[r] = sum(diffprop * fcprec(_,r));  
            }
        oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
        newpriorbit = 0.5 * sum(quadprop * diffprop);      
        
        // Likelihood ratio
        lpold = pow((phioffset(j,_) - phinew(j,_)),2);
        lpnew = pow((phioffset(j,_) - propphi),2);
        oldlikebit = 0.5 * sum(lpold / nu2);
        newlikebit = 0.5 * sum(lpnew / nu2);                             

        // Accept or reject the value
        acceptance = exp(oldpriorbit - newpriorbit + oldlikebit - newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            phinew(j,_) = propphi;
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



//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List multinomialbetaupdateRW(NumericMatrix X, const int nsites, const int J, const int p, 
                             const int col, NumericMatrix beta, NumericMatrix offset, NumericMatrix y,                            
                             NumericVector prior_meanbeta, NumericVector prior_varbeta, 
                             const int nblock, double beta_tune, List block_list, NumericVector zeros)
{
    // Compute the acceptance probability for beta
    //Create new objects
    NumericMatrix lp_current(nsites, J), lp_proposal(nsites, J);
    NumericVector p_current(nsites), p_proposal(nsites);
    int accept=0;
    double oldlikebit=0, newlikebit=0, priorbit=0;
    double acceptance;
    
    
    // Create a beta old and new matrix
    NumericMatrix beta_new(p, (J-1));
    for(int j = 0; j < (J-1); j++)
    {
        beta_new(_,j) = beta(_,j);
    }
    
    
    // Update each block in turn
    for(int r=0; r<nblock; r++)
    {
        // Determine the block to update
        IntegerVector idx = block_list[r];
        int len = block_list[(nblock+r)];
        
        // Propose a value
        for(int g=0; g<len; g++)
        {
            beta_new(idx[g],(col-1)) = rnorm(1, beta(idx[g],(col-1)), beta_tune)[0];
        }
        
        // Compute the linear predictors
        lp_current(_, 0) = zeros;
        lp_proposal(_, 0) = zeros;   
        
        for(int g=1; g<J; g++)
        {
            lp_current(_,g) = linpredcompute(X, nsites, p, beta(_, (g-1)), offset(_, (g-1)));    
            lp_proposal(_,g) = linpredcompute(X, nsites, p, beta_new(_, (g-1)), offset(_, (g-1)));    
        }
        
        // Compute the probabilities and the likelihood component of the MH step
        oldlikebit = 0;
        newlikebit=0;
        for(int j = 0; j < nsites; j++)     
        {
            p_current = exp(lp_current(j, _)) / sum(exp(lp_current(j, _)));
            p_proposal = exp(lp_proposal(j, _)) / sum(exp(lp_proposal(j, _))); 
            oldlikebit = oldlikebit + sum(y(j, _) * log(p_current));
            newlikebit = newlikebit + sum(y(j, _) * log(p_proposal));
        }
        
        // Compute the prior component of the MH step  
        for(int g = 0; g < len; g++)     
        {
            priorbit = priorbit + 0.5 * pow((beta(idx[g],(col-1))-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new(idx[g],(col-1))-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
        }
        
        // Accept or reject the proposal      
        acceptance = exp(newlikebit - oldlikebit + priorbit);
        if(runif(1)[0] <= acceptance) 
        {
            for(int g=0; g<len; g++)
            {
                beta(idx[g], (col-1)) = beta_new(idx[g], (col-1));  
            }
            accept = accept + 1;
        }
        else
        { 
            for(int g=0; g<len; g++)
            {
                beta_new(idx[g], (col-1)) = beta(idx[g], (col-1));  
            }   
        }
    }
    
    
    // Compute the acceptance probability and return the value
    List out(2);
    out[0] = beta_new;
    out[1] = accept;
    return out;    
}



//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List multinomialmcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                          const int nsites,  const int nvar, NumericMatrix phi, 
                          NumericMatrix Y, NumericMatrix phioffset, 
                          NumericVector denoffset,  NumericMatrix Sigmainv, 
                          double rho, double phi_tune, NumericMatrix innovations)
{
    // Update the spatially correlated random effects 
    //Create new objects
    NumericMatrix fcprec((nvar-1), (nvar-1));
    int rowstart=0, rowend=0, accept=0;
    NumericVector sumphi((nvar-1)), fcmean((nvar-1)), propphi((nvar-1));
    NumericVector diffcurrent(nvar), diffprop(nvar);        
    NumericVector quadcurrent(nvar), quadprop(nvar);  
    NumericVector lpold(nvar), lpnew(nvar), pold(nvar), pnew(nvar);
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;
    NumericMatrix phinew(nsites, nvar);
    phinew = clone(phi);
    
    //  Update each random effect in turn
    for(int j = 0; j < nsites; j++)
    {
    // Calculate the prior precision and mean
        for(int r=0; r<(nvar-1); r++)
        {
        fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
        }
        
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    sumphi = rep(0,(nvar-1));
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew((Wtriplet(l,1) - 1),_);
        fcmean = rho * sumphi / denoffset[j]; 
        
        // Propose a possible value
        for(int r=0; r<(nvar-1); r++)
        {
        propphi[r] = phinew(j,r) + innovations(j, r);
        }

    // Compute the prior ratio
    diffcurrent = phinew(j,_) - fcmean;
    diffprop = propphi - fcmean;
        for(int r=0; r<(nvar-1); r++)
        {
        quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
        quadprop[r] = sum(diffprop * fcprec(_,r));  
        }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);      
        
    // Likelihood ratio
    // compute the linear predictor
    lpold[0] = 0;
    lpnew[0] = 0;   
        for(int g=1; g<nvar; g++)
        {
        lpold[g] =  phinew(j, (g-1))  + phioffset(j, (g-1));
        lpnew[g] =  propphi[(g-1)]  + phioffset(j, (g-1));
        }
    
    // Compute the probabilities and the likelihood component of the MH step
    pold = exp(lpold) / sum(exp(lpold));
    pnew = exp(lpnew) / sum(exp(lpnew));; 
    oldlikebit = sum(Y(j, _) * log(pold));
    newlikebit = sum(Y(j, _) * log(pnew));

    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
        phinew(j,_) = propphi;
        accept = accept + 1;
        }
        else
        { 
        }
    }     
    
    
    // Return the results
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}



//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
NumericVector gaussiancarmultilevelupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                                          NumericVector Wtripletsum, NumericVector n_individual,
                                          const int nsites, NumericVector phi, double tau2, 
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
        fcprecision = (n_individual[j]/nu2) + (1/priorvar);
        fcsd = pow((1/fcprecision),0.5);
        fcmean = (priormean / priorvar + offset[j]) / fcprecision;
        phinew[j] = rnorm(1, fcmean, fcsd)[0];      
    }
    
    
    return phinew;
}



//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List binomialcarmultilevelupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                                 NumericVector Wtripletsum, List ind_area_list, NumericVector n_individual,
                                 const int nsites, NumericVector phi, double tau2, 
                                 const NumericVector y, const NumericVector failures, const double phi_tune, 
                                 double rho, NumericVector offset)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0, rowstart=0, rowend=0, n_current=0, datapoint=0;
    double acceptance, sumphi;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit, likebittotal;
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
        // Prior part
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
        
        // Likelihood part
        // Determine the set of data points relating to area j
        n_current = n_individual[j];
        NumericVector individuals(n_current);
        individuals = ind_area_list[j];
        
        // Compute the data likelihood
        likebittotal = 0;
        for(int r = 0; r < n_current; r++)
        {
            datapoint = individuals[r] - 1;
            pold = exp(offset[datapoint] + phinew[j]) / (1 + exp(offset[datapoint] + phinew[j]));
            pnew = exp(offset[datapoint] + propphi) / (1 + exp(offset[datapoint] + propphi)); 
            oldlikebit = y[datapoint] * log(pold) + failures[datapoint] * log((1-pold));
            newlikebit = y[datapoint] * log(pnew) + failures[datapoint] * log((1-pnew));
            likebittotal = likebittotal + newlikebit - oldlikebit;   
        }
        
        // Compute the acceptance probability and accept or reject the proposal    
        acceptance = exp(oldpriorbit - newpriorbit + likebittotal);
        if(runif(1)[0] <= acceptance) 
        {
            phinew[j] = propphi;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    
    // Return the results
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}




//////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List poissoncarmultilevelupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                                NumericVector Wtripletsum, List ind_area_list, NumericVector n_individual,
                                const int nsites, NumericVector phi, 
                                double tau2, const NumericVector y, const double phi_tune, 
                                double rho, NumericVector offset)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0,rowstart=0, rowend=0, n_current=0, datapoint=0;
    double acceptance, sumphi;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit, likebittotal;
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
        // Prior part
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
        
        // Likelihood part
        // Determine the set of data points relating to area j
        n_current = n_individual[j];
        NumericVector individuals(n_current);
        individuals = ind_area_list[j];
        
        // Compute the data likelihood
        likebittotal = 0;
        for(int r = 0; r < n_current; r++)
        {
            datapoint = individuals[r] - 1;
            lpold = offset[datapoint] + phinew[j];
            lpnew = offset[datapoint] + propphi; 
            oldlikebit = y[datapoint] * lpold - exp(lpold);
            newlikebit = y[datapoint] * lpnew - exp(lpnew);
            likebittotal = likebittotal + newlikebit - oldlikebit;   
        }
        
        // Compute the acceptance probability and accept or reject the proposal    
        acceptance = exp(oldpriorbit - newpriorbit + likebittotal);
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



