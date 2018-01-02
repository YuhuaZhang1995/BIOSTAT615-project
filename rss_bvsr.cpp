#include <RcppArmadillo.h>
#include <cstdlib>
#include <math.h>
#include <random>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec cal_post_bvsr(
	arma::vec& betahat,  //p by 1
	arma::vec& se, //p by 1
	arma::uvec& rank,  //the index of snps included in the current model, 1 by n_gamma
	double psi,  //the current value of the prior variance of beta
	double logpi) 
{
	//in current model, I want to implement the identity matrix type
	double logpost=0;
	int n_gamma=rank.size();
	int p=betahat.size();
	arma::vec Psiz(n_gamma);
	Psiz.fill(psi);

	arma::vec RSBSMz = betahat.elem(rank);
	arma::vec tmp(n_gamma);
	tmp.ones();
	arma::vec S2iz= tmp/((se.elem(rank))%(se.elem(rank)));
	arma::vec OmegaInv = S2iz + tmp/Psiz;
	arma::vec muz=RSBSMz / OmegaInv;

	arma::vec betapost_err(n_gamma);
	betapost_err.randn();
	arma::rowvec betapost;
	betapost=(muz + betapost_err).t();

	double logdet= -0.5*sum(log(OmegaInv))-0.5*sum(log(Psiz));
	//Rcout<<logdet<<endl;
	double logquad=0.5*dot(RSBSMz,muz);
	//Rcout<<logquad<<endl;
	logpost=logpost+logdet+logquad;

	double loglik=logpost;
	//Rcout<<loglik<<endl;
	logpost=logpost+(p-n_gamma)*log(1-exp(logpi))+(n_gamma-1)*logpi;
	//Rcout<<logpost<<endl;
	betapost.insert_cols(n_gamma,1);
	betapost(n_gamma)=loglik;
	betapost.insert_cols(n_gamma+1,1);
	betapost(n_gamma+1)=loglik;
	return betapost;
}

// [[Rcpp::export]]
double calc_beta_variance(arma::vec& xxyy, double logpi, double h){
	double xxyysum=sum(xxyy);
	double pival = exp(logpi);
	double psi = h/(pival*xxyysum);
	return psi;
}

// [[Rcpp::export]]
double propose_h(double h_old,int rep){
	double h=h_old;
	srand(unsigned(time(NULL)));
	for(int i=0;i<rep;i++){
		double dice_roll=(double)rand()/(RAND_MAX+1);
		//Rcout<<dice_roll<<endl;
		h=h+(dice_roll-0.5)*2;
		while(true){
			if(h<0){
				h=0-h;
			}
			if(h>1){
				h=2-h;
			}
			if(h>=0 and h<=1){break;}
		}
	}
	return h;
}

// [[Rcpp::export]]
arma::vec propose_logpi(double logpi_old, int rep, int p){
	double pi_logratio=0;
	double logpi_new=0;
	arma::vec tmp(2);

	srand(unsigned(time(NULL)));
	for(int i=0;i<rep;i++){
		double dice_roll=(double)rand()/(RAND_MAX+1);
		logpi_new=logpi_old+(dice_roll-0.5)*0.1;
		while(true){
			if(logpi_new<log(1/p)){
				logpi_new=2*log(1/p)-logpi_new;
			}
			if(logpi_new>log(1-1e-8)){
				logpi_new=2*log(1)-logpi_new;
			}
			if(logpi_new >=log(1/p) and logpi_new<=log(1-1e-8)){break;}
		}
		pi_logratio=pi_logratio+(logpi_new-logpi_old);
		logpi_old=logpi_new;
	}
	tmp(0)=pi_logratio;
	tmp(1)=logpi_new;
	return tmp;
}

// [[Rcpp::export]]
arma::vec pmf_ugmix(int p, double geomean){
	double gp=1-exp(-1/geomean);
	arma::vec unif_term(p);
	unif_term.fill((double)0.3/p);
	arma::vec geom_seq(p);
	geom_seq(0)=1;
	for(int i=1;i<p;i++){
		geom_seq(i)=geom_seq(i-1)*(1-gp);
	}
	arma::vec geom_term(p);
	geom_term=0.7*gp/(1-std::pow((1-gp),p))*geom_seq;

	arma::vec p_gamma=unif_term+geom_term;
	p_gamma=p_gamma/sum(p_gamma);
	//Rcout<<p_gamma<<endl;
	return p_gamma;
}

// [[Rcpp::export]]
arma::uvec initiate_model(int p, arma::vec abz, arma::vec& gamma_start){
	gamma_start.zeros();
	//arma::uvec snp_rank;
	double q_genome = R::qnorm(1-(0.025/p),0,1,true,false);
	//Rcout<<q_genome<<endl;
	arma::uvec in_loci=find(abz>q_genome);
	//Rcout<<"in_loci"<<in_loci<<endl;
	int ngamma_start=sum(in_loci);
	arma::uvec snp_rank=sort_index(abz,"descend");
	//Rcout<<"snp_rank"<<snp_rank<<endl;
	int baseline=min(10,p-1);
	arma::uvec tmp(baseline);
	for(int i=0;i<baseline;i++){
		tmp(i)=snp_rank(i);
	}
	if(ngamma_start<baseline){
		ngamma_start=baseline;
		gamma_start.elem(tmp).ones();
		//Rcout<<"gamma"<<gamma_start<<endl;
	}
	
	else{
		gamma_start.elem(in_loci).ones();
		//Rcout<<"gamma"<<gamma_start<<endl;
	}
	snp_rank.insert_rows(p,1);
	snp_rank(p)=ngamma_start;
	return snp_rank;
}

//Remove the unusual conditions
// [[Rcpp::export]]
double propose_gamma(
	arma::uvec& rank_new,
	arma::vec& p_gamma, //the pmf of the mixture of discrete uniform and geometric rv, p by 1
	int rep, 
	arma::vec gamma_old //0 and 1 vector, indicate the SNPs that are currently in the model
	){
	int p=p_gamma.size();
	double gamma_logratio=0;

	int ngamma_old=rank_new.size();
	int ngamma_new=ngamma_old;
	//rank_new=rank_old;
	int r_add;
	//
	srand(unsigned(time(NULL)));
	for(int i=0;i<rep;i++){
		//rare case 1: no snp in the model
		if (ngamma_new==0){
			//Rcout<<"sublabel8"<<endl;
			double dice_roll=R::runif(0,1);
			double sum_tmp=0;
			for(int j=0;j<p;j++){
				sum_tmp+=p_gamma(j);
				if (dice_roll<=sum_tmp){
					r_add=j;
					break;
				}
			}
			gamma_old(r_add)=1;
			rank_new.set_size(1);
			rank_new(0)=r_add;
			ngamma_new +=1;
			gamma_logratio -=log((double)p_gamma(r_add));
			continue;
		}

		if(ngamma_new==p){
			int col_id=rand()%ngamma_new;
			int r_remove= rank_new(col_id);

			gamma_old(r_remove)=0;
			rank_new.shed_row(col_id);

			gamma_logratio 	= gamma_logratio+log((double)ngamma_new);
			ngamma_new -=1;
			//Rcout<<"sublabel1"<<endl;
			continue;
		}
		//Rcout<<"sublabel2"<<endl;

		double gamma_flag=(double)rand()/(RAND_MAX+1);
		//Rcout<<gamma_flag<<endl;
		if ( gamma_flag < 0.4 ){
			//sample a SNP that is not in the current model, return the rank r_add
			while(true){
				double dice_roll=R::runif(0,1);
				double sum_tmp=0;
				for(int j=0;j<p;j++){
					sum_tmp+=p_gamma(j);
					if (dice_roll<=sum_tmp){
						r_add=j;
						break;
					}
				}
				//Rcout<<r_add;
				if (gamma_old(r_add)==0){break;}
				//Rcout<<"sublabel3"<<endl;
			}
			//Rcout<<"sublabel4"<<endl;
			double prob_total=1;
			prob_total=prob_total-sum(p_gamma.elem(rank_new));

			gamma_old(r_add)=1;
			rank_new.insert_rows(ngamma_new,1);
			rank_new(ngamma_new)=r_add;
			ngamma_new += 1;
			gamma_logratio=gamma_logratio-log(p_gamma(r_add)/prob_total)-log(ngamma_new);
		}
			//Rcout<<"label4"<<endl;
		else if (gamma_flag >=0.4 and gamma_flag< 0.8){
			//Rcout<<"label2"<<endl;
			int col_id=rand()%ngamma_new;
			int r_remove=rank_new(col_id);

			double prob_total=1;
			prob_total=prob_total-sum(p_gamma.elem(rank_new));
			prob_total=prob_total+p_gamma(r_remove);

			gamma_old(r_remove)=0;
			rank_new.shed_row(col_id);

			gamma_logratio 	= gamma_logratio + log(p_gamma(r_remove) / prob_total) + log(ngamma_new);
			ngamma_new 	= ngamma_new - 1;
			//Rcout<<"sublabel5"<<endl;
		}
		
		else{
			//Rcout<<"sublabel6"<<endl;
			int col_id=rand()%ngamma_new;
			int r_remove=rank_new(col_id);

			while(true){
				double dice_roll=R::runif(0,1);
				double sum_tmp=0;
				for(int j=0;j<p;j++){
					sum_tmp+=p_gamma(j);
					if (dice_roll<=sum_tmp){
						r_add=j;
						break;
					}
				}
				if (gamma_old(r_add)==0){break;}
			}
			//Rcout<<"label3"<<endl;
			double prob_total=1;
			prob_total=prob_total-sum(p_gamma.elem(rank_new));
			gamma_logratio 	= gamma_logratio + log( p_gamma(r_remove) / (prob_total + p_gamma(r_remove) - p_gamma(r_add)));
			gamma_logratio 	= gamma_logratio - log( p_gamma(r_add) / prob_total );

			gamma_old(r_remove)=0;
			gamma_old(r_add)=1;
			rank_new.shed_row(col_id);
			rank_new.insert_rows(ngamma_new-1,1);
			rank_new(ngamma_new-1)=r_add;
			//Rcout<<"sublabel7"<<endl;
		}
	}
	return gamma_logratio;
}

// [[Rcpp::export]]
List rss_bvsr(
	arma::vec betahat,
	arma::vec se,
	arma::vec Nsnp,
	int Ndraw, int Nburn, int nthin){

	//pre-compute
	int p=betahat.size();

	//declare a tmp vector stores all ones;
	arma::vec tmp(p); 
	tmp.ones();

	arma::vec xxyy = tmp/(Nsnp%se%se);
	arma::vec q = betahat/(se%se);
	//arma::vec Si=tmp/se;
	//arma::mat SiRiS=...;

	//preallocate memory to store the output
	//int Nsam=Ndraw-Nburn;
	//arma::mat betasam(p,Nsam);
	arma::vec out_put(p);
	arma::vec logpi((Ndraw-Nburn)/nthin);
	arma::vec h_store((Ndraw-Nburn)/nthin);
	out_put.zeros();

	//initialize model parameter
	arma::vec abz=abs(betahat/se);
	arma::vec gamma_start(p);
	arma::uvec snp_rank = initiate_model(p, abz, gamma_start);
	double ngamma_start=snp_rank(p);
	snp_rank.shed_row(p);

	arma::vec gamma_old=gamma_start;
	//Rcout<<gamma_old<<endl;
	arma::uvec rank_old= find(gamma_old==1);
	//Rcout<<rank_old<<endl;
	double logpi_old=log((double)ngamma_start/p);
	double h_old=R::runif(0,1);

	//rank-based gamma proposal distribution
	arma::vec p_gamma=pmf_ugmix(p,2000);
	p_gamma=p_gamma.elem(snp_rank);

	//Rcout<<"label1"<<endl;

	//calculate the posterior for te first iteration
	double psi_old=calc_beta_variance(xxyy,logpi_old,h_old);
	arma::rowvec betapost_old=cal_post_bvsr(betahat, se, rank_old, psi_old, logpi_old); 
	double logpost_old=betapost_old(rank_old.size()+1);
	betapost_old.shed_col(rank_old.size()+1);
	betapost_old.shed_col(rank_old.size());
	//Rcout<<"label3"<<endl;

	arma::vec beta_old(p);
	beta_old.zeros();
	beta_old.elem(rank_old)=betapost_old;

	//store the return value from logpi
	arma::vec tmp_pi(2);

	for(int i=0;i<Ndraw;i++){
		int rep=0;
		double dice_roll=R::runif(0,1);
		if(dice_roll<0.33){
			rep=rand()%19+2;
		}
		else{
			rep=1;
		}
		double logMHratio=0;
		//propose h
		double h_new=propose_h(h_old,rep);
		logMHratio=logMHratio;

		//Rcout<<"label10"<<endl;
		//propose gamma
		arma::uvec rank_new;
		rank_new.set_size(size(rank_old));
		rank_new=rank_old;
		//Rcout<<"label11"<<endl;
		double gamma_logratio=propose_gamma(rank_new, p_gamma, rep, gamma_old);
		logMHratio=logMHratio+gamma_logratio;

		//propose log(pi)
		tmp_pi=propose_logpi(logpi_old, rep, p);
		double pi_logratio=tmp_pi(0);
		double logpi_new=tmp_pi(1);
		logMHratio=logMHratio+pi_logratio;
		//Rcout<<"label4"<<endl;

		//calculate posterior
		double psi_new=calc_beta_variance(xxyy,logpi_new,h_new);
		arma::rowvec betapost_new=cal_post_bvsr(betahat, se, rank_new, psi_new, logpi_new); 
		//Rcout<<"label7"<<endl;

		double logpost_new=betapost_new(rank_new.size()+1);
		betapost_new.shed_col(rank_new.size()+1);
		betapost_new.shed_col(rank_new.size());
		
		//Rcout<<"label6"<<endl;
		arma::vec beta_new(p);
		beta_new.zeros();
		beta_new.elem(rank_new)=betapost_new;
		//Rcout<<"label5"<<endl;

		logMHratio=logMHratio+logpost_new-logpost_old;

		double dice_r=R::runif(0,1);
		if(logMHratio>0 or log(dice_r)<logMHratio){
			logpost_old=logpost_new;
			rank_old.resize(rank_new.size());
			rank_old=rank_new;
			beta_old=beta_new;
			h_old=h_new;
			logpi_old=logpi_new;
		}
		//Rcout<<"label9"<<endl;
		gamma_old.zeros();
		gamma_old.elem(rank_old).ones();

		if(i>=Nburn and (i-Nburn)%nthin==0){
			for(int j=0;j<p;j++){
				if(beta_old(j)!=0){
					out_put(j)+=1;
				}
			}
			logpi((i-Nburn)/nthin)=logpi_old;
			h_store((i-Nburn)/nthin)=h_old;
		}
		//Rcout<<betasam<<endl;
	}
	//arma::vec out_put(p);
	//Rcout<<"label8"<<endl;
	//out_put=mean(betasam,1);
	tmp.fill(sum(out_put));
	out_put=out_put/tmp;
	return List::create(Named("output")=out_put,Named("h")=h_store,Named("logpi")=logpi);
}