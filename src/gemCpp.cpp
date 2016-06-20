
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export]]

List gemCpp (List dList, CharacterVector met){
	
	List datalist(dList);
	int k = datalist.size();  
	std::string method = as<std::string>(met);

	arma::mat temp = as<arma::mat>(datalist[0]);
	arma::mat beta(k,temp.n_cols);
    arma::mat xFrame = temp.cols(1, (temp.n_cols-1));
	arma::vec varY(k);
	arma::vec N(k);
    for(int i=0; i<k; i++) {     
         arma::mat ll = as<arma::mat>(datalist[i]);
		 N(i) = ll.n_rows;
		 arma::mat y = ll.col(0);
		 varY(i) = var(y.col(0));
		 arma::mat X = join_rows(arma::ones(ll.n_rows,1), ll.cols(1, (ll.n_cols-1)));

		 if (i > 0)
		 	{
				xFrame = join_cols(xFrame, ll.cols(1, (ll.n_cols-1)));
		 	}
         arma::mat coef = arma::solve(X, y); 
		 beta.row(i) = coef.t();
    }
	arma::mat Beta = beta.cols(1,(beta.n_cols-1));
	arma::vec pro = N/sum(N);
    arma::vec pro1 = (N-1)/sum(N);
	arma::mat xMean = mean(xFrame);
	arma::mat co = cov(xFrame);
    arma::vec eigval;
    arma::mat eigvec;
	eig_sym(eigval, eigvec, co);
	arma::mat sqrtco = eigvec * diagmat(sqrt(eigval)) * eigvec.t();

	arma::mat Beta_bar = mean(Beta);
	arma::mat B = arma::zeros(Beta.n_cols,Beta.n_cols);
	arma::mat D = arma::zeros(Beta.n_cols,Beta.n_cols);
    arma::mat A = arma::zeros(Beta.n_cols,Beta.n_cols);
	for(int i=0; i<k; i++) {     
         B = B + pro(i) * (Beta.row(i)-Beta_bar).t() * (Beta.row(i)-Beta_bar);
         D = D + pro1(i) * Beta.row(i).t() * Beta.row(i);
		 arma::mat temp;  temp.eye(size(A));
		 A = A + temp * (varY(i) * (N(i)-1)/sum(N));
	}
	A = A - sqrtco * D * sqrtco;
	arma::mat astar;
	if (method == "nu"){
		eig_sym(eigval, eigvec, sqrtco * B * sqrtco);	
		astar = inv_sympd(sqrtco) * eigvec.tail_cols(1);
		if (astar(0) < 0){
			astar = astar * (-1);
		}
	}else if (method == "de"){
		eig_sym(eigval, eigvec, sqrtco * D * sqrtco);	
		astar = inv_sympd(sqrtco) * eigvec.tail_cols(1);
		if (astar(0) < 0){
			astar = astar * (-1);
		}
	}else if (method == "F"){
		eig_sym(eigval, eigvec, inv(A) * sqrtco * B * sqrtco);	
		astar = inv_sympd(sqrtco) * eigvec.tail_cols(1);
		astar = astar.each_row()/(sqrt(astar.t() * co * astar));
		if (astar(0) < 0){
			astar = astar * (-1);
		}
	}

	arma::mat gamma(k,2);
    for(int i=0; i<k; i++) {     
         arma::mat ll = as<arma::mat>(datalist[i]);
		 arma::mat y = ll.col(0);
		 arma::mat X = join_rows(arma::ones(ll.n_rows,1), ll.cols(1, (ll.n_cols-1)) * astar);
         arma::mat coef = arma::solve(X, y); 
		 gamma.row(i) = coef.t();
    }

	return List::create(Named("alpha")  =  astar,
						Named("GammaByTrt")   = gamma,
						Named("betaByTrt")   = beta);
}


