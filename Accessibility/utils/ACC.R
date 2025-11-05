library(Rcpp)
library(RcppParallel)

e2sfca <- function(p, n, D, d0) {
    # Create binary mask where distances <= d0
    mask <- D <= d0
    
    # Apply Gaussian decay to distances (only where mask is TRUE)
    decay <-  mask * exp(- (D^2) / (396^2))
    
    # Step 1: Compute supply-to-demand ratio (Rj)
    step1 <- sweep(decay, 1, p, "*")   # Multiply each row of decay by p
    denom <- colSums(step1)
    Rj <- n / denom
    Rj[!is.finite(Rj)] <- 0           # Avoid Inf or NaN
    
    # Step 2: Compute accessibility scores
    step2 <- sweep(decay, 2, Rj, "*")  # Multiply each column of decay by Rj
    access <- rowSums(step2)
    
    return(access)
}


# cppFunction('
# NumericVector e2sfca_cpp(NumericVector p, NumericVector n,
#                          NumericMatrix D, double d0) {
#   int m = D.nrow();  // number of demand locations
#   int k = D.ncol();  // number of supply locations
# 
#   NumericMatrix D0(m, k);
#   NumericMatrix K(m, k);  // weight matrix
#   double sigma = 396.0;
# 
#   // Step 1: compute Gaussian kernel weights and D0 mask
#   for (int i = 0; i < m; i++) {
#     for (int j = 0; j < k; j++) {
#       if (D(i, j) <= d0) {
#         D0(i, j) = 1.0;
#         double d = D(i, j);
#         K(i, j) = std::exp(-(d * d) / (sigma * sigma));
#       } else {
#         D0(i, j) = 0.0;
#         K(i, j) = 0.0;
#       }
#     }
#   }
# 
#   // Step 2: compute Rj
#   NumericVector Rj(k);
#   for (int j = 0; j < k; j++) {
#     double sum_k = 0.0;
#     for (int i = 0; i < m; i++) {
#       sum_k += K(i, j) * D0(i, j) * p[i];
#     }
#     if (sum_k > 0) {
#       Rj[j] = n[j] / sum_k;
#     } else {
#       Rj[j] = 0.0;
#     }
#   }
# 
#   // Step 3: compute accessibility scores for each demand location
#   NumericVector Ai(m);
#   for (int i = 0; i < m; i++) {
#     double sum_a = 0.0;
#     for (int j = 0; j < k; j++) {
#       if (D0(i, j) == 1.0) {
#         sum_a += Rj[j] * K(i, j);
#       }
#     }
#     Ai[i] = sum_a;
#   }
# 
#   return Ai;
# }
# ')

# # Rcpp functioncall
# accessibility <- e2sfca_cpp(p, n, D, d0)
# print(accessibility)


# even faster
cppFunction('
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// setting up the par worker
struct Step3Worker : public Worker {
  const RMatrix<double> K;
  const RVector<double> Rj;
  RVector<double> Ai;

  Step3Worker(const NumericMatrix& K,
              const NumericVector& Rj,
              NumericVector& Ai)
    : K(K), Rj(Rj), Ai(Ai) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      double acc = 0.0;
      for (std::size_t j = 0; j < Rj.length(); j++) {
        acc += Rj[j] * K(i, j);
      }
      Ai[i] = acc;
    }
  }
};

// [[Rcpp::export]]
NumericVector e2sfca_cpp_parallel(NumericVector p, NumericVector n,
                                        NumericMatrix D, double d0) {
  int m = D.nrow();
  int k = D.ncol();
  double sigma_sq = 396.0 * 396.0;

  NumericMatrix K(m, k);

  // Step 1: compute kernel matrix only
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < k; j++) {
      double d = D(i, j);
      if (d <= d0) {
        K(i, j) = std::exp(- (d * d) / sigma_sq);
      } else {
        K(i, j) = 0.0;
      }
    }
  }

  // Step 2: compute Rj
  NumericVector Rj(k);
  for (int j = 0; j < k; j++) {
    double denom = 0.0;
    for (int i = 0; i < m; i++) {
      denom += K(i, j) * p[i];
    }
    Rj[j] = (denom > 0.0) ? n[j] / denom : 0.0;
  }

  // Step 3: compute Ai (accessibility, in parallel)
  NumericVector Ai(m);
  Step3Worker worker(K, Rj, Ai);
  parallelFor(0, m, worker);

  return Ai;
}
', depends = "RcppParallel")



accdist <- function(d_spatvect,s_spatvect, roads = NA){
    
    library(Rfast)
    
    
    
    road_dist <- function(orig,dest,roads){
        
        library(dplyr)
        library(magrittr)
        library(tidygraph)
        library(sf)
        library(sfnetworks)
        library(cppRouting)
        
        roads_sf <- st_as_sf(roads); orig_sf <- st_as_sf(orig); dest_sf <- st_as_sf(dest)
        
        if (!(st_crs(roads_sf)==st_crs(orig_sf) && st_crs(orig_sf) == st_crs(dest_sf))) stop("Check the CRS of inputs")
        
        
        roads_n <- roads_sf[!st_is_empty(roads_sf), ] |> st_make_valid() |> st_cast("LINESTRING") |>  as_sfnetwork(directed = FALSE)
        
        
        # foreach O&D poinds, calculate shortest distance to the road network
        origin <- st_as_sf(orig_sf) |> st_centroid()  %>%
            mutate(d_r = st_distance(., st_union(roads_sf))) |> select(d_r)
        destination <- st_as_sf(dest_sf) |> st_centroid()  %>%
            mutate(d_r = st_distance(., st_union(roads_sf))) |>  select(d_r)
        
        
        # calculate blend operatnon, takes about 3 minutes
        road_blend <- st_network_blend(roads_n,st_as_sfc(bind_rows(origin,destination))) #|> activate("nodes") |> mutate(node_I = row_number())
        
        
        
        
        
        net <- road_blend|> #as_sfnetwork(roads_n, directed = FALSE) |>
            activate("edges") |>
            mutate(cost = units::drop_units(edge_length()))             
        
        # faster O-D distance calculation using `cppRouting`
        edges <- net |> activate("edges") |> as_tibble() |>
            dplyr::transmute(from = .data$from, to = .data$to, cost = as.numeric(cost)) |> st_drop_geometry()
        
        
        nodes <- net |> activate("nodes") |> st_as_sf() |>
            dplyr::mutate(node_ID = dplyr::row_number()) |>
            dplyr::transmute(node_ID,
                             X = st_coordinates(geometry)[,1],
                             Y = st_coordinates(geometry)[,2])
        
        
        
        # 1. distance calculation on the graph
        RcppParallel::setThreadOptions(numThreads ="auto")
        g  <- makegraph(as.data.frame(edges), directed = FALSE, coords = as.data.frame(st_drop_geometry(nodes)))
        gc <- cpp_contract(g)
        D  <- get_distance_matrix(gc,
                                  from = st_nearest_feature(x=origin,y=nodes), to = st_nearest_feature(x=destination,y=nodes), 
                                  algorithm = "phast")
        D[is.na(D)] <- Inf
        
        
        
        # checking the order of cols, and rows of the D matrix
        #(st_nearest_feature(x=supply,y=nodes) == rownames(D)) |> table()
        #(st_nearest_feature(x=demand,y=nodes) == colnames(D)) |> table()
        
        # > dim(D)
        # [1] 34343 10709
        
        # should vectorize
        # for(j in 1:dim(D)[2]){
        #     
        #     for (i in 1:dim(D)[1]) {
        #     
        #         test = D[i,j]
        #             
        #     }
        #     
        # }
        
        # the vectorized version
        D <- sweep(D, 1, as.numeric(origin$d_r), "+")
        D <- sweep(D, 2, as.numeric(destination$d_r), "+")  
        
        coord_dest <- st_coordinates(destination) |>  as.matrix()
        coord_orig <- st_coordinates(origin) |>  as.matrix()
        eD <- Rfast::dista(coord_orig, coord_dest, type = "euclidean",parallel = T)
        
        
        
        if (!all(dim(eD) == dim(D))) stop("debug")
    
        
        A <- outer(as.numeric(origin$d_r), as.numeric(destination$d_r), "+")
        if (!all(dim(A) == dim(D))) stop("dim mismatch")
        
        D <- D + ((eD - D) * as.numeric(eD < A))
        
        return(D)   
        
    }
    
    
    
    
    if (is.na(roads)[1]) {
        
        
        cat("calculating distance based on the euclidean\n")
        coord_dem <- centroids(d_spatvect) |> crds() |> as.matrix()
        coord_sup <- centroids(s_spatvect) |> crds() |> as.matrix()
        distmat <- Rfast::dista(coord_dem, coord_sup, type = "euclidean",parallel = T)
        dim(distmat)
        
    } else {
        
        cat("calculating distance based on the road provided\n")
        distmat <- road_dist(orig = d_spatvect, dest = s_spatvect, roads = roads)
        dim(distmat)
    }
    
    
    return(distmat)
    
}



acc <- function(d_spatvect,d_attr,s_spatvect,s_attr, distmat,
                d0 = 1609, acolname = "tsfca"){
    
    library(terra)

    
    if (length(d_spatvect[[d_attr]][[1]])!=dim(distmat)[1] || length(s_spatvect[[s_attr]][[1]])!=dim(distmat)[2]) stop("dim mismatch")
    
    
    accres <- e2sfca_cpp_parallel(p = replace(d_spatvect[[d_attr]][[1]], is.na(d_spatvect[[d_attr]][[1]]), 0),
                                  n = replace(s_spatvect[[s_attr]][[1]], is.na(s_spatvect[[s_attr]][[1]]), 0),
                                  D = distmat,
                                  d0 = d0)
    
    ret <- d_spatvect
    
    ret[[acolname]] <- accres
    
    return(ret)
    
}
