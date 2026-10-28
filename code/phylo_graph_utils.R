#' Read a phylogenetic network in extended Newick format as an igraph DAG
#'
#' @param text Character with a network in extended Newick format,
#' with inheritance probabilities (gamma) after a double colon
#' (e.g., "#H1:0.5::0.3"), as written by SNaQ/PhyloNetworks.
#'
#' @return A directed `igraph` object (a DAG) with edge attributes
#' 'length' (branch length) and 'gamma' (inheritance probability), and
#' vertex attribute 'is_tip'.
#'
read_network_dag <- function(text) {

    if(grepl("#H", text)) {
        net <- ape::read.evonet(text = text)
    } else { # tree (no reticulations)
        net <- ape::read.tree(text = text)
        net$reticulation <- matrix(integer(0), ncol = 2)
    }
    ntip <- length(net$tip.label)
    node_names <- c(net$tip.label, paste0("node", seq_len(net$Nnode)))

    # Tree edges (including major hybrid edges)
    edges <- data.frame(
        from = node_names[net$edge[, 1]],
        to = node_names[net$edge[, 2]],
        length = net$edge.length,
        gamma = 1
    )

    # Get length and gamma of hybrid edges from the Newick string
    hyb <- regmatches(text, gregexpr("\\)?#H[0-9]+:[^,);]*", text))[[1]]
    hyb_df <- data.frame(
        label = sub(":.*", "", sub("^\\)", "", hyb)),
        defining = startsWith(hyb, ")"), # occurrence with a subtree
        length = as.numeric(sapply(strsplit(hyb, ":"), `[`, 2)),
        gamma = as.numeric(sapply(strsplit(hyb, "::"), `[`, 2))
    )

    # Update major hybrid edges (parent edge of each hybrid node)
    hnode_labels <- net$node.label[net$reticulation[, 2] - ntip]
    hnodes <- node_names[net$reticulation[, 2]]
    major_gamma <- vapply(hnode_labels, function(h) {
        hyb_df$gamma[hyb_df$label == h & hyb_df$defining]
    }, numeric(1))
    idx <- match(edges$to, hnodes)
    edges$gamma[!is.na(idx)] <- major_gamma[idx[!is.na(idx)]]

    # Add minor hybrid edges
    minor_edges <- lapply(seq_len(nrow(net$reticulation)), function(i) {
        minor <- hyb_df[hyb_df$label == hnode_labels[i] & !hyb_df$defining, ]
        data.frame(
            from = node_names[net$reticulation[i, 1]],
            to = node_names[net$reticulation[i, 2]],
            length = minor$length,
            gamma = minor$gamma
        )
    })
    edges <- rbind(edges, do.call(rbind, minor_edges))

    # Create DAG
    dag <- igraph::graph_from_data_frame(
        edges, directed = TRUE,
        vertices = data.frame(
            name = node_names, is_tip = seq_along(node_names) <= ntip
        )
    )

    return(dag)
}


#' Compute the Brownian motion variance-covariance matrix on a phylogenetic DAG
#'
#' Implements the recursive algorithm of Bastide et al. (2018, Syst Biol):
#' nodes are visited in topological order (parents before children); for tree
#' nodes, covariances are inherited from the parent; for hybrid nodes,
#' the trait is a weighted average (weights = gamma) of the parents' traits.
#' Works for trees (DAGs without hybrid nodes) and networks.
#'
#' @param dag A directed `igraph` object as returned by `read_network_dag()`.
#'
#' @return A tips x tips variance-covariance matrix.
#'
network_vcv <- function(dag) {

    ord <- igraph::as_ids(igraph::topo_sort(dag, mode = "out"))
    n <- igraph::vcount(dag)
    nodes <- igraph::V(dag)$name
    vcv <- matrix(0, n, n, dimnames = list(nodes, nodes))

    # Update covariances of node `v` given those of previously visited nodes
    update_vcv <- function(vcv, i) {
        v <- ord[i]
        visited <- ord[seq_len(i - 1)]
        in_edges <- igraph::incident(dag, v, mode = "in")
        if(length(in_edges) == 0) return(vcv) # root

        parents <- igraph::tail_of(dag, in_edges)$name
        g <- in_edges$gamma
        l <- in_edges$length

        # Covariance with previously visited nodes
        vcv[v, visited] <- vcv[visited, v] <- colSums(
            g * vcv[parents, visited, drop = FALSE]
        )

        # Variance
        vcv[v, v] <- sum(g^2 * (diag(vcv)[parents] + l)) +
            sum(outer(g, g) * vcv[parents, parents] * (1 - diag(length(g))))

        return(vcv)
    }
    vcv <- Reduce(update_vcv, seq_along(ord), init = vcv)

    tips <- nodes[igraph::V(dag)$is_tip]
    return(vcv[tips, tips])
}


#' Convert a variance-covariance matrix to a weighted species graph
#'
#' @param vcv A species x species variance-covariance matrix.
#'
#' @return An undirected `igraph` object with edge weights representing
#' covariances (i.e., shared evolutionary history), and vertex attribute
#' 'variance' storing the diagonal of `vcv` (i.e., root-to-tip distances).
#'
vcv2graph <- function(vcv) {
    g <- igraph::graph_from_adjacency_matrix(
        vcv, mode = "undirected", weighted = TRUE, diag = FALSE
    )
    igraph::V(g)$variance <- diag(vcv)[igraph::V(g)$name]
    return(g)
}


#' Extract a variance-covariance matrix from a species graph
#'
#' @param graph An undirected `igraph` object as returned by `vcv2graph()`.
#'
#' @return A species x species variance-covariance matrix.
#'
graph2vcv <- function(graph) {
    C <- as.matrix(igraph::as_adjacency_matrix(graph, attr = "weight"))
    diag(C) <- igraph::V(graph)$variance
    return(C)
}


#' Fit a phylogenetic GLS model (intercept only) using a species graph
#'
#' Models a trait as evolving under Brownian motion (BM), with the
#' covariance among species obtained from a `colGraph` of a
#' `GraphExperiment` object. Optionally, Pagel's lambda is estimated
#' by maximum likelihood.
#'
#' @param ge A `GraphExperiment` object with traits in rows and
#' species in columns.
#' @param trait Character, name of the trait (row) to model.
#' @param graph Character, name of the graph in `colGraphs` with the
#' covariance structure (as returned by `vcv2graph()`).
#' @param model Character, either "BM" or "lambda".
#' @param assay Character, name of the assay to use.
#'
#' @return A data frame with the estimated parameters, log-likelihood, and AIC.
#'
graph_gls <- function(ge, trait, graph, model = "BM", assay = 1) {

    vals <- SummarizedExperiment::assay(ge, assay)[trait, ]
    vals <- vals[!is.na(vals)]
    C <- graph2vcv(GraphExperiment::colGraph(ge, graph))[names(vals), names(vals)]
    n <- length(vals)
    dm <- matrix(1, n, 1) # design matrix (intercept only)

    # Profile log-likelihood for a given lambda
    loglik <- function(lambda) {
        Cl <- C * lambda
        diag(Cl) <- diag(C)
        Ci <- solve(Cl)
        beta <- solve(t(dm) %*% Ci %*% dm, t(dm) %*% Ci %*% vals)
        r <- vals - dm %*% beta
        sigma2 <- as.numeric(t(r) %*% Ci %*% r) / n
        ll <- -n / 2 * log(2 * pi * sigma2) -
            as.numeric(determinant(Cl)$modulus) / 2 - n / 2
        list(ll = ll, beta = as.numeric(beta), sigma2 = sigma2)
    }

    lambda <- 1
    if(model == "lambda") {
        lambda <- optimize(function(l) loglik(l)$ll, c(0, 1), maximum = TRUE)$maximum
    }
    fit <- loglik(lambda)
    k <- 2 + (model == "lambda") # number of parameters

    res <- data.frame(
        trait = trait, graph = graph, model = model, n = n,
        root_state = fit$beta, sigma2 = fit$sigma2, lambda = lambda,
        logLik = fit$ll, AIC = -2 * fit$ll + 2 * k
    )
    return(res)
}
