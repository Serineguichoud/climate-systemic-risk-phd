########################Install Packages#########################################################
install.packages("igraph")
install.packages(c("readxl", "tidyverse"))
install.packages("ggplot2")
install.packages("ggnetwork")
install.packages("NetworkRiskMeasures")
library(readxl)
library(tidyverse)
library(igraph)
library(NetworkRiskMeasures)

########################Graph SETUP#########################################################

md_mat2 <- as.matrix(MatrixAll[-1])

head(md_mat2)

dim(md_mat2)

# rownames and colnames for the matrix
rownames(md_mat2) <- colnames(md_mat2) <- MatrixAll$Flows

#Once we have our network, we can visualize it either using igraph or ggplot2 along with the ggnetwork package. Below we give an example with ggplot2 – it’s useful to remember that we have an assets matrix, so a → b means that node a has an asset with node b:

library(ggplot2)
library(ggnetwork)
library(igraph)

# converting our network to an igraph object
gmd2 <- graph_from_adjacency_matrix(md_mat2, weighted = T)

head(gmd2)

# adding other node attributes to the network
V(gmd2)$buffer <- Raw_data_with_all_damages$buffer
V(gmd2)$weights <- Raw_data_with_all_damages$weights/sum(Raw_data_with_all_damages$weights)
V(gmd2)$exports  <- Raw_data_with_all_damages$assets
V(gmd2)$imports <- Raw_data_with_all_damages$liabilities

# ploting with ggplot and ggnetwork
netdf2 <- ggnetwork(gmd2)

ggplot(netdf2, aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_edges(arrow = arrow(length = unit(6, "pt"), type = "closed"), 
             color = "grey50", curvature = 0.1, alpha = 0.5) + 
  geom_nodes(aes(size = weights)) + 
  ggtitle("Global Value Chains as a graph") + 
  theme_blank()

# network densityh
edge_density(gmd2)

# assortativity
assortativity_degree(gmd2)


########################IMPACT ANALYSIS#########################################################
Raw_data_with_all_damages$degree <- igraph::degree(gmd2)
Raw_data_with_all_damages$btw    <- igraph::betweenness(gmd2)
Raw_data_with_all_damages$close  <- igraph::closeness(gmd2)
Raw_data_with_all_damages$eigen  <- igraph::eigen_centrality(gmd2)$vector
Raw_data_with_all_damages$imps <- impact_susceptibility(exposures = gmd2, buffer = Raw_data_with_all_damages$buffer)
Raw_data_with_all_damages$impd <- impact_diffusion(exposures = gmd2, buffer = Raw_data_with_all_damages$buffer, weights = Raw_data_with_all_damages$weights)$total

head(Raw_data_with_all_damages)

# Value Chain simulation with debtrank method
contdr2 <-  contagion(exposures = gmd2, buffer = Raw_data_with_all_damages$buffer, weights = Raw_data_with_all_damages$weights, 
                      shock = "all", method = "debtrank", verbose = F)

plot(contdr2)

contdr2_summary <- summary(contdr2)
Raw_data_with_all_damages$Debtrank <- contdr2_summary$summary_table$additional_stress

Raw_data_with_all_damages$Debtrank

# Traditional default cascades simulation
contthr <-  contagion(exposures = md_mat2, buffer = Raw_data_with_all_damages$buffer, weights = Raw_data_with_all_damages$weights, 
                      shock = "all", method = "threshold", verbose = F)
summary(contthr)

plot(contthr)

contthr_summary <- summary(contthr)
Raw_data_with_all_damages$cascade <- contthr_summary$summary_table$additional_stress

head(Raw_data_with_all_damages)

# save network particularities as excel
setwd("C:\\Users\\sguichoud\\Desktop\\Article 1")
write_xlsx(Raw_data_with_all_damages, "Results All final.xlsx")

rankings <- Raw_data_with_all_damages[1]
rankings <- cbind(rankings, lapply(Raw_data_with_all_damages[c("DebtRank","cascade","degree","eigen","impd","assets", "liabilities", "buffer")], 
                                   function(x) as.numeric(factor(-1*x))))
rankings <- rankings[order(rankings$DebtRank), ]
head(rankings, 10)

setwd("C:\\Users\\sguichoud\\Desktop\\Article 1")
write_xlsx(rankings, "Results All final as ranking.xlsx")










# Fonction de contagion inspirée de Debtrank mais appelée "climate extremes"
contagion_climate_extremes <- function(exposures, buffer, shock, weights = NULL, 
                                       method = "climate_extremes", verbose = TRUE, 
                                       max.iter = 1000, abs.tol = 1e-8) {
  
  # Vérification des dimensions des entrées
  n <- nrow(exposures)
  if (length(shock) != n) {
    stop("La longueur de 'shock' doit être égale au nombre de noeuds (lignes de 'exposures').")
  }
  
  # Si aucun poids n'est donné, on assigne des poids égaux
  if (is.null(weights)) {
    weights <- rep(1, n)
  }
  weights <- weights / sum(weights)  # Normalisation des poids
  
  # Initialisation de la matrice d'impact (v)
  v <- exposures * buffer  # Matrice d'impact initiale (exposures ajustées par le buffer)
  
  # Initialisation des variables pour la propagation des chocs
  cap <- rep(0, n)  # Variable pour la capacité de propagation des chocs
  s <- shock  # Le choc initial
  w <- weights  # Poids des noeuds
  
  # Boucle de contagion avec propagation des chocs
  for (iter in 1:max.iter) {
    # Calcul de la propagation : pmin entre capacité et choc
    cap <- pmin(cap + s, 1)  # Assure que la capacité est toujours entre 0 et 1
    s <- cap * v %*% w  # Propagation du choc en fonction de la matrice d'impact et des poids
    
    # Vérification de la convergence
    if (max(abs(s - cap)) < abs.tol) {
      if (verbose) cat("Convergence atteinte après", iter, "iterations.\n")
      break
    }
  }
  
  # Vérification que les résultats sont bien des vecteurs numériques
  s <- as.numeric(s)
  cap <- as.numeric(cap)
  
  # Résultat final
  results <- list(final_shock = s, propagation_capacity = cap)
  class(results) <- "contagion"
  
  return(results)
}

# Appel de la fonction avec les données d'exemple
contagion_results <- contagion_climate_extremes(
  exposures = md_mat2, 
  buffer = Raw_data_with_all_damages$buffer, 
  shock = Raw_data_with_all_damages$shock, 
  weights = Raw_data_with_all_damages$weights, 
  verbose = TRUE
)

# Affichage des résultats
print(contagion_results)

# Optionnel : analyse de l'objet contagion
plot(contagion_results)
summary(contagion_results)

# Inspecter la structure des résultats
str(contagion_results)

# Vérifier que les éléments sont numériques
print(is.numeric(contagion_results$final_shock))  # Doit retourner TRUE
print(is.numeric(contagion_results$propagation_capacity))  # Doit retourner TRUE

# Plot des chocs finaux
plot(contagion_results$final_shock, type = "b", col = "blue", 
     main = "Propagation des chocs finaux", xlab = "Noeuds", ylab = "Chocs")

# Plot de la capacité de propagation
plot(contagion_results$propagation_capacity, type = "b", col = "red", 
     main = "Capacité de propagation des chocs", xlab = "Noeuds", ylab = "Capacité")



# Fonction de contagion simplifiée (Climate extremes)
contagion_climate_extremes <- function(exposures, 
                                       buffer, 
                                       shock, 
                                       weights = NULL, 
                                       method = "debtrank", 
                                       verbose = TRUE, 
                                       max.it = 1000, 
                                       abs.tol = 1e-6) {
  
  # Vérification des dimensions
  n_nodes <- nrow(exposures)
  
  if (length(shock) != n_nodes) {
    stop("La longueur de 'shock' doit correspondre au nombre de nœuds.")
  }
  
  # Si les poids ne sont pas fournis, on les initialise à 1 pour chaque nœud
  if (is.null(weights)) {
    weights <- rep(1, n_nodes)
    warning("Aucun poids fourni : utilisation de poids égaux pour tous les nœuds.")
  }
  
  # Normalisation des poids
  weights <- weights / sum(weights)
  
  # Application des chocs
  simulations <- list()
  for (i in seq_along(shock)) {
    if (verbose) cat("\nApplication du choc sur le nœud", i, "avec un choc de", shock[i] * 100, "%\n")
    
    # Initialisation du vecteur de choc spécifique
    shock_vector <- rep(0, n_nodes)
    shock_vector[i] <- shock[i]
    
    # Simulation de la propagation
    simulations[[i]] <- contagion_engine(
      exposures = exposures, 
      buffer = buffer, 
      shock_vector = shock_vector, 
      weights = weights, 
      method = method, 
      max.it = max.it, 
      abs.tol = abs.tol, 
      verbose = verbose
    )
  }
  
  # Résultats
  results <- list(simulations = simulations)
  return(results)
}

# Fonction pour effectuer le calcul de propagation (simplifié ici)
contagion_engine <- function(exposures, 
                             buffer, 
                             shock_vector, 
                             weights, 
                             method, 
                             max.it, 
                             abs.tol, 
                             verbose) {
  # Implémentation simplifiée de l'engin de contagion
  
  # Calcul de l'impact initial en fonction des chocs
  impact <- exposures %*% shock_vector
  
  # Propagation du choc (simplification : propagation itérative)
  for (i in 1:max.it) {
    new_impact <- impact * (1 - buffer)  # Réduction par le buffer à chaque itération
    
    # Condition d'arrêt
    if (max(abs(new_impact - impact)) < abs.tol) {
      break
    }
    
    impact <- new_impact
  }
  
  return(list(impact = impact))
}

# Utilisation de la fonction sur vos données

# Supposons que vos données sont dans les objets suivants
# md_mat2 = Matrice d'expositions
# Raw_data_with_all_damages$buffer = Matrice de buffer
# Raw_data_with_all_damages$shock = Vecteur de chocs

contagion_results <- contagion_climate_extremes(
  exposures = md_mat2, 
  buffer = Raw_data_with_all_damages$buffer, 
  shock = Raw_data_with_all_damages$shock, 
  weights = Raw_data_with_all_damages$weights, 
  method = "debtrank", 
  verbose = TRUE
)

# Résumé et graphique des résultats
print(summary(contagion_results))

# Fonction de plot pour visualiser les résultats
plot_contagion_results <- function(results) {
  # Extraction des impacts simulés
  impacts <- sapply(results$simulations, function(sim) sim$impact)
  
  # Affichage des résultats (par exemple, les impacts pour chaque nœud)
  matplot(impacts, type = "l", lty = 1, col = 1:ncol(impacts), 
          main = "Propagation des chocs (Contagion Climate Extremes)", 
          xlab = "Itérations", ylab = "Impact")
}

# Affichage du graphique des résultats
plot_contagion_results(contagion_results)
summary(contagion_results)

# Supposons que 'md_mat2' soit ta matrice d'exposition, et que 'contagion_results' contienne les résultats de la contagion.

# Exemple d'exposition pour chaque nœud (si ta matrice est md_mat2)
exposition_par_noeud <- rowSums(md_mat2)  # La somme des expositions pour chaque nœud (par ligne)

# Supposons que tu veuilles aussi inclure les résultats de contagion (par exemple 'stress_values')
stress_values <- sapply(contagion_results$simulations, function(x) sum(x$st))  # Ou toute autre extraction selon ta structure de résultats

# Créer un tableau avec l'exposition et les résultats de contagion pour chaque nœud
tableau_exposition <- data.frame(
  Noeud = rownames(md_mat2),  # Noms des nœuds (si disponibles dans ta matrice)
  Exposition = exposition_par_noeud,  # Exposition par nœud
  Stress = stress_values  # Résultats de contagion (par exemple, le stress ou l'impact)
)

# Afficher le tableau
print(tableau_exposition)

# Si tu veux sauvegarder ce tableau en fichier CSV
write.csv(tableau_exposition, "exposition_par_noeud.csv", row.names = FALSE)

