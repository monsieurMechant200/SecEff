
# SecEff – Cross Section Analysis in R

**SecEff** est une bibliothèque R dédiée à l’analyse, la modélisation et la visualisation des sections efficaces (*cross sections*) en physique nucléaire et des particules, avec une extension conceptuelle innovante vers l’analyse de performance en marketing et sciences comportementales.

Le package propose un cadre mathématique rigoureux inspiré de la physique fondamentale, appliqué aussi bien à des phénomènes physiques qu’à des processus décisionnels et d’interaction en sciences sociales.

---

## Présentation générale

En physique, la section efficace mesure la probabilité qu’une interaction se produise entre des particules.
**SecEff** transpose ce formalisme probabiliste vers d’autres domaines, notamment le marketing, afin de modéliser l’efficacité d’une exposition, d’une campagne ou d’un message, tout en conservant la cohérence mathématique du concept original.

Le package est structuré de manière modulaire afin de faciliter son usage aussi bien en recherche académique qu’en analyse appliquée.

---

## Fonctionnalités principales

**Module Physique (MATH_TOOLS)**

* Calculs de sections efficaces différentielles et totales
* Modèle de diffusion de Rutherford
* Résonances de Breit–Wigner (simples et multiples)
* Propagation statistique des incertitudes
* Loi d’atténuation exponentielle
* Conversions densité atomique – section efficace

**Module Marketing (MARKETING_TOOLS)**

* Calcul d’efficacité de campagnes publicitaires
* Modélisation de la fatigue publicitaire
* Analyse de conversion inspirée des sections efficaces
* Optimisation de l’exposition média

**Module Utilitaires (UTILS)**

* Conversions d’unités (barn ↔ m²)
* Normalisation de spectres
* Génération de grilles d’énergie logarithmiques
* Interpolation de données expérimentales

**Module Visualisation (VIZ_TOOLS)**

* Graphiques de sections efficaces (échelles linéaires et logarithmiques)
* Distributions angulaires (polaires et cartésiennes)
* Cartes de chaleur énergie–angle
* Comparaison multi-courbes
* Tableaux de bord complets avec résidus
* Export d’animations image par image

---

## Installation

Installation directe depuis GitHub :

```r
remotes::install_github("monsieurMechant200/SecEff")
```

---

## Prise en main rapide

**Diffusion de Rutherford**

```r
library(SecEff)

theta <- seq(10, 170, by = 1)
Z1 <- 2
Z2 <- 79
E  <- 5

sigma <- sigma_rutherford(theta, Z1, Z2, E)

plot_angular_distribution(
  theta,
  sigma,
  type = "polar",
  main = "Alpha Scattering on Gold"
)
```

**Résonances de Breit–Wigner multiples**

```r
E <- grille_energie(0.1, 10, n = 500, log = TRUE)

params <- data.frame(
  E0 = c(1.5, 3.0, 5.5),
  Gamma = c(0.2, 0.3, 0.4),
  sigma0 = c(5, 8, 3)
)

sigma <- sigma_multiresonance(E, params)

plot_resonance(E, sigma, title = "Multiple Resonances")
```

**Application marketing**

```r
impressions <- 1e6
conversions <- 5000
density <- 0.7

eff <- sigma_marketing(conversions, impressions, density)

exposures <- seq(1, 20, by = 1)
efficiency <- sigma_fatigue(eff, mu = 0.15, exposures)

plot_sigma(
  exposures,
  efficiency,
  xlab = "Number of Exposures",
  ylab = "Relative Efficiency",
  main = "Advertising Fatigue Curve"
)
```

---

## Concepts théoriques

**Section efficace microscopique**

\[\sigma = \frac{N_{\text{reactions}}}{N_{\text{projectiles}} \cdot n_{\text{targets}} \cdot t}\]


Unité : barn
(1 barn = 10⁻²⁸ m²)

**Diffusion de Rutherford**

\[\frac{d\sigma}{d\Omega}=\left(\frac{Z_1 Z_2 e^2}{4E}\right)^2\frac{1}{\sin^4\left(\frac{\theta}{2}\right)}\]


**Formule de Breit–Wigner**

\[\sigma(E)=\sigma_0\frac{\Gamma^2 / 4}{(E - E_0)^2 + \Gamma^2 / 4}\]


---

## Flux de travail typiques

**Analyse expérimentale complète**

* Import des données expérimentales
* Calcul des sections efficaces
* Ajustement théorique
* Visualisation avec résidus

**Optimisation marketing**

* Modélisation de l’exposition
* Estimation de l’efficacité
* Analyse de la fatigue
* Recherche de densité optimale

---

## Citation

Si vous utilisez ce package dans vos travaux, veuillez citer :

```bibtex
@software{secteff2025,
  author = {Aat Ndongo David Meilleur and Ngneumeu Yvan},
  title  = {SectEff: Cross Section Analysis with R},
  year   = {2025}
}
```

---

## Contact

Email : [meilleurd2001@email.com](mailto:meilleurd2001@email.com)

Développé avec rigueur scientifique pour la communauté R.
