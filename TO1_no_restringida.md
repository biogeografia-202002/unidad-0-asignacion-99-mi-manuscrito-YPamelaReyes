Técnicas de ordenación. <br> Parte 1: Ordenación no restringida. <br> PCA, CA y PCoA
================
JR
21 de noviembre, 2020

``` r
knitr::opts_chunk$set(fig.width=12, fig.height=8)
```

Preámbulo
---------

### Cargar paquetes

``` r
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-6

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.2.1     ✔ purrr   0.3.3
    ## ✔ tibble  2.1.3     ✔ dplyr   0.8.3
    ## ✔ tidyr   1.0.0     ✔ stringr 1.4.0
    ## ✔ readr   1.3.1     ✔ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(sf)
```

    ## Linking to GEOS 3.8.0, GDAL 3.0.4, PROJ 7.0.0

``` r
library(mapview)
source('biodata/funciones.R')
```

### Cargar datos

``` r
load('biodata/Euphorbiaceae.Rdata')
load('biodata/matriz_ambiental.Rdata')
mi_fam <- mc_ephrb
(colnames(mi_fam) <- make.cepnames(colnames(mi_fam)))
```

    ##  [1] "Acaldive" "Acalmacr" "Adeltril" "Alchcost" "Alchlati" "Crotbill"
    ##  [7] "Hieralch" "Huracrep" "Sapibroa" "Sapiglan"

``` r
(df_equivalencias <- data.frame(
  nombre_original = colnames(mc_ephrb),
  colnames(mi_fam)))
```

    ##            nombre_original colnames.mi_fam.
    ## 1    Acalypha diversifolia         Acaldive
    ## 2    Acalypha macrostachya         Acalmacr
    ## 3           Adelia triloba         Adeltril
    ## 4  Alchornea costaricensis         Alchcost
    ## 5      Alchornea latifolia         Alchlati
    ## 6     Croton billbergianus         Crotbill
    ## 7  Hieronyma alchorneoides         Hieralch
    ## 8           Hura crepitans         Huracrep
    ## 9         Sapium broadleaf         Sapibroa
    ## 10      Sapium glandulosum         Sapiglan

``` r
bci_env_grid %>% tibble
```

    ## # A tibble: 50 x 1
    ##     .$id $categoria_de_e… $geologia $habitat $quebrada $heterogeneidad…
    ##    <dbl> <fct>            <fct>     <fct>    <fct>                <dbl>
    ##  1     1 c3               Tb        OldSlope Yes                 0.627 
    ##  2     2 c3               Tb        OldLow   Yes                 0.394 
    ##  3     3 c3               Tb        OldLow   No                  0     
    ##  4     4 c3               Tb        OldLow   No                  0     
    ##  5     5 c3               Tb        OldSlope No                  0.461 
    ##  6     6 c3               Tb        OldLow   No                  0.0768
    ##  7     7 c3               Tb        OldLow   Yes                 0.381 
    ##  8     8 c3               Tb        OldLow   Yes                 0.211 
    ##  9     9 c3               Tb        OldLow   No                  0     
    ## 10    10 c3               Tb        OldLow   No                  0     
    ## # … with 40 more rows, and 33 more variables: $UTM.EW <dbl>,
    ## #   $UTM.NS <dbl>, $geomorf_llanura_pct <dbl>, $geomorf_pico_pct <dbl>,
    ## #   $geomorf_interfluvio_pct <dbl>, $geomorf_hombrera_pct <dbl>,
    ## #   $`geomorf_espolón/gajo_pct` <dbl>, $geomorf_vertiente_pct <dbl>,
    ## #   $geomorf_vaguada_pct <dbl>, $geomorf_piedemonte_pct <dbl>,
    ## #   $geomorf_valle_pct <dbl>, $geomorf_sima_pct <dbl>, $Al <dbl>,
    ## #   $B <dbl>, $Ca <dbl>, $Cu <dbl>, $Fe <dbl>, $K <dbl>, $Mg <dbl>,
    ## #   $Mn <dbl>, $P <dbl>, $Zn <dbl>, $N <dbl>, $N.min. <dbl>, $pH <dbl>,
    ## #   $elevacion_media <dbl>, $pendiente_media <dbl>,
    ## #   $orientacion_media <dbl>, $curvatura_perfil_media <dbl>,
    ## #   $curvatura_tangencial_media <dbl>, $geometry <POLYGON [m]>,
    ## #   $abundancia_global <dbl>, $riqueza_global <int>

``` r
# grupos_upgma_k2 <- readRDS('grupos_upgma_k2.RDS')
# table(grupos_upgma_k2)
# grupos_ward_k3 <- readRDS('grupos_ward_k3.RDS')
# table(grupos_ward_k3)
```

Ordenación
----------

La ordenación se basa en los mismos principios que la medición de asociación (similaridad) y el agrupamiento: un objeto se caracteriza por sus propiedades en un espacio n-dimensional, donde cada dimensión es una variable, un descriptor. Un simple diagrama de dispersion nos serviría para representar el caso especial de objetos descritos por sólo dos variables, pero no es lo común. Sin embargo, no podremos encontrar patrones consistentes analizando nuestros datos por pares de variables (e.g. paneles de correlación).

A diferencia del análisis de agrupamiento, o como complemento de éste, el análisis de ordenación abarca un conjunto de técnicas que procuran reducir la dimensionalidad de los datos. Se intenta representar, en ejes ortogonales (comúnmente dos), la complejidad de todo un conjunto. Todas las técnicas de ordenación representan las principales tendencias de variación de los datos en un espacio de dimensiones reducidas, ordenando de manera convencional los ejes con grado descendente de varianza explicada en cada eje sucesivo (e.g. en n dimensiones, el eje 1 explica la mayor varianza, el eje n explica la mínima).

El análisis de ordenación puede ser no restringido (o simple) y restringido (o 'canónico'). En el primer caso, las tendencias detectadas en el conjunto de datos de interés, no están restringidas por otro conjunto, por ejemplo, si buscamos tendencias en una matriz de comunidad, y sólo en ella. En el segundo caso, las tendencias detectadas en un conjunto de datos se asocian a otro conjunto, por ejemplo, si buscamos tendencias en una matriz de comunidad pero restringiéndolas a una matriz ambiental. En este script me concentraré en la ordenación no restringida o simple.

Las principales técnicas de ordenación no restringida son análisis de componentes principales o PCA (siglas de *principal components analysis*), análisis de correspondencia o CA (*correspondence analysis*), análisis de correspondencia múltiple o MCA (*multiple correspondence analysis*), análisis de coordenadas principales o PCoA (*principal coordinates analysis*) y escalamiento multidimensional no métrico o NMDS (*non-metric multidimensional scaling*). Salvo el NMDS, todas estas técnicas se basan en la extracción de los vectores propios (*eigenvectors*) de una matriz de asociación. Explicaré PCA, CA y PCoA a continuación.

### Análisis de componentes principales (PCA)

Es el método tradicional basado en vectores propios que comúnmente se aplica a datos cuantitativos no procesados, que preserva la distancia euclídea; también se aplica a datos de especies, previa transformación de los datos. Por esta razón, es más común aplicarlo a variables ambientales (matriz ambiental), pero también se aplica a datos transformados de composición (matriz de comunidad). Un requisito fundamental para garantizar la eficiencia de este método, es que las variables deben tener algún grado de correlación entre sí, un supuesto a veces imposible de lograr con datos no procesados de matrices de comunidad, pero que sí es bastante común en matrices ambientales. Primero explicaré su uso con un subconjunto de variables ambientales (suelo), para luego aplicarlo a una matriz de comunidad.

#### PCA aplicado a datos ambientales

Para aplicar PCA a datos ambientales, es necesario que todas las variables sean numéricas y "comparables" en cuanto a escalas de medición. Esto se consigue "escalándolas" (convirtiéndolas en puntuaciones z). A partir de la matriz escalada, se generará una matriz de correlaciones.

Dado que se requiere que las variables de entrada sean exclusivamente numéricas, el primer paso que realizaré será obtener un conjunto de columnas numéricas y, de éstas, seleccionaré sólo las de suelo.

¡IMPORTANTE! Haré esta demostración sólo con las variables de suelo, **pero puedes (y debes) ordenar los sitios en función de otras variables, por ejemplo, las geomorfológicas combinadas con la de hábitat y la de heterogeneidad ambiental**. Por tal razón, haz también el PCA para las demás variables de la matriz ambiental.

A partir de los datos de suelo, la función `rda`, de `vegan` realizará los siguientes pasos: escalar las variables originales, calcular matriz de correlaciones y obtener vectores propios para el PCA.

``` r
env_suelo <- bci_env_grid %>%
  st_drop_geometry %>%
  dplyr::select(matches('^[A-T,Z]|^pH$', ignore.case = F))
env_suelo %>% tibble
```

    ## # A tibble: 50 x 1
    ##     .$Al    $B   $Ca   $Cu   $Fe    $K   $Mg   $Mn    $P   $Zn    $N
    ##    <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ##  1  901. 0.794 1680.  6.20  135. 142.   279.  267.  1.95  2.97  18.5
    ##  2  954. 0.670 1503.  6.03  142. 137.   280.  320.  2.25  2.53  21.6
    ##  3 1114. 0.595 1182.  6.80  157.  98.7  230.  445.  1.95  2.25  20.2
    ##  4 1024. 0.568 1558.  6.63  153.  98.4  229.  408.  2.63  2.44  20.8
    ##  5 1002. 0.399 1242.  6.44  149.  94.1  203.  251.  1.86  2.14  16.9
    ##  6 1091. 0.731 1442.  6.50  174. 132.   277.  477.  1.62  2.63  20.3
    ##  7 1184. 0.340 1111.  5.55  138. 117.   242.  301.  2.13  2.16  20.1
    ##  8 1256. 0.322 1029.  6.28  147. 104.   185.  204.  3.11  2.07  21.5
    ##  9 1122. 0.464 1230.  7.18  153. 110.   207.  415.  1.99  2.33  21.4
    ## 10 1172. 0.314 1127.  6.89  133. 105.   172.  330.  1.69  2.05  18.3
    ## # … with 40 more rows, and 2 more variables: $N.min. <dbl>, $pH <dbl>

``` r
env_suelo_pca <- rda(env_suelo, scale = TRUE)
env_suelo_pca
```

    ## Call: rda(X = env_suelo, scale = TRUE)
    ## 
    ##               Inertia Rank
    ## Total              13     
    ## Unconstrained      13   13
    ## Inertia is correlations 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##   PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8   PC9  PC10  PC11  PC12 
    ## 7.871 1.698 1.378 0.899 0.334 0.298 0.149 0.114 0.101 0.075 0.038 0.029 
    ##  PC13 
    ## 0.017

``` r
summary(env_suelo_pca)
```

    ## 
    ## Call:
    ## rda(X = env_suelo, scale = TRUE) 
    ## 
    ## Partitioning of correlations:
    ##               Inertia Proportion
    ## Total              13          1
    ## Unconstrained      13          1
    ## 
    ## Eigenvalues, and their contribution to the correlations 
    ## 
    ## Importance of components:
    ##                          PC1    PC2    PC3     PC4     PC5     PC6    PC7
    ## Eigenvalue            7.8705 1.6985 1.3777 0.89872 0.33400 0.29828 0.1495
    ## Proportion Explained  0.6054 0.1307 0.1060 0.06913 0.02569 0.02294 0.0115
    ## Cumulative Proportion 0.6054 0.7361 0.8421 0.91119 0.93688 0.95982 0.9713
    ##                            PC8      PC9     PC10     PC11     PC12
    ## Eigenvalue            0.113523 0.100540 0.075125 0.038299 0.028503
    ## Proportion Explained  0.008733 0.007734 0.005779 0.002946 0.002193
    ## Cumulative Proportion 0.980052 0.987786 0.993565 0.996511 0.998703
    ##                           PC13
    ## Eigenvalue            0.016858
    ## Proportion Explained  0.001297
    ## Cumulative Proportion 1.000000
    ## 
    ## Scaling 2 for species and site scores
    ## * Species are scaled proportional to eigenvalues
    ## * Sites are unscaled: weighted dispersion equal on all dimensions
    ## * General scaling constant of scores:  5.023829 
    ## 
    ## 
    ## Species scores
    ## 
    ##            PC1       PC2      PC3      PC4      PC5      PC6
    ## Al      0.7651 -0.303596  1.07416 -0.04457  0.07698 -0.04464
    ## B      -1.2880  0.285829  0.12078 -0.10042 -0.25949 -0.08173
    ## Ca     -1.3475 -0.154936 -0.05566  0.09869 -0.21974  0.10661
    ## Cu     -1.1654 -0.297329  0.28359 -0.50726  0.13822  0.18506
    ## Fe     -1.1094 -0.492148  0.40677  0.19816  0.32966  0.18518
    ## K      -1.3645  0.008666  0.04784  0.03985  0.01269  0.07345
    ## Mg     -1.2703 -0.170100 -0.07264  0.25627 -0.35788  0.20705
    ## Mn     -0.9548 -0.317673  0.44289 -0.74565 -0.12314 -0.31689
    ## P       0.2877  0.748145  0.98749  0.47477 -0.23540 -0.03654
    ## Zn     -1.3040  0.004603 -0.02793  0.32706  0.01696 -0.03240
    ## N      -0.5124  1.155952  0.14154 -0.29652  0.26625  0.33702
    ## N.min. -1.1324 -0.163100  0.03454  0.59856  0.34289 -0.30849
    ## pH     -0.9569  0.846451 -0.23774 -0.15686  0.10741 -0.36153
    ## 
    ## 
    ## Site scores (weighted sums of species scores)
    ## 
    ##         PC1      PC2       PC3      PC4       PC5       PC6
    ## 1   0.50347 -0.39861 -0.798154 -0.32512 -1.705624  1.097611
    ## 2   0.47396 -0.28681 -0.518234 -0.23993 -0.940674  0.556514
    ## 3   0.62571 -0.61177 -0.017990 -1.17056 -0.689271  0.325703
    ## 4   0.48722 -0.36275 -0.122325 -0.53899 -0.492944 -0.192532
    ## 5   0.70142 -0.57473 -0.607379 -0.09754  0.095730  0.096622
    ## 6   0.32387 -0.57941 -0.127142 -0.96126 -0.468442 -0.308392
    ## 7   0.70592 -0.29567 -0.231908 -0.18306 -0.001098 -0.372107
    ## 8   0.89397 -0.16407  0.193673 -0.10215  0.098856  0.724943
    ## 9   0.57971 -0.40666 -0.071601 -1.10761 -0.086292 -0.006721
    ## 10  0.72551 -0.49205 -0.281459 -0.67883  0.394079 -0.616348
    ## 11 -0.21317 -0.50177  0.077290 -0.97703 -1.759653 -1.347073
    ## 12  0.72031 -0.71165  0.452589 -0.03515 -0.284845 -0.876996
    ## 13  1.02666 -0.35257  0.902772  1.42743  0.056272  0.769304
    ## 14  0.71660 -0.80220  0.308626 -0.07447  0.179499  0.319996
    ## 15  0.38786 -1.25976  0.502979 -1.12783  0.489936 -0.004906
    ## 16  0.06081 -0.65839 -0.217405 -0.15093 -0.009161 -0.418417
    ## 17  0.63013 -0.10688  0.461756  0.69800 -0.855514 -0.800082
    ## 18  0.96721 -0.48010  0.768611  1.68547  0.280984  0.804510
    ## 19  0.65612 -0.39287  0.642368  0.81392  0.334820 -0.130497
    ## 20  0.46141 -0.62193 -0.007429 -0.28657  0.574129 -0.059174
    ## 21 -0.33476  0.07816 -0.494215 -0.82582  0.082001  0.466656
    ## 22  0.37375  0.39755  0.354272  0.65129  0.466086 -0.169757
    ## 23  0.68109 -0.03137  0.521714  1.37446  0.819763  1.100024
    ## 24  0.12448  0.29626  0.059695 -0.01573  0.236925  0.723204
    ## 25  0.08692  0.66652 -0.374643 -0.06148 -0.269650  1.201049
    ## 26 -0.06428  0.72854 -1.263340 -0.34008  0.931914 -0.269882
    ## 27  0.35229  0.90701 -0.793591  0.27022  0.822847 -0.974726
    ## 28  0.69848  1.15175 -1.855664  0.48677  0.412453 -0.297122
    ## 29  0.19286  1.08572  0.087991  0.11127  0.042980  0.495350
    ## 30  0.34849  0.97558 -0.862374  0.01778  0.495742  0.840577
    ## 31 -0.27158  0.72562 -0.688703 -0.30711  0.513849 -0.289476
    ## 32 -0.23882  1.06756  0.284612  0.39185 -0.583222 -1.160850
    ## 33  0.45683  1.33454 -1.222664  0.21290  0.043306 -0.906784
    ## 34 -0.06262  1.32157  0.707260 -0.38230 -0.554752  0.096821
    ## 35 -0.12287  0.56386  0.842747 -0.37166  0.215049 -0.170818
    ## 36 -0.77718  0.81121 -0.017228 -0.75987  0.240578  0.996530
    ## 37 -0.36979  0.78099  0.600582 -0.68959  0.433235  0.023059
    ## 38 -0.46379  1.17932  1.466360 -0.34391 -0.858632  0.095262
    ## 39 -0.52517  0.80160  1.526464  0.17949 -1.074136 -1.051545
    ## 40 -0.21649  0.62426  1.199057  0.71392 -0.023249 -0.665797
    ## 41 -1.67247  0.07797  0.320315 -0.42982 -0.542412  1.926398
    ## 42 -0.79902 -0.40933  0.170367 -0.68347  0.962695  0.340097
    ## 43 -0.86973 -0.46973  0.977941 -0.73908  1.841802  0.388190
    ## 44 -0.87106 -0.51244  0.510545  0.11302  1.520391 -0.888169
    ## 45 -1.21483 -0.12869  0.318452  0.59274 -0.206188 -0.656346
    ## 46 -1.15473 -0.77489 -0.740832  1.60533 -1.220076  0.792066
    ## 47 -0.53188 -1.10833 -0.737304  1.09916 -0.007989 -1.090101
    ## 48 -1.25160 -1.07535 -0.391085  0.26679  0.818315 -0.567806
    ## 49 -1.60947 -0.72664 -0.790057  0.20147 -0.148172  0.030082
    ## 50 -1.32772 -0.27816 -1.026309  1.09364 -0.622240  0.081856

Para agilizar la producción de scripts analíticos de referencia, trasladaré las explicaciones de cada resultado a los vídeos regulares que alojo en el repositorio de la asignatura. En ellos explicaré cómo interpretar éste y otros resultados.

En el vídeo asociado, explico el significado de:

-   Inercia, *Inertia*
-   Valores propios, autovalores, *Eigenvalues*
-   Escalamiento, *Scaling*
-   Puntuaciones de "especies", *Species scores*
-   Puntuaciones de "sitios", *Site scores*

``` r
screeplot(env_suelo_pca, bstick = TRUE)
```

![](TO1_no_restringida_files/figure-markdown_github/unnamed-chunk-5-1.png)

Usando función `cleanplot.pca`

``` r
par(mfrow = c(1, 2))
cleanplot.pca(env_suelo_pca, scaling = 1, mar.percent = 0.08, cex.char1 = 1.5)
cleanplot.pca(env_suelo_pca, scaling = 2, mar.percent = 0.04, cex.char1 = 1.5)
```

![](TO1_no_restringida_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
par(mfrow = c(1, 1))
```

Comparar distribución de los sitios en biplots con distribución real en el mapa:

### Generar mapa de cuadros sin simbología

``` r
mapa_cuadros <- mapView(
  bci_env_grid,
  col.regions = 'grey80',
  alpha.regions = 0.3,
  map.types = 'OpenTopoMap',
  legend = F, zoom = 14,
  zcol = 'id') %>% addStaticLabels() %>%
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_cuadros
```

![](TO1_no_restringida_files/figure-markdown_github/unnamed-chunk-7-1.png)

Comparar con resultados de un análisis de agrupamiento del mismo conjunto de datos. Primero agrupo mis sitios basado en la misma matriz ambiental fuente del PCA (`env_suelo`), escalándola.

``` r
(env_agrupamiento <- hclust(dist(scale(env_suelo)), 'ward.D'))
```

    ## 
    ## Call:
    ## hclust(d = dist(scale(env_suelo)), method = "ward.D")
    ## 
    ## Cluster method   : ward.D 
    ## Distance         : euclidean 
    ## Number of objects: 50

``` r
(env_grupos <- cutree(env_agrupamiento, k = 3))
```

    ##  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 
    ##  1  1  1  1  1  1  1  1  1  1  2  1  1  1  1  2  1  1  1  1  2  1  1  2  2 
    ## 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 
    ##  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  3  3  3  3  3  3  3  3  3  3

``` r
(mi_cluster <- factor(env_grupos))
```

    ##  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 
    ##  1  1  1  1  1  1  1  1  1  1  2  1  1  1  1  2  1  1  1  1  2  1  1  2  2 
    ## 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 
    ##  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  3  3  3  3  3  3  3  3  3  3 
    ## Levels: 1 2 3

``` r
(mi_cluster_l <- levels(mi_cluster))
```

    ## [1] "1" "2" "3"

``` r
(mi_cluster_l_seq <- 1:length(mi_cluster_l))
```

    ## [1] 1 2 3

Observa que estoy generando un agrupamiento basado en los datos de suelo. No estoy comparando un agrupamiento externo o anterior (e.g. como los creados en los scripts "aa\_analisis\_de\_agrupamiento\*"). Sin embargo, dicha comparación es deseable y posible.

Luego calculo las puntuaciones de los sitios para usarlas luego como coordenadas de los puntos que añadires al gráfico:

``` r
(puntuaciones <- scores(env_suelo_pca, display = 'wa', scaling = 1))
```

    ##            PC1         PC2
    ## 1   0.39174510 -0.14407979
    ## 2   0.36878719 -0.10367009
    ## 3   0.48685902 -0.22112748
    ## 4   0.37909774 -0.13111921
    ## 5   0.54576828 -0.20773969
    ## 6   0.25199854 -0.20942963
    ## 7   0.54926902 -0.10687106
    ## 8   0.69558882 -0.05930492
    ## 9   0.45106847 -0.14698808
    ## 10  0.56450991 -0.17785248
    ## 11 -0.16586882 -0.18136857
    ## 12  0.56046370 -0.25723110
    ## 13  0.79883126 -0.12743768
    ## 14  0.55758047 -0.28995916
    ## 15  0.30178903 -0.45534568
    ## 16  0.04731375 -0.23797984
    ## 17  0.49029611 -0.03863153
    ## 18  0.75257803 -0.17353606
    ## 19  0.51051696 -0.14200354
    ## 20  0.35902048 -0.22479886
    ## 21 -0.26047083  0.02825202
    ## 22  0.29081497  0.14369533
    ## 23  0.52994887 -0.01133905
    ## 24  0.09685586  0.10708545
    ## 25  0.06763353  0.24091762
    ## 26 -0.05001680  0.26333581
    ## 27  0.27411641  0.32784310
    ## 28  0.54347823  0.41630466
    ## 29  0.15006273  0.39243810
    ## 30  0.27115559  0.35262894
    ## 31 -0.21131378  0.26227860
    ## 32 -0.18582397  0.38587422
    ## 33  0.35545742  0.48237789
    ## 34 -0.04872498  0.47768948
    ## 35 -0.09560728  0.20380867
    ## 36 -0.60471858  0.29321598
    ## 37 -0.28773331  0.28229391
    ## 38 -0.36087173  0.42627111
    ## 39 -0.40862742  0.28974135
    ## 40 -0.16844745  0.22564349
    ## 41 -1.30133161  0.02818187
    ## 42 -0.62171296 -0.14795567
    ## 43 -0.67672756 -0.16978708
    ## 44 -0.67776192 -0.18522268
    ## 45 -0.94524751 -0.04651670
    ## 46 -0.89848343 -0.28008841
    ## 47 -0.41385309 -0.40061278
    ## 48 -0.97385861 -0.38869065
    ## 49 -1.25231536 -0.26264706
    ## 50 -1.03308844 -0.10054309
    ## attr(,"const")
    ## [1] 5.023829

Luego creo el gráfico base, coloco los puntos sobre el gráfico usando las puntuaciones, les coloco rótulos y, finalmente, coloco leyenda:

``` r
grafico_base <- plot(
  env_suelo_pca,
  display = "wa",
  scaling = 1,
  type = "n",
  main = "PCA y grupos"
)
abline(v = 0, lty = "dotted")
abline(h = 0, lty = "dotted")
for (i in mi_cluster_l_seq) {
  points(puntuaciones[mi_cluster == i, ],
         pch = (14 + i),
         cex = 2,
         col = i + 1)
}
text(puntuaciones, row.names(env_suelo), cex = 1, pos = 3)
legend(
  "topright", # Otras alternativas: "bottomleft", "bottomright" y "topleft"
  paste("Grupo", c(mi_cluster_l_seq)),
  pch = 14 + c(mi_cluster_l_seq),
  col = 1 + c(mi_cluster_l_seq),
  pt.cex = 2
)
```

![](TO1_no_restringida_files/figure-markdown_github/unnamed-chunk-10-1.png)

Es razonable que el análisis cluster y el biplot muestren patrones consistentes, puesto que se basan en la misma matriz ambiental.

Si hago lo mismo, pero usando mi análisis de agrupamiento anterior (*scripts* "aa\_analisis\_de\_agrupamiento\_\*"), no obtengo resultados consistentes, al menos en mi caso.

``` r
# (mi_cluster_anterior <- grupos_upgma_k2)
# (mi_cluster_anterior <- grupos_ward_k3)
# (mi_cluster_anterior_l <- levels(mi_cluster_anterior))
# (mi_cluster_anterior_l_seq <- 1:length(mi_cluster_anterior_l))
# grafico_base <- plot(
#   env_suelo_pca,
#   display = "wa",
#   scaling = 1,
#   type = "n",
#   main = "PCA y grupos"
# )
# abline(v = 0, lty = "dotted")
# abline(h = 0, lty = "dotted")
# for (i in mi_cluster_anterior_l_seq) {
#   points(puntuaciones[mi_cluster_anterior == i, ],
#          pch = (14 + i),
#          cex = 2,
#          col = i + 1)
# }
# text(puntuaciones, row.names(env_suelo), cex = 1, pos = 3)
# legend(
#   "topright", # Otras alternativas: "bottomleft", "bottomright" y "topleft"
#   paste("Grupo", c(mi_cluster_anterior_l_seq)),
#   pch = 14 + c(mi_cluster_anterior_l_seq),
#   col = 1 + c(mi_cluster_anterior_l_seq),
#   pt.cex = 2
# )
```

Esto podría significar que las tendencias/patrones de mi matriz de comunidad (cuadros de 1 Ha de BCI según composición), no se asocian/no son consistentes con variables de suelo según el PCA. Es probable que, usando una combinación diferente de variables ambientales, se puedan extraer patrones. No recomiendo identificar variables ambientales de forma meramente heurística, porque sería equivalente a pescar; recomiendo construir una matriz ambiental de variables seleccionadas a partir de patrones de dependencia identificados en scripts anteriores. Concretamente, en el script [aa\_analisis\_de\_agrupamiento\_3\_variables\_ambientales\_segun\_grupos.R](aa_analisis_de_agrupamiento_3_variables_ambientales_segun_grupos.R) identifiqué posibles variables asociadas según los distintos agrupamientos realizados. Si fuese tu caso, te recomiendo construir una matriz ambiental con dichas variables y probar la consistencia de los métodos de ordenación y agrupamiento.

#### PCA aplicado a datos de comunidad transformados

``` r
mi_fam_hel <- decostand(mi_fam, method = 'hellinger')
mi_fam_hel %>% tibble
```

    ## # A tibble: 50 x 1
    ##    .$Acaldive $Acalmacr $Adeltril $Alchcost $Alchlati $Crotbill $Hieralch
    ##         <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>
    ##  1      0.811     0.310     0         0.310     0         0.388     0    
    ##  2      0.655     0         0         0.378     0         0.378     0.535
    ##  3      0.302     0.174     0.739     0.522     0         0.174     0    
    ##  4      0.240     0.186     0.322     0.577     0         0.661     0.186
    ##  5      0.544     0         0.136     0.360     0         0.694     0.192
    ##  6      0.811     0         0         0.169     0.169     0.478     0    
    ##  7      0.648     0.141     0.141     0.283     0         0.632     0.2  
    ##  8      0.683     0         0         0.447     0         0.516     0    
    ##  9      0.378     0         0.309     0.378     0         0.724     0.309
    ## 10      0.764     0         0         0.577     0         0.289     0    
    ## # … with 40 more rows, and 3 more variables: $Huracrep <dbl>,
    ## #   $Sapibroa <dbl>, $Sapiglan <dbl>

``` r
mi_fam_hel_pca <- rda(mi_fam_hel)
summary(mi_fam_hel_pca)
```

    ## 
    ## Call:
    ## rda(X = mi_fam_hel) 
    ## 
    ## Partitioning of variance:
    ##               Inertia Proportion
    ## Total          0.1935          1
    ## Unconstrained  0.1935          1
    ## 
    ## Eigenvalues, and their contribution to the variance 
    ## 
    ## Importance of components:
    ##                           PC1     PC2     PC3     PC4     PC5     PC6
    ## Eigenvalue            0.06052 0.03923 0.02728 0.02056 0.01663 0.01558
    ## Proportion Explained  0.31282 0.20274 0.14102 0.10627 0.08594 0.08050
    ## Cumulative Proportion 0.31282 0.51557 0.65659 0.76286 0.84881 0.92931
    ##                           PC7      PC8       PC9     PC10
    ## Eigenvalue            0.01034 0.002344 0.0006353 0.000361
    ## Proportion Explained  0.05343 0.012113 0.0032838 0.001866
    ## Cumulative Proportion 0.98274 0.994850 0.9981342 1.000000
    ## 
    ## Scaling 2 for species and site scores
    ## * Species are scaled proportional to eigenvalues
    ## * Sites are unscaled: weighted dispersion equal on all dimensions
    ## * General scaling constant of scores:  1.754715 
    ## 
    ## 
    ## Species scores
    ## 
    ##               PC1       PC2      PC3      PC4       PC5      PC6
    ## Acaldive -0.67554 -0.055951  0.02452  0.09044 -0.218205  0.00584
    ## Acalmacr  0.10396  0.052367 -0.26322 -0.05346  0.011515  0.42354
    ## Adeltril  0.33626  0.618466 -0.07393  0.13412 -0.188020 -0.06296
    ## Alchcost  0.01744  0.017740 -0.31418 -0.09966  0.286459 -0.16408
    ## Alchlati -0.01248 -0.012637  0.01281  0.01826 -0.005994  0.01080
    ## Crotbill  0.59077 -0.403792  0.15683  0.10408 -0.104958  0.02736
    ## Hieralch  0.09301  0.040561  0.07512 -0.52141 -0.177416 -0.04569
    ## Huracrep -0.11532  0.243672  0.42152 -0.05810  0.208296  0.16904
    ## Sapibroa -0.02481 -0.003075  0.01586  0.01541  0.003395  0.01494
    ## Sapiglan  0.10453  0.106418  0.22764  0.04031  0.116372 -0.07622
    ## 
    ## 
    ## Site scores (weighted sums of species scores)
    ## 
    ##         PC1        PC2        PC3      PC4        PC5        PC6
    ## 1  -0.21580 -0.2319770 -0.3256287  0.24364 -9.235e-02  0.4524647
    ## 2  -0.09306 -0.2006705 -0.1073717 -0.62473 -2.574e-01 -0.2252215
    ## 3   0.27608  0.7082160 -0.5378199  0.31631  1.886e-01 -0.1966840
    ## 4   0.47374 -0.0292612 -0.4005631 -0.05283  2.018e-01 -0.0358374
    ## 5   0.21312 -0.2502419  0.0320147  0.07526 -9.383e-02 -0.2117736
    ## 6  -0.22726 -0.2301982  0.2333463  0.33257 -1.092e-01  0.1967319
    ## 7   0.08132 -0.1819581  0.0226051  0.02615 -2.146e-01  0.2295254
    ## 8  -0.10940 -0.2247441  0.0523049  0.21196  3.008e-01  0.0195098
    ## 9   0.39900 -0.1084700 -0.0897381 -0.10476 -2.342e-01 -0.2303127
    ## 10 -0.27236 -0.1817581 -0.3702026  0.16802  2.623e-01 -0.2659596
    ## 11  0.08843 -0.2071340 -0.2591244  0.29772  1.144e-01 -0.2416894
    ## 12  0.01217  0.0004932  0.0520161 -0.16007  6.996e-02 -0.0824702
    ## 13  0.18639 -0.2784362  0.2716981  0.10876 -2.899e-01  0.0302823
    ## 14  0.16275 -0.2008679 -0.2410217 -0.12028 -1.106e-01  0.0205786
    ## 15  0.32978 -0.0040433 -0.3271610  0.34444  1.631e-02  0.3386392
    ## 16  0.20019 -0.1127758 -0.1801318 -0.45562  3.028e-01  0.9984571
    ## 17  0.60812 -0.3798766  0.0092549  0.18366  7.647e-01 -0.2860754
    ## 18  0.32509 -0.3857639  0.2315527  0.06552 -4.193e-01  0.0006059
    ## 19  0.03922 -0.3298774 -0.0557423 -0.36153 -2.236e-01 -0.1738540
    ## 20 -0.07456  0.1483311 -0.5452460 -0.10638 -1.791e-01  0.6578724
    ## 21 -0.25562 -0.0807240 -0.2902487 -0.57542 -1.474e-02 -0.3210961
    ## 22 -0.05533  0.2774237 -0.2734233 -0.21381  3.904e-02  0.0996508
    ## 23  0.05566 -0.2613643  0.2978174 -0.30932 -1.712e-01 -0.0185157
    ## 24 -0.31713 -0.2052215  0.0783357  0.37925 -2.342e-01 -0.0436043
    ## 25  0.06875 -0.1860648  0.0347027  0.00886 -5.403e-01 -0.0782217
    ## 26 -0.14506  0.1431738  0.5198868 -0.34047  6.527e-01  0.1343694
    ## 27 -0.36849  0.0237401  0.1957500  0.20817  3.642e-01 -0.0178085
    ## 28 -0.37896  0.0235203  0.1744239 -0.26747  1.870e-02 -0.0355608
    ## 29  0.23103  0.4798564  0.2670760 -0.11252 -2.273e-01  0.2953892
    ## 30  0.33850  0.2927709  0.2222307  0.13682 -2.656e-01 -0.0764641
    ## 31  0.09404  0.1961144 -0.0005275 -0.08782 -1.181e-01 -0.1380764
    ## 32 -0.09579  0.3849601  0.0105086 -0.21798 -4.219e-02 -0.1125037
    ## 33  0.07074  0.1709820 -0.1375102 -0.04133  1.151e-01  0.1836580
    ## 34 -0.20737  0.1102076  0.0161092  0.32332 -1.856e-02 -0.0015148
    ## 35  0.05367  0.4271124 -0.0990375 -0.06814 -3.556e-01 -0.3765018
    ## 36 -0.41233  0.3020175 -0.0827888 -0.17382  6.385e-02 -0.1899600
    ## 37 -0.06627 -0.0138276  0.1132722 -0.03111 -2.215e-01  0.2819099
    ## 38  0.03235  0.1157753  0.1048838  0.16199 -8.932e-02 -0.1084828
    ## 39  0.19193  0.1674042  0.0542184  0.38693  2.168e-02 -0.0315564
    ## 40  0.03663  0.5283553  0.1440704 -0.02864 -8.567e-02 -0.1871996
    ## 41 -0.14306 -0.1194528 -0.1132442 -0.29012  3.397e-01 -0.1739024
    ## 42 -0.19990 -0.2128261  0.0130492  0.05168  3.947e-02 -0.0279455
    ## 43 -0.39084 -0.1416651 -0.3133487  0.20289  9.468e-02 -0.1990058
    ## 44 -0.35241  0.2137153  0.0622977  0.33714 -1.504e-01  0.0666108
    ## 45 -0.40417 -0.0501034  0.2583829  0.25101  5.531e-02  0.2434319
    ## 46  0.36961  0.0808489  0.4797336 -0.08566  2.988e-01  0.0097413
    ## 47  0.16893  0.1185054  0.4549411 -0.13532  2.227e-01  0.0224662
    ## 48 -0.13755  0.2223881  0.4120056  0.11662 -6.303e-05  0.1133010
    ## 49 -0.24051 -0.1381524 -0.1422344 -0.18886  1.334e-01 -0.1693746
    ## 50  0.05598 -0.1884559  0.0736258  0.21532  7.764e-02 -0.1380233

``` r
screeplot(
  mi_fam_hel_pca,
  bstick = TRUE,
  npcs = length(mi_fam_hel_pca$CA$eig)
)
```

![](TO1_no_restringida_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
mi_fam_hel_pca_sc1 <- scores(mi_fam_hel_pca,
                             display = "species", scaling = 1)
mi_fam_hel_pca_sc2 <- scores(mi_fam_hel_pca,
                             display = "species", scaling = 2)
par(mfrow = c(1, 2))
cleanplot.pca(mi_fam_hel_pca, scaling = 1, mar.percent = 0.06, cex.char1 = 0.7)
cleanplot.pca(mi_fam_hel_pca, scaling = 2, mar.percent = 0.06, cex.char1 = 0.7)
```

![](TO1_no_restringida_files/figure-markdown_github/unnamed-chunk-12-2.png)

``` r
par(mfrow = c(1, 1))
```

Si intentáramos realizar el PCA a datos de comunidad no transformados, no recogeríamos apropiadamente las tendencias y patrones, debido a la presencia de doble-ceros y valores extremos.

Las especies que contribuyen mucho a los ejes 1 y 2 del PCA (aquellas cuyos vectores sobresalen el círculo de contribución uniforme), podrían coincidir con aquellas que podríamos considerar como especies indicadoras o que muestran preferencias por determinados hábitats.

Evaluaré el ajuste del PCA de datos de comunidad a datos ambientales, mediante la función `envfit`

``` r
biplot(
  mi_fam_hel_pca,
  main = "PCA, escalamiento 2, ajuste a variables ambientales")
(mi_fam_hel_pca_envfit <- envfit(mi_fam_hel_pca, env_suelo, scaling = 2))
```

    ## 
    ## ***VECTORS
    ## 
    ##             PC1      PC2     r2 Pr(>r)    
    ## Al      0.88198 -0.47128 0.1214  0.046 *  
    ## B      -0.72678  0.68687 0.1685  0.013 *  
    ## Ca     -0.90760  0.41984 0.0787  0.143    
    ## Cu     -0.95327  0.30212 0.1205  0.050 *  
    ## Fe     -0.99739 -0.07225 0.0702  0.176    
    ## K      -0.92817  0.37215 0.1633  0.021 *  
    ## Mg     -0.97066  0.24046 0.0305  0.495    
    ## Mn     -0.71610  0.69800 0.0821  0.131    
    ## P       0.74036  0.67221 0.0482  0.307    
    ## Zn     -0.73896  0.67374 0.1138  0.049 *  
    ## N      -0.61960  0.78492 0.3026  0.001 ***
    ## N.min. -0.84074  0.54144 0.0696  0.178    
    ## pH     -0.73946  0.67320 0.3859  0.001 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 999

``` r
plot(mi_fam_hel_pca_envfit, p.max = 0.05 , col = 3)
```

![](TO1_no_restringida_files/figure-markdown_github/unnamed-chunk-13-1.png)

Comento los resultados en el vídeo asociado. También probaré ajuste con todas las numéricas de la matriz ambiental, excluyendo por supuesto la columna `id`:

NOTA: te recomiendo probar otros métodos de selección de variables, como por ejemplo, usando la función `step` para seleccionar fórmulas de modelos basadas en AIC.

``` r
env_num <- bci_env_grid %>%
  select_if(is.numeric) %>%
  select(-id) %>%
  st_drop_geometry
(mi_fam_hel_pca_envfit_num <- envfit(mi_fam_hel_pca, env_num, scaling = 2))
```

    ## 
    ## ***VECTORS
    ## 
    ##                                 PC1      PC2     r2 Pr(>r)    
    ## heterogeneidad_ambiental   -0.89596  0.44414 0.0071  0.868    
    ## UTM.EW                     -0.69425  0.71973 0.1945  0.008 ** 
    ## UTM.NS                      0.80718  0.59031 0.0281  0.490    
    ## geomorf_llanura_pct        -0.92772  0.37328 0.0104  0.797    
    ## geomorf_pico_pct           -0.23623 -0.97170 0.0199  0.684    
    ## geomorf_interfluvio_pct     0.60776 -0.79412 0.0038  0.916    
    ## geomorf_hombrera_pct       -0.97028 -0.24199 0.0146  0.711    
    ## geomorf_espolón/gajo_pct   -0.99998 -0.00559 0.0018  0.959    
    ## geomorf_vertiente_pct       0.83320 -0.55297 0.0160  0.694    
    ## geomorf_vaguada_pct         0.29169  0.95651 0.0167  0.676    
    ## geomorf_piedemonte_pct     -0.26523  0.96419 0.0321  0.439    
    ## geomorf_valle_pct           0.98956  0.14409 0.0285  0.492    
    ## geomorf_sima_pct            0.96498 -0.26233 0.0567  0.205    
    ## Al                          0.88198 -0.47128 0.1214  0.050 *  
    ## B                          -0.72678  0.68687 0.1685  0.015 *  
    ## Ca                         -0.90760  0.41984 0.0787  0.152    
    ## Cu                         -0.95327  0.30212 0.1205  0.054 .  
    ## Fe                         -0.99739 -0.07225 0.0702  0.200    
    ## K                          -0.92817  0.37215 0.1633  0.017 *  
    ## Mg                         -0.97066  0.24046 0.0305  0.499    
    ## Mn                         -0.71610  0.69800 0.0821  0.130    
    ## P                           0.74036  0.67221 0.0482  0.330    
    ## Zn                         -0.73896  0.67374 0.1138  0.061 .  
    ## N                          -0.61960  0.78492 0.3026  0.001 ***
    ## N.min.                     -0.84074  0.54144 0.0696  0.201    
    ## pH                         -0.73946  0.67320 0.3859  0.001 ***
    ## elevacion_media             0.07222  0.99739 0.0125  0.752    
    ## pendiente_media             0.10920 -0.99402 0.0230  0.580    
    ## orientacion_media          -0.73757 -0.67527 0.0829  0.131    
    ## curvatura_perfil_media      0.40111 -0.91603 0.0031  0.939    
    ## curvatura_tangencial_media -0.91149 -0.41131 0.0613  0.236    
    ## abundancia_global          -0.74394  0.66825 0.0297  0.457    
    ## riqueza_global              0.22565 -0.97421 0.1587  0.020 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 999

``` r
biplot(
  mi_fam_hel_pca,
  main = "PCA, escalamiento 2, ajuste a variables ambientales")
plot(mi_fam_hel_pca_envfit_num, p.max = 0.05 , col = 3)
```

![](TO1_no_restringida_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
biplot(
  mi_fam_hel_pca,
  main = "PCA, escalamiento 2, ajuste a variables ambientales")
plot(mi_fam_hel_pca_envfit_num, p.max = 0.1 , col = 3)
```

![](TO1_no_restringida_files/figure-markdown_github/unnamed-chunk-14-2.png)

Comento los resultados en el vídeo asociado.

¿Cuándo o a qué datos aplicar PCA?

-   PCA no es especialmente sensible a datos muy desviados de la normalidad.
-   Como toda técnica, PCA tiene limitaciones.
-   Las variables deben ser dimensionalmente homogéneas (unidades comparables o adimensionales).
-   No usar en matriz transpuestas (no hace sentido la covarianza entre objetos).
-   Es posible usar PCA con dato de presencia/ausencia, en cuyo caso, la matriz de comunidad debe transformarse a Hellinger, cuerdas o log-chord.
-   Las relaciones entre variables se miden por ángulos, no por proximidad de las puntas de los vectores.

### Análisis de correspondencia (CA)

``` r
mi_fam_ca <- cca(mi_fam)
summary(mi_fam_ca)
```

    ## 
    ## Call:
    ## cca(X = mi_fam) 
    ## 
    ## Partitioning of scaled Chi-square:
    ##               Inertia Proportion
    ## Total          0.7826          1
    ## Unconstrained  0.7826          1
    ## 
    ## Eigenvalues, and their contribution to the scaled Chi-square 
    ## 
    ## Importance of components:
    ##                          CA1    CA2    CA3     CA4     CA5     CA6     CA7
    ## Eigenvalue            0.1899 0.1695 0.1613 0.08067 0.05744 0.04670 0.03260
    ## Proportion Explained  0.2427 0.2166 0.2062 0.10309 0.07339 0.05967 0.04166
    ## Cumulative Proportion 0.2427 0.4593 0.6655 0.76856 0.84196 0.90163 0.94329
    ##                           CA8     CA9
    ## Eigenvalue            0.02772 0.01665
    ## Proportion Explained  0.03542 0.02128
    ## Cumulative Proportion 0.97872 1.00000
    ## 
    ## Scaling 2 for species and site scores
    ## * Species are scaled proportional to eigenvalues
    ## * Sites are unscaled: weighted dispersion equal on all dimensions
    ## 
    ## 
    ## Species scores
    ## 
    ##               CA1     CA2      CA3      CA4      CA5      CA6
    ## Acaldive -0.41822 -0.1200  0.05675  0.03412 -0.09566  0.07627
    ## Acalmacr  1.27929 -1.6292  1.77693 -0.33351 -0.09258  0.00203
    ## Adeltril  0.80759 -0.6186 -1.13984  0.03480 -0.27775  0.17875
    ## Alchcost  0.11171 -0.1383 -0.05666  0.48799  0.29982 -0.28897
    ## Alchlati -1.28896  0.1684  0.32770 -0.56237 -2.11252  0.18981
    ## Crotbill  0.35049  0.5382  0.13659 -0.05343 -0.11484 -0.04342
    ## Hieralch  0.15879  0.1740  0.01099 -0.20653  0.72675  0.60922
    ## Huracrep -0.28744 -0.2889 -0.40925 -0.93540  0.34248 -0.49877
    ## Sapibroa -2.03063 -0.6759 -0.01491 -1.74229 -1.13817 -1.81387
    ## Sapiglan  0.08226  0.2021 -0.46154 -0.64361  0.15532 -0.15116
    ## 
    ## 
    ## Site scores (weighted averages of species scores)
    ## 
    ##         CA1       CA2       CA3        CA4      CA5      CA6
    ## 1  -0.46758 -0.986691  1.381249  0.3618978 -1.05045  0.34466
    ## 2  -0.35723  0.326968  0.240982  0.2193256  3.36152  3.41050
    ## 3   2.55304 -2.436296 -3.644517  1.5363078 -1.39300  0.42372
    ## 4   1.57670  0.400622 -0.075774  1.5650853  0.65144 -1.52765
    ## 5   0.43812  1.227909  0.232059  0.2082098 -0.30029 -0.33175
    ## 6  -1.28896  0.168382  0.327704 -0.5623727 -2.11252  0.18981
    ## 7   0.08331  0.649323  0.489237 -0.0117518 -0.58528  0.20465
    ## 8  -0.51882  0.239616  0.150569  0.4575298  0.13109 -1.43539
    ## 9   1.22082  1.195694 -0.222808  0.3748638  0.20496  0.46924
    ## 10 -0.93477 -0.420238  0.158683  2.2078051  0.60183 -1.18740
    ## 11  0.14325  0.547964  0.146229  1.5904260  0.03186 -1.35981
    ## 12 -0.16822  0.156542 -0.148729  0.2907156  1.08275 -0.59010
    ## 13  0.37122  1.483452  0.412749 -0.4230666 -0.96468  0.26190
    ## 14  0.42968  0.589879  0.451307  0.8381807  0.47076 -0.12501
    ## 15  1.65504 -0.446369  0.765368  0.2204333 -1.02478 -0.62585
    ## 16  2.41911 -2.545362  3.607609 -1.8257848  1.13315 -0.28415
    ## 17  1.35471  1.929500  0.083748  0.1010131  0.39337 -2.53280
    ## 18  0.88224  2.052800  0.559976 -0.5312466 -0.94544  0.45838
    ## 19 -0.01963  0.818168  0.371249  0.3474879  1.47862  1.59574
    ## 20  1.27799 -2.954569  2.658906 -0.2848562 -0.43826  0.69705
    ## 21 -0.71551 -0.260348  0.124074  1.2904170  3.02327  1.74466
    ## 22  0.11990 -1.013219 -0.350985  0.9672137  1.08780 -0.01008
    ## 23 -0.02511  0.982252  0.281518 -0.6845147  1.15564  1.60155
    ## 24 -1.52844 -0.157518  0.240952  0.0276402 -1.33327  0.91006
    ## 25 -0.09657  0.773327  0.210679  0.0008766 -0.57069  1.31507
    ## 26 -0.86327 -0.542234 -0.994560 -3.9560466  3.35362 -3.59763
    ## 27 -1.40778 -0.484155 -0.223880 -0.6205911  0.18110 -0.88315
    ## 28 -1.42810 -0.445052  0.004666 -0.3526016  0.84677  0.98824
    ## 29  1.20039 -0.979696 -1.936304 -1.8627769  0.16381  1.00161
    ## 30  1.46077 -0.086218 -1.851492 -0.6087086 -1.21358  0.81354
    ## 31  0.36326 -0.226074 -0.906608  0.3419582  0.39264  0.29082
    ## 32 -0.13432 -0.884090 -1.187333 -0.1070408  1.12951  0.38267
    ## 33  0.34907 -0.555982 -0.301556  0.4669224  0.54493 -0.97837
    ## 34 -0.94043 -0.488850 -0.297998  0.2454752 -0.60121 -0.23871
    ## 35  0.64231 -0.952152 -1.983618  0.0444437 -0.06822  2.09752
    ## 36 -1.23712 -0.740689 -0.233132  0.5977449  0.56490  0.47740
    ## 37 -0.48601 -0.098663  0.332084 -0.4151768 -0.26975  0.69931
    ## 38 -0.06669 -0.036728 -0.559148  0.1240530 -0.53492 -0.09729
    ## 39  0.77879 -0.001562 -1.120771  0.0022577 -0.89674 -0.74764
    ## 40  0.58503 -1.080692 -2.173028 -0.6815031 -0.03587  0.85475
    ## 41 -0.47557 -0.169176  0.027728  1.5211738  1.95141 -1.21259
    ## 42 -0.91525  0.036775  0.250892  0.6522084 -0.23992 -0.25683
    ## 43 -1.41964 -0.562414  0.220457  1.5987748 -0.18319 -0.17841
    ## 44 -1.40465 -0.822217 -0.362447 -0.1284044 -1.00566  0.42290
    ## 45 -2.03063 -0.675917 -0.014908 -1.7422870 -1.13817 -1.81387
    ## 46  0.77829  0.927812 -0.731363 -2.3757015  1.13228 -1.75154
    ## 47  0.15994  0.397376 -0.708559 -2.2965214  1.50095 -1.47232
    ## 48 -0.86075 -0.456278 -0.722730 -1.9421954  0.17154 -0.94705
    ## 49 -0.89280 -0.261529  0.156683  1.4041531  0.95283 -0.27661
    ## 50 -0.09285  0.722173  0.175545  0.4180559 -0.41930 -0.74205

``` r
summary(mi_fam_ca, scaling = 1)
```

    ## 
    ## Call:
    ## cca(X = mi_fam) 
    ## 
    ## Partitioning of scaled Chi-square:
    ##               Inertia Proportion
    ## Total          0.7826          1
    ## Unconstrained  0.7826          1
    ## 
    ## Eigenvalues, and their contribution to the scaled Chi-square 
    ## 
    ## Importance of components:
    ##                          CA1    CA2    CA3     CA4     CA5     CA6     CA7
    ## Eigenvalue            0.1899 0.1695 0.1613 0.08067 0.05744 0.04670 0.03260
    ## Proportion Explained  0.2427 0.2166 0.2062 0.10309 0.07339 0.05967 0.04166
    ## Cumulative Proportion 0.2427 0.4593 0.6655 0.76856 0.84196 0.90163 0.94329
    ##                           CA8     CA9
    ## Eigenvalue            0.02772 0.01665
    ## Proportion Explained  0.03542 0.02128
    ## Cumulative Proportion 0.97872 1.00000
    ## 
    ## Scaling 1 for species and site scores
    ## * Sites are scaled proportional to eigenvalues
    ## * Species are unscaled: weighted dispersion equal on all dimensions
    ## 
    ## 
    ## Species scores
    ## 
    ##              CA1     CA2      CA3     CA4     CA5       CA6
    ## Acaldive -0.9597 -0.2914  0.14129  0.1201 -0.3992  0.352964
    ## Acalmacr  2.9356 -3.9570  4.42378 -1.1742 -0.3863  0.009395
    ## Adeltril  1.8532 -1.5025 -2.83770  0.1225 -1.1590  0.827165
    ## Alchcost  0.2563 -0.3360 -0.14105  1.7181  1.2510 -1.337237
    ## Alchlati -2.9578  0.4090  0.81584 -1.9800 -8.8148  0.878362
    ## Crotbill  0.8043  1.3073  0.34005 -0.1881 -0.4792 -0.200943
    ## Hieralch  0.3644  0.4226  0.02736 -0.7271  3.0325  2.819160
    ## Huracrep -0.6596 -0.7018 -1.01886 -3.2933  1.4291 -2.308058
    ## Sapibroa -4.6598 -1.6417 -0.03712 -6.1341 -4.7491 -8.393738
    ## Sapiglan  0.1888  0.4910 -1.14903 -2.2660  0.6481 -0.699491
    ## 
    ## 
    ## Site scores (weighted averages of species scores)
    ## 
    ##          CA1        CA2       CA3        CA4       CA5       CA6
    ## 1  -0.203764 -0.4062474  0.554815  0.1027907 -0.251747  0.074479
    ## 2  -0.155672  0.1346217  0.096797  0.0622956  0.805611  0.737005
    ## 3   1.112563 -1.0030891 -1.463918  0.4363612 -0.333843  0.091566
    ## 4   0.687093  0.1649468 -0.030437  0.4445350  0.156122 -0.330123
    ## 5   0.190924  0.5055634  0.093213  0.0591383 -0.071966 -0.071690
    ## 6  -0.561704  0.0693276  0.131631 -0.1597321 -0.506279  0.041018
    ## 7   0.036304  0.2673437  0.196515 -0.0033379 -0.140266  0.044225
    ## 8  -0.226091  0.0986563  0.060480  0.1299533  0.031418 -0.310186
    ## 9   0.532008  0.4922997 -0.089497  0.1064735  0.049120  0.101403
    ## 10 -0.407353 -0.1730232  0.063739  0.6270882  0.144234 -0.256595
    ## 11  0.062427  0.2256117  0.058737  0.4517325  0.007636 -0.293854
    ## 12 -0.073308  0.0644526 -0.059741  0.0825727  0.259488 -0.127519
    ## 13  0.161770  0.6107773  0.165792 -0.1201646 -0.231191  0.056596
    ## 14  0.187245  0.2428693  0.181279  0.2380705  0.112820 -0.027014
    ## 15  0.721233 -0.1837821  0.307431  0.0626102 -0.245596 -0.135246
    ## 16  1.054200 -1.0479943  1.449093 -0.5185821  0.271567 -0.061405
    ## 17  0.590355  0.7944273  0.033640  0.0286910  0.094274 -0.547335
    ## 18  0.384464  0.8451932  0.224929 -0.1508913 -0.226581  0.099056
    ## 19 -0.008553  0.3368621  0.149122  0.0986978  0.354361  0.344837
    ## 20  0.556922 -1.2164759  1.068021 -0.0809084 -0.105032  0.150631
    ## 21 -0.311804 -0.1071921  0.049837  0.3665203  0.724546  0.377019
    ## 22  0.052249 -0.4171697 -0.140982  0.2747201  0.260699 -0.002178
    ## 23 -0.010943  0.4044199  0.113079 -0.1944244  0.276956  0.346093
    ## 24 -0.666062 -0.0648544  0.096785  0.0078507 -0.319526  0.196663
    ## 25 -0.042082  0.3183997  0.084625  0.0002490 -0.136769  0.284184
    ## 26 -0.376194 -0.2232522 -0.399492 -1.1236455  0.803719 -0.777442
    ## 27 -0.613483 -0.1993398 -0.089928 -0.1762680  0.043403 -0.190847
    ## 28 -0.622337 -0.1832399  0.001874 -0.1001503  0.202934  0.213557
    ## 29  0.523103 -0.4033672 -0.777768 -0.5290890  0.039259  0.216446
    ## 30  0.636574 -0.0354981 -0.743701 -0.1728930 -0.290843  0.175805
    ## 31  0.158301 -0.0930806 -0.364163  0.0971272  0.094098  0.062846
    ## 32 -0.058536 -0.3640036 -0.476924 -0.0304031  0.270695  0.082694
    ## 33  0.152118 -0.2289129 -0.121128  0.1326211  0.130596 -0.211425
    ## 34 -0.409819 -0.2012728 -0.119699  0.0697229 -0.144084 -0.051584
    ## 35  0.279906 -0.3920268 -0.796773  0.0126235 -0.016349  0.453270
    ## 36 -0.539113 -0.3049619 -0.093644  0.1697789  0.135381  0.103165
    ## 37 -0.211793 -0.0406222  0.133390 -0.1179237 -0.064648  0.151120
    ## 38 -0.029062 -0.0151221 -0.224597  0.0352351 -0.128198 -0.021024
    ## 39  0.339379 -0.0006431 -0.450188  0.0006413 -0.214910 -0.161565
    ## 40  0.254944 -0.4449500 -0.872855 -0.1935690 -0.008597  0.184710
    ## 41 -0.207242 -0.0696544  0.011138  0.4320627  0.467668 -0.262039
    ## 42 -0.398846  0.0151411  0.100777  0.1852483 -0.057499 -0.055500
    ## 43 -0.618648 -0.2315609  0.088553  0.4541039 -0.043902 -0.038554
    ## 44 -0.612118 -0.3385289 -0.145586 -0.0364710 -0.241013  0.091388
    ## 45 -0.884906 -0.2782931 -0.005988 -0.4948660 -0.272769 -0.391975
    ## 46  0.339164  0.3820052 -0.293771 -0.6747763  0.271360 -0.378504
    ## 47  0.069700  0.1636103 -0.284612 -0.6522866  0.359713 -0.318166
    ## 48 -0.375096 -0.1878619 -0.290304 -0.5516465  0.041111 -0.204656
    ## 49 -0.389064 -0.1076787  0.062936  0.3988250  0.228353 -0.059776
    ## 50 -0.040462  0.2973382  0.070513  0.1187414 -0.100487 -0.160355

Screeplot

``` r
screeplot(mi_fam_ca, bstick = TRUE, npcs = length(mi_fam_ca$CA$eig))
```

![](TO1_no_restringida_files/figure-markdown_github/unnamed-chunk-16-1.png)

Biplots

``` r
par(mfrow = c(1, 2))
plot(mi_fam_ca,
     scaling = 1,
     main = "Análisis de correspondencia, escalamiento 1"
)
plot(mi_fam_ca,
     scaling = 2, # Por defecto scaling=2, lo escribo sólo para fines didáticos
     main = "Análisis de correspondencia, escalamiento 2")
```

![](TO1_no_restringida_files/figure-markdown_github/unnamed-chunk-17-1.png)

``` r
par(mfrow = c(1, 1))
```

Excluyendo especie *Thevetia ahouai*, abreviada como *Thevahou*.

``` r
# mi_fam_ca <- cca(mi_fam[, -grep('Thevahou', colnames(mi_fam))])
# summary(mi_fam_ca)
# summary(mi_fam_ca, scaling = 1)
# screeplot(mi_fam_ca, bstick = TRUE, npcs = length(mi_fam_ca$CA$eig))
# par(mfrow = c(1, 2))
# plot(mi_fam_ca,
#      scaling = 1,
#      main = "CA, escalamiento 1, sin Thevetia ahouai"
# )
# plot(mi_fam_ca,
#      scaling = 2,
#      main = "CA, escalamiento 2, sin Thevetia ahouai")
# par(mfrow = c(1, 1))
```

Análisis de coordenadas principales (PCoA)

Las técnicas de ordenación anteriores preservan la distancia euclídea entre los objetos. Si necesitaras ordenar objetos usando una métrica diferente, por ejemplo, la de Gower para datos mixtos, entonces PCA y CA serían inútiles y, en su lugar, podrías usar PCoA. Paso el resto de la explicación al vídeo asociado.

La función que realiza el PCoA en `{vegan}` es `cmdscale` (de *Classical (Metric) Multidimensional Scaling*), y se le suministra una matriz de distancias.

``` r
mi_fam_d_bray <- vegdist(mi_fam, method = 'bray') # En realidad, 'bray' es la opción por defecto
mi_fam_d_bray_pcoa <- cmdscale(
  mi_fam_d_bray,
  k = (nrow(mi_fam) - 1),
  add = T,
  eig = TRUE)
round(mi_fam_d_bray_pcoa$eig, 2)
```

    ##  [1] 4.51 3.02 2.00 1.71 1.51 1.30 0.92 0.80 0.73 0.68 0.57 0.50 0.48 0.42
    ## [15] 0.36 0.33 0.32 0.28 0.25 0.24 0.21 0.19 0.19 0.19 0.17 0.16 0.16 0.15
    ## [29] 0.14 0.14 0.13 0.13 0.12 0.12 0.12 0.11 0.11 0.10 0.10 0.10 0.09 0.09
    ## [43] 0.09 0.07 0.07 0.06 0.05 0.02 0.00 0.00

``` r
round(sum(mi_fam_d_bray_pcoa$eig[mi_fam_d_bray_pcoa$eig<0]),2)
```

    ## [1] 0

``` r
round(sum(mi_fam_d_bray_pcoa$eig[mi_fam_d_bray_pcoa$eig>=0]),2)
```

    ## [1] 24.31

``` r
ordiplot(scores(mi_fam_d_bray_pcoa, choices = c(1, 2)),
         type = "t",
         main = "PCoA con promedios ponderados de especies")
```

    ## species scores not available

``` r
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
mi_fam_d_bray_pcoa_wa <- wascores(mi_fam_d_bray_pcoa$points[, 1:2], mi_fam)
text(
  mi_fam_d_bray_pcoa_wa,
  rownames(mi_fam_d_bray_pcoa_wa),
  cex = 0.7, col = "red")
(mi_fam_d_bray_pcoa_env <- envfit(mi_fam_d_bray_pcoa, env_num))
```

    ## 
    ## ***VECTORS
    ## 
    ##                                Dim1     Dim2     r2 Pr(>r)   
    ## heterogeneidad_ambiental    0.93607 -0.35183 0.1057  0.077 . 
    ## UTM.EW                      0.68684  0.72681 0.1484  0.020 * 
    ## UTM.NS                      0.31898 -0.94776 0.0352  0.450   
    ## geomorf_llanura_pct         0.13865  0.99034 0.0150  0.709   
    ## geomorf_pico_pct            0.90944 -0.41584 0.0064  0.908   
    ## geomorf_interfluvio_pct    -0.40096  0.91609 0.0143  0.743   
    ## geomorf_hombrera_pct        0.31807  0.94807 0.0094  0.783   
    ## geomorf_espolón/gajo_pct   -0.66742  0.74468 0.0062  0.870   
    ## geomorf_vertiente_pct       0.17595 -0.98440 0.0328  0.449   
    ## geomorf_vaguada_pct        -0.94297  0.33289 0.0086  0.813   
    ## geomorf_piedemonte_pct     -0.70424  0.70996 0.0314  0.454   
    ## geomorf_valle_pct          -0.94766  0.31927 0.0719  0.180   
    ## geomorf_sima_pct           -0.73138 -0.68197 0.0619  0.243   
    ## Al                         -0.34992 -0.93678 0.0598  0.237   
    ## B                           0.67238  0.74021 0.1231  0.040 * 
    ## Ca                          0.88503  0.46554 0.1139  0.048 * 
    ## Cu                          0.82847  0.56003 0.1528  0.014 * 
    ## Fe                          0.97790  0.20907 0.2315  0.002 **
    ## K                           0.78737  0.61648 0.2054  0.006 **
    ## Mg                          0.97144  0.23727 0.0861  0.119   
    ## Mn                          0.59066  0.80692 0.0496  0.298   
    ## P                           0.17959 -0.98374 0.0109  0.765   
    ## Zn                          0.80897  0.58784 0.1029  0.072 . 
    ## N                           0.25926  0.96581 0.1340  0.032 * 
    ## N.min.                      0.81242  0.58308 0.0838  0.125   
    ## pH                          0.24541  0.96942 0.2419  0.003 **
    ## elevacion_media             0.90791 -0.41916 0.0306  0.500   
    ## pendiente_media             0.83084 -0.55652 0.0721  0.147   
    ## orientacion_media           0.96026  0.27911 0.0756  0.150   
    ## curvatura_perfil_media     -0.05313 -0.99859 0.0132  0.733   
    ## curvatura_tangencial_media  0.77988  0.62593 0.0914  0.118   
    ## abundancia_global           0.84914  0.52816 0.0314  0.515   
    ## riqueza_global              0.69266 -0.72126 0.2399  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 999

``` r
plot(mi_fam_d_bray_pcoa_env, p.max = 0.05, col = 3)
```

![](TO1_no_restringida_files/figure-markdown_github/unnamed-chunk-19-1.png)
