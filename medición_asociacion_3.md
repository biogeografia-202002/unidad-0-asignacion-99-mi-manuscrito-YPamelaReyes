Medición de asociación. Modo R aplicado a mi familia asignada
================
JR
3 de noviembre, 2020

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
library(adespatial)
```

    ## Registered S3 methods overwritten by 'adegraphics':
    ##   method         from
    ##   biplot.dudi    ade4
    ##   kplot.foucart  ade4
    ##   kplot.mcoa     ade4
    ##   kplot.mfa      ade4
    ##   kplot.pta      ade4
    ##   kplot.sepan    ade4
    ##   kplot.statis   ade4
    ##   scatter.coa    ade4
    ##   scatter.dudi   ade4
    ##   scatter.nipals ade4
    ##   scatter.pco    ade4
    ##   score.acm      ade4
    ##   score.mix      ade4
    ##   score.pca      ade4
    ##   screeplot.dudi ade4

    ## Registered S3 method overwritten by 'spdep':
    ##   method   from
    ##   plot.mst ape

    ## Registered S3 methods overwritten by 'adespatial':
    ##   method             from       
    ##   plot.multispati    adegraphics
    ##   print.multispati   ade4       
    ##   summary.multispati ade4

``` r
library(broom)
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
library(gclus)
```

    ## Loading required package: cluster

    ## Registered S3 method overwritten by 'gclus':
    ##   method         from 
    ##   reorder.hclust vegan

``` r
source('biodata/funciones.R')
```

### Cargar datos

``` r
load('biodata/matriz_ambiental.Rdata')
load('biodata/Euphorbiaceae.Rdata')
```

Modo R: matrices de dependencia entre variables (índice de correlación)
-----------------------------------------------------------------------

### Modo R para datos cuantitativos de especies (abundancia)

En este caso, las variables usaré los valores de abundancias de especies como variables. Es decir, **compararé el grado de asociación entre especies, NO entre sitios**.

Aunque se podría usar el índice de correlación como métrica de la dependencia (tal como mostré en el script `aed_5_correlaciones_variables_ambientales.R`), tanto los doble-ceros (ausencias de una misma especie en dos lugares), como los *outliers* de abundancia, contribuirían a aumentar de manera ficticia el índice de correlación de Pearson.

Por tal razón, es recomendable aplicar la transformación *Chi* a la matriz de comunidad transpuesta. Al utilizar una matriz transpuesto, lograrás comparar especies, NO sitios (recuerda que en modo R comparamos descriptores, no objetos). Explico el procedimiento a continuación, paso a paso.

Primero, sustituyo el caracter de espacio por un <enter> en los nombres de las especies (caracter \\n), para facilitar la lectura de los nombres en la diagonal del mapa de calor. Luego transpongo la matriz.

``` r
mi_fam_t <- mc_ephrb %>% 
  rename_all(gsub, pattern = ' ', replacement = '\n') %>% 
  t()
mi_fam_t %>% tibble
```

    ## # A tibble: 10 x 1
    ##    .[,"1"] [,"2"] [,"3"] [,"4"] [,"5"] [,"6"] [,"7"] [,"8"] [,"9"] [,"10"]
    ##      <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>
    ##  1      48      3      3      5     16     23     21      7      3       7
    ##  2       7      0      1      3      0      0      1      0      0       0
    ##  3       0      0     18      9      1      0      1      0      2       0
    ##  4       7      1      9     29      7      1      4      3      3       4
    ##  5       0      0      0      0      0      1      0      0      0       0
    ##  6      11      1      1     38     26      8     20      4     11       1
    ##  7       0      2      0      3      2      0      2      0      2       0
    ##  8       0      0      0      0      0      2      1      1      0       0
    ##  9       0      0      0      0      0      0      0      0      0       0
    ## 10       0      0      1      0      2      0      0      0      0       0
    ## # … with 40 more variables: [,"11"] <dbl>, [,"12"] <dbl>, [,"13"] <dbl>,
    ## #   [,"14"] <dbl>, [,"15"] <dbl>, [,"16"] <dbl>, [,"17"] <dbl>,
    ## #   [,"18"] <dbl>, [,"19"] <dbl>, [,"20"] <dbl>, [,"21"] <dbl>,
    ## #   [,"22"] <dbl>, [,"23"] <dbl>, [,"24"] <dbl>, [,"25"] <dbl>,
    ## #   [,"26"] <dbl>, [,"27"] <dbl>, [,"28"] <dbl>, [,"29"] <dbl>,
    ## #   [,"30"] <dbl>, [,"31"] <dbl>, [,"32"] <dbl>, [,"33"] <dbl>,
    ## #   [,"34"] <dbl>, [,"35"] <dbl>, [,"36"] <dbl>, [,"37"] <dbl>,
    ## #   [,"38"] <dbl>, [,"39"] <dbl>, [,"40"] <dbl>, [,"41"] <dbl>,
    ## #   [,"42"] <dbl>, [,"43"] <dbl>, [,"44"] <dbl>, [,"45"] <dbl>,
    ## #   [,"46"] <dbl>, [,"47"] <dbl>, [,"48"] <dbl>, [,"49"] <dbl>,
    ## #   [,"50"] <dbl>

Segundo, transformo la matriz transpuesta usando estandarización *Chi*.

``` r
mi_fam_t_chi <- decostand(mi_fam_t, "chi.square")
mi_fam_t_chi %>% tibble
```

    ## # A tibble: 10 x 1
    ##    .[,"1"] [,"2"] [,"3"] [,"4"] [,"5"] [,"6"] [,"7"] [,"8"] [,"9"] [,"10"]
    ##      <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>
    ##  1   0.270 0.0545 0.0251 0.0258 0.105  0.187  0.143  0.0869 0.0315  0.0972
    ##  2   0.775 0      0.165  0.304  0      0      0.134  0      0       0     
    ##  3   0     0      1.08   0.332  0.0468 0      0.0487 0      0.150   0     
    ##  4   0.128 0.0589 0.244  0.484  0.148  0.0263 0.0881 0.121  0.102   0.180 
    ##  5   0     0      0      0      0      8.32   0      0      0       0     
    ##  6   0.100 0.0295 0.0136 0.318  0.276  0.105  0.221  0.0805 0.187   0.0225
    ##  7   0     0.315  0      0.134  0.113  0      0.118  0      0.182   0     
    ##  8   0     0      0      0      0      0.175  0.0732 0.134  0       0     
    ##  9   0     0      0      0      0      0      0      0      0       0     
    ## 10   0     0      0.214  0      0.335  0      0      0      0       0     
    ## # … with 40 more variables: [,"11"] <dbl>, [,"12"] <dbl>, [,"13"] <dbl>,
    ## #   [,"14"] <dbl>, [,"15"] <dbl>, [,"16"] <dbl>, [,"17"] <dbl>,
    ## #   [,"18"] <dbl>, [,"19"] <dbl>, [,"20"] <dbl>, [,"21"] <dbl>,
    ## #   [,"22"] <dbl>, [,"23"] <dbl>, [,"24"] <dbl>, [,"25"] <dbl>,
    ## #   [,"26"] <dbl>, [,"27"] <dbl>, [,"28"] <dbl>, [,"29"] <dbl>,
    ## #   [,"30"] <dbl>, [,"31"] <dbl>, [,"32"] <dbl>, [,"33"] <dbl>,
    ## #   [,"34"] <dbl>, [,"35"] <dbl>, [,"36"] <dbl>, [,"37"] <dbl>,
    ## #   [,"38"] <dbl>, [,"39"] <dbl>, [,"40"] <dbl>, [,"41"] <dbl>,
    ## #   [,"42"] <dbl>, [,"43"] <dbl>, [,"44"] <dbl>, [,"45"] <dbl>,
    ## #   [,"46"] <dbl>, [,"47"] <dbl>, [,"48"] <dbl>, [,"49"] <dbl>,
    ## #   [,"50"] <dbl>

Tercero, calculo la distancia euclídea.

``` r
mi_fam_t_chi_d <- dist(mi_fam_t_chi)
mi_fam_t_chi_d %>% tidy
```

    ## # A tibble: 45 x 3
    ##    item1                      item2                    distance
    ##    <fct>                      <fct>                       <dbl>
    ##  1 "Acalypha\nmacrostachya"   "Acalypha\ndiversifolia"    2.87 
    ##  2 "Adelia\ntriloba"          "Acalypha\ndiversifolia"    1.80 
    ##  3 "Alchornea\ncostaricensis" "Acalypha\ndiversifolia"    0.893
    ##  4 "Alchornea\nlatifolia"     "Acalypha\ndiversifolia"    8.20 
    ##  5 "Croton\nbillbergianus"    "Acalypha\ndiversifolia"    1.03 
    ##  6 "Hieronyma\nalchorneoides" "Acalypha\ndiversifolia"    1.22 
    ##  7 "Hura\ncrepitans"          "Acalypha\ndiversifolia"    1.34 
    ##  8 "Sapium\nbroadleaf"        "Acalypha\ndiversifolia"    6.40 
    ##  9 "Sapium\nglandulosum"      "Acalypha\ndiversifolia"    1.54 
    ## 10 "Adelia\ntriloba"          "Acalypha\nmacrostachya"    3.16 
    ## # … with 35 more rows

Finalmente, creo el "mapa de calor".

``` r
coldiss(mi_fam_t_chi_d, diag = TRUE)
```

![](medición_asociacion_3_files/figure-markdown_github/unnamed-chunk-7-1.png)

En el mapa de calor **ordenado** (el de la derecha), se identifica al menos un patrón de dependencia entre las especies relacionadas en la diagonal desde *Chrysophyllum cainito* hasta *Trichilia pallida* (cuadros de color rosa centrales). También se observan las especies que no parecen asociarse con otras, situadas en los extremos de la diagonal, y relacionadas con otras por medio de valores pequeños de distancia (cuadros azules), como *Rauvolfia littoralis* y *Pouteria fossicola* y *Cedrela odorata*.

### Modo R para datos binarios (presencia/ausencia)

Arriba usé la distancia de Jaccard para evaluar asociación entre sitios. Dicha métrica también se puede usar para evaluar la distancia entre especies, usando como fuente la matriz de comunidad transpuesta convertida a binaria (presencia/ausencia)

``` r
mi_fam_t_jac <- vegdist(mi_fam_t, "jaccard", binary = TRUE)
mi_fam_t_jac %>% tidy
```

    ## # A tibble: 45 x 3
    ##    item1                      item2                    distance
    ##    <fct>                      <fct>                       <dbl>
    ##  1 "Acalypha\nmacrostachya"   "Acalypha\ndiversifolia"   0.755 
    ##  2 "Adelia\ntriloba"          "Acalypha\ndiversifolia"   0.367 
    ##  3 "Alchornea\ncostaricensis" "Acalypha\ndiversifolia"   0.02  
    ##  4 "Alchornea\nlatifolia"     "Acalypha\ndiversifolia"   0.980 
    ##  5 "Croton\nbillbergianus"    "Acalypha\ndiversifolia"   0.0400
    ##  6 "Hieronyma\nalchorneoides" "Acalypha\ndiversifolia"   0.286 
    ##  7 "Hura\ncrepitans"          "Acalypha\ndiversifolia"   0.367 
    ##  8 "Sapium\nbroadleaf"        "Acalypha\ndiversifolia"   0.980 
    ##  9 "Sapium\nglandulosum"      "Acalypha\ndiversifolia"   0.58  
    ## 10 "Adelia\ntriloba"          "Acalypha\nmacrostachya"   0.697 
    ## # … with 35 more rows

``` r
coldiss(mi_fam_t_jac, diag = TRUE)
```

![](medición_asociacion_3_files/figure-markdown_github/unnamed-chunk-8-1.png)

### Modo R para datos cuantitativos, NO de abundancia de especies (variables ambientales)

En modo R evalúas asociación entre descriptores, es decir, entre variables. La métrica comúnmente usada es el índice de correlación de Pearson. Sin embargo, si los datos no presentan distribución normal, puedes emplear métricas más flexibles, como el índice *rho* de Spearman o **tau** de Kendall.

En este ejemplo, mostraré la correlación entre variables de suelo y la abundancia y riqueza globales y de mi familia asignada. Haré lo propio con variables geomorfológicas. Ya tuviste ocasión de usar estos métodos en el [*script* de análisis exploratorio de datos, sección correlación](aed_5_correlaciones_variables_ambientales.md). En este *script*, añadirás al análisis el índice *rho* de Spearman.

``` r
env_num <- bci_env_grid %>%
  dplyr::select_if(is.numeric) %>%
  dplyr::select(-id, -matches('^U.*')) %>% 
  st_drop_geometry %>% 
  mutate(
    riqueza_mifam = specnumber(mc_ephrb),
    abundancia_mifam = rowSums(mc_ephrb)) %>% 
  rename_all(gsub, pattern = '_pct$', replacement = '') %>% 
  rename_all(gsub, pattern = '_| ', replacement = '\n')
env_num %>% tibble
```

    ## # A tibble: 50 x 1
    ##    .$`heterogeneid… $`geomorf\nllan… $`geomorf\npico` $`geomorf\ninte…
    ##               <dbl>            <dbl>            <dbl>            <dbl>
    ##  1           0.627             10.0              0                0.83
    ##  2           0.394             34.8              0                0.36
    ##  3           0                  0                0                0   
    ##  4           0                  0                0                0.16
    ##  5           0.461              2.58             0                0   
    ##  6           0.0768             0                0.17             3.01
    ##  7           0.381              0                0.53             2.87
    ##  8           0.211              0                0                0   
    ##  9           0                  0                0                0   
    ## 10           0                  1.03             0                0   
    ## # … with 40 more rows, and 29 more variables: $`geomorf\nhombrera` <dbl>,
    ## #   $`geomorf\nespolón/gajo` <dbl>, $`geomorf\nvertiente` <dbl>,
    ## #   $`geomorf\nvaguada` <dbl>, $`geomorf\npiedemonte` <dbl>,
    ## #   $`geomorf\nvalle` <dbl>, $`geomorf\nsima` <dbl>, $Al <dbl>, $B <dbl>,
    ## #   $Ca <dbl>, $Cu <dbl>, $Fe <dbl>, $K <dbl>, $Mg <dbl>, $Mn <dbl>,
    ## #   $P <dbl>, $Zn <dbl>, $N <dbl>, $N.min. <dbl>, $pH <dbl>,
    ## #   $`elevacion\nmedia` <dbl>, $`pendiente\nmedia` <dbl>,
    ## #   $`orientacion\nmedia` <dbl>, $`curvatura\nperfil\nmedia` <dbl>,
    ## #   $`curvatura\ntangencial\nmedia` <dbl>, $`abundancia\nglobal` <dbl>,
    ## #   $`riqueza\nglobal` <int>, $`riqueza\nmifam` <int>,
    ## #   $`abundancia\nmifam` <dbl>

``` r
p_cor_suelo_ar <- env_num %>%
  dplyr::select(matches('^[A-T,Z]|abundancia|riqueza|^pH$', ignore.case = F)) %>%
  ezCorM(r_size_lims = c(4,8), label_size = 3, method = 'pearson')
```

    ## -------------------------------------------------------------------------

    ## You have loaded plyr after dplyr - this is likely to cause problems.
    ## If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
    ## library(plyr); library(dplyr)

    ## -------------------------------------------------------------------------

    ## 
    ## Attaching package: 'plyr'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize

    ## The following object is masked from 'package:purrr':
    ## 
    ##     compact

    ## 
    ## Attaching package: 'reshape2'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     smiths

``` r
p_cor_suelo_ar
```

![](medición_asociacion_3_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
p_cor_suelo_ar_spearman <- env_num %>%
  dplyr::select(matches('^[A-T,Z]|abundancia|riqueza|^pH$', ignore.case = F)) %>%
  ezCorM(r_size_lims = c(4,8), label_size = 3, method = 'spearman')
p_cor_suelo_ar_spearman
```

![](medición_asociacion_3_files/figure-markdown_github/unnamed-chunk-9-2.png)

``` r
png(
  filename = 'matriz_correlacion_suelo_abun_riq_spearman.png',
  width = 1920, height = 1080, res = 125
)
p_cor_suelo_ar_spearman
dev.off() #NO OLVIDAR ESTA IMPORTANTE SENTENCIA
```

    ## png 
    ##   2

``` r
p_cor_geomorf_ar <- env_num %>%
  dplyr::select(-matches('^[A-T,Z]|pH', ignore.case = F)) %>%
  ezCorM(r_size_lims = c(4,8), label_size = 3, method = 'pearson')
p_cor_geomorf_ar
```

![](medición_asociacion_3_files/figure-markdown_github/unnamed-chunk-9-3.png)

``` r
p_cor_geomorf_ar_spearman <- env_num %>%
  dplyr::select(-matches('^[A-T,Z]|pH', ignore.case = F)) %>%
  ezCorM(r_size_lims = c(4,8), label_size = 3, method = 'spearman')
p_cor_geomorf_ar_spearman
```

![](medición_asociacion_3_files/figure-markdown_github/unnamed-chunk-9-4.png)

``` r
png(
  filename = 'matriz_correlacion_geomorf_abun_riq_spearman.png',
  width = 1920, height = 1080, res = 110
)
p_cor_geomorf_ar_spearman
dev.off() #NO OLVIDAR ESTA IMPORTANTE SENTENCIA
```

    ## png 
    ##   2
