# load ----
require(shiny)
library(shinycssloaders)
library(tidyverse)
library(DT)
library(scico)


options(DT.options = list(pageLength = 10, searching = FALSE),
        shiny.useragg = TRUE)
dt_output = function(title, id) {
  fluidRow(column(6, hr(), DTOutput(id)
  ))
}
theme_set(theme_light() +
            theme(panel.grid.major = ggplot2::element_blank(),
                  panel.grid.minor = ggplot2::element_blank(),
                  strip.background = ggplot2::element_rect(fill = NA, colour = NA),
                  strip.text.x = element_text(colour = "black"),
                  strip.text.y = element_text(colour = "black"),
                  panel.border = element_rect(fill = NA),
                  legend.key.size = grid::unit(0.9, "lines"),
                  legend.key = ggplot2::element_rect(colour = NA, fill = NA),
                  legend.background = ggplot2::element_rect(colour = NA, fill = NA)))

# biological functions ----

#' von Bertalanffy growth curve
#'
#' @param age fish age
#' @param linf asymptotic average length
#' @param k growth rate coefficient
#' @param t0 age when the average length is zero
#' @param vonb_sd provides error about estimates

vonb <- function(age, linf = 80.2, k = 0.222, t0 = -1.9, vonb_sd = 0){
  rnorm(1, linf * (1 - exp(-1.0 * k * (age - t0))), vonb_sd)
}

#' Weight at age/length
#' @description log transformation formula of Le Cren (1951) for describing the weight at length relationship
#' @param length fish length
#' @param a intercept
#' @param b slope

waa <- function(length, a, b){
  exp(log(a) + log(length) * b)
}

#' Biological data function
#' @description Takes user inputs and estimates the proportion mature at age, length, and age/length along with associated levels of skip spawning at age.
#'
#' @param age fish age
#' @param lw dataframe of length and weight by age
#' @param b0 maturity at age slope
#' @param b1 maturity at age intercept
#' @param b0_l maturity at length slope
#' @param b1_l maturity at length intercept
#' @param minskip minimum age that skip spawning occurs
#' @param maxskip maximum age that skip spawning occurs
#' @param skip the amount of skip spawning
#' @param shigh the amount of skip spawning at the oldest age
#' @param sims the number of times to simulate each age

biological <- function(age = 50, lw,
                       b0 = 0.84, b1 = 6.60, b0_l = 0.4, b1_l = 65, minskip = 5,
                       maxskip = 25, skip = 5, shigh = 5,
                       sims = 500){

  data.frame(age = minskip:maxskip) |>
    mutate(s = c(skip, rep(NA, length(age) - 2 ), shigh),
           s = approx(age, s, age)$y / 100) -> skipping

  lw %>%
    left_join(skipping) %>%
    replace_na(list(s = 0))  %>%
    mutate(prob_mat = 1 / (1 + exp(-b0 * (age - b1))),
           prob_mat_l = 1 / (1 + exp(-b0_l * (length - b1_l))),
           mature = purrr::pmap_dbl(list(1, 1, prob_mat), rbinom),
           mature_l = purrr::pmap_dbl(list(1, 1, prob_mat_l), rbinom),
           mature_al = purrr::pmap_dbl(list(mature, mature_l), min),
           skip = ifelse(purrr::pmap_dbl(list(1, 1, s), rbinom) == 1, 0, mature),
           skip_l = ifelse(purrr::pmap_dbl(list(1, 1, s), rbinom) == 1, 0, mature_l),
           skip_al = ifelse(purrr::pmap_dbl(list(1, 1, s), rbinom) == 1, 0, mature_al)) %>%
    dplyr::select(-c(sims, s))

}

#' Population estimation
#' @description produces population for 50 years based upon the max age and population dynamics inputs
#'
#' @param age number of ages to examine, user defined
#' @param meanpop mean population size - log scale, user defined
#' @param sigr recruitment variation, user defined
#' @param failure can set a number of years with recruitment failures: cohort population * runif(1,0,0.1)
#' @param type default is NULL, otherwise a cyclical population is generated
#' @param m  natural mortality, user defined
#' @param fish fishing mortality, user defined
#' @param s50 age at 50% selectivity, user defined

population <- function(age = 50, meanpop = 6.5, sigr = 1, failure = NULL,
                       type = NULL, m = 0.098, fish = 0.117, s50 = 3.7, bin = 40){

  slx = 1 / (1 + exp(-log(19) * (1:age - s50)/((s50 + 2) - s50))) # selectivity at age

  # numbers
  if(!is.null(type)){
    strt = sample(1:10, 1)
    n = matrix(c(unlist(replicate(50 + age, c(rlnorm(1, meanpop, sigr),
                                              rlnorm(9, 3.5, 0.7)),
                                  simplify = FALSE))[strt:(49+age+strt)],
                 rep(NA, (50 + age) * (age - 1))),
               nrow = age, byrow = TRUE)

  } else if(is.null(type)){
    n = matrix(c(rlnorm(50 + age, meanpop, sigr),
                 rep(NA, (50 + age)* (age - 1))),
               nrow = age, byrow = TRUE)

    if(!is.null(failure)){
      fail = sample(ncol(n), failure)

      for(i in 1:failure){
        n[1,fail[i]] = n[1,fail[i]] * runif(1, 0, 0.1)
      }

    }
  }

  b0 = n1 = n

  for(i in 2:age){
    for(j in 1:(50 + age)){
      n[i,j] = n[i-1, j] * exp(-m - fish * slx[i])
      b0[i,j] = b0[i-1, j] * exp(-m)
    }
  }

  # binned population
  if(bin < age){

    slx = 1 / (1 + exp(-log(19) * (1:bin - s50)/((s50 + 2) - s50)))

    for(i in 2:(bin-1)){
      for(j in 1:(50 + age)){
        n1[i,j] = n1[i-1, j] * exp(-m - fish * slx[i-1])
      }
    }

    for(j in 1:(50 + age)){
      n1[bin,j] = n1[bin-1, j] * exp(-m - fish * slx[bin-1]) / (1 - exp(-m - fish * slx[bin]))
      for(i in (bin+1):age){
        n1[i,j] = -1
      }
    }
  } else {
    n1 = n
  }

  data.frame(n, id = "n") %>%
    bind_rows(data.frame(b0, id = "b0")) %>%
    bind_rows(data.frame(n1, id = "plus")) %>%
    dplyr::mutate(age = rep(1:age, 3)) %>%
    tidyr::pivot_longer(cols = dplyr::starts_with("X"),
                        names_to = "year") %>%
    dplyr::mutate(year = as.numeric(substring(year,2)),
                  year = year + age) %>%
    filter(value > -1) %>%
    dplyr::filter(year %in% max(age):(max(age) + 49)) %>%
    dplyr::mutate(year = year - max(age) + 1)
}

model_a <- function(bio, n = 100, mat_bin = 40, a_k = 10){

  bio |>
    group_by(age) |>
    sample_n(n) |>
    ungroup() -> dat2

  m1 = glm(skip ~ age, data = dat2, family = binomial)
  m2 = mgcv::gam(skip ~ s(age, k = a_k), data = dat2, gamma = 1.4, family = binomial)

  bio %>%
    group_by(age) %>%
    summarise(true = sum(skip) / n()) %>%
    ungroup() -> dat3

  p1 = predict(m1, dat3, type = "response", se = T)
  p2 = predict(m2, dat3, type = "response", se = T)

  dat3 %>%
    mutate(lo = NA,
           hi = NA,
           glm = p1$fit,
           lo_glm = glm - 1.96 * p1$se.fit,
           hi_glm = glm + 1.96 * p1$se.fit,
           gam = p2$fit,
           lo_gam = gam - 1.96 * p2$se.fit,
           hi_gam = gam + 1.96 * p2$se.fit) -> dat3

  if(mat_bin < max(bio$age)){
    dat2 %>%
      mutate(age = ifelse(age > mat_bin, mat_bin, age)) -> dat4

    m3 = glm(skip ~ age, data = dat4, family = binomial)
    m4 = mgcv::gam(skip ~ s(age, k = a_k), data = dat4, gamma = 1.4, family = binomial)

    p3 = predict(m3, dat3, type = "response", se = T)
    p4 = predict(m4, dat3, type = "response", se = T)

    dat3 %>%
      mutate(glm_bin = p3$fit,
             lo_glm_bin = glm_bin - 1.96 * p3$se.fit,
             hi_glm_bin = glm_bin + 1.96 * p3$se.fit,
             gam_bin = p4$fit,
             lo_gam_bin = gam_bin - 1.96 * p4$se.fit,
             hi_gam_bin = gam_bin + 1.96 * p4$se.fit) %>%
      pivot_longer(c(true, glm, gam, glm_bin, gam_bin),
                   names_to = "Model", values_to = "fit") %>%
      mutate(lo = case_when(Model == "glm" ~ lo_glm,
                            Model == "gam" ~ lo_gam,
                            Model == "glm_bin" ~ lo_glm_bin,
                            Model == "gam_bin" ~ lo_gam_bin),
             hi = case_when(Model == "glm" ~ hi_glm,
                            Model == "gam" ~ hi_gam,
                            Model == "glm_bin" ~ hi_glm_bin,
                            Model == "gam_bin" ~ hi_gam_bin)) %>%
      dplyr::select(-c(ends_with("glm"), ends_with("gam"), ends_with("bin")))

  } else {

    dat3 %>%
      pivot_longer(c(true, glm, gam),
                   names_to = "Model", values_to = "fit") %>%
      mutate(lo = case_when(Model == "glm" ~ lo_glm,
                            Model == "gam" ~ lo_gam),
             hi = case_when(Model == "glm" ~ hi_glm,
                            Model == "gam" ~ hi_gam)) %>%
      dplyr::select(-c(ends_with("glm"), ends_with("gam"), ends_with("bin")))
  }
}

model_l <- function(bio, n = 100, mat_bin = 50, l_k = 10){

  bio |>
    group_by(age) |>
    sample_n(n) |>
    ungroup() -> dat2

  m1 = glm(skip_l ~ length, data = dat2, family = binomial)
  m2 = mgcv::gam(skip_l ~ s(length, k = l_k), data = dat2, gamma = 1.4, family = binomial)

  bio %>%
    group_by(length) %>%
    summarise(true = sum(skip_l) / n()) %>%
    ungroup() -> dat3

  p1 = predict(m1, dat3, type = "response", se = T)
  p2 = predict(m2, dat3, type = "response", se = T)

  dat3 %>%
    mutate(lo = NA,
           hi = NA,
           glm = p1$fit,
           lo_glm = glm - 1.96 * p1$se.fit,
           hi_glm = glm + 1.96 * p1$se.fit,
           gam = p2$fit,
           lo_gam = gam - 1.96 * p2$se.fit,
           hi_gam = gam + 1.96 * p2$se.fit) -> dat3

  if(mat_bin < max(bio$age)){
    dat2 %>%
      mutate(age = ifelse(age > mat_bin, mat_bin, age)) -> dat4

    m3 = glm(skip_l ~ length, data = dat4, family = binomial)
    m4 = mgcv::gam(skip_l ~ s(length, k = l_k), data = dat4, gamma = 1.4, family = binomial)

    p3 = predict(m3, dat3, type = "response", se = T)
    p4 = predict(m4, dat3, type = "response", se = T)

    dat3 %>%
      mutate(glm_bin = p3$fit,
             lo_glm_bin = glm_bin - 1.96 * p3$se.fit,
             hi_glm_bin = glm_bin + 1.96 * p3$se.fit,
             gam_bin = p4$fit,
             lo_gam_bin = gam_bin - 1.96 * p4$se.fit,
             hi_gam_bin = gam_bin + 1.96 * p4$se.fit) %>%
      pivot_longer(c(true, glm, gam, glm_bin, gam_bin),
                   names_to = "Model", values_to = "fit") %>%
      mutate(lo = case_when(Model == "glm" ~ lo_glm,
                            Model == "gam" ~ lo_gam,
                            Model == "glm_bin" ~ lo_glm_bin,
                            Model == "gam_bin" ~ lo_gam_bin),
             hi = case_when(Model == "glm" ~ hi_glm,
                            Model == "gam" ~ hi_gam,
                            Model == "glm_bin" ~ hi_glm_bin,
                            Model == "gam_bin" ~ hi_gam_bin)) %>%
      dplyr::select(-c(ends_with("glm"), ends_with("gam"), ends_with("bin")))

  } else {

    dat3 %>%
      pivot_longer(c(true, glm, gam),
                   names_to = "Model", values_to = "fit") %>%
      mutate(lo = case_when(Model == "glm" ~ lo_glm,
                            Model == "gam" ~ lo_gam),
             hi = case_when(Model == "glm" ~ hi_glm,
                            Model == "gam" ~ hi_gam)) %>%
      dplyr::select(-c(ends_with("glm"), ends_with("gam"), ends_with("bin")))
  }
}

model_al <- function(bio, n = 100, mat_bin = 40, aal_k = 10, al_k = 10){

  bio |>
    group_by(age) |>
    sample_n(n) |>
    ungroup() -> dat2

  m1 = glm(skip_al ~ age * length, data = dat2, family = binomial)
  m2 = mgcv::gam(skip_al ~ s(age, k = aal_k) + s(length, k = al_k), data = dat2, gamma = 1.4, family = binomial)

  bio %>%
    group_by(age, length) %>%
    summarise(true = sum(skip_al) / n()) %>%
    ungroup() -> dat3

  p1 = predict(m1, dat3, type = "response", se = T)
  p2 = predict(m2, dat3, type = "response", se = T)

  dat3 %>%
    mutate(lo = NA,
           hi = NA,
           glm = p1$fit,
           lo_glm = glm - 1.96 * p1$se.fit,
           hi_glm = glm + 1.96 * p1$se.fit,
           gam = p2$fit,
           lo_gam = gam - 1.96 * p2$se.fit,
           hi_gam = gam + 1.96 * p2$se.fit) -> dat3

  if(mat_bin < max(bio$age)){
    dat2 %>%
      mutate(age = ifelse(age > mat_bin, mat_bin, age)) -> dat4

    m3 = glm(skip_al ~ age * length, data = dat4, family = binomial)
    m4 = mgcv::gam(skip_al ~ s(age, k = aal_k) + s(length, k = al_k), data = dat4, gamma = 1.4, family = binomial)

    p3 = predict(m3, dat3, type = "response", se = T)
    p4 = predict(m4, dat3, type = "response", se = T)

    dat3 %>%
      mutate(glm_bin = p3$fit,
             lo_glm_bin = glm_bin - 1.96 * p3$se.fit,
             hi_glm_bin = glm_bin + 1.96 * p3$se.fit,
             gam_bin = p4$fit,
             lo_gam_bin = gam_bin - 1.96 * p4$se.fit,
             hi_gam_bin = gam_bin + 1.96 * p4$se.fit) %>%
      pivot_longer(c(true, glm, gam, glm_bin, gam_bin),
                   names_to = "Model", values_to = "fit") %>%
      mutate(lo = case_when(Model == "glm" ~ lo_glm,
                            Model == "gam" ~ lo_gam,
                            Model == "glm_bin" ~ lo_glm_bin,
                            Model == "gam_bin" ~ lo_gam_bin),
             hi = case_when(Model == "glm" ~ hi_glm,
                            Model == "gam" ~ hi_gam,
                            Model == "glm_bin" ~ hi_glm_bin,
                            Model == "gam_bin" ~ hi_gam_bin)) %>%
      dplyr::select(-c(ends_with("glm"), ends_with("gam"), ends_with("bin")))

  } else {

    dat3 %>%
      pivot_longer(c(true, glm, gam),
                   names_to = "Model", values_to = "fit") %>%
      mutate(lo = case_when(Model == "glm" ~ lo_glm,
                            Model == "gam" ~ lo_gam),
             hi = case_when(Model == "glm" ~ hi_glm,
                            Model == "gam" ~ hi_gam)) %>%
      dplyr::select(-c(ends_with("glm"), ends_with("gam"), ends_with("bin")))
  }
}
ssb <- function(mods, pop, lw){

  mods %>%
    left_join(lw) %>%
    left_join(pop) %>%
    mutate(ssb = value/2 * weight * fit) %>%
    group_by(year, Model, id)  %>%
    summarise(ssb = sum(ssb, na.rm = T) / 1000) %>%
    ungroup()
}

# plotting ----
plot_al <- function(data){
  data |>
    ggplot(aes(age, length)) +
    geom_point() +
    expand_limits(y = 0) +
    xlab("Age") +
    ylab("Length")
}
plot_wl <- function(data){
  data |>
    ggplot(aes(length, weight)) +
    geom_line() +
    expand_limits(y = 0) +
    xlab("Length") +
    ylab("Weight (kg)")
}
plot_skip <- function(data){
  data |>
    ggplot(aes(age, s)) +
    geom_point() +
    geom_line() +
    expand_limits(y = 0) +
    ylab("% skip spawning")
}
plot_pop <- function(pop){
  pop %>%
    dplyr::filter(id == "plus") %>%
    ggplot2::ggplot(ggplot2::aes(year, age, size = value)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::scale_size_area() +
    xlab("Year") +
    ylab("Age")
}
plot_slx <- function(s50, age){

  data.frame(age = 1:age) %>%
    mutate(slx = 1 / (1 + exp(-log(19) * (age - s50)/((s50 + 2) - s50)))) %>%
    ggplot(aes(age, slx)) +
    geom_line() +
    expand_limits(x = 0, y = 0) +
    xlab("Age") +
    ylab("Selectivity")
}
plot_ssb <- function(sb, age=50, bin=50){

  sb %>%
    filter(id == "b0", Model == "true") %>%
    group_by(year) %>%
    summarise(b40 = ssb *  0.4) -> blim

  sb %>%
    filter(id == "n", Model == "true") -> true

  if(age>bin){
    sb %>%
      filter(Model != "true", id!= "b0" ,
             (Model %in% c("gam", "glm") & id == "n")|
               (Model %in% c("gam", "glm") & id == "plus")) %>%
      mutate(group = case_when(id == "n" ~ "full",
                               T ~ "plus"),
             Model = gsub("\\_.*", "", Model)) %>%
      ggplot(aes(year, ssb)) +
      geom_line(data = true, lty = 3) +
      geom_line(aes(color = Model)) +
      # geom_line(data = blim, aes(y = b40), lty = 3) +
      scale_y_continuous(label = scales::comma) +
      facet_wrap(~group) +
      theme(legend.position = "bottom") +
      scale_color_scico_d(palette = "roma", begin = 0.2)  +
      ylab("SSB (KT)")
  } else if("gam_bin" %in% sb$Model){
    sb %>%
      filter(Model != "true", id!= "b0" ,
             (Model %in% c("gam", "glm") & id == "n")|
               (Model %in% c("gam_bin", "glm_bin") & id == "plus")) %>%
      mutate(group = case_when(id == "n" ~ "full",
                               T ~ "plus"),
             Model = gsub("\\_.*", "", Model)) %>%
      ggplot(aes(year, ssb)) +
      geom_line(data = true, lty = 3) +
      geom_line(aes(color = Model)) +
      # geom_line(data = blim, aes(y = b40), lty = 3) +
      scale_y_continuous(label = scales::comma) +
      facet_wrap(~group) +
      theme(legend.position = "bottom") +
      scale_color_scico_d(palette = "roma", begin = 0.2)  +
      ylab("SSB (KT)")
  } else {
    sb %>%
      filter(Model != "true", id!= "b0" ,
             (Model %in% c("gam", "glm") & id == "n")|
               (Model %in% c("gam_bin", "glm_bin") & id == "plus")) %>%
      mutate(group = case_when(id == "n" ~ "full",
                               T ~ "plus"),
             Model = gsub("\\_.*", "", Model)) %>%
      ggplot(aes(year, ssb)) +
      geom_line(data = true, lty = 3) +
      geom_line(aes(color = Model)) +
      # geom_line(data = blim, aes(y = b40), lty = 3) +
      scale_y_continuous(label = scales::comma) +
      facet_wrap(~group) +
      theme(legend.position = "bottom") +
      scale_color_scico_d(palette = "roma", begin = 0.2)  +
      ylab("SSB (KT)")
  }
}
plot_pd <- function(sb, age, bin){


  sb  %>%
    filter(id != "b0", !(Model=="true" & id =="plus")) %>%
    mutate(group = case_when(id == "n" ~ "full",
                             T ~ "plus"),
           Model = gsub("\\_.*", "", Model)) %>%
    group_by(year) %>%
    mutate(true = ifelse(Model=="true", ssb, NA),
           true = mean(true, na.rm = T),
           pd = (true - ssb) / ((true + ssb) / 2) * 100) %>%
    filter(Model != "true")  %>%
    ggplot(aes(year, pd, fill = Model)) +
    geom_area(alpha = 0.4, position = "identity") +
    theme(legend.position = c(0.9, 0.85)) +
    scale_fill_manual(values = c("#B48A2C", "#8CDED9")) +
    ylab("% difference from true") +
    geom_hline(yintercept = 0, lty = 3) -> p

  if(age>bin | "gam_bin" %in% sb$Model){
    p + facet_wrap(~group)
  } else {
     p
  }

}
plot_mature <- function(data, type = "age"){

  if(type=="age" & "glm_bin" %in% data$Model){
    data %>%
      mutate(Model = case_when(Model == "glm_bin" ~ "glm_plus",
                               Model == "gam_bin" ~ "gam_plus",
                               TRUE ~ Model)) %>%
      ggplot(aes(age, fit)) +
      geom_line(aes(linetype = Model, color = Model)) +
      geom_point(aes(color = Model, shape = Model)) +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = Model), alpha = 0.1, show.legend=FALSE) +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = Model), alpha = 0.1, show.legend=FALSE) +
      scale_linetype_manual("Model", values = c(1, 1, 1, 1, 0)) +
      scale_shape_manual("Model", values = c(NA, NA, NA, NA, 19)) +
      theme(legend.position = c(0.8, 0.2)) +
      scale_color_manual("Model", values = c("#7E1900", "#C1A53A", "#D1ECC9", "#479BC5", "#1A3399")) +
      scale_fill_manual("Model", values = c("#7E1900", "#C1A53A", "#D1ECC9", "#479BC5", "#1A3399")) +
      xlab("Age") +
      ylab("Maturity") +
      expand_limits(x = 0, y = 0)
  } else if(type=="age" & !("glm_bin" %in% data$Model)){
    data %>%
      ggplot(aes(age, fit)) +
    geom_line(aes(linetype = Model, color = Model)) +
      geom_point(aes(color = Model, shape = Model)) +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = Model), alpha = 0.1, show.legend=FALSE) +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = Model), alpha = 0.1, show.legend=FALSE) +
      scale_linetype_manual("Model", values = c(1, 1, 0)) +
      scale_shape_manual("Model", values = c(NA, NA, 19)) +
      theme(legend.position = c(0.8, 0.2)) +
      scale_color_manual("Model", values = c("#7E1900",  "#D1ECC9", "#479BC5", "#1A3399")) +
      scale_fill_manual("Model", values = c("#7E1900",  "#D1ECC9", "#479BC5", "#1A3399")) +
      xlab("Age") +
      ylab("Maturity") +
      expand_limits(x = 0, y = 0)
  } else if(type=="length" & "glm_bin" %in% data$Model){
    data %>%
      mutate(Model = case_when(Model == "glm_bin" ~ "glm_plus",
                               Model == "gam_bin" ~ "gam_plus",
                               TRUE ~ Model)) %>%
      ggplot(aes(length, fit)) +
    geom_line(aes(linetype = Model, color = Model)) +
      geom_point(aes(color = Model, shape = Model)) +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = Model), alpha = 0.1, show.legend=FALSE) +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = Model), alpha = 0.1, show.legend=FALSE) +
      scale_linetype_manual("Model", values = c(1, 1, 1, 1, 0)) +
      scale_shape_manual("Model", values = c(NA, NA, NA, NA, 19)) +
      theme(legend.position = c(0.8, 0.2)) +
      scale_color_manual("Model", values = c("#7E1900", "#C1A53A", "#D1ECC9", "#479BC5", "#1A3399")) +
      scale_fill_manual("Model", values = c("#7E1900", "#C1A53A", "#D1ECC9", "#479BC5", "#1A3399")) +
      xlab("Age") +
      ylab("Maturity") +
      expand_limits(x = 0, y = 0)
  } else if(type=="length" & !("glm_bin" %in% data$Model)){
    data %>%
      ggplot(aes(length, fit)) +
    geom_line(aes(linetype = Model, color = Model)) +
      geom_point(aes(color = Model, shape = Model)) +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = Model), alpha = 0.1, show.legend=FALSE) +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = Model), alpha = 0.1, show.legend=FALSE) +
      scale_linetype_manual("Model", values = c(1, 1, 0)) +
      scale_shape_manual("Model", values = c(NA, NA, 19)) +
      theme(legend.position = c(0.8, 0.2)) +
      scale_color_manual("Model", values = c("#7E1900",  "#D1ECC9", "#479BC5", "#1A3399")) +
      scale_fill_manual("Model", values = c("#7E1900",  "#D1ECC9", "#479BC5", "#1A3399")) +
      xlab("Age") +
      ylab("Maturity") +
      expand_limits(x = 0, y = 0)
  } else if(type=="al"& "glm_bin" %in% data$Model){
      data %>%
        pivot_wider(names_from = Model, values_from = fit) %>%
        group_by(age, length) %>%
        summarise(true = sum(true, na.rm = T),
                  glm = sum(true, na.rm = T) - sum(glm, na.rm = T),
                  gam = sum(true, na.rm = T) - sum(gam, na.rm = T),
                  glm_plus = sum(true, na.rm = T) - sum(glm_bin, na.rm = T),
                  gam_plus = sum(true, na.rm = T) - sum(gam_bin, na.rm = T)) %>%
        pivot_longer(c(glm, glm_plus, gam, gam_plus), names_to = "Model", values_to = "Difference") %>%
        mutate(cols = ifelse(Difference <=0, "1", "0")) %>%
        ungroup() %>%
        ggplot(aes(age, length)) +
        geom_point(aes(size = Difference, fill = cols), pch = 21, alpha = 0.7) +
        scale_size_area(guide = "none") +
        facet_wrap(~Model) +
        scale_fill_manual("Residuals",
                          values = c("#8CDED9", "#7E1900"),
                          labels = c("positive", "negative")) +
        xlab("Age") +
        ylab("Length (cm)\n") +
        expand_limits(x = 0, y = 0) +
        theme(legend.position = c(0.85, 0.2))
   } else {
            data %>%
              pivot_wider(names_from = Model, values_from = fit) %>%
              group_by(age, length) %>%
              summarise(true = sum(true, na.rm = T),
                        glm = sum(true, na.rm = T) - sum(glm, na.rm = T),
                        gam = sum(true, na.rm = T) - sum(gam, na.rm = T)) %>%
              pivot_longer(c(glm, gam), names_to = "Model", values_to = "Difference") %>%
              mutate(cols = ifelse(Difference <=0, "1", "0")) %>%
              ungroup() %>%
              ggplot(aes(age, length)) +
              geom_point(aes(size = Difference, fill = cols), pch = 21, alpha = 0.7) +
              scale_size_area(guide = "none") +
              facet_wrap(~Model) +
              scale_fill_manual("Residuals",
                                values = c("#8CDED9", "#7E1900"),
                                labels = c("positive", "negative")) +
              xlab("Age") +
              ylab("Length (cm)\n") +
              expand_limits(x = 0, y = 0) +
              theme(legend.position = c(0.85, 0.2))
          }

}
