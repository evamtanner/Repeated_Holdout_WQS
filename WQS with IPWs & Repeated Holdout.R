######################################################################################################
# WQS with Inverse Probability Weighting & Repeated Holdout Validation
## Edited gwqs Function to Incorporate Regression Weights
## w_name = name of IPW variable
######################################################################################################

gwqs.ipw = 
  
  function (formula, mix_name, data, q = 4, validation = 0.6, 
            valid_var = NULL, b = 100, b1_pos = TRUE, b1_constr = FALSE, 
            family = "gaussian", w_name, seed = NULL, wqs2 = FALSE, plots = FALSE, 
            tables = FALSE) 
  {
    gWQS:::.check.function(formula, mix_name, data, q, validation, 
                           valid_var, b, b1_pos, family, seed, wqs2, plots, tables)
    y_name = all.vars(formula)[1]
    covrts = as.matrix(model.matrix(formula, data)[, -1])
    covrts = covrts[match(rownames(data), rownames(covrts)), 
                    ]
    covrts = as.data.frame(covrts)
    if (dim(covrts)[2] == 1) 
      names(covrts) = all.vars(formula)[-1]
    cov_name = names(covrts)
    data = as.data.frame(data)
    if (is.null(valid_var)) 
      data_f = as.data.frame(suppressWarnings(cbind(data[, 
                                                         c(y_name, mix_name, w_name), drop = FALSE], covrts)))
    else data_f = as.data.frame(cbind(data[, c(y_name, mix_name, w_name), 
                                           drop = FALSE], covrts, data[, valid_var, drop = FALSE]))
    data_f = na.omit(data_f)
    suppressWarnings(RNGversion("3.5.0"))
    set.seed(seed)
    if (is.null(q)) {
      q_name = mix_name
    }
    else {
      data_f = gWQS:::quantile_f(data_f, mix_name, q)
      q_name = paste(mix_name, "q", sep = "_")
    }
    if (is.null(valid_var)) {
      splt = gWQS:::split_f(data_f, validation, seed)
      data_t = splt$data_t
      data_v = splt$data_v
    }
    else {
      unique_valid_var = unique(unlist(data_f[, valid_var, 
                                              drop = FALSE]))
      if (identical(unique_valid_var[order(unique_valid_var)], 
                    c(0, 1))) {
        data_t = data_f[data_f[, valid_var] == 0, ]
        data_v = data_f[data_f[, valid_var] == 1, ]
      }
      else stop("valid_var values must be 0 and 1")
    }
    par_model = gWQS:::par.modl.est(data_t, y_name, q_name, cov_name, 
                                    b, b1_pos, b1_constr, family, seed)
    wght_matrix = par_model$wght_matrix
    b1 = par_model$b1
    conv = par_model$conv
    p_val = par_model$p_val
    index_b = par_model$index_b
    wb1pm <- as.data.frame(cbind(wght_matrix, b1, p_val))
    names(wb1pm) = c(mix_name, "b1", "p_val")
    if (b1_pos) {
      mean_weight = colMeans(wb1pm[wb1pm$b1 > 0 & conv != 
                                     2, mix_name, drop = FALSE])
      if (dim(wb1pm[wb1pm$b1 > 0 & conv != 2, mix_name, drop = FALSE])[1] == 
          0) 
        stop("There are no positive b1 in the bootstrapped models")
    }
    else {
      mean_weight = colMeans(wb1pm[wb1pm$b1 < 0 & conv != 
                                     2, mix_name, drop = FALSE])
      if (dim(wb1pm[wb1pm$b1 < 0 & conv != 2, mix_name, drop = FALSE])[1] == 
          0) 
        stop("There are no negative b1 in the bootstrapped models")
    }
    wqs_model = gWQS:::model.fit(data_v[, q_name, drop = FALSE], data_v[, 
                                                                        y_name, drop = FALSE], mean_weight, family, data_v[, 
                                                                                                                           cov_name, drop = FALSE], wqs2)
    if (dim(covrts)[2] == 0) 
      y_adj = data_v[, y_name, drop = FALSE]
    else {
      y = as.matrix(data_v[, y_name, drop = FALSE])
      x = as.matrix(data_v[, cov_name, drop = FALSE])
      w = as.matrix(data_v[, w_name, drop = FALSE])
      if (family == "gaussian") {
        fit = glm(y ~ x, family = gaussian(link = "identity"), weights = w)
        y_adj = mean(as.matrix(data_v[, y_name, drop = FALSE])) + 
          fit$residuals
      }
      else if (family == "binomial") {
        fit = glm(y ~ x, family = binomial(link = "logit"), weights = w)
        y_adj = resid(fit, type = "pearson")
      }
    }
    data_plot = data.frame(mix_name, mean_weight)
    data_plot = data_plot[order(data_plot$mean_weight), ]
    pos = match(data_plot$mix_name, sort(mix_name))
    data_plot$mix_name = factor(data_plot$mix_name, levels(data_plot$mix_name)[pos])
    y_adj_wqs_df = as.data.frame(cbind(y_adj, wqs_model$wqs))
    names(y_adj_wqs_df) = c("y_adj", "wqs")
    if (plots == TRUE) 
      gWQS:::plots(data_plot, y_adj_wqs_df, q, mix_name, mean_weight)
    data_plot = data_plot[order(data_plot$mean_weight, decreasing = TRUE), 
                          ]
    y_adj = as.numeric(unlist(y_adj))
    wqs_index = as.numeric(unlist(wqs_model$wqs))
    if (tables == TRUE) 
      gWQS:::tables(data_plot, wqs_model$m_f, wqs_model$m_f2, wqs_model$aov)
    results = list(wqs_model$m_f, conv, wb1pm, y_adj, wqs_index, 
                   index_b, data_t, data_v, data_plot, wqs_model$m_f2, 
                   wqs_model$aov)
    names(results) = c("fit", "conv", "wb1pm", "y_adj", "wqs", 
                       "index_b", "data_t", "data_v", "final_weights", "fit_2", 
                       "aov")
    return(results)
  }


######################################################################################################
# Example
######################################################################################################

X26 = c('MEP','MBP','MBzP','DEHP','DINP','MOiNCH','MHiDP','MCiNP',
        'DPP','TCP','BPS','BPA','BPF','Triclosan','PBA','OHPH',
        'PCB','HCB','Nonachlor','DDT_DDE','PFNA','PFDA','PFUnDA','PFHxS','PFOA','PFOS')

ipw = 'w_nomiss'

set.seed(213684854)

IQ_rep = lapply(1:100,function(i){
  gwqs.ipw(IQ ~ mom_edu + mom_IQ + mom_smoke_3 + mom_weight + mom_age + parity_bin + Male,
           mix_name = X26, data = data, 
           q = 10, validation = 0.6, b = 100, b1_pos = F, b1_constr = T, 
           family = "gaussian", w_name = ipw, wqs2 = F, plots = F, tables = F)
})
