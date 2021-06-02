get_ss = function(model) {
  require(car) 
  
  a = Anova(model)
  out = data.frame(
    PREDICTOR = rownames(a),
    SUMSQ = a$`Sum Sq`,
    pVAR = round(a$`Sum Sq` / sum(a$`Sum Sq`) *100, 2),
    F_val = round(a$`F value`, 4),
    p_val = signif(a$`Pr(>F)`, 3)
  )
  return(out)
}

