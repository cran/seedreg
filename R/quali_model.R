#' Analysis: generalized linear models for factor qualitative
#'
#' Performs the deviance analysis for the generalized linear model using binomial or quasibinomial family. The function also returns multiple comparison test with tukey adjustment
#' @param trat Numerical or complex vector with treatments
#' @param resp Numerical vector containing the response of the experiment.
#' @param n Number of seeds per repetition
#' @param family a description of the error distribution and link function to be used in the model. For glm this can be a character string naming a family function, a family function or the result of a call to a family function.
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @return The function returns analysis by glm (binomial or quasibinomial family), post-hoc and column graph
#' @export
#' @importFrom emmeans emmeans
#' @importFrom emmeans regrid
#' @importFrom multcomp cld
#' @examples
#' library(seedreg)
#' data("aristolochia")
#' attach(aristolochia)
#' quali_model(trat, resp, n=25, family="quasibinomial")

quali_model=function(trat,
                     resp,
                     n=50,
                     family="binomial",
                     ylab="Germination (%)",
                     xlab=expression("Temperature ("^"o"*"C)"),
                     theme=theme_classic()){
  requireNamespace("ggplot2")
  requireNamespace("drc")
  requireNamespace("multcomp")
  requireNamespace("emmeans")

  nseeds=resp*n/100
  trat=as.factor(trat)
  mod=glm(cbind(nseeds,n-nseeds)~trat,
          family = family)
  hnp::hnp(mod, print.on=T)
  anova=car::Anova(mod)
  aaa=cld(regrid(emmeans(mod,~trat)),Letters = letters,
          reversed = T,sort = F,adjust="tukey")
  graph=ggplot(aaa,aes(y=aaa$prob*100,x=as.factor(trat)))+
  geom_col(color="black",fill="gray")+
  theme+
  geom_errorbar(aes(ymin=aaa$asymp.LCL*100,
                    ymax=aaa$asymp.UCL*100),width=0.2)+
  geom_label(aes(y=aaa$asymp.UCL*100+5,
                 label=paste(round(aaa$prob*100,1),aaa$.group)))+
  theme(axis.text = element_text(size=12,color="black"))+
  labs(x=xlab, y=ylab)
  cat("\n=======================================================\n")
  print(anova)
  cat("\n=======================================================\n")
  print(aaa)

  cat("\n=======================================================\n")
  list(graph)[[1]]
  }
