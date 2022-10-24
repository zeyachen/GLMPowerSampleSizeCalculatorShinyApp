# Load R packages
library(shiny)
library(shinythemes)
library(ggplot2)

# Define UI
ui <- fluidPage(theme = shinytheme("cosmo"),
                navbarPage(
                  # theme = "cerulean", 
                  "Power Calculator for GLM",
                  tabPanel("Calculator",
                           sidebarPanel(
                             tags$h3("Input:"),
                             radioButtons("powersample", "What would you like to compute?", choices = list("Power"=TRUE, "Sample Size"= FALSE), selected = TRUE),
                             radioButtons("scorelrt", "What test you like to use?", choices = list("Score"=TRUE, "LRT"= FALSE), selected = TRUE),
                             numericInput("zn", "Number of covariates(z)", value = 1,min = 1, max = 2),
                             numericInput("xn", "Number of covariates(x)", value = 1,min = 0, max = 2),
                             radioButtons("wgttf", "Equal weights among different categories?", choices = list("YES"=TRUE, "NO"= FALSE), selected = TRUE),
                             conditionalPanel(condition = "input.wgttf == 'FALSE'",textInput("weight","give a list of weights",value = "0.25,0.25,0.25,0.25")),
                             numericInput("psi1", "Type in value of psi1(z1)", value = log(2)),
                             conditionalPanel(condition = "input.xn > '0'",numericInput("lambda1", "Type in value of lambda1(x1)", value = log(2))),
                             conditionalPanel(condition = "input.zn > '1'",numericInput("psi2", "Type in value of psi2(z2)", value = log(2))),
                             conditionalPanel(condition = "input.xn > '1'",numericInput("lambda2", "Type in value of lambda2(x2)", value = log(2))),
                             numericInput("size", "Type in value of alpha significance", value = 0.05, min = 0, max = 1),
                             numericInput("mu_hat", "Type in value of mu hat", value = 0.15),
                             conditionalPanel(condition = "input.powersample == 'FALSE'",numericInput("nominal_pow", "Type in value of nominal power", value = 0.8, min = 0, max = 1)),
                             conditionalPanel(condition = "input.powersample == 'TRUE'",numericInput("n", "Type in value of sample size", value = 537)),
                             selectInput("family", "Choose a family", choices = c("binomial","gaussian","Gamma","inverse.gaussian","poisson")),
                             conditionalPanel(condition = "input.plot == 'TRUE'",radioButtons("delim", "Use Default Bound", choices = list("Yes"=TRUE, "No"= FALSE), selected = TRUE)),
                             conditionalPanel(condition = "input.delim == 'FALSE'",numericInput("Upper", "Type in Upper Bound", value = NA)),
                             conditionalPanel(condition = "input.delim == 'FALSE'",numericInput("Lower", "Type in Lower Bound", value = NA)),
                             conditionalPanel(condition = "input.delim == 'FALSE'",numericInput("INTV", "Type in Number of Estimates", value = NA, min = 2, max = 15)),
                             numericInput("diffpsi", "Type in difference between psi1", value = 0.1)
                           ), # sidebarPanel
                           mainPanel(
                             h1("Calculator Output"),
                             
                             h4("Design Matrix"),
                             uiOutput('matrix'),
                             
                             conditionalPanel(condition = "input.powersample == 'TRUE'",{h4("Power")}),
                             conditionalPanel(condition = "input.powersample == 'TRUE'",{textOutput("txtout")}),
                             
                             conditionalPanel(condition = "input.powersample == 'FALSE'",{h4("Sample Size")}),
                             conditionalPanel(condition = "input.powersample == 'FALSE'",{textOutput("txtout1")}),
                             
                             h4("Plot"),
                             conditionalPanel(condition = "input.powersample == 'TRUE'",{plotOutput("PPlot")}),
                             conditionalPanel(condition = "input.powersample == 'FALSE'",{plotOutput("SPlot")})
                            
                             
                           ) # mainPanel
                           
                  ), # Calculater, tabPanel
                  tabPanel("README",
                           withMathJax(),
                           h2('$$g(E[Y]) = g(\\mu) = \\eta = Z^{T}\\psi + X^{T}\\lambda$$'),
                           h4(HTML("<ul><li>\\(\\mu\\) is the mean of Y, \\(\\eta\\) is the linear predictor,
                                    </li><li>where Z is a vector of p covariates, X is a vector of q covariates, </li>
                                    <li>and \\(\\psi\\)  and \\(\\lambda\\) denote their respective regression coefficients.</li></ul>")),
                           h4('Our calculater returns power or sample size based on score test with \\(H_0:\\psi = \\psi_0\\). \\(\\lambda\\) is/are treated as nuisance parameter'),
                           h4('Here are the steps to use the Calculator.'),
                           h4(HTML("<ul>
                                    <li> 1. Start by selecting power or sample size you would like to compute based on Score or Likelihood Ratio test. </li>
                                    <li> 2. Z is/are the variable/s of interest. Minimum 1 z and maximum 2 each. Input numbers of Z and X, ie \\(z_{1}\\), \\(z_{2}\\)... </li>
                                    <li> 3. Set the weights of design matrix. If choose unequal, input weights by comma separated form and make sure they adds up to 1.</li>
                                    <li> 4. Input the coefficient from your glm with \\(\\psi_i\\) correponds to \\(z_{i}\\), \\(\\lambda_i\\) correponds to \\(x_i\\).</li>
                                    <li> 5. Input \\(\\hat{\\psi_i}\\) from your glm, desired \\(\\alpha_{0}\\), and sample size or nominal power.</li>
                                    <li> 6. Select corresponding family.</li>
                                    </ul>")),
                           h4('Here are instructions for plotting.'),
                           h4(HTML("<ul>
                                    <li> Our plot is done by linear interpolation of several estimates. Each polyline corresponds to a \\(\\psi_1\\) value.  </li>
                                    <li> You could use default setup to decide the number of estimates, default has 5 etimates with 3rd one being the input power or sample size and  1st, 2nd, 4th, 5th being the 0.8, 0.9, 1.1, 1.2 times of the input.  </li>
                                    <li> Or you can set it up mannually by specifying the upper bound, lower bound, and number of estimates. Calculator will generate a list of inputs evenly spread out in the range. Number of estimates is capped at 15. </li>
                                    <li> \\(\\delta\\):Difference of \\(\\psi_1\\) value is used to produce multiple polyline, \\(\\psi_1\\) values are input\\(\\psi_1\\) - 2\\(\\delta\\), input\\(\\psi_1\\) - \\(\\delta\\), input\\(\\psi_1\\), input\\(\\psi_1\\) + \\(\\delta\\), input\\(\\psi_1\\) + 2\\(\\delta\\).</li>
                                    </ul>"))
                           
                  ), # README, tabPanel  
                  tabPanel("Contact",
                           h4('Questions and requests should be sent to Zeya Chen (zeya.chen@mail.utoronto.ca) or Osvaldo Espin-Garcia (osvaldo.espin-garcia@uhnresearch.ca)'),
                           h4('Author sincerely appreicates help and codes from supervisor Dr. Espin-Garcia.')
                  ) # Manual, tabPanel         
                ) # navbarPage
                
              
) # fluidPage


# Define server function  
server <- function(input, output) {

  source("powerJTC_v1.R")
  #dt <- data.frame(z=rep(c(0,1),each=2),x=rep(c(0,1),times=2))
  #wgts <- reactive(ifelse(input$wgttf,rep(1/2^(input$xn + 1),2^(input$xn + 1)),c(strsplit(input$weight,","))))
  #res <- reactive(powerglm(Xvar=as.formula(ifelse(input$xn == 0,"~",paste0("~",paste0("x",1:input$xn,collapse="+")))), Zvar=as.formula(ifelse(input$zn == 0,"~",paste0("~",paste0("z",1:input$zn,collapse="+")))), lambda=input$lambda, psi=input$psi, psi0=0, mu_hat=input$mu_hat, family=input$family, dt, wgts=wgts(), N=ifelse(input$powersample,input$n,1), signif_alpha=input$size))

  
  powerglm <- function(Xvar, Zvar, lambda, psi, psi0, mu_hat, family, data, wgts, N, signif_alpha, disp=1, sorl){
    if (sorl)
    {power.glm(Xvar, Zvar, lambda, psi, psi0, mu_hat, family, data, wgts, N, signif_alpha, disp=1)}
    else 
    {power.lrt.glm(Xvar, Zvar, lambda, psi, psi0, mu_hat, family, data, wgts, N, signif_alpha, disp=1)}  
  }
    
  

    
    
  dt <- reactive({
    matrix <- data.frame(as.factor(rep(c(0,1),each=2^(input$xn + input$zn - 1))))
    #ifelse(input$zn!=0, z=rep(c(0,1),each=2), NA),x=rep(c(0,1),times=2)
    names(matrix) <- c("z1")
    if(input$zn == 2){
      matrix$z2 <- as.factor(rep(c(0,1),time=2^(input$xn + input$zn - 1)))
      if(input$xn != 0){
        matrix$x1 <- as.factor(rep(c(0,0,1,1),time=2^(input$xn + input$zn - 2)))
        if(input$xn == 2){
          matrix$x2 <- as.factor(rep(c(0,0,0,0,1,1,1,1),time=2))
        }
      }
    }
    else{
      if(input$xn != 0){
        matrix$x1 <- as.factor(rep(c(0,1),time=2^(input$xn + input$zn - 1)))
        if(input$xn == 2){
          matrix$x2 <- as.factor(c(0,0,1,1,0,0,1,1))
        }
      }
    }
    matrix
  })
  
  wgts <- reactive({
    if(input$wgttf){
      rep(1/2^(input$xn + input$zn),2^(input$xn + input$zn))
    }
    else{
      as.numeric(unlist(strsplit(input$weight,",")))
    }
  })
  
  xvar <- reactive({
    #as.formula(ifelse(input$xn == 0,"~0",paste0("~",paste0("x",1:input$xn,collapse="+"))))
    if(input$xn == 0){
      formula("~0")
    }
    else if(input$xn == 1){
      formula("~x1")
    }
    else {
      formula("~x1+x2")
    }
  })
  
  zvar <- reactive({
    #as.formula(ifelse(input$zn == 0,"~0",paste0("~",paste0("z",1:input$zn,collapse="+"))))
    if(input$zn == 0){
      formula("~0")
    }
    else if(input$zn == 1){
      formula("~z1")
    }
    else {
      formula("~z1+z2")
    }
  })
  
  la <- reactive({
    if (input$xn == 1){
      input$lambda1
    }
    else if (input$xn == 0) {
      NULL
    }
    else{
      c(input$lambda1,input$lambda2)
    }
  })
  
  ps <- reactive({
    if (input$zn == 1){
      input$psi1
    }
    else{
      c(input$psi1,input$psi2)
    }
  })
  
  ps0 <- reactive({
    if (input$zn == 1){
      0
    }
    else{
      c(0,0)
    }
  })
  
  res <- reactive({powerglm(Xvar=xvar(), Zvar=zvar(), lambda=la(), psi=ps(), psi0=ps0(), mu_hat=input$mu_hat, family=input$family, dt(), wgts(), N=ifelse(input$powersample,input$n,1), signif_alpha=input$size, sorl=input$scorelrt)})
  
  output$matrix <- renderTable({
    matrix <- dt()
    matrix$weights <- wgts()
    matrix
  })
  
  output$txtout <- renderText({
    ifelse(input$powersample,res()$asy.power,NA)
  })
  output$txtout1 <- renderText({
    n.sol <- reactive(ceiling(uniroot(function(n){ pchisq(qchisq(1-input$size,df=1),df=1,ncp=n*res()$ncp,lower.tail = FALSE) - input$nominal_pow }, c(0, 1e100), tol = 1e-10)$root))
    ifelse(input$powersample,NA,n.sol())
  })
  
  output$PPlot <- renderPlot({
    if(input$delim){
      sl <- c(ceiling(input$n*0.8),ceiling(input$n*0.9),input$n,ceiling(input$n*1.1),ceiling(input$n*1.2))
    }
    else{
      sl <- ceiling(seq(input$Lower,input$Upper,length.out = input$INTV))
    }
    psil <- c(input$psi1 - 2*input$diffpsi,input$psi1 - input$diffpsi,input$psi1,input$psi1 + input$diffpsi,input$psi1 + 2*input$diffpsi)
    psi1_value <- c()
    pl <- c()
    for (j in psil) {
      for(i in sl){
        if (input$zn == 1){
          psival = j 
        }
        else{
          psival =c(j,input$psi2)
        }
        nres <- powerglm(Xvar=xvar(), Zvar=zvar(), lambda=la(), psi=psival, psi0=ps0(), mu_hat=input$mu_hat, family=input$family, data=dt(), wgts=wgts(), N = i, signif_alpha=input$size, sorl=input$scorelrt)
        pl <- c(pl,nres$asy.power)
        psi1_value <- c(psi1_value,j)
      }
    }
    sl <- rep(sl,length(psil))
    psi1_value <- factor(psi1_value,labels = round(psil,2))
    df <- data.frame(cbind(sl, pl, psi1_value)) 
    #plot(sl,pl,type="b",main = "Power Plot",xlab = "Sample Size",ylab = "Estimated Power")
    ggplot(df,aes(sl, pl, color = factor(psi1_value,labels = round(psil,2)), group=psi1_value)) +
      geom_line() +
      geom_point() +
      scale_y_continuous(labels = scales::percent) +
      #scale_x_log10() +
      labs(x = "Sample Size (N)", y = "Estimated Power (%)", title = "Estimated Power Plot")+
      theme_classic()+
      guides(color = guide_legend(title = "psi1"))
  })
  
  output$SPlot <- renderPlot({
    if(input$delim){
      pl <- c(input$nominal_pow*0.8,input$nominal_pow*0.9,input$nominal_pow,input$nominal_pow*1.1,input$nominal_pow*1.2)
    }
    else{
      pl <- seq(input$Lower,input$Upper,length.out = input$INTV)
    }
    sl <- c()
    psil <- c(input$psi1 - 2*input$diffpsi,input$psi1 - input$diffpsi,input$psi1,input$psi1 + input$diffpsi,input$psi1 + 2*input$diffpsi)
    psi1_value <- c()
    
    for (j in psil) {
      for(i in pl){
        if (input$zn == 1){
          psival = j 
        }
        else{
          psival =c(j,input$psi2)
        }
        nres <- powerglm(Xvar=xvar(), Zvar=zvar(), lambda=la(), psi=psival, psi0=ps0(), mu_hat=input$mu_hat, family=input$family, data=dt(), wgts=wgts(), N = 1, signif_alpha=input$size, sorl=input$scorelrt)
        sl <- c(sl,ceiling(uniroot(function(n){ pchisq(qchisq(1-input$size,df=1),df=1,ncp=n*nres$ncp,lower.tail = FALSE) - i }, c(0, 1e100), tol = 1e-10)$root))
        psi1_value <- c(psi1_value,j)
      }
    }
    #plot(pl,sl,type="b",main = "Sample Size Plot",ylab = "Estimated Sample Size",xlab = "Nominal Power")
    pl <- rep(pl,length(psil))
    psi1_value <- factor(psi1_value,labels = round(psil,2))
    df <- data.frame(cbind(pl, sl, psi1_value)) 
    #plot(sl,pl,type="b",main = "Power Plot",xlab = "Sample Size",ylab = "Estimated Power")
    ggplot(df,aes(pl,sl, color = factor(psi1_value,labels = round(psil,2)), group=psi1_value)) +
      geom_line() +
      geom_point() +
      scale_x_continuous(labels = scales::percent) +
      #scale_y_log10() +
      labs(y = "Estimated Sample Size (N)", x = "Nominal Power (%)", title = "Estimated Sample Size Plot") +
      theme_classic() +
      guides(color = guide_legend(title = "psi1"))
  })
  
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)
