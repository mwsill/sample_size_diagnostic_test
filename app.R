#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(markdown)
library(binom)

samplesizeNormapprox <- function(alpha, beta, prop, nullp){ 
  ceiling(((qnorm((1-alpha))*sqrt(nullp*(1-nullp))) +
             (qnorm((1-beta))*sqrt(prop*(1-prop))))^2/((prop-nullp)^2))
}

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Sample size for phase 2 studies: Retrospective validation of a binary test"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(width=3,
        # Simple integer interval
        sliderInput("alpha", "alpha:", 
                    min=0, max=1, value=.05, step= 0.01),
        
        # Decimal interval with step valu1e
        sliderInput("beta", "beta:", 
                    min=0, max=1, value=.2, step= 0.01),
        
        # Specification of range within an interval
        sliderInput("tpf1", "Sensitivity H1 (tpf1):", 
                    min=0, max=1, value=.80, step= 0.01),
        
        # Specification of range within an interval
        sliderInput("tpf0", "Sensitivity H0 (tpf0)", 
                    min=0, max=1, value=.63, step= 0.01),
        
        # Specification of range within an interval
        sliderInput("fpf1", "1-Specificity H1 (fpf1):", 
                    min=0, max=1, value=.1, step= 0.01),
        
        # Specification of range within an interval
        sliderInput("fpf0", "1-Specificity H0 (fpf0):", 
                    min=0, max=1, value=.16, step= 0.01),
        
        # Specification of range within an interval
        sliderInput("b", "SimulationSteps (b):", 
                    min=1000, max=10000, value=1000, step= 1000),
        
        selectInput("ci", "binomial CI for simulation:",
                    list("Wilson" = "wilson",
                         "exact Pearson-Klopper" = "exact", 
                         "asymptotic Wald"  = "asymptotic", 
                         "Agresti-Coull" = "ac")),
        plotOutput("cregion",height = "250px")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        "This tool can be used to estimates sample size for a phase 2
        study where the aim is to test that the performance of medical
        test is sufficiently good to undergo further development.
        The user needs to define values for the false positive fraction (FPF; specificity)
        and true positive fraction (TPF; sensitivity) that are minimally acceptable in order to design the study.
        The study will test the null hypothesis:
        ",
        withMathJax("$$H_0 : \\left[ \\text{TPF} \\leq \\text{TPF}_0 \\text{ or } \\text{FPF} \\geq \\text{FPF}_0 \\right] $$"),
        "From a study that rejects the null hypothesis it will be concluded that TPF and FPF meet the minimal criteria.
        Sample sizes are calculated by using the formula based on asymptotic normal distribution theory as described in Pepe (section 8.2, 2003) or by simulations
        using different exact or approximate confidence intervals for the difference of binomial
        proportions as described in Agresti and Coull (1998).",
       
        #includeMarkdown("ref.md"),
        
         plotOutput("plot",height = "800px")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$plot <- renderPlot({
     
     layout(matrix(c(1,2,3,4),2))
     
     #print(input$beta)
     #print(input$alpha)
     
     betas <- seq(0.5,0.01,by=-.01)
     beta_star <- round(1-sqrt(1-input$beta),3)
     betas <- sort(c(betas,beta_star),decreasing=TRUE)
     alpha_star <- round(1-sqrt(1-input$alpha),3)
     power <- 1-betas
     n <- samplesizeNormapprox(alpha_star,betas,input$tpf1,input$tpf0) 
     plot(power~n,type="b",main="Achieved power for the one-sided test (TPF)",lwd=2,pch=16,xlab="N diseased group",ylab="Power")
     n <- seq(min(n),max(n),by=1)
     simbeta1 <- function(n) sum(binom.confint(rbinom(input$b,n,input$tpf1),n,1-(alpha_star*2),methods=input$ci)$lower > input$tpf0)/input$b
     betasex <- simplify2array(lapply(n,simbeta1))
     powerex <- 1-betasex
     points(betasex~n,type="b",col="blue",lwd=2,pch=16)
     abline(h=sqrt(1-input$beta),col="red",lwd=2)
     nnorm <- samplesizeNormapprox(alpha_star,beta_star,input$tpf1,input$tpf0)
     abline(v=nnorm ,col="black",lwd=2)
     nexact <- n[which(betasex >= sqrt(1-input$beta))[1]]
     bexact <- which(betasex >= sqrt(1-input$beta))[1]
     abline(v=n[bexact],col="blue",lwd=2)
     #abline(h=betasex[bexact],col="blue",lwd=2)
     legend("bottomright",cex=1.5,legend=c(paste0("normal approximation N = ",nnorm),paste0("simulation binom.confint N = ",nexact),
                                           paste0("1-beta_star = ",1-beta_star)),text.col=c("black","blue","red"),pch=16,lwd=2,
                                           col=c("black","blue","red"),bty="n")
     
     
  
     simbeta1 <- function(n) sum(binom.confint(rbinom(input$b,n,input$tpf1),n,1-(alpha_star*2),methods=input$ci)$lower <= input$tpf0)/input$b
     alphasex <- simplify2array(lapply(n,simbeta1))
     plot(alphasex~n,type="b",main="Simulated alpha",lwd=1,pch=16,xlab="N diseased group",ylab="alpha",col="blue")
     abline(h=alpha_star,col="red",lwd=2) 
     nalpha<-which(alphasex<=alpha_star)[1]
     abline(v=n[nalpha],col="blue",lwd=2,lty=2) 
     abline(v=nnorm ,col="black",lwd=2)
     abline(v=n[bexact],col="blue",lwd=2)
     text(x=0.05,y=n[nalpha],labels=paste(n[nalpha]),col="blue")
     legend("topright",cex=1.5,legend=c(paste0("normal approximation N = ",nnorm),paste0("simulation binom.confint N = ",nexact),
                                           paste0("alpha_star = ",alpha_star)),text.col=c("black","blue","red"),pch=16,lwd=2,
                                           col=c("black","blue","red"),bty="n")
     
     
     
     n <- samplesizeNormapprox(alpha_star,betas,input$fpf1,input$fpf0) 
     plot(power~n,type="b",main="Achieved power for the one-sided test (FPF)",lwd=2,pch=16,xlab="N non-diseased group",ylab="Power")
     n <- seq(min(n),max(n),by=1)
     simbeta2 <- function(n) sum(binom.confint(rbinom(input$b,n,input$fpf1),n,1-((1-sqrt(1-input$alpha))*2),methods=input$ci)$upper < input$fpf0)/input$b
     betasex <- simplify2array(lapply(n,simbeta2))
     powerex <- 1-betasex
     points(betasex~n,type="b",col="blue",lwd=2,pch=16)
     abline(h=sqrt(1-input$beta),col="red",lwd=2)
     nnorm <- samplesizeNormapprox(alpha_star,beta_star,input$fpf1,input$fpf0)
     abline(v=nnorm ,col="black",lwd=2)
     nexact <- n[which(betasex >= sqrt(1-input$beta))[1]]
     bexact <-which(betasex >= sqrt(1-input$beta))[1]
     abline(v=n[bexact],col="blue",lwd=2)
     legend("bottomright",cex=1.5,legend=c(paste0("normal approximation N = ",nnorm),paste0("simulation binom.confint N = ",nexact),
                                           paste0("1-beta_star = ",1-beta_star)),text.col=c("black","blue","red"),pch=16,lwd=2,
            col=c("black","blue","red"),bty="n")
     
     n <- samplesizeNormapprox(alpha_star,betas,input$fpf1,input$fpf0) 
     n <- seq(min(n),max(n),by=1)
     simbeta2 <- function(n) sum(binom.confint(rbinom(input$b,n,input$fpf1),n,1-((1-sqrt(1-input$alpha))*2),methods=input$ci)$upper >= input$fpf0)/input$b
     alphasex <- simplify2array(lapply(n,simbeta2))
     plot(alphasex~n,type="b",main="Simulated alpha",lwd=2,pch=16,xlab="N non-diseased group",ylab="alpha",col="blue")
     abline(h=alpha_star,col="red",lwd=2) 
     nalpha<-which(alphasex<=alpha_star)[1]
     abline(v=n[nalpha],col="blue",lwd=2,lty=2) 
     abline(v=nnorm ,col="black",lwd=2)
     abline(v=n[bexact],col="blue",lwd=2)
     text(x=0.05,y=n[nalpha],labels=paste(n[nalpha]),col="blue")
     legend("topright",cex=1.5,legend=c(paste0("normal approximation N = ",nnorm),paste0("simulation binom.confint N = ",nexact),
                                        paste0("alpha_star = ",alpha_star)),text.col=c("black","blue","red"),pch=16,lwd=2,
            col=c("black","blue","red"),bty="n")
     
     
   })
   
   output$cregion <- renderPlot({
     plot(0,0,ylim=c(0,1),xlim=c(0,1),type="n",axes=FALSE,ylab="Sensitivity",xlab="1-Specificity")
     abline(0,1,pch=2)
     abline(h=input$tpf0,col="red",lty=2,lwd=2)
     abline(v=input$fpf0,col="red",lty=2,lwd=2)
     points(x=input$fpf1,y=input$tpf1,cex=2,pch=19)
     axis(1)
     axis(2)
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

