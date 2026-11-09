library(shiny)
library(ggplot2)
library(plotly)
library(viridis)
library(dplyr)
library(scales)

# Source parameter file
source("Parameters.R")

region_coef <- list(
  "Western Europe" = list(
    intercept = 37.94,
    log_age = -5.09,
    log_bmi = -6.34
  ),
  "South Asia" = list(
    intercept = 42.68,
    log_age = -7.03,
    log_bmi = -6.92
  ),
  "Eastern Europe" = list(
    # Refit 2026-07-26 after the contributing investigators corrected five
    # erroneous BMI records against source documents; n = 3287 (no exclusions).
    # Raw intercept 36.8935; the value below includes the prevalence offset to
    # set prevalence to 5.2%.
    # Must match Table 2. Superseded: intercept 26.55, log_age -5.09, log_bmi -3.32.
    intercept = 31.20,
    log_age = -4.55,
    log_bmi = -5.31
  ),
  "Northern Europe" = list(
    intercept = 24.71,
    log_age = -4.45,
    log_bmi = -3.01
  )
)

# Age at diagnosis actually spanned by each contributing dataset. The input
# range accepted by the tool stays 15-60 for every region; these bounds drive
# an annotation, not a restriction, so the user can see where an estimate
# extrapolates beyond the data the model was fitted to. Must match Table 1.
region_age_range <- list(
  "Western Europe"  = list(min = 18, max = 50),  # Exeter
  "South Asia"      = list(min = 15, max = 40),  # Mohan Clinics
  "Eastern Europe"  = list(min = 15, max = 50),  # Ukraine T1D Exomes (age at onset 15-50 in the analysis set)
  "Northern Europe" = list(min = 18, max = 60)   # Scania ANDIS
)

# Upper age at diagnosis of the South Asia cohort from which the C-peptide
# gamma distributions were derived. Above this the fitted distributions
# understate residual secretion, and so understate the evidence for T1D at
# values around 0.3 nmol/L (see Limitation 2).
CPEPTIDE_AGE_LIMIT <- 40




calc_cpeptide_lr <- function(cps, cps_type = "stimulated") {
  # Return 0 if invalid input or not fasting
  if (is.na(cps) || !is.finite(cps) || cps <= 0) return(0)
  if (cps_type != "fasting") return(0)
  
  cpep_params <- params$cpeptide$fasting
  
  # Calculate gamma densities with error handling
  density_t1d <- tryCatch({
    dgamma(cps, shape = cpep_params$t1d$alpha, rate = cpep_params$t1d$lambda)
  }, error = function(e) 0)
  
  density_t2d <- tryCatch({
    dgamma(cps, shape = cpep_params$t2d$alpha, rate = cpep_params$t2d$lambda)
  }, error = function(e) 0)
  
  # Check for valid densities
  if (density_t2d == 0 || density_t1d == 0 || 
      !is.finite(density_t1d) || !is.finite(density_t2d)) return(0)
  
  # Calculate likelihood ratio
  lr <- density_t1d / density_t2d
  if (lr <= 0 || !is.finite(lr)) return(0)
  
  # Return log likelihood ratio
  log_lr <- tryCatch(log(lr), error = function(e) 0)
  return(log_lr)
}

calc_t1d_prob <- function(age, bmi, male = NA, fh = NA, gad_status = NA, ia2_status = NA,
                          zinc_status = NA, iaa_status = NA, ethnicity = "Western Europe", cps = NA, cps_type = "stimulated") {
  
  # Input validation
  if (is.na(age) || is.na(bmi) || !is.finite(age) || !is.finite(bmi) || age <= 0 || bmi <= 0) {
    return(NA_real_)
  }
  if (age < 15 || age > 60) return(NA_real_)
  
  # BMI range validation by ethnicity
  if (ethnicity == "Western Europe") {
    if (bmi < 17.5 || bmi > 70) return(NA_real_)
  } else {
    if (bmi < 15 || bmi > 60) return(NA_real_)
  }
  
  # Check if ethnicity is implemented
  if (!ethnicity %in% names(region_coef)) return(NA_real_)
  
  # Calculate base logit
  coef <- region_coef[[ethnicity]]
  hp_logit <- coef$intercept + (coef$log_age * log(age)) + (coef$log_bmi * log(bmi))
  
  # Add pooled coefficients (only male sex and T2D family history)
  if (!is.na(male) && is.finite(male)) {
    hp_logit <- hp_logit + (params$pooled$male_sex * male)
  }
  if (!is.na(fh) && is.finite(fh)) {
    hp_logit <- hp_logit + (params$pooled$fh * fh)
  }
  
  # Add antibody likelihood ratios
  if (!is.na(gad_status) && is.finite(gad_status)) {
    hp_logit <- hp_logit + ifelse(gad_status == 1, params$lr$logLR_pos[1], ifelse(gad_status == 0, params$lr$logLR_neg[1], 0))
  }
  if (!is.na(ia2_status) && is.finite(ia2_status)) {
    hp_logit <- hp_logit + ifelse(ia2_status == 1, params$lr$logLR_pos[2], ifelse(ia2_status == 0, params$lr$logLR_neg[2], 0))
  }
  if (!is.na(zinc_status) && is.finite(zinc_status)) {
    hp_logit <- hp_logit + ifelse(zinc_status == 1, params$lr$logLR_pos[3], ifelse(zinc_status == 0, params$lr$logLR_neg[3], 0))
  }
  if (!is.na(iaa_status) && is.finite(iaa_status)) {
    hp_logit <- hp_logit + ifelse(iaa_status == 1, params$lr$logLR_pos[4], ifelse(iaa_status == 0, params$lr$logLR_neg[4], 0))
  }
  
  # Add C-peptide likelihood ratio
  if (!is.na(cps) && is.finite(cps) && cps > 0) {
    cps_lr_log <- calc_cpeptide_lr(cps, cps_type)
    if (is.finite(cps_lr_log) && cps_lr_log != 0) {
      hp_logit <- hp_logit + cps_lr_log
    }
  }
  
  # Final validation of logit
  if (is.na(hp_logit) || !is.finite(hp_logit)) return(NA_real_)
  
  # Convert logit to probability with numerical stability
  if (hp_logit > 700) return(1.0)
  if (hp_logit < -700) return(0.0)
  
  # Use built-in plogis function for numerical stability
  prob <- tryCatch({
    plogis(hp_logit)
  }, error = function(e) {
    # Fallback calculation
    if (hp_logit > 500) return(1.0)
    if (hp_logit < -500) return(0.0)
    
    exp_logit <- exp(hp_logit)
    if (!is.finite(exp_logit)) return(1.0)
    
    exp_logit / (1 + exp_logit)
  })
  
  if (is.na(prob) || !is.finite(prob)) return(NA_real_)
  
  # Ensure probability is in [0,1]
  return(pmax(0, pmin(1, prob)))
}

generate_nomogram_data <- function(ethnicity = "Western Europe", male = NA, fh = NA, 
                                   gad_status = NA, ia2_status = NA, zinc_status = NA, iaa_status = NA,
                                   cps = NA, cps_type = "stimulated") {
  
  # Define ranges based on ethnicity
  age_range <- seq(15, 60, by = 0.5)
  bmi_range <- if (ethnicity == "Western Europe") {
    seq(17.5, 50, by = 0.5)
  } else {
    seq(15, 60, by = 0.5)
  }
  
  # Create grid
  grid <- expand.grid(age = age_range, bmi = bmi_range)
  
  # Calculate probabilities with error handling
  grid$probability <- tryCatch({
    mapply(calc_t1d_prob, grid$age, grid$bmi,
           MoreArgs = list(male = male, fh = fh, gad_status = gad_status, 
                           ia2_status = ia2_status, zinc_status = zinc_status, iaa_status = iaa_status,
                           ethnicity = ethnicity, cps = cps, cps_type = cps_type),
           SIMPLIFY = TRUE)
  }, error = function(e) {
    rep(NA_real_, nrow(grid))
  })
  
  # Return only valid data
  valid_data <- grid[!is.na(grid$probability) & is.finite(grid$probability) & 
                       grid$probability >= 0 & grid$probability <= 1, ]
  
  if (nrow(valid_data) == 0) {
    return(data.frame(age = numeric(0), bmi = numeric(0), probability = numeric(0)))
  }
  
  return(valid_data)
}

js <- "
Shiny.addCustomMessageHandler('resetHtWt', function(msg){
  ['height_val','weight_val'].forEach(function(id){
    var el = document.getElementById(id);
    if(el) el.value = '';
  });
});

$(document).ready(function(){
  var style = document.createElement('style');
  style.type = 'text/css';
  style.innerHTML = `
    #age::-webkit-outer-spin-button, #age::-webkit-inner-spin-button,
    #bmi::-webkit-outer-spin-button, #bmi::-webkit-inner-spin-button,
    #height_val::-webkit-outer-spin-button, #height_val::-webkit-inner-spin-button,
    #weight_val::-webkit-outer-spin-button, #weight_val::-webkit-inner-spin-button,
    #cps::-webkit-outer-spin-button, #cps::-webkit-inner-spin-button {
      -webkit-appearance: none;
      margin: 0;
    }
    
    #age[type=number], #bmi[type=number], #height_val[type=number],
    #weight_val[type=number], #cps[type=number] {
      -moz-appearance: textfield;
    }
  `;
  
  if (document.head) document.head.appendChild(style);
    
  function disableIncrement(selector) {
    $(document).on('keydown', selector, function(e) {
      if (e.which === 38 || e.which === 40) {
        e.preventDefault();
        return false;
      }
    });
    
    $(document).on('wheel', selector, function(e) {
      e.preventDefault();
      if (this.blur) this.blur();
      return false;
    });
  }
  
  ['#age', '#bmi', '#height_val', '#weight_val', '#cps'].forEach(disableIncrement);
  
  // Disable unimplemented region options
  setTimeout(function() {
    $('#ethnicity option').each(function() {
      if ($(this).text().includes('(not implemented)')) {
        $(this).prop('disabled', true).css('color', '#999');
      }
    });
    
    // Disable stimulated C-peptide option
    $('#cps_type input[value=\"stimulated\"]').prop('disabled', true);
    $('#cps_type label:contains(\"Stimulated\")').css('color', '#999');
    
    // Disable T1D family history option
    $('#fh input[value=\"na_t1d\"]').prop('disabled', true);
    $('#fh label:contains(\"T1D (not implemented)\")').css('color', '#999');
    
    // DKA field is now enabled - removed disabling code
    
    // Disable Acanthosis Nigricans field
    // $('#an input').prop('disabled', true);
    // $('#an label').css('color', '#999');
    // $('#an .control-label').css('color', '#999');
  }, 100);
});
"

ui <- fluidPage(
  title = "Prototype Tool to Estimate the Probability of Type 1 vs Type 2 Diabetes",
  tags$head(tags$script(HTML(js))),
  titlePanel(
    div(
      h3("Prototype Tool to Estimate the Probability of Type 1 vs Type 2 Diabetes", style = "margin-bottom:2px; margin-top:5px;"),
      h6("For adults diagnosed with stage-3 overt hyperglycemia within the past year.", style = "color:gray; margin-bottom:10px;")
    )
  ),
  fluidRow(
    column(6,
           wellPanel(
             h5("Required inputs (Region, Age, & BMI) – others optional", style = "margin-bottom:8px;"),
             selectInput("ethnicity", "Region:",
                         choices = c("Western Europe", "South Asia", "Eastern Europe", "Northern Europe", "East Asia (not implemented)", "Middle East/North Africa (not implemented)", "Latin America (not implemented)", "Sub-Saharan Africa (not implemented)"),
                         selected = "Western Europe"),
             conditionalPanel(
               condition = "input.ethnicity == 'Western Europe'",
               tags$p("Valid Age: 15–60 yrs   |   Valid BMI: 17.5–70 kg/m²", style = "color:gray; margin-bottom:1px; font-size:12px;")
             ),
             conditionalPanel(
               condition = "input.ethnicity != 'Western Europe'",
               tags$p("Valid Age: 15–60 yrs   |   Valid BMI: 15–60 kg/m²", style = "color:gray; margin-bottom:1px; font-size:12px;")
             ),
             numericInput("age", "Age at diagnosis (yrs):", value = NA, min = 15, max = 60, step = 1),
             tags$p("Enter BMI directly, or use height/weight below.", style = "font-size:12px; margin-bottom:1px; color:gray;"),
             numericInput("bmi", "BMI (kg/m²):", value = NA, min = 15, max = 70, step = 0.1),
             tags$hr(style = "margin:2px 0;"),
             tags$div(
               radioButtons("units", "Unit system:",
                            choices = c("Metric (cm / kg)" = "metric", "Imperial (in / lb)" = "imperial"),
                            selected = "metric", inline = TRUE),
               uiOutput("htwt_inputs"),
               actionButton("use_hw", "Use Height + Weight", size = "sm"),
               style = "margin-bottom:4px;"
             ),
             tags$hr(style = "margin:2px 0;"),
             radioButtons("sex", "Sex:",
                          choices = c("Female" = "0", "Male" = "1", "Not specified" = "na"),
                          selected = "na", inline = TRUE),
             radioButtons("fh", "History of diabetes in a parent:",
                          choices = c("T1D (not implemented)" = "na_t1d", "T2D or diabetes not otherwise specified" = "1", "No known history of diabetes in a parent" = "na_none"),
                          selected = "na_none", inline = TRUE),
             tags$div(
               style = "display: flex; gap: 20px; align-items: flex-start;",
               tags$div(
                 style = "flex: 1;",
                 radioButtons("dka", "DKA (with antibody positive, rules in T1D):",
                              choices = c("Yes" = "1", "No" = "0"),
                              selected = "0", inline = TRUE)
               )
               # tags$div(
               #   style = "flex: 1;",
               #   radioButtons("an", "Acanthosis Nigricans (not implemented):",
               #                choices = c("Yes" = "1", "No" = "0"),
               #                selected = "0", inline = TRUE)
               # )
             ),
             tags$hr(style = "margin:2px 0;"),
             conditionalPanel(
               condition = "input.age && input.bmi && input.ethnicity",
               tags$div(
                 tags$strong("Pre-lab T1D probability: ", style = "font-size:14px; color:#2c5282;"),
                 tags$span(textOutput("prelab_prediction", inline = TRUE), style = "font-weight:bold; font-size:18px; color:#1a365d;"),
                 style = "background-color: #e6f3ff; padding: 8px; border-radius: 5px; margin-bottom: 6px; border-left: 4px solid #3182ce;"
               )
             ),
             tags$hr(style = "margin:2px 0;"),
             h6("Lab Results", style = "margin-bottom:4px; margin-top:4px;"),
             radioButtons("gad", "GADA status:",
                          choices = c("Positive" = "1", "Negative" = "0", "Not tested" = "na"),
                          selected = "na", inline = TRUE),
             radioButtons("ia2", "IA-2 status:",
                          choices = c("Positive" = "1", "Negative" = "0", "Not tested" = "na"),
                          selected = "na", inline = TRUE),
             radioButtons("zinc", "ZNT8 status:",
                          choices = c("Positive" = "1", "Negative" = "0", "Not tested" = "na"),
                          selected = "na", inline = TRUE),
             radioButtons("iaa", "Anti-insulin autoantibody (IAA), pre-insulin:",
                          choices = c("Positive" = "1", "Negative" = "0", "Not tested" = "na"),
                          selected = "na", inline = TRUE),
             numericInput("cps", "Fasting C-peptide (nmol/L) within 12 months of diagnosis/onset:", value = NA, min = 0, max = 10, step = 0.1),
             conditionalPanel(
               condition = "input.cps && !isNaN(parseFloat(input.cps)) && parseFloat(input.cps) >= 0",
               radioButtons("cps_type", "C-peptide type:",
                            choices = c("Stimulated (not implemented)" = "stimulated", "Fasting" = "fasting"),
                            selected = "fasting", inline = TRUE)
             )
           )
    ),
    column(6,
           wellPanel(
             actionButton("reset", "Reset all"),
             br(),
             tabsetPanel(
               tabPanel("Model predictions",
                        tags$div(
                          textOutput("prediction"),
                          style = "background-color: #fff5f5; padding: 12px; border-radius: 8px; margin-bottom: 12px;
                                   border: 2px solid #fc8181; font-weight: bold; font-size: 16px; color: #742a2a;
                                   box-shadow: 0 2px 4px rgba(0,0,0,0.1);"
                        ),
                        uiOutput("extrapolation_notes"),
                        checkboxInput("show_patient_point", "Show patient on nomogram", TRUE),
                        plotlyOutput("nomogram_plot", height = "500px")),
               tabPanel("Model information",
                        "The model uses region-specific intercepts and coefficients for log(age), and log(BMI). ",
                        "Male sex and history of diabetes in a parent (T2D) use pooled likelihood ratios across all regions. ",
                        "GADA, IA-2, ZNT8, and C-peptide likelihood ratios also use pooled likelihood ratios across all regions. ",
                        "Not implemented: regions East Asia, Middle East/North Africa, Latin America, and Sub-Saharan Africa; T1D history of diabetes in a parent.",
                        tags$br(), tags$br(),
                        tags$h6("Autoantibody Likelihood Ratios"),
                        textOutput("antibody_lr_display"),
                        tags$br(),
                        tags$h6("C-peptide Implementation"),
                        "C-peptide likelihood ratios are calculated using gamma distributions for fasting C-peptide only. ",
                        textOutput("cpeptide_params_display"),
                        "Stimulated C-peptide parameters are not yet implemented."),
               tabPanel("Diabetes definitions",
                        "T1D = at least 3 years of follow-up with sustained insulin treatment by 3 years, and random (non-fasting) or stimulated C-peptide < 0.2 nmol/L (0.6 ng/mL) measured ≥ 3 years after diagnosis;",
                        tags$br(),
                        "T2D = residual C-peptide (≥ 0.2 nmol/L) measured ≥ 3 years after diagnosis, or no sustained insulin requirement.")
             )
           )
    )
  )
)

server <- function(input, output, session) {
  
  # Input validation with notifications
  observeEvent(input$age, {
    if (!is.na(input$age) && (input$age < 15 || input$age > 60)) {
      showNotification("Age must be between 15 and 60 years", type = "warning")
    }
  })
  
  observeEvent(input$bmi, {
    if (!is.na(input$bmi)) {
      if (input$ethnicity == "Western Europe") {
        if (input$bmi < 17.5 || input$bmi > 70) {
          showNotification("BMI must be between 17.5 and 70 for Western Europe", type = "warning")
        }
      } else if (input$ethnicity %in% names(region_coef)) {
        if (input$bmi < 15 || input$bmi > 60) {
          showNotification("BMI must be between 15 and 60 for this region", type = "warning")
        }
      }
    }
  })
  
  observeEvent(input$cps, {
    if (!is.na(input$cps) && (input$cps < 0 || input$cps > 10)) {
      showNotification("C-peptide must be between 0 and 10 nmol/L", type = "warning")
    }
  })
  
  # Handle ethnicity changes
  observeEvent(input$ethnicity, {
    # Check for unimplemented regions and reset to Western Europe
    if (grepl("\\(not implemented\\)", input$ethnicity)) {
      showNotification("This region is not yet implemented. Please select an available region.", type = "error")
      updateSelectInput(session, "ethnicity", selected = "Western Europe")
      return()
    }
    
    current_bmi <- input$bmi
    if (input$ethnicity == "Western Europe") {
      updateNumericInput(session, "bmi", min = 17.5, max = 70, value = current_bmi)
    } else {
      updateNumericInput(session, "bmi", min = 15, max = 60, value = current_bmi)
      if (!is.na(current_bmi) && input$ethnicity %in% names(region_coef) && 
          (current_bmi < 15 || current_bmi > 60)) {
        updateNumericInput(session, "bmi", value = NA)
        showNotification("BMI cleared due to region change - please re-enter", type = "warning")
      }
    }
  })
  
  # Process inputs reactive
  processed_inputs <- reactive({
    list(
      male = if (input$sex == "na") NA else as.numeric(input$sex),
      fh = if (input$fh == "1") 1 else NA,  # T2D or Diabetes NOS history in a parent
      gad_status = if (input$gad == "na") NA else as.numeric(input$gad),
      ia2_status = if (input$ia2 == "na") NA else as.numeric(input$ia2),
      zinc_status = if (input$zinc == "na") NA else as.numeric(input$zinc),
      iaa_status = if (input$iaa == "na") NA else as.numeric(input$iaa),
      dka_status = if (input$dka == "na") NA else as.numeric(input$dka),
      cps = if (is.na(input$cps)) NA else as.numeric(input$cps),
      cps_type = if (is.na(input$cps)) "fasting" else input$cps_type,
      ethnicity = input$ethnicity
    )
  })
  
  # Nomogram data caching
  nomogram_cache <- reactiveValues(
    last_key = NULL,
    data = data.frame(age = numeric(0), bmi = numeric(0), probability = numeric(0))
  )
  
  nomogram_data <- reactive({
    req(input$ethnicity)
    p <- processed_inputs()
    
    # Create cache key
    cache_key <- paste(p$ethnicity, p$male, p$fh, p$gad_status, p$ia2_status, 
                       p$zinc_status, p$iaa_status, p$cps, p$cps_type, sep = "_")
    
    # Generate data if cache miss
    if (is.null(nomogram_cache$last_key) || cache_key != nomogram_cache$last_key) {
      nomogram_cache$last_key <- cache_key
      nomogram_cache$data <- tryCatch({
        generate_nomogram_data(p$ethnicity, p$male, p$fh, p$gad_status, 
                               p$ia2_status, p$zinc_status, p$iaa_status, p$cps, p$cps_type)
      }, error = function(e) {
        data.frame(age = numeric(0), bmi = numeric(0), probability = numeric(0))
      })
    }
    
    return(nomogram_cache$data)
  })
  
  # Individual prediction
  prediction_prob <- reactive({
    req(!is.na(input$age), !is.na(input$bmi))
    p <- processed_inputs()
    tryCatch({
      calc_t1d_prob(input$age, input$bmi, p$male, p$fh, p$gad_status, 
                    p$ia2_status, p$zinc_status, p$iaa_status, p$ethnicity, p$cps, p$cps_type)
    }, error = function(e) NA_real_)
  })
  
  # Pre-lab prediction (only Region, Age, BMI, Sex, Family history)
  prelab_prediction_prob <- reactive({
    req(!is.na(input$age), !is.na(input$bmi), input$ethnicity)
    p <- processed_inputs()
    tryCatch({
      # Call calc_t1d_prob with only non-lab inputs (no antibodies or C-peptide)
      calc_t1d_prob(input$age, input$bmi, p$male, p$fh, 
                    gad_status = NA, ia2_status = NA, zinc_status = NA, iaa_status = NA,
                    p$ethnicity, cps = NA, cps_type = "fasting")
    }, error = function(e) NA_real_)
  })
  
  # Height/weight input UI
  output$htwt_inputs <- renderUI({
    if (input$units == "metric") {
      tags$div(
        style = "display: flex; gap: 10px;",
        tags$div(
          numericInput("height_val", "Height (cm):", value = NA, min = 50, max = 250, step = 1),
          style = "flex: 1; margin-bottom: 5px;"
        ),
        tags$div(
          numericInput("weight_val", "Weight (kg):", value = NA, min = 20, max = 250, step = 0.1),
          style = "flex: 1; margin-bottom: 5px;"
        )
      )
    } else {
      tags$div(
        style = "display: flex; gap: 10px;",
        tags$div(
          numericInput("height_val", "Height (in):", value = NA, min = 20, max = 100, step = 0.1),
          style = "flex: 1; margin-bottom: 5px;"
        ),
        tags$div(
          numericInput("weight_val", "Weight (lb):", value = NA, min = 50, max = 500, step = 0.1),
          style = "flex: 1; margin-bottom: 5px;"
        )
      )
    }
  })
  
  # Reset functionality
  observeEvent(input$reset, {
    updateNumericInput(session, "age", value = NA)
    updateNumericInput(session, "bmi", value = NA)
    updateRadioButtons(session, "units", selected = "metric")
    session$sendCustomMessage("resetHtWt", list())
    updateNumericInput(session, "cps", value = NA)
    updateRadioButtons(session, "cps_type", selected = "fasting")
    updateRadioButtons(session, "sex", selected = "na")
    updateRadioButtons(session, "dka", selected = "0")
    updateRadioButtons(session, "fh", selected = "na_none")
    updateRadioButtons(session, "gad", selected = "na")
    updateRadioButtons(session, "ia2", selected = "na")
    updateRadioButtons(session, "zinc", selected = "na")
    updateRadioButtons(session, "iaa", selected = "na")
    updateSelectInput(session, "ethnicity", selected = "Western Europe")
    
    # Clear cache
    nomogram_cache$last_key <- NULL
    nomogram_cache$data <- data.frame(age = numeric(0), bmi = numeric(0), probability = numeric(0))
  })
  
  # BMI calculation from height/weight
  observeEvent(input$use_hw, {
    req(input$height_val, input$weight_val)
    
    if (is.na(input$height_val) || is.na(input$weight_val) || 
        input$height_val <= 0 || input$weight_val <= 0) {
      showNotification("Please enter valid height and weight values", type = "error")
      return()
    }
    
    # Calculate BMI based on units
    if (input$units == "metric") {
      if (input$height_val < 50 || input$height_val > 250) {
        showNotification("Height must be 50–250 cm", type = "error")
        return()
      }
      if (input$weight_val < 20 || input$weight_val > 250) {
        showNotification("Weight must be 20–250 kg", type = "error")
        return()
      }
      h_meters <- input$height_val / 100
      w_kg <- input$weight_val
    } else {
      if (input$height_val < 20 || input$height_val > 100) {
        showNotification("Height must be 20–100 in", type = "error")
        return()
      }
      if (input$weight_val < 50 || input$weight_val > 500) {
        showNotification("Weight must be 50–500 lb", type = "error")
        return()
      }
      h_meters <- input$height_val * 0.0254
      w_kg <- input$weight_val * 0.453592
    }
    
    if (h_meters <= 0) {
      showNotification("Invalid height calculation", type = "error")
      return()
    }
    
    calculated_bmi <- round(w_kg / (h_meters * h_meters), 1)
    
    # Validate calculated BMI
    if (input$ethnicity == "Western Europe") {
      if (calculated_bmi < 17.5 || calculated_bmi > 70) {
        showNotification("Calculated BMI is outside valid range for Western Europe (17.5-70)", type = "warning")
      }
    } else if (input$ethnicity %in% names(region_coef)) {
      if (calculated_bmi < 15 || calculated_bmi > 60) {
        showNotification("Calculated BMI is outside valid range for this region (15-60)", type = "warning")
      }
    }
    
    updateNumericInput(session, "bmi", value = calculated_bmi)
  })
  
  # Pre-lab prediction output
  output$prelab_prediction <- renderText({
    if (is.na(input$age) || is.na(input$bmi)) {
      return("Please enter both age and BMI.")
    }
    if (input$age < 15 || input$age > 60) {
      return("Age must be between 15 and 60.")
    }
    
    # BMI validation by ethnicity
    if (input$ethnicity == "Western Europe") {
      if (input$bmi < 17.5 || input$bmi > 70) {
        return("BMI must be between 17.5 and 70 for Western Europe.")
      }
    } else if (input$ethnicity %in% names(region_coef)) {
      if (input$bmi < 15 || input$bmi > 60) {
        return("BMI must be between 15 and 60 for this region.")
      }
    }
    
    # Check for unimplemented regions
    if (input$ethnicity == "Middle East/North Africa") {
      return(paste(input$ethnicity, "region is not yet implemented. Please select a different region."))
    }
    
    # Calculate pre-lab prediction
    prob <- prelab_prediction_prob()
    if (is.na(prob)) {
      return("Unable to calculate prediction. Please check your inputs.")
    }
    
    # Format probability with special rounding rules
    prob_pct <- prob * 100
    if (prob_pct < 0.5) {
      prob_display <- "<1"
    } else if (prob_pct >= 99.5) {
      prob_pct <- min(prob_pct, 99.9)
      prob_display <- sprintf("%.1f", prob_pct)
    } else {
      prob_display <- sprintf("%.0f", prob_pct)
    }

    # Add DKA disclaimer if DKA is positive
    dka_disclaimer <- ""
    if (!is.na(input$dka) && input$dka == "1") {
      dka_disclaimer <- " (not accounting for DKA)"
    }

    return(paste0(prob_display, "%", dka_disclaimer))
  })
  
  # C-peptide parameters display for Model Information tab
  output$cpeptide_params_display <- renderText({
    t1d_alpha <- params$cpeptide$fasting$t1d$alpha
    t1d_lambda <- params$cpeptide$fasting$t1d$lambda
    t2d_alpha <- params$cpeptide$fasting$t2d$alpha
    t2d_lambda <- params$cpeptide$fasting$t2d$lambda
    
    paste0("Current fasting C-peptide parameters: T1D α=", t1d_alpha, ", λ=", t1d_lambda, 
           "; T2D α=", t2d_alpha, ", λ=", t2d_lambda, ". ")
  })
  
  # Autoantibody likelihood ratios display for Model Information tab
  output$antibody_lr_display <- renderText({
    lr_data <- params$lr
    
    # Create formatted string for each antibody
    lr_strings <- mapply(function(biomarker, lr_pos, lr_neg) {
      paste0(biomarker, ": LR+ = ", lr_pos, ", LR− = ", lr_neg)
    }, lr_data$biomarker, lr_data$LR_pos, lr_data$LR_neg)
    
    paste0("Current autoantibody likelihood ratios: ", paste(lr_strings, collapse = "; "), ".")
  })
  
  # Extrapolation annotations
  #
  # The tool still computes and displays an estimate for every accepted input.
  # These notes flag the two situations where the inputs sit outside the data
  # the underlying estimates were derived from: (1) an age at diagnosis beyond
  # the range of the selected region's dataset, and (2) a C-peptide value
  # entered for a person older than the South Asia cohort that supplied the
  # gamma distributions. Both can fire at once; they say different things.
  output$extrapolation_notes <- renderUI({
    if (is.na(input$age) || is.na(input$bmi)) return(NULL)
    if (!input$ethnicity %in% names(region_coef)) return(NULL)

    notes <- list()

    rng <- region_age_range[[input$ethnicity]]
    if (!is.null(rng) && (input$age < rng$min || input$age > rng$max)) {
      notes <- c(notes, list(sprintf(
        paste0("You entered %s and age at diagnosis %s years. This model was ",
               "derived from a %s dataset spanning ages at diagnosis of %d to ",
               "%d years, so the estimate above extrapolates the effect of age ",
               "beyond the range of the data."),
        input$ethnicity, format(input$age), input$ethnicity, rng$min, rng$max
      )))
    }

    if (!is.na(input$cps) && is.finite(input$cps) && input$cps > 0 &&
        input$age > CPEPTIDE_AGE_LIMIT) {
      notes <- c(notes, list(sprintf(
        paste0("You entered a C-peptide value for a person older than %d years. ",
               "The C-peptide likelihood ratios were derived from a cohort ",
               "diagnosed between ages 15 and %d years and do not condition on ",
               "age at diagnosis. In older-onset diabetes they understate the ",
               "evidence for type 1 diabetes at C-peptide values around ",
               "0.3 nmol/L."),
        CPEPTIDE_AGE_LIMIT, CPEPTIDE_AGE_LIMIT
      )))
    }

    if (length(notes) == 0) return(NULL)

    tags$div(
      tags$strong("Note:"),
      tags$ul(
        lapply(notes, tags$li),
        style = "margin: 4px 0 0 0; padding-left: 18px;"
      ),
      style = paste0("background-color: #fffaf0; border-left: 4px solid #dd6b20; ",
                     "padding: 8px 10px; border-radius: 5px; margin-bottom: 12px; ",
                     "font-size: 13px; color: #7b341e;")
    )
  })

  # Main prediction output
  output$prediction <- renderText({
    if (is.na(input$age) || is.na(input$bmi)) {
      return("Please enter both age and BMI.")
    }
    if (input$age < 15 || input$age > 60) {
      return("Age must be between 15 and 60.")
    }
    
    # BMI validation by ethnicity
    if (input$ethnicity == "Western Europe") {
      if (input$bmi < 17.5 || input$bmi > 70) {
        return("BMI must be between 17.5 and 70 for Western Europe.")
      }
    } else if (input$ethnicity %in% names(region_coef)) {
      if (input$bmi < 15 || input$bmi > 60) {
        return("BMI must be between 15 and 60 for this region.")
      }
    }
    
    # C-peptide validation
    if (!is.na(input$cps)) {
      if (input$cps < 0 || input$cps > 10) {
        return("C-peptide must be between 0 and 10 nmol/L.")
      }
    }
    
    # Check for unimplemented regions
    if (input$ethnicity == "Middle East/North Africa") {
      return(paste(input$ethnicity, "region is not yet implemented. Please select a different region."))
    }
    
    # Calculate prediction
    prob <- prediction_prob()
    if (is.na(prob)) {
      return("Unable to calculate prediction. Please check your inputs.")
    }
    
    # Format probability with special rounding rules
    prob_pct <- prob * 100
    if (prob_pct < 0.5) {
      prob_display <- "<1"
    } else if (prob_pct >= 99.5) {
      prob_pct <- min(prob_pct, 99.9)
      prob_display <- sprintf("%.1f", prob_pct)
    } else {
      prob_display <- sprintf("%.0f", prob_pct)
    }

    # Check for special cases
    dka_positive <- !is.na(input$dka) && input$dka == "1"
    any_antibody_positive <- (!is.na(input$gad) && input$gad == "1") || 
                            (!is.na(input$ia2) && input$ia2 == "1") || 
                            (!is.na(input$zinc) && input$zinc == "1") ||
                            (!is.na(input$iaa) && input$iaa == "1")
    
    if (dka_positive && any_antibody_positive) {
      # Calculate probability without DKA for comparison
      prob_without_dka <- tryCatch({
        p <- processed_inputs()
        calc_t1d_prob(input$age, input$bmi, p$male, p$fh, p$gad_status, 
                      p$ia2_status, p$zinc_status, p$iaa_status, p$ethnicity, p$cps, p$cps_type)
      }, error = function(e) NA_real_)
      
      if (!is.na(prob_without_dka)) {
        prob_without_dka_pct <- prob_without_dka * 100
        if (prob_without_dka_pct < 0.5) {
          prob_without_display <- "<1"
        } else if (prob_without_dka_pct >= 99.5) {
          prob_without_dka_pct <- min(prob_without_dka_pct, 99.9)
          prob_without_display <- sprintf("%.1f", prob_without_dka_pct)
        } else {
          prob_without_display <- sprintf("%.0f", prob_without_dka_pct)
        }
        prediction_text <- paste0("DKA & antibody positive, prob(T1D) > 99% (Without DKA, model predicts ", prob_without_display, "%)")
      } else {
        prediction_text <- paste0("DKA & antibody positive, prob(T1D) > 99% (Model predicts ", prob_display, "%)")
      }
    } else if (dka_positive && !any_antibody_positive) {
      # DKA positive but no antibodies positive
      prediction_text <- paste0("Predicted probability of Type 1 diabetes: ", prob_display, "% (not accounting for DKA)")
    } else if (!is.na(input$gad) && !is.na(input$ia2) && input$gad == "1" && input$ia2 == "1") {
      prediction_text <- paste0("GADA & IA-2 positive, prob(T1D) > 99% (Model predicts ", prob_display, "%)")
    } else if (!is.na(input$gad) && !is.na(input$zinc) && input$gad == "1" && input$zinc == "1") {
      prediction_text <- paste0("GADA & ZnT8 positive, prob(T1D) > 99% (Model predicts ", prob_display, "%)")
    } else if (!is.na(input$gad) && !is.na(input$iaa) && input$gad == "1" && input$iaa == "1") {
      prediction_text <- paste0("GADA & IAA positive, prob(T1D) > 99% (Model predicts ", prob_display, "%)")
    } else if (!is.na(input$ia2) && !is.na(input$zinc) && input$ia2 == "1" && input$zinc == "1") {
      prediction_text <- paste0("IA-2 & ZnT8 positive, prob(T1D) > 99% (Model predicts ", prob_display, "%)")
    } else if (!is.na(input$ia2) && !is.na(input$iaa) && input$ia2 == "1" && input$iaa == "1") {
      prediction_text <- paste0("IA-2 & IAA positive, prob(T1D) > 99% (Model predicts ", prob_display, "%)")
    } else if (!is.na(input$zinc) && !is.na(input$iaa) && input$zinc == "1" && input$iaa == "1") {
      prediction_text <- paste0("ZnT8 & IAA positive, prob(T1D) > 99% (Model predicts ", prob_display, "%)")
    } else {
      prediction_text <- paste0("Predicted probability of Type 1 diabetes: ", prob_display, "%")
    }
    
    return(prediction_text)
  })
  
  # Nomogram plot
  output$nomogram_plot <- renderPlotly({
    grid <- nomogram_data()
    
    if (nrow(grid) == 0) {
      p <- ggplot() + 
        geom_text(aes(x = 0, y = 0, label = "Unable to generate nomogram"), size = 5) +
        theme_void()
      return(ggplotly(p))
    }
    
    # Prepare contour data
    xs <- sort(unique(grid$age))
    ys <- sort(unique(grid$bmi))
    
    if (length(xs) == 0 || length(ys) == 0) {
      p <- ggplot() + 
        geom_text(aes(x = 0, y = 0, label = "Insufficient data for nomogram"), size = 5) +
        theme_void()
      return(ggplotly(p))
    }
    
    # Contour levels and labels
    levels <- c(0.10, 0.25, 0.33, 0.50, 0.67, 0.75, 0.90)
    lbl_df <- data.frame(x = numeric(0), y = numeric(0), label = character(0))
    
    # Generate contour labels with error handling
    tryCatch({
      zm <- matrix(grid$probability, nrow = length(ys), ncol = length(xs))
      
      if (!any(is.na(zm)) && !any(!is.finite(zm)) && 
          min(zm) < max(zm) && length(unique(as.vector(zm))) > 1) {
        
        cl <- contourLines(xs, ys, zm, levels = levels)
        if (length(cl) > 0) {
          lbl_df <- do.call(rbind, lapply(cl, function(l) {
            if (length(l$x) > 0) {
              m <- max(1, round(length(l$x) / 2))
              level_pct <- l$level * 100
              label_text <- if (level_pct < 0.5) {
                "<1%"
              } else if (level_pct >= 99.5) {
                paste0(sprintf("%.1f", min(level_pct, 99.9)), "%")
              } else {
                paste0(round(level_pct), "%")
              }
              data.frame(x = l$x[m], y = l$y[m],
                         label = label_text,
                         stringsAsFactors = FALSE)
            } else {
              data.frame(x = numeric(0), y = numeric(0), label = character(0),
                         stringsAsFactors = FALSE)
            }
          }))
        }
      }
    }, error = function(e) {})
    
    # Set axis limits based on ethnicity
    y_min <- if (input$ethnicity == "Western Europe") 17.5 else 15
    y_max <- if (input$ethnicity == "Western Europe") 50 else 60
    
    # Create base plot
    gp <- ggplot(grid, aes(age, bmi)) +
      geom_tile(aes(fill = probability)) +
      scale_fill_gradient(name = "T1D\nProbability", labels = scales::percent,
                          low = "white", high = "darkred") +
      labs(title = "Risk surface with iso-risk contours",
           x = "Age (years)", y = "BMI (kg/m²)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5)) +
      coord_cartesian(ylim = c(y_min, y_max))
    
    # Add contours and labels if available
    if (nrow(lbl_df) > 0) {
      tryCatch({
        gp <- gp + 
          geom_contour(aes(z = probability), breaks = levels, colour = "white", linewidth = 0.5) +
          geom_text(data = lbl_df, aes(x, y, label = label),
                    colour = "black", fontface = "bold",
                    size = 3, inherit.aes = FALSE)
      }, error = function(e) {})
    }
    
    # Add patient point if requested and valid
    if (input$show_patient_point && !is.na(input$age) && !is.na(input$bmi) &&
        is.finite(input$age) && is.finite(input$bmi) && 
        input$age >= 15 && input$age <= 60) {

      within_bounds <- if (input$ethnicity == "Western Europe") {
        input$bmi >= 17.5 && input$bmi <= 70
      } else if (input$ethnicity %in% names(region_coef)) {
        input$bmi >= 15 && input$bmi <= 60
      } else {
        FALSE
      }

      if (within_bounds) {
        gp <- gp + annotate("point", x = input$age, y = input$bmi,
                            colour = "yellow", size = 3, shape = 16)
      }
    }

    # Convert to plotly with error handling
    tryCatch({
      suppressMessages(suppressWarnings({
        ggplotly(gp, tooltip = c("x", "y", "fill"), height = 500) %>%
          layout(hovermode = "closest")
      }))
    }, error = function(e) {
      tryCatch({
        suppressMessages(suppressWarnings(ggplotly(gp, height = 500)))
      }, error = function(e2) {
        p <- ggplot() +
          geom_text(aes(x = 0, y = 0, label = "Plot rendering error"), size = 5) +
          theme_void()
        ggplotly(p)
      })
    })
  })
}

shinyApp(ui, server)
