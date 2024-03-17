# Load necessary packages
library(shiny)
library(dplyr)
library(data.tree)
library(collapsibleTree)
library(shinycssloaders)
library(ggtree)
library(tidyverse)

setwd("/Users/Elijah/Documents/BINP29_Files/Population_Genetics_Project")

# Load data
dat <- readxl::read_xlsx("AADR_Annotation.xlsx")
col_Names <- c("S/N", "Genetic_ID", "Master_ID", "Skeletal_code", "Skeletal_element", 
               "Publication_Year", "Publication","Date_Determination_Method", "Date_mean_BP", 
               "Date_standard_deviation_BP", "Full_Date","Age_at_Death", "Group_ID", "Locality", 
               "Political_Entity", "Latitude", "Longitude","Pulldown_Strategy", "Data_source", 
               "No_Libraries", "1240k_coverage", "Autosomal_SNPs_hit", "Autosomal_SNPs_hit_HO", 
               "Molecular_Sex", "Family_ID_and_Position", "Y_haplogroup_terminal", "Y_haplogroup_ISOGG", 
               "mtDNA_coverage", "mtDNA_haplogroup", "mtDNA_match_to_consensus", "Damage_rate", 
               "Sex_ratio","Library_type", "Libraries", "Assessment", "Assessment_Warnings")
colnames(dat) <- col_Names


# Filter the dataframe to remove missing entries and special characters (*, :, (, ), trailing spaces, ~, and -) from the entries
selected_data <- dat %>%
  select(Date_mean_BP, Political_Entity, Y_haplogroup_ISOGG) %>%
  filter(!Y_haplogroup_ISOGG %in% c("n/a (female)", "..", "n/a  (sex unknown)", "n/a (Sex Unknown)", 
                                    "n/a (not published)", "n/a", "na")) %>%
  mutate(Y_haplogroup_ISOGG = gsub("\\(.*|~.*|-.*|;.*|\\*.*|\\s+", "", Y_haplogroup_ISOGG))


# Read in ancestral tree file
tree_file <- readr::read_tsv("chrY_hGrpTree_isogg2016.txt")
colnames(tree_file) <- c("child", "parent")
tree <- FromDataFrameNetwork(tree_file)

# Define UI function
ui <- fluidPage(
  titlePanel("The Haplogroup Storyteller"),
  
  sidebarLayout(
    sidebarPanel(
      # Prompt user to provide input
      textInput("haplogroup", "Enter Your Haplogroup (ISOGG Format):"),
      # Submit BUTTON
      actionButton("submit_button", "Submit", style = "background-color: #007bff; color: #fff; border: none; margin-top: 10px;"),
      # Reset BUTTON
      actionButton("reset_button", "Reset", style = "background-color: #dc3545; color: #fff; border: none; margin-top: 10px;"),
      br(),
      # Diaplay reports
      uiOutput("reports")
    ),
    mainPanel(
      # Title display
      uiOutput("haplogroup_title"),
      # Haplogroup message display
      uiOutput("haplogroup_message"),
      # Display the network graph
      collapsibleTreeOutput('ancestry_tree', height='700px') %>% withSpinner(color = "green")
    )
  )
)

generateAncestryTree <- function(haplogroup, selected_data, tree) {
  if (!is.null(haplogroup)) {
    node <- FindNode(tree, haplogroup)
    if (!is.null(node)) {
      parent_node <- ifelse(is.null(node$parent), "", node$parent$name)
      child_nodes <- ifelse(is.null(node$children), character(0), sapply(node$children, function(x) x$name))
      tree_data <- data.frame(name = c(haplogroup, parent_node, child_nodes))
      links <- data.frame(source = c(haplogroup, parent_node, rep(parent_node, length(child_nodes))),
                          target = c(parent_node, haplogroup, child_nodes))
      #return(dendroNetwork(tree_data, linkData = links))
    }
  }
  return(NULL)
}

server <- function(input, output, session) {
  
  # Filter data based on user input
  filtered_data <- eventReactive(input$submit_button, {
    subset(selected_data, Y_haplogroup_ISOGG == input$haplogroup)
  })
  
  # Find the parent haplogroup
  parent_haplogroup <- eventReactive(input$submit_button, {
    node <- FindNode(tree, input$haplogroup)
    if (!is.null(node)) {
      parent_node <- ifelse(is.null(node$parent), "", node$parent$name)
      return(parent_node)
    } else {
      return(NULL)
    }
  })
  
  # Find the child haplogroup 
  child_haplogroup <- eventReactive(input$submit_button, {
    node <- FindNode(tree, input$haplogroup)
    if (!is.null(node)) {
      child_node <- ifelse(is.null(node$children), character(0), sapply(node$children, function(x) x$name))
      return(child_node)
    } else {
      return(NULL)
    }
  })
  
  # Find the most ancient occurrence of the parent haplogroup
  oldest_date <- eventReactive(input$submit_button, {
    if (!is.null(parent_haplogroup())) {
      values <- selected_data$Date_mean_BP[selected_data$Y_haplogroup_ISOGG == parent_haplogroup()]
      if (length(values) > 0) {
        max(values, na.rm = TRUE) 
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  })
  # Find the most recent occurrence of the parent haplogroup
  newest_date <- eventReactive(input$submit_button, {
    if (!is.null(parent_haplogroup())) {
      values <- selected_data$Date_mean_BP[selected_data$Y_haplogroup_ISOGG == parent_haplogroup()]
      if (length(values) > 0) {
        min(values, na.rm = TRUE)
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  })
  
  # Find the countries of the descendants of the next paternal haplogroup in the lineage
  countries_of_DNA_tested_descendants <- eventReactive(input$submit_button, {
    country_message <- NULL  # Initialize to NULL
    
    # Check if the child haplogroup is not null and is found in the database
    if (!is.null(child_haplogroup()) && child_haplogroup() %in% selected_data$Y_haplogroup_ISOGG) {
      filtered_countries <- selected_data$Political_Entity[selected_data$Y_haplogroup_ISOGG == child_haplogroup()]
      countries <- unique(filtered_countries)
      
      # Find the list of countries in the database
      if (length(countries) > 1) {
        country_message <- paste(paste(countries[-length(countries)], collapse = ", "), "and", countries[length(countries)])
      } else if (length(countries) == 1) {
        country_message <- countries
      }
      # No need to handle length(countries) == 0 explicitly, it will return NULL by default
    }
    
    return(country_message)  # Return list of countries or NULL
  })
  

  # Find the countries of the descendants if the query haplogroup is the last paternal line in the lineage
  DNA_tested_descendants_of_parent_haplogroup <- eventReactive(input$submit_button, {
    #Check if the child haplogroup is not null and is found in the database
    if (!is.null(parent_haplogroup()) && parent_haplogroup() %in% (selected_data$Y_haplogroup_ISOGG)) {
      filtered_countries <- selected_data$Political_Entity[selected_data$Y_haplogroup_ISOGG == parent_haplogroup()]
      countries <- unique(filtered_countries)
      
      # Find the list of countries in the database
      if (length(countries) > 1) {
        country_message <- paste(paste(countries[-length(countries)], collapse = ", "), "and", countries[length(countries)])
      } else if (length(countries) == 1) {
        country_message <- countries
      } else if (length(countries) == 0) {
        return(NULL)
      }
    } else {
      return(NULL)
    }
    return(country_message) # return list of countries
  })
  
  
  # Count the number of countries associated with the child haplogroup
  num_countries <- eventReactive(input$submit_button, {
    if (!is.null(child_haplogroup())) {
      # Filter the selected_data to include only entries matching the child_haplogroup
      filtered_countries <- selected_data$Political_Entity[selected_data$Y_haplogroup_ISOGG == child_haplogroup()]
      # Calculate the total number of occurrences of countries where the haplogroup is found
      total_occurrences <- length(filtered_countries)
      return(total_occurrences)
    } else {
      return(NULL)
    }
  })
  
  
  # Count the number of countries associated with the parent haplogroup if it's the last paternal line
  num_countries2 <- eventReactive(input$submit_button, {
    if (!is.null(parent_haplogroup())) {
      # Filter the selected_data to include only entries matching the child_haplogroup
      filtered_countries <- selected_data$Political_Entity[selected_data$Y_haplogroup_ISOGG == parent_haplogroup()]
      # Calculate the total number of occurrences of countries where the haplogroup is found
      total_occurrences <- length(filtered_countries)
      return(total_occurrences)
    } else {
      return(NULL)
    }
  })
  
  # Generate reports
  output$reports <- renderUI({
    # Check if the input haplogroup exists in the database
    
    if (is.null(input$haplogroup) || input$haplogroup == "") {
      # Tell the user to enter a valid haplogroup ID
      return(HTML( "<span style='color: red;'>", "Please provide a haplogroup ID", "</span>"))
      
    } else if (!is.null(input$haplogroup) && !(input$haplogroup %in% filtered_data()$Y_haplogroup_ISOGG)) {
      # Display initial report message
      return(
        tagList(
          br(),
          div(
            "Sorry, no data was found for your haplogroup. We're constantly improving the database as new data becomes available. You may type in a different haplogroup ID, or try again some other time. Thank you!",
            style = "font-size: 16px; margin-bottom: 10px; background-color: #f8f9fa; 
              border-radius: 10px; padding: 10px; box-shadow: 0 4px 8px 0 rgba(0,0,0,0.2);"
          )
        )
      )
    } else {
      # Generate usual reports if the input haplogroup exists in the database
      reports <- list()
      #reports[[1]] <- 
      if (!is.null(parent_haplogroup()) && !is.null(oldest_date())) {
        reports[[1]] <- HTML(paste0("Your haplogroup's paternal line was formed when it branched off from the ancestor ", "<span style='color: navy;'><strong>",
                                    parent_haplogroup(), "</strong></span>", " around the year ", "<span style='color: navy;'><strong>", oldest_date(), " BCE.", "</strong></span>"))
      } else if (is.null(parent_haplogroup()) || is.null(oldest_date())) {
        reports[[1]] <- HTML(paste0("The exact date your haplogroup's paternal line was formed is not known, but available data shows that it branched off from the ancestor ", "<span style='color: navy;'><strong>",
                                    parent_haplogroup(), "</strong></span>"))
      }
      if (!is.null(parent_haplogroup()) & !is.null(newest_date())) {
        reports[[2]] <- HTML(paste0("The man who is the most recent common ancestor of this lineage is estimated to have been around the year ", "<span style='color: navy;'><strong>",
                                    newest_date(), " BCE.", "</strong></span>"))
      }
      # Compare the database with the parent haplogroup in the tree and find the number of descendants, and their countries
      # Report 3
      if (!is.null(child_haplogroup()) && !is.null(countries_of_DNA_tested_descendants()) && !is.null(num_countries())) {
        reports[[3]] <- HTML(paste0("He is the ancestor of at least ", "<span style='color: navy;'><strong>", num_countries(), "</strong></span>", " descendant lineages known as ", "<span style='color: navy;'><strong>", child_haplogroup(), "</strong></span>", "."))
      } else if (!is.null(parent_haplogroup()) && !is.null(child_haplogroup())) {
        reports[[3]] <- HTML(paste0("He is the most recent paternal line ancestor of all members of this group, with ", "<span style='color: navy;'><strong>", child_haplogroup(), "</strong></span>", " being the most recent known descendant group."))
      }
      
      # Report 4
      if (!is.null(child_haplogroup()) && !is.null(countries_of_DNA_tested_descendants()) && !is.null(num_countries())) {
        reports[[4]] <- HTML(paste0("There are ", "<span style='color: navy;'><strong>", num_countries(), "</strong></span>", " DNA-tested descendants, and they specified that their earliest known origins are from: ", "<span style='color: navy;'><strong>",
                                    countries_of_DNA_tested_descendants(), "</strong></span>", "."))
        # Fifth report
        reports[[5]] <- HTML(paste0("But the story does not end here. To find out which of the ", "<span style='color: navy;'><strong>", child_haplogroup(), "</strong></span>", " branches and descendants' country or countries ", " you're directly linked to, consider doing a comprehensive DNA test!"))
      } else if (!is.null(DNA_tested_descendants_of_parent_haplogroup()) && !is.null(num_countries2())) {
        reports[[4]] <- HTML(paste0("While we don't have records for any DNA-tested descendants from ", "<span style='color: navy;'><strong>", child_haplogroup(), "</strong></span>", ", there are at least ", "<span style='color: navy;'><strong>", num_countries2(), "</strong></span>", " DNA-tested descendants of ", 
                                    "<span style='color: navy;'><strong>", parent_haplogroup(), "</strong></span>", " from: ", "<span style='color: navy;'><strong>", DNA_tested_descendants_of_parent_haplogroup(), "</strong></span>", "."))
        # Fifth report becomes
        reports[[5]] <- HTML(paste0("But the story does not end here. To find out which of the ", "<span style='color: navy;'><strong>", parent_haplogroup(), "</strong></span>", " branches and descendants' country or countries ", " you're directly linked to, consider doing a comprehensive DNA test!"))
      } else {
        reports[[4]] <- "There are no DNA-tested descendants associated with this haplogroup in the database."
      }
      
      lapply(reports, function(report) {
        tagList(
          br(),
          div(report, style = "font-size: 18px; margin-bottom: 5px; background-color: #f8f9fa; animation: mymove 5s infinite; 
              border-radius: 10px; padding: 10px; box-shadow: 0 4px 8px 0 rgba(0,0,0,0.2);")
        )
      })
    }
  })
  # Show user the message header along with the input haplogroup
  output$haplogroup_title <- eventReactive(input$submit_button, {
    if (is.null(input$submit_button) || input$haplogroup == "") {
      # Tell the user to enter a valid haplogroup ID
      return("")
    } else if (!(input$submit_button %in% filtered_data()$Y_haplogroup_ISOGG)) {
      # Display initial report message
      return(HTML(paste("<span style='font-size: 25px;'><strong>", "Your Haplogroup Story:", "<span style='color: blue; font-weight: bold;'>", input$haplogroup, "</span>", "</span></strong>")))
    }
  })
  # Show user the additional message behind how haplogroup ancestry is obtained
  output$haplogroup_message <- eventReactive(input$submit_button, {
    if (is.null(input$submit_button) || input$haplogroup == "") {
      # Tell the user to enter a valid haplogroup ID
      return("")
    } else if (!(input$submit_button %in% filtered_data()$Y_haplogroup_ISOGG)) {
      # Display initial report message
      HTML("<span style='font-size: 16px; text-align: justify;'>","The Y chromosome is passed from father to son remaining mostly unaltered across generations, except for small traceable changes in the DNA. By tracing these changes, we constructed a family tree of humankind where all male lineages trace back to a single common ancestor who lived hundreds of thousands of years ago. This human tree allows us to explore lineages through time and place and to uncover the modern history of your direct paternal surname line and the ancient history of our shared ancestors.",
           "</span>")
    }
  })
 
  # Generate ancestry tree
  output$ancestry_tree <- renderCollapsibleTree({
    generateAncestryTree(input$haplogroup, selected_data, tree)
  })  
  
  # Reset button
  observeEvent(input$reset_button, {
    # Reset text input
    updateTextInput(session, "haplogroup", value = "")
    # Clear outputs in the main panel
    output$haplogroup_title <- renderText(NULL)
    output$haplogroup_message <- renderText(NULL)
    #output$ancestry_tree <- renderText(NULL)
  })
  
}


shinyApp(ui = ui, server = server)

