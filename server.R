library(scales)

shinyServer(function(input, output, session) {

# Data Sourcing -----------------------------------------------------------
    
    # Changing the file to load here will allow this app to run with any
    # proper Phyloseq data object
    # The contents of the rds file should be exactly 1 phyloseq object
    # recommendations for the phyloseq object:
    #     the tax table should have row names corresponding to the ASV sequence
    #     the tax table should include a column "SVHash" that contains a unique hash for each ASV sequence
    #         see the function "seqToHash64Base62" in accessory_scripts.R for an example that you can use
    # if you have difficulty, examine the phyloseq object included in niv_data.rds and make your own object match it
  
    inputPhy <- reactive({
        # this sample data includes only the NIV Stations, not the controls
        return(readRDS(file.path("niv_data.rds")))
    })

# Samples Tab -------------------------------------------------------------

    observeEvent(
        {
            inputPhy()
        },
        {
            #browser()
            updateSelectInput(
                session,
                "samplesFeatureSelect",
                choices = sample_variables(inputPhy()),
                selected = sample_variables(inputPhy())
            )
        }
    )
    
    ignore_strains <- reactive({
      # if you have ASVs that are universally uninteresting for your analysis, you can put them here
      # alternatively, you could build a simple user input to exclude certain ASVs
      # place any strains to ignore in all analysis in this list
      strains_to_ignore <- c("fakehash")
      
      strains_to_ignore
    })

    featureFilteredPhy <- reactive({
        req(inputPhy(), input$samplesFeatureSelect)
        xtry <- try({
            x <- inputPhy()
            
            # caution: this filtering step is slow
            if(input$filterSamples) {
              sample_data(x) <- sample_data(x)[, as.character(input$samplesFeatureSelect), drop = FALSE]
            }
            
        })

        if(inherits(xtry, "try-error")) {
            return(NULL)
        } else {
          
            # if the user wants to remove strains, do it here because everything else in the app depends on this
            # note that this is the most downstream phyloseq object that the differential analysis uses
            if (input$ignoreStrains == TRUE) {
              
              withProgress(message = "Removing strains",
                           detail = "Please wait...",
                           value = 0.25,
                           {
                            ignore_strain_hashes <- ignore_strains()
              
                            all_taxa_hashes <- as.vector(tax_table(x)[,'SVhash'])
                            all_taxa_seqs <- as.vector(row.names(tax_table(x)))
                            all_taxa_df <- data.frame(all_taxa_seqs)

                            row.names(all_taxa_df) <- all_taxa_hashes
              
                            remaining_taxa_hashes <- setdiff(all_taxa_hashes, ignore_strain_hashes)
                            remaining_taxa_seqs <- as.vector(all_taxa_df[remaining_taxa_hashes, ])
              
                            pmps <- prune_taxa(remaining_taxa_seqs, x)
                           })
              
              return(pmps)
            } else {
              return(x)
            }
            
        }
    })
    
    featureFilteredPhy_ordered <- reactive({
      #this step can be slow with large datasets
      sampleData = sample_data(featureFilteredPhy())
      sampleData_ordered_cols <- sampleData[ , order(names(sampleData))]
      setcolorder(sampleData_ordered_cols, c('sample_name', 'run', 'date'))

      sampleData_ordered_cols
    })
  
    # if you wish, for whatever reason, to have some samples pre-selected, the code
    # this commented out below demonstrates how
    output$samplesDT <- DT::renderDataTable(
        data.frame(featureFilteredPhy_ordered()),
        class = 'compact stripe',
        server = FALSE,
        rownames = FALSE,
        filter = list(position = "top"),
        extensions = c('Buttons'),
        options = list(
            pageLength = 20,
            lengthMenu = c(10, 20, 40, 200, 1000),
            dom = 'Blfrtip',
            scrollX = TRUE,
            searchHighlight = TRUE
        ),
        selection = list(
            target = 'row' #,
            #selected = 1:(min(3, nsamples(featureFilteredPhy()))),
            #selected = NULL
        )
    )

    samplesDT_proxy <- DT::dataTableProxy("samplesDT", deferUntilFlush = TRUE)

    # these observeEvent statements allow the user to easily select a large subset of samples
    observeEvent(
        input$samplesDT_filtered,
        DT::selectRows(samplesDT_proxy, input$samplesDT_rows_all)
    )

    observeEvent(
        input$samplesDT_visible,
        DT::selectRows(samplesDT_proxy, input$samplesDT_rows_current)
    )

    observeEvent(
        input$samplesDT_reset,
        DT::selectRows(samplesDT_proxy, NULL)
    )
    
    sampleFilteredPhy <- eventReactive(
      {
        input$samplesDT_apply
      },
      {
      req(featureFilteredPhy(), input$samplesDT_rows_selected)
      
      xtry <- try(
        x <- prune_samples(sample_names(featureFilteredPhy())[input$samplesDT_rows_selected], featureFilteredPhy())
      )
      
      if(inherits(xtry, "try-error")) {
        return(NULL)
      } else {
        return(x)
      }
      
    })
    
    # Because the user can select all filtered samples with one click, it is useful feedback for them to see the exact #
    # This number may differ from the number of samples that appears in the left pane because the left pane shows
    # the number of samples AFTER all subsequent filtering.
    
    number_of_samples_selected <- eventReactive(
      {
        input$samplesDT_apply
      },  
      {
      length(input$samplesDT_rows_selected)
      })
    
    output$n_samples_selected <- renderText({
      ns <- as.character(number_of_samples_selected())
      if(is.null(ns)) {
        ns <- '0'
      }
      paste0(ns, ' total samples selected from table')
    })
    


# Taxonomy Tab ------------------------------------------------------------

    observeEvent(
        {
            sampleFilteredPhy()
        },
        {
            x <- rank_names(sampleFilteredPhy())

            updateSelectInput(
                session,
                "taxonomyRankSelect",
                choices = x,
                selected = x[1]
            )
        }
    )

    observeEvent(
        {
            input$taxonomyRankSelect
        },
        {
            x <- sort(sub("^$", "<unclassified>", get_taxa_unique(sampleFilteredPhy(), taxonomic.rank = input$taxonomyRankSelect)))

            updateSelectInput(
                session,
                "taxonomyTaxaSelect",
                choices = x,
                selected = x
            )
        }
    )
    
    # filter the phyloseq object (remove ASVs) based on user-selected taxonomy restrictions
    taxaFilteredPhy <- reactive({
        req(sampleFilteredPhy(), input$taxonomyRankSelect, input$taxonomyTaxaSelect)

        xtry <- try(
            x <- prune_taxa(as.logical(tax_table(sampleFilteredPhy())[, input$taxonomyRankSelect] %in% sub("<unclassified>", "", as.character(input$taxonomyTaxaSelect))), sampleFilteredPhy())
        )

        if(inherits(xtry, "try-error")) {
            return(NULL)
        } else {
            return(x)
        }
    })


# Counts Tab --------------------------------------------------------------
    
    # filter the phyloseq object (remove samples) based on user-specified minimum number of reads per sample
    sampleCountFilteredPhy <- reactive({
        req(taxaFilteredPhy(), input$countFilteringSample)

        xtry <- try(
            x <- prune_samples(sample_sums(taxaFilteredPhy()) >= input$countFilteringSample, taxaFilteredPhy())
        )

        if(inherits(xtry, "try-error")) {
            return(NULL)
        } else {
            return(x)
        }
    })
    
    # filter the phyloseq object (remove ASVs) based on user-specified minimum number of reads per sample/samples
    seqCountFilteredPhy <- reactive({
        req(sampleCountFilteredPhy(), input$countFilteringA, input$countFilteringK)

        xtry <- try(
            x <- filter_taxa(sampleCountFilteredPhy(), filterfun(kOverA(input$countFilteringK, input$countFilteringA)), prune = TRUE)
        )

        if(inherits(xtry, "try-error")) {
            return(NULL)
        } else {
            return(x)
        }
    })
    
    # Transform the counts per ASV per sample based on user input.
    # This helps to normalize for differing coverage across samples.
    # Note that this is NOT rarefaction but rather dividing and rounding ASV counts.
    tPhy <- reactive({
        req(seqCountFilteredPhy(), input$transformationMethod)

        xtry <- try(
            x <- switch(
                EXPR = input$transformationMethod,
                none = seqCountFilteredPhy(),
                prevalence = transform_sample_counts(seqCountFilteredPhy(), function(x) {round(x / sum(x) * 10000)}),
                binary = transform_sample_counts(seqCountFilteredPhy(), function(x) {ifelse(x > 0, 1, 0)}),
                seqCountFilteredPhy()
            )
        )

        if(inherits(xtry, "try-error")) {
            return(NULL)
        } else {
            return(x)
        }
    })

    # this allows connecting the hash to the phylogeny of the taxa
    hash_dictionary <- reactive({
      d <- data.frame(tax_table(tPhy()), stringsAsFactors = FALSE)
      d[, 'sequence'] <- row.names(d)
      row.names(d) <- d[, 'SVhash']
      
      #form a new column with the combined phylogeny by joining the columns with semicolon
      p <- data.frame(row.names(d), paste(d$Kingdom, d$Phylum, d$Class, d$Order, d$Family, d$Genus, d$Species, sep = ';'))

      row.names(p) <- row.names(d)
      colnames(p) <- c('SVHash', 'phylogeny')
      
      p
    })
    
    # can be used to give brief details about the row that the user has selected
    selected_hash_table <- reactive({
      selected_col <- input$countsCountsDT_columns_selected
      s <- hash_dictionary()[selected_col,1:2]

      s
    })
    
    # This code allows the user to either include or exclude the phylogeny information from the counts table.
    # This is useful because the phylogeny is lengthy and can get in the way in the app, but may be useful, 
    # especially when the user wants to download the data.
    
    counts_without_phylogeny <- reactive({
      x <- otu_table(tPhy())@.Data
      rownames(x) <- tax_table(tPhy())[, "SVhash"]
      x
    })
    
    counts_with_phylogeny <- reactive({
      x <- otu_table(tPhy())@.Data
      rownames(x) <- tax_table(tPhy())[, "SVhash"]
      
      h <- data.frame(hash_dictionary(), stringsAsFactors = FALSE)
      hv = as.vector(h[,2])
      hvr <- gsub(";", "; ", hv)

      p <- cbind(hvr, x)
      colnames(p)[1] = 'phylogeny'
      
      p
    })
    
    counts_table_to_render <- reactive({
      if(input$counts_phylogeny_row == TRUE) {
        cttr <- counts_with_phylogeny()
      } else {
        cttr <- counts_without_phylogeny()
      }
      cttr
    })
    
    # this download function downloads ALL pages of the table, whereas the other buttons only download the visible rows
    output$countsTableDownload <- downloadHandler(
      filename = function(){"counts.csv"},
      content = function(fname){
        write.csv(counts_table_to_render(), fname)
      }
    )
    
    output$countsDownloadWarning <- renderText(
      "WARNING: The 3 buttons below (Copy, CSV, and Print) will export ONLY the visible rows on the current page of the table."
    )
    
    # Ideally, we would fix the phylogeny data (when visible) to not scroll horizontally,
    # but it is not working (although it doesn't break anything either).
    nFixedCols <- reactive({
      if(input$counts_phylogeny_row == TRUE) {
        nfc <- 2
      } else {
        nfc <- 1
      }
      nfc
    })
    
    output$countsCountsDT <- DT::renderDataTable(
      counts_table_to_render(),
      class = 'compact stripe',
      server = TRUE,
      rownames = TRUE,
      filter = list(position = "top"),
      extensions = c('Buttons', 'FixedColumns'),
      options = list(
        pageLength = 20,
        lengthMenu = c(10, 20, 40, 200, 1000),
        dom = 'Blfrtip',
        scrollX = TRUE,
        fixedColumns = list(leftColumns = nFixedCols()),
        searchHighlight = TRUE
      ),
      selection = list(
        target = 'row',
        mode = 'single'
      )
    )
  
    output$countsHashDT <- DT::renderDataTable(
        data.frame(hash = unname(tax_table(tPhy())[, "SVhash"]), sequence = rownames(tax_table(tPhy())), stringsAsFactors = FALSE),
        class = 'compact stripe',
        server = FALSE,
        rownames = FALSE,
        filter = list(position = "top"),
        extensions = c('Buttons', 'FixedColumns'),
        options = list(
            pageLength = 20,
            lengthMenu = c(10, 20, 40, 200, 1000),
            dom = 'Blfrtip',
            scrollX = TRUE,
            fixedColumns = TRUE,
            searchHighlight = TRUE
        ),
        selection = 'none'
    )
    
    output$selected_hash <- DT::renderDataTable(
      data.frame(selected_hash_table()),
      options = list(
        lengthChange = FALSE,
        searching = FALSE,
        ordering = FALSE,
        paging = FALSE,
        info = FALSE
      ),
      rownames = FALSE
    )

# Physeq Object -----------------------------------------------------------
    
    taxa_names_display <- reactiveValues()
    
    # Only compute the tree if the user wants to because it is slow
    # The tree becomes part of the phyloseq object and is necessary for tree-based distance metrics
    physeq <- reactive({
      req(tPhy())
      
      if (input$build_tree == TRUE) {
          withProgress(message = "Building tree",
                       detail = "This may take some time...",
                       value = 0.25,
                       {
                         tree <- compute_phylogenetic_tree(tPhy())
                         incProgress(0.4)
                         
                         taxa_names_display$node_names <- taxa_names(tree)
                         
                         t_corrected_taxa_names <- tree
                         taxa_names(t_corrected_taxa_names) <- taxa_names(otu_table(tPhy()))
                         
                         phy_tax <- tax_table(tPhy())
                         phy_sample_data <- sample_data(tPhy())
                         phy_otu <- otu_table(tPhy())
                         incProgress(0.2)
                         
                         phy_with_tree <- phyloseq(tax_table(phy_tax), sample_data(phy_sample_data), 
                                                   otu_table(phy_otu), phy_tree(t_corrected_taxa_names))
                       })
          
          return(phy_with_tree)
        } else {
          return(tPhy())
        }
    })

    output$samplesBox <- renderValueBox({
        valueBox(
            nsamples(physeq()),
            "Samples",
            icon = icon("list"),
            color = "light-blue"
        )
    })

    output$seqBox <- renderValueBox({
        valueBox(
            ntaxa(physeq()),
            "Sequence Variants",
            icon = icon("list"),
            color = "light-blue"
        )
    })

# Phylogenetic Tree Tab ---------------------------------------------------
    
    # It is tricky to create simple code for displaying a tree that works for all cases.
    # This is an imperfect compromise. Ideal trees will need to be generated and displayed with custom code
    # outside of a general-purpose Shiny app.
    
    tree_width <- reactive(
      800
    )
    
    tree_height <- reactive(
      10 * ntaxa(physeq())
    )
    
    output$ggtree_graphic <- renderPlot(
      {
      t <- phy_tree(physeq())
      
      # convert the names on the tree for ease of intepretation
      taxa_names(t) <- taxa_names_display$node_names 
      
      ft <- fortify(t) %>% dplyr::as_data_frame()   #unclear if this is useful
      p <- ggtree(ft) + geom_tiplab(size=4, align=T, linesize=0.3, hjust=0) + geom_treescale(width=0.5, color='white')
      # the geom_treescale is a hack - it draws a scale bar at the bottom with the side effect of managing the width of the plot
      # the scale bar is hidden by making it white
      
      p
    },
    width = function(){tree_width()},
    height = function(){tree_height()}
    )
    
    output$treeBox <- renderValueBox({
        valueBox(
            if(is.null(physeq()@phy_tree)) {"X"} else {"NJ"},
            "Phylogenetic Tree",
            icon = icon(if(is.null(physeq()@phy_tree)) {"thumbs-down"} else {"thumbs-up"}),
            color = if(is.null(physeq()@phy_tree)) {"red"} else {"green"}
        )
    })

# Ordination Tab ----------------------------------------------------------

    get_formula_ord <- reactive({
        if(is.null(av(input$ordinationConstraint))){
            return(NULL)
        } else {
            formstring <- paste("~", paste(input$ordinationConstraint, collapse = "+"))
            return(as.formula(formstring))
        }
    })

    observe({
        options <- "none"
        
        if(input$ordinationPlotType %in% c("samples", "split"))
            options <- c(options, sample_variables(physeq()))

        if(input$ordinationPlotType %in% c("taxa", "split"))
            options <- c(options, rank_names(physeq()))

        options <- as.list(options[order(options)])
        
        ordinationConstraintChoices <- c("none", sample_variables(physeq()))
        ordinationConstraintChoices <- as.list(ordinationConstraintChoices[order(ordinationConstraintChoices)])
        
        updateSelectInput(
            session,
            "ordinationConstraint",
            choices = ordinationConstraintChoices,
            selected = "none"
        )

        updateSelectInput(
            session,
            "ordinationPlotColor",
            choices = options,
            selected = "none"
        )

        updateSelectInput(
            session,
            "ordinationPlotShape",
            choices = options,
            selected = "none"
        )

        updateSelectInput(
            session,
            "ordinationPlotLabel",
            choices = options,
            selected = "none"
        )

        updateSelectInput(
            session,
            "ordinationPlotText",
            choices = options,
            selected = "sample_name"
        )

    })

    observeEvent(
        {
            input$ordinationMethod
        },
        {
            if(input$ordinationMethod %in% c("DCA", "CCA", "RDA", "DPCoA")) {
                updateSelectInput(
                    session,
                    "ordinationDistance",
                    choices = list("none" = "bray"),  # because distance is not used for these ordination methods
                    selected = "bray"
                )
            } else {
                updateSelectInput(
                    session,
                    "ordinationDistance",
                    choices = distlist,
                    selected = "bray"
                )
            }
        }
    )

    output$ordinationPlot <- renderPlotly({
        p <- plot_ordination(
            physeq(),
            ordinate(
                physeq(),
                method = input$ordinationMethod,
                distance = input$ordinationDistance,
                formula = get_formula_ord()
            ),
            type = input$ordinationPlotType,
            color = av(input$ordinationPlotColor),
            shape = av(input$ordinationPlotShape),
            label = av(input$ordinationPlotLabel)
        )
        
        # we have 18 reasonable shapes to work with, so they need to be repeated if there are more than 18 distinct values
        # of the variable that the user has chosen to define the shape
        if (input$ordinationPlotShape != 'none') {
          n_shapes <- length(levels(factor(p$data[, input$ordinationPlotShape])))
          good_shapes = c(1:18)
          n_reps <- ceiling(n_shapes / length(good_shapes))
          shapes_to_use <- rep(good_shapes, n_reps)
          p <- p + scale_shape_manual(values = shapes_to_use)
        }
        
        ggplotly(p)
    })

# Barplot Tab -------------------------------------------------------------

    observeEvent(
        {
            physeq()
        },

        {
            updateSelectInput(
                session,
                "barplotPlotX",
                choices = c("Sample", sample_variables(physeq()), rank_names(physeq()), "OTU"),
                selected = "Sample"
            )

            updateSelectInput(
                session,
                "barplotPlotColor",
                choices = c("none", "Sample", sample_variables(physeq()), rank_names(physeq()), "OTU"),
                selected = "none"
            )

            updateSelectInput(
                session,
                "barplotPlotGridRow",
                choices = c("Sample", sample_variables(physeq()), rank_names(physeq()), "OTU"),
                selected = NULL
            )

            updateSelectInput(
                session,
                "barplotPlotGridColumn",
                choices = c("Sample", sample_variables(physeq()), rank_names(physeq()), "OTU"),
                selected = NULL
            )

            updateSelectInput(
                session,
                "barplotPlotText",
                choices = c("none", "Sample", sample_variables(physeq()), rank_names(physeq()), "OTU"),
                selected = "none"
            )

        }
    )

    output$barplotPlot <- renderPlotly({
        p <- plot_bar(
            physeq(),
            x = input$barplotPlotX,
            y = input$barplotPlotY,
            fill = av(input$barplotPlotColor),
            facet_grid = get_facet_grid(input$barplotPlotGridRow, input$barplotPlotGridColumn)
        )

        if(!is.null(av(input$barplotPlotText))) {
            p <- p + geom_bar(mapping = aes_string(text = av(input$barplotPlotText)), stat = "identity", position = "stack", color = "black")
            print('writing text to barplot')
        }
        
        # If there are too many ASVs, the black lines that separate them can dominate the graph.
        # We allow users to remove these lines for increased clarity in such cases.
        if(input$remove_black_lines == TRUE) {
          p <- p + geom_bar(stat="identity", position="stack")
        }
        
        ggplotly(p)

    })

# Heatmap Tab -------------------------------------------------------------

    # We note that the phyloseq plot_heatmap function can fail unexpectedly if it is unable to cluster the samples for
    # whatever reason. Often, simply increasing the number of samples/ASVs will resolve this.
    observeEvent(
        {
            physeq()
        },

        {
            updateSelectInput(
                session,
                "heatmapPlotSample",
                choices = c(sample_variables(physeq())),
                selected = sample_variables(physeq())[[1]]
            )

            updateSelectInput(
                session,
                "heatmapPlotTaxa",
                choices = c(rank_names(physeq())),
                selected = rank_names(physeq())[[1]]
            )
        }
    )

    observeEvent(
        {
            input$heatmapPlotMethod
        },
        {
            if(input$heatmapPlotMethod %in% c("DCA", "CCA", "RDA", "DPCoA")) {
                updateSelectInput(
                    session,
                    "heatmapPlotDistance",
                    choices = list("none" = "bray"),  # because distance is not used for these ordination methods
                    selected = "bray"
                )
            } else {
                updateSelectInput(
                    session,
                    "heatmapPlotDistance",
                    choices = distlist,
                    selected = "bray"
                )
            }
        }
    )

    output$heatmapPlot <- renderPlotly({
        ggplotly(
            plot_heatmap(
                physeq(),
                method = input$heatmapPlotMethod,
                distance = input$heatmapPlotDistance,
                sample.label = input$heatmapPlotSample,
                taxa.label = input$heatmapPlotTaxa,
                max.label = 100,
                low = "#000033",
                high = "#66CCFF",
                trans = log_trans(4)
            )
        )
    })
    
# Replicate selection Tab -------------------------------------------------------------
    
    # Grouping samples into groups of replicates for differential abundance analysis is tricky.
    # We have implemented the most flexible possible way here at the expense of clumsy code and 
    # time-consuming user-GUI interactions. Allowing the user to export/import replicate assignments
    # as text files is a stop-gap measure to make this bearable.
    
    # Possibly the clearest way to implement replicate assignment in the future would be via sample metadata.  
    # For example, all samples with the same date, location, and sample_type could be grouped. This, however, 
    # requires either perfect planning ahead of time or ongoing modification to the sample_data.
    
    featureFilteredPhy_ordered_subset <- reactive({
      req(sampleFilteredPhy(), featureFilteredPhy_ordered)
      selected_sample_row_names = row.names(sample_data(sampleFilteredPhy()))
      
      # remove samples from selection table if they have already been assigned to a replicate
      selected_sample_row_names_minus_assigned = selected_sample_row_names
      
      n_groups_used <- length(user_replicates[['replicate_list']])
      if (n_groups_used > 0) {
        for (group_number in 1:n_groups_used) {
          n_reps_in_group <- length(user_replicates[['replicate_list']][[group_number]])
          for (rep_number in 1:n_reps_in_group) {
            if (user_replicates[['replicate_list']][[group_number]][[rep_number]] %in% selected_sample_row_names_minus_assigned) {
              selected_sample_row_names_minus_assigned <- selected_sample_row_names_minus_assigned[!selected_sample_row_names_minus_assigned %in% user_replicates[['replicate_list']][[group_number]][[rep_number]]]
            }
          }
        }
      }

      if (length(selected_sample_row_names_minus_assigned) > 0) {
        t <- featureFilteredPhy_ordered()[selected_sample_row_names_minus_assigned, ]
      } else {
        t <- data.frame()
      }

      t
    })
    
    observeEvent(
      {
        featureFilteredPhy_ordered_subset()
      },
      {
        output$replicateSelectionDT <- DT::renderDataTable(
          data.frame(featureFilteredPhy_ordered_subset()),
          class = 'compact stripe',
          server = FALSE,
          rownames = FALSE,
          filter = list(position = "top"),
          options = list(
            pageLength = 20,
            lengthMenu = c(10, 20, 40, 200, 1000),
            scrollX = TRUE,
            searchHighlight = TRUE
          ),
          selection = list(
            target = 'row'
          )
        )
      }
    )
    
    current_row_replicate_selection <- reactive({
      r <- input$replicateSelectionDT_rows_selected
      r
    })
    
    replicate_subset_table <- reactive({
      if (length(current_row_replicate_selection()) > 0) {
        replicate_subset <- featureFilteredPhy_ordered_subset()[current_row_replicate_selection(), ]
      } else {
        replicate_subset <- data.frame()
      }
      
      replicate_subset
    })
    
    observeEvent(
      {
        replicate_subset_table()
      },
      {
        output$currentReplicateTableDT <- DT::renderDataTable(
          replicate_subset_table(),
          class = 'compact stripe',
          server = FALSE,
          rownames = FALSE,
          options = list(
            pageLength = 20,
            scrollX = TRUE,
            searching = FALSE
          ),
          selection = 'none'
        )
        
      }
    )
    
    user_replicates <- reactiveValues()
    user_replicates[['replicate_list']] <- list()
    user_replicates[['table']] <- data.frame()

    observeEvent(
      {
        input$storeReplicate
      },
      {
        
        # validate replicate name
        if (input$replicateName %in% names(user_replicates[['replicate_list']])) {
          output$replicateNameError <- renderText('ERROR: must input a unique Replicate Name!')
        } else if (nchar(input$replicateName) == 0) {
          output$replicateNameError <- renderText('ERROR: must input a Replicate Name!')
        } else {
          output$replicateNameError <- renderText('Stored Replicate!')
          user_replicates[['replicate_list']][[input$replicateName]] <- row.names(featureFilteredPhy_ordered_subset()[current_row_replicate_selection(), ])
          updateTextInput(session, 'replicateName', value = "")
          
          # display all replicate assignments in a table
          replicate_names <- names(user_replicates[['replicate_list']])
          n_replicates <- length(replicate_names)
          replicate_samples <- vector()
          for (rn in 1:n_replicates) {
            replicate_samples[[rn]] <- paste(user_replicates[['replicate_list']][[replicate_names[[rn]]]], collapse=';')
          }
          
          user_replicates[['table']] <- data.table(replicate_names, replicate_samples)
          
          }
      }
    )
    
    output$replicateAssignmentTableDT <- DT::renderDataTable(
      user_replicates[['table']],
      class = 'compact stripe',
      server = FALSE,
      rownames = FALSE,
      options = list(
        pageLength = 20,
        scrollX = TRUE,
        searching = FALSE
      ),
      selection = 'none'
    )

    # export replicates file
    output$exportReplicates <- downloadHandler(
      filename = function() {
        paste0('replicates_', format(Sys.time(), "%Y%m%d_%H%M%S"))
      },
      content = function(file) {
        replicate_names <- names(user_replicates[['replicate_list']])
        n_replicates <- length(replicate_names)
        replicate_samples <- vector()
        for (rn in 1:n_replicates) {
          replicate_samples[[rn]] <- paste(user_replicates[['replicate_list']][[replicate_names[[rn]]]], collapse=';')
        }
        user_replicate_table = data.table(replicate_names, replicate_samples)
        write.csv(user_replicate_table, file, row.names = FALSE)
      }
    )
    
    # import replicates file
    replicates_from_file <- reactive({
      replicate_file_data <- read.csv(input$replicate_input_file$datapath, header = TRUE, stringsAsFactors = FALSE)
      
      replicate_file_data
    })
    
    observeEvent(
      {
        input$import_replicates_button
      },
      {
        user_replicates[['table']] <- replicates_from_file()
        
        n_replicates <- nrow(user_replicates[['table']])
        
        user_replicates[['replicate_list']] <- list()
        for(r in 1:n_replicates) {
          user_replicates[['replicate_list']][[user_replicates[['table']][[r, 'replicate_names']]]] <- vector()
          samples_string_split <- strsplit(user_replicates[['table']][r, 'replicate_samples'], ';')[[1]]
          n_samples <- length(samples_string_split)
          
          for (s in 1:n_samples) {
            user_replicates[['replicate_list']][[user_replicates[['table']][[r, 'replicate_names']]]][[s]] <- samples_string_split[[s]]
          }
        }
        updateTabsetPanel(session, 'replicate_tabs', selected = "Replicate assignments")
      }
    )
    
    observeEvent(input$finalizeReplicates,
                 {
                   updateTabsetPanel(session, 'replicate_tabs', selected = "Samples with replicate information")  
                   updateSelectInput(session, 'diff_numerator', choices = names(user_replicates[['replicate_list']]), 
                                     selected = names(user_replicates[['replicate_list']][1]))
                   updateSelectInput(session, 'diff_denominator', choices = names(user_replicates[['replicate_list']]), 
                                     selected = names(user_replicates[['replicate_list']][length(names(user_replicates[['replicate_list']]))]))
                   
                   updateSelectInput(session, 'pairwise_permanova_1', choices = names(user_replicates[['replicate_list']]), 
                                     selected = names(user_replicates[['replicate_list']][1]))
                   updateSelectInput(session, 'pairwise_permanova_2', choices = names(user_replicates[['replicate_list']]), 
                                     selected = names(user_replicates[['replicate_list']][2]))
                 }
                 )
    
    # THE PHYLOSEQ OBJECT IS LARGE - COPYING EACH TIME THE USER ADDS A REPLICATE BY HAND WILL BE SLOW
    # THUS, USE A BUTTON TO "FINALIZE REPLICATE ASSIGNMENTS"
    featureAndReplicateFilteredPhy <- eventReactive(
      {
        input$finalizeReplicates
      },
      {
        # for differential analysis of taxa, starting with the (mostly) unfiltered dataset makes the most sense
        # this gives the permanova and DESeq2 algorithms all of the data to start with, and does not confuse it with 
        # any normalization or removal of samples
        # Note that if this is to be changed, significant rework will be required to get everything to work ideally
        x <- featureFilteredPhy()   
        s <- data.frame(sample_data(x))
        
        # add user_replicate column to s
        # all replicates need to be defined at once to avoid confusion
        s$user_replicate <- NA

        r <- user_replicates[['replicate_list']]
        n_rep_groups <- length(r)
        
        all_user_rep_names <- vector()
        for (i in 1:n_rep_groups) {
          n_current_reps <- length(r[[i]])
          for (j in 1:n_current_reps) {
            all_user_rep_names <- c(all_user_rep_names, r[[i]][[j]])
          }
        }
        # clear selection then select
        DT::selectRows(samplesDT_proxy, NULL)
        
        # must select by row number, not row name
        # get row numbers from the names from featureFilteredPhy_ordered()
        n_user_reps <- length(all_user_rep_names)
        all_user_rep_row_numbers <- vector()
        for (i in 1:n_user_reps) {
          rl <- which(row.names(featureFilteredPhy_ordered()) == all_user_rep_names[[i]])
          all_user_rep_row_numbers <- c(all_user_rep_row_numbers, rl)
        }
        
        DT::selectRows(samplesDT_proxy, all_user_rep_row_numbers)
        
        for (i in 1:n_rep_groups) {
          n_rep_members <- length(r[[i]])
          for (j in 1:n_rep_members) {
            s[r[[i]][[j]], "user_replicate"] <- names(r)[[i]]
          }
        }

        # add ss back into x as its sample data
        ss <- sample_data(s)
        sample_data(x) <- ss
        
        xrf <- subset_samples(x, complete.cases(user_replicate))
        xrf <- prune_taxa(taxa_sums(xrf) > 0, xrf)

        xrf
      })
    
    output$phyWithReplicatesSamplesDT <- DT::renderDataTable(
      data.frame(sample_data(featureAndReplicateFilteredPhy())),
      class = 'compact stripe',
      server = FALSE,
      rownames = FALSE,
      filter = list(position = "top"),
      extensions = c('Buttons'),
      options = list(
        pageLength = 20,
        lengthMenu = c(10, 20, 40, 200, 1000),
        dom = 'Blfrtip',
        scrollX = TRUE,
        searchHighlight = TRUE
      )
    )
    
    # Permanova Tab -------------------------------------------------------------
    observeEvent(
      input$runPermanova,
      {
        set.seed(43)  # for April 3
        
        # This example of permanova:
        # http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
        # uses scaled data for permanova; they scale to the smallest library size, which is inconsistent with the latest advice
        # we transform to 10000 reads instead, but don't claim that this is best practice
        # At worst, it seems like this would reduce statistical power rather than give false positives
        scaled_frfp <- transform_sample_counts(featureAndReplicateFilteredPhy(), function(x) {round(x / sum(x) * 10000)})
        
        if (input$build_tree == TRUE) {
          phy_tree(scaled_frfp) <- phy_tree(physeq())
        }
        
        
        xrf_dist <- phyloseq::distance(scaled_frfp, method = input$permanovaDistanceMetric)
        sample_df <- data.frame(sample_data(scaled_frfp))
        
        # for Weighted Unifrac, there are NaNs appearing in the distance - this seems to cause adonis to give an error
        
        # Adonis test
        # Is there an opportunity to enhance the design formula with user-specified details?
        a <- adonis(xrf_dist ~ user_replicate, data = sample_df)
        beta <- betadisper(xrf_dist, sample_df$user_replicate)
        perm <- permutest(beta)
        
        output$permanovaTextResults <- renderPrint(a)
        output$betaPermResults <- renderPrint(perm)
        
        # add an ordination option here, colored by replicate to visualize the centriods and dispersion
        output$replicateOrdinationPlot <- renderPlotly({
          p <- plot_ordination(
            scaled_frfp,
            ordinate(
              scaled_frfp,
              method = input$replicateOrdinationMethod,
              distance = input$permanovaDistanceMetric,
              formula = NULL
            ),
            type = "samples",
            color = "user_replicate"
          )
          
          if(!is.null(av(input$ordinationPlotText)))
            p <- p + geom_point(aes_string(text = av(input$ordinationPlotText)))
          
          ggplotly(p)
        })
        
      }
    )
    
    # Pairwise Permanova Tab ------------------------------------------------
    observeEvent(
      input$runPairwisePermanova,
      {
        set.seed(59)  # for May 9 but could be changed to anything
        
        # same logic as the overall permanova
        scaled_frfp <- transform_sample_counts(featureAndReplicateFilteredPhy(), function(x) {round(x / sum(x) * 10000)})
        if (input$build_tree == TRUE) {
          phy_tree(scaled_frfp) <- phy_tree(physeq())
        }
        
        # prune samples from scaled_frfp so that it contains only the samples in the two groups selected
        sd <- data.frame(sample_data(scaled_frfp))
        keep_row <- sd[, 'user_replicate'] == input$pairwise_permanova_1 | sd[, 'user_replicate'] == input$pairwise_permanova_2
        scaled_frfp_subset <- prune_samples(keep_row, scaled_frfp)
        
        xrf_dist <- phyloseq::distance(scaled_frfp_subset, method = input$pairwisePermanovaDistanceMetric)
        sample_df <- data.frame(sample_data(scaled_frfp_subset))
        
        # Adonis test
        a <- adonis(xrf_dist ~ user_replicate, data = sample_df)
        beta <- betadisper(xrf_dist, sample_df$user_replicate)
        perm <- permutest(beta)
        
        output$pairwisePermanovaTextResults <- renderPrint(a)
        output$pairwiseBetaPermResults <- renderPrint(perm)
        
        output$pairwiseReplicateOrdinationPlot <- renderPlotly({
          p <- plot_ordination(
            scaled_frfp_subset,
            ordinate(
              scaled_frfp_subset,
              method = input$pairwiseReplicateOrdinationMethod,
              distance = input$pairwisePermanovaDistanceMetric,
              formula = NULL
            ),
            type = "samples",
            color = "user_replicate"
          )
      })
      }
    )
    
    # Differential Abundance Tab ---------------------------------------------
    
    # replace amplicon sequences with hashes for taxa names
    frfp <- reactive({
      f <- featureAndReplicateFilteredPhy()
      taxa_names(f) <- as.vector(tax_table(f)[, 'SVhash'])
      
      f
    })
    
    dds_phy_filtered <- eventReactive(
      {
        input$runDifferentialAbundance
      },
      {
        # TODO: provide an option to specify data columns to control with in the design formula?
        # i.e. ~ treatment + user_replicate to control for the effect of treatment
        
        # TODO: implement a way to ensure that the control level is the first level of the user_replicate column?
        dds_pf <- phyloseq_to_deseq2(frfp(), ~ user_replicate)

        # code from https://github.com/joey711/phyloseq/issues/387
        # this avoids an error that happens if every ASV has a zero in at least one sample
        # calculate geometric means prior to estimate size factors
        gm_mean <- function(x, na.rm=TRUE){
          exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
        }
        
        geoMeans <- apply(counts(dds_pf), 1, gm_mean)
        dds_pf <- estimateSizeFactors(dds_pf, geoMeans = geoMeans)
        dds_pf <- DESeq(dds_pf, test = "Wald", fitType="local")
        
        dds_pf
      }
    )
        
    dds_phy_filtered_results <- reactive({
      r <- results(dds_phy_filtered(), contrast=c('user_replicate', input$diff_numerator, input$diff_denominator), cooksCutoff = FALSE)
      
      r
    })
    
    significant_results_df <- reactive({
      alpha <- input$differential_abundance_p
      sr_table <- dds_phy_filtered_results()[which(dds_phy_filtered_results()$padj < alpha), ]

      # account for the possibility of nothing in sr_table being less than alpha - no significant results
      if (nrow(sr_table) > 0) {
        sr_table <- cbind(as(sr_table, "data.frame"), 
                          as(tax_table(frfp())[rownames(sr_table), ], "matrix"))
        
        sr_df <- data.frame(sr_table)
        setcolorder(sr_df, c('SVhash', 'log2FoldChange', 'padj', 'Species', 'Genus', 'Family', 'Order', 
                             'Class', 'Phylum', 'Kingdom', 'baseMean', 'lfcSE', 'stat', 'pvalue'))
      }
      else {
        sr_df <- data.frame("No significantly different ASVs")
      }
      
      sr_df
    })
    
    output$differential_abundance_table <- DT::renderDataTable(
      significant_results_df(),
      class = 'compact stripe',
      server = FALSE,
      rownames = FALSE,
      filter = list(position = "top"),
      extensions = c('Buttons'),
      options = list(
        pageLength = 20,
        lengthMenu = c(10, 20, 40, 200, 1000),
        dom = 'Blfrtip',
        scrollX = TRUE,
        searchHighlight = TRUE
      ),
      selection = list(
        mode = 'single',
        target = 'row',
        selected = 1)
    )
    
    differential_abundance_DT_proxy <- DT::dataTableProxy("differential_abundance_table", deferUntilFlush = TRUE)
    
    output$differential_abundance_volcano <- renderPlotly({
      v <- significant_results_df()[, c('log2FoldChange', 'padj', 'SVhash', 'Family', 'Genus', 'Species')]
      v[, 'neglog10_padj'] <- -log10(v[, 'padj'])
      rownames(v) <- significant_results_df()[['SVhash']]
      
      volcano_plot <- ggplot(v, 
        aes(x=log2FoldChange, y=neglog10_padj, hash=SVhash, padj=padj, Family=Family, Genus=Genus, Species=Species)) + 
        geom_point() + 
        xlab('log2(fold change)') + ylab('-log10(p adjusted)') + 
        geom_hline(aes(yintercept = -log10(input$differential_abundance_p)), color='red')
      
      volcano_plotly <- ggplotly(volcano_plot, tooltip = c('SVhash', 'Family', 'Genus', 'Species', 'log2FoldChange', 'padj'), source = "V" )
      
      volcano_plotly
    })
    
    volcano_plot_selected_row <- reactive({
      d <- event_data("plotly_click", source = "V")
      
      #plotly pointNumber is 0-indexed
      row_number <- d$pointNumber + 1
      
      row_number
    })
    
    observeEvent(
      volcano_plot_selected_row(),
      DT::selectRows(differential_abundance_DT_proxy, volcano_plot_selected_row())
    )
    
    output$volcano_plot_point_table <- DT::renderDataTable(
      significant_results_df()[volcano_plot_selected_row(), ],
      server = FALSE,
      rownames = FALSE,
      options = list(
        scrollX = TRUE,
        dom = 't'
      ),
      selection = 'none'
    )
    
    differential_abundance_selected_hash <- reactive({
      r <- input$differential_abundance_table_rows_selected
      h <- as.character(significant_results_df()[[r, 'SVhash']])

      h
    })
    
    genus_taxonomy_plot_data <- reactive({
      # code modified from: https://joey711.github.io/phyloseq-extensions/DESeq2.html
      
      req(significant_results_df())

      # Phylum order
      sig_t <- significant_results_df()
      x = tapply(sig_t$log2FoldChange, sig_t$Phylum, function(x) max(x))
      x = sort(x, TRUE)
      sig_t$Phylum = factor(as.character(sig_t$Phylum), levels=names(x))
      
      # Genus order
      x = tapply(sig_t$log2FoldChange, sig_t$Genus, function(x) max(x))
      x = sort(x, TRUE)
      sig_t$Genus = factor(as.character(sig_t$Genus), levels=names(x))
      
      sig_t
    })
    
    output$genus_taxonomy_plot <- renderPlotly({
      # code modified from: https://joey711.github.io/phyloseq-extensions/DESeq2.html
      
      theme_set(theme_bw())
      scale_fill_discrete <- function(palname = "Set1", ...) {
        scale_fill_brewer(palette = palname, ...)
      }
      
      pt <- ggplot(genus_taxonomy_plot_data(), aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=1) + 
        theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
      
      pt_plotly <- ggplotly(pt, tooltip = c('SVhash', 'Family', 'Genus', 'Species', 'log2FoldChange', 'padj'), source = "G")
      
      pt_plotly
    })
    
    genus_plot_selected_row <- reactive({
      # For this plot, the pointNumber reported by the event_data does not match the row of the plot data
      # This is a known limitation of plotly events (matching the row number isn't always guaranteed)
      # We need to find the point that matches the x and y coords
      
      req(event_data("plotly_click", source = "G"))

      d <- event_data("plotly_click", source = "G")
      
      g <- genus_taxonomy_plot_data()
      current_genus <- levels(g$Genus)[d$x]
      current_approx_fold_change <- d$y
      current_hash <- as.character(filter(g, Genus == current_genus, near(log2FoldChange, current_approx_fold_change, 0.1))$SVhash)
      row_number <- which(row.names(g) == current_hash)

      row_number
    })
    
    observeEvent(
      genus_plot_selected_row(),
      DT::selectRows(differential_abundance_DT_proxy, genus_plot_selected_row())
    )
    
    output$genus_plot_point_table <- DT::renderDataTable(
      significant_results_df()[genus_plot_selected_row(), ],
      server = FALSE,
      rownames = FALSE,
      options = list(
        scrollX = TRUE,
        dom = 't'
      ),
      selection = 'none'
    )
    
    output$counts_plot <- renderPlot({
      counts_plot_data <- plotCounts(dds_phy_filtered(), differential_abundance_selected_hash(), intgroup = 'user_replicate', 
                                     normalized = as.logical(input$diff_abundance_counts_normalized),
                                     replaced = as.logical(input$diff_abundance_counts_outlier),
                                     returnData = TRUE)
      counts_plot <- ggplot(counts_plot_data, aes(x=user_replicate, y=count, color=user_replicate)) + 
        geom_boxplot(notch = FALSE, outlier.size = 0) + geom_jitter(shape=15, position = position_jitter(width=0.2, height=0)) + 
        stat_summary(fun.y=mean, geom="point", shape=23, size=4) + ylab('deseq2 adjusted counts') +
        theme(legend.position = "none")
      
      if (as.logical(input$diff_abundance_counts_transform) == TRUE) {
        counts_plot <- counts_plot + scale_y_continuous(trans = 'log10')
      }

      counts_plot
    })
    
    output$raw_counts_plot <- renderPlot({
      f <- frfp()
      t <- otu_table(f)
      t_df <- data.frame(t@.Data)
      m <- melt(t_df[differential_abundance_selected_hash(), ])
      colnames(m) <- c('sample_name', 'raw_counts')
      row.names(m) <- as.character(m[, 'sample_name'])
      
      ur <- user_replicates[['replicate_list']]
      ur_m <- melt(ur)
      colnames(ur_m) <- c('sample_name', 'user_replicate')
      row.names(ur_m) <- as.character(ur_m[, 'sample_name'])
      
      m_reps <- merge(m, ur_m, by=0)[, c('Row.names', 'raw_counts', 'user_replicate')]
      row.names(m_reps) <- as.character(m_reps[, 'Row.names'])
      
      raw_counts_plot <- ggplot(m_reps, aes(x=user_replicate, y=raw_counts, color=user_replicate)) + 
        geom_boxplot(notch = FALSE, outlier.size = 0) + geom_jitter(shape=15, position = position_jitter(width=0.2, height=0)) + 
        stat_summary(fun.y=mean, geom="point", shape=23, size=4) + ylab('phyloseq counts') +
        theme(legend.position = "none")
      
      if (as.logical(input$diff_abundance_counts_transform) == TRUE) {
        raw_counts_plot <- raw_counts_plot + scale_y_continuous(trans = 'log10')
      }

      raw_counts_plot
    })
    
})
