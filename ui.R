library(shinycssloaders)

function(request) {
    dashboardPage(
        dashboardHeader(
            title = "Microbiome barcoding"
        ),

        dashboardSidebar(
            sidebarMenu(
                id = "menu",
                menuItem(
                    "Data Selection",
                    icon = icon("list", lib = "font-awesome"),
                    startExpanded = TRUE,
                    menuSubItem("Samples", tabName = "samplesTab", icon = icon("bug")),
                    menuSubItem("Taxonomy", tabName = "taxonomyTab", icon = icon("book")),
                    menuSubItem("Counts", tabName = "countsTab", icon = icon("table"))

                ),
                menuItem("Phylogenetic Tree", tabName = "treeTab", icon = icon("sitemap")),
                menuItem(
                    "Diversity",
                    icon = icon("leaf"),
                    startExpanded = TRUE,
                    menuSubItem("Alpha", tabName = "diversityTab", icon = icon("pie-chart")),
                    menuSubItem("Beta/Ordination", tabName = "ordinationTab", icon = icon("exchange"))
                ),
                menuItem(
                    "Plots",
                    icon = icon("area-chart"),
                    startExpanded = TRUE,
                    menuSubItem("Barplot", tabName = "barplotTab", icon = icon("bar-chart")),
                    menuSubItem("Heatmap", tabName = "heatmapTab", icon = icon("th-large")),
                    menuSubItem("Scatter", tabName = "scatterTab", icon = icon("check-square-o"))
                ),
                menuItem(
                  "Sample group eval",
                  icon = icon("globe"),
                  startExpanded = TRUE,
                  menuSubItem("Replicate Selection", tabName = "replicateTab", icon = icon("check-double")),
                  menuSubItem("Permanova", tabName = "permanovaTab", icon = icon("dna")),
                  menuSubItem("Pairwise Permanova", tabName = "pairwisePermanovaTab", icon = icon("dna")),
                  menuSubItem("Differential Abundance", tabName = "differentialTab", icon = icon("search"))
                ),
                
                menuItem("Mining", tabName = "miningTab", icon = icon("magic")),
                valueBoxOutput("samplesBox", width = 12),
                valueBoxOutput("seqBox", width = 12),
                valueBoxOutput("treeBox", width = 12) # ,
                # bookmarkButton()
            )
        ),

        dashboardBody(
            tabItems(
                tabItem(tabName = "samplesTab",
                        fluidRow(
                            box(
                                title = "Sample Feature Filter",
                                width = 12,
                                selectInput(
                                    "samplesFeatureSelect",
                                    label = "Retain Features",
                                    choices = NULL,
                                    multiple = T
                                ),
                                actionButton("filterSamples", "Filter Samples")
                            )
                        ),
                        fluidRow(
                            box(
                                title = "Strain Filter (changing this resets much of the app - double check ALL selection/options and click Apply Selection afterward!!)",
                                width = 12,
                                checkboxInput("ignoreStrains", "Remove some strains from analysis")
                            )
                        ),
                        fluidRow(
                            box(
                                width = 12,
                                actionButton("samplesDT_filtered", label = "Select Filtered"),
                                actionButton("samplesDT_visible", label = "Select Visible"),
                                actionButton("samplesDT_reset", label = "Reset Selection"),
                                withSpinner(DT::dataTableOutput("samplesDT")),
                                verbatimTextOutput("n_samples_selected"),
                                checkboxInput("build_tree", "Build Phylogenetic Tree"),
                                actionButton("samplesDT_apply", "Apply Selection")
                            )
                        )
                ),

                tabItem(tabName = "taxonomyTab",
                        fluidRow(
                            box(
                                title = "Taxonomy Filter",
                                width = 12,
                                column(
                                    width = 1,
                                    selectInput(
                                        "taxonomyRankSelect",
                                        label = "Rank",
                                        choices = NULL,
                                        multiple = F
                                    )
                                ),
                                column(
                                    width = 11,
                                    selectInput(
                                        "taxonomyTaxaSelect",
                                        label = "Retain Taxa",
                                        choices = NULL,
                                        multiple = T
                                    )
                                )
                            )
                        )
                ),

                tabItem(tabName = "countsTab",
                        fluidRow(
                            box(
                                title = "Sample Filter",
                                width = 3,
                                height = "140px",
                                numericInput(
                                    "countFilteringSample",
                                    label = "Min counts",
                                    value = 10000,
                                    min = 0,
                                    step = 1000
                                )
                            ),
                            box(
                                title = "Sequence Variant Filter",
                                width = 6,
                                height = "140px",
                                column(
                                    width = 6,
                                    numericInput(
                                        "countFilteringA",
                                        label = "Exceed A counts (default of 50 may be too large for your purposes)...",
                                        value = 50,
                                        min = 0,
                                        step = 1
                                    )
                                ),
                                column(
                                    width = 6,
                                    numericInput(
                                        "countFilteringK",
                                        label = "... in at least K samples",
                                        value = 1,
                                        min = 0,
                                        step = 1
                                    )
                                )
                            ),
                            box(
                                title = "Transformation",
                                width = 3,
                                height = "140px",
                                selectInput(
                                    "transformationMethod",
                                    label = "Method",
                                    choices = list("None" = "none", "Normalized to 10,000 counts" = "prevalence", "Presence/Absence" = "binary"),
                                    selected = "prevalence"
                                )
                            )
                        ),

                        fluidRow(
                            tabBox(
                                id = "countsTabBox",
                                title = NULL,
                                width = 12,
                                tabPanel(
                                    id = "countsTabPanelCounts",
                                    title = "Counts",
                                    checkboxInput('counts_phylogeny_row', "Show phylogeny details", value=FALSE),
                                    downloadButton('countsTableDownload', 'Download entire counts table'),
                                    verbatimTextOutput("countsDownloadWarning"),
                                    withSpinner(DT::dataTableOutput("countsCountsDT")),
                                    DT::dataTableOutput("selected_hash")
                                ),
                                tabPanel(
                                    id = "countsTabPanelHash",
                                    title = "Hash Dictionary",
                                    withSpinner(DT::dataTableOutput("countsHashDT"))
                                )
                            )
                            
                        )
                ),

                tabItem(tabName = "treeTab",
                    fluidRow(
                        box(
                            width = 12,
                            height = '10000px',
                            withSpinner(plotOutput("ggtree_graphic", inline=TRUE))
                        )
                    )
                ),

                tabItem(tabName = "diversityTab",
                        fluidRow(
                            box(
                                width = 12,
                                helpText("Coming soon!")
                            )
                        )
                ),

                tabItem(tabName = "ordinationTab",
                        fluidRow(
                            box(
                                width = 9,
                                height = "640px",
                                withSpinner(plotlyOutput("ordinationPlot", width = "auto", height = "625px"))
                            ),
                            box(
                                width = 3,
                                height = "640px",
                                selectInput(
                                    "ordinationPlotType",
                                    label = "Type",
                                    choices = list("Samples" = "samples", "Sequence Variants" = "taxa", "Samples and Sequence Variants" = "split"),
                                    selected = "split"
                                ),
                                selectInput(
                                    "ordinationMethod",
                                    label = "Method",
                                    choices = ordlist,
                                    selected = "PCoA"
                                ),
                                selectInput(
                                    "ordinationDistance",
                                    label = "Distance",
                                    choices = distlist,
                                    selected = "bray"
                                ),
                                selectInput(
                                    "ordinationConstraint",
                                    label = "Constraint",
                                    choices = "none",
                                    selected = "none"
                                ),
                                selectInput(
                                    "ordinationPlotColor",
                                    label = "Color",
                                    choices = "none",
                                    selected = "none"
                                ),
                                selectInput(
                                    "ordinationPlotShape",
                                    label = "Shape",
                                    choices = "none",
                                    selected = "none"
                                ),
                                selectInput(
                                    "ordinationPlotLabel",
                                    label = "Label",
                                    choices = "none",
                                    selected = "none"
                                )#,
                                # selectInput(
                                #     "ordinationPlotText",
                                #     label = "Text",
                                #     choices = "none",
                                #     selected = "none"
                                # )
                            )
                        )
                ),

                tabItem(tabName = "barplotTab",
                        fluidRow(
                            box(
                                width = 9,
                                height = "640px",
                                withSpinner(plotlyOutput("barplotPlot", width = "auto", height = "625px"))
                            ),
                            box(
                                width = 3,
                                height = "640px",
                                selectInput(
                                    "barplotPlotX",
                                    label = "X",
                                    choices = c("Sample"),
                                    selected = "Sample"
                                ),
                                selectInput(
                                    "barplotPlotY",
                                    label = "Y",
                                    choices = c("Abundance"),
                                    selected = "Abundance"
                                ),
                                selectInput(
                                    "barplotPlotColor",
                                    label = "Color",
                                    choices = c("none", "OTU"),
                                    selected = "none"
                                ),
                                selectInput(
                                    "barplotPlotGridRow",
                                    label = "Grid Row",
                                    choices = "",
                                    multiple = TRUE
                                ),
                                selectInput(
                                    "barplotPlotGridColumn",
                                    label = "Grid Column",
                                    choices = "",
                                    multiple = TRUE
                                ),
                                selectInput(
                                    "barplotPlotText",
                                    label = "Text",
                                    choices = "none",
                                    selected = "none"
                                ),
                                checkboxInput('remove_black_lines', 'Remove lines separating ASVs')
                            )
                        )
                ),

                tabItem(tabName = "heatmapTab",
                        fluidRow(
                            box(
                                width = 9,
                                height = "600px",
                                withSpinner(plotlyOutput("heatmapPlot", width = "auto", height = "625px"))
                            ),
                            box(
                                width = 3,
                                height = "600px",
                                selectInput(
                                    "heatmapPlotMethod",
                                    label = "Method",
                                    choices = ordlist,
                                    selected = "NMDS"
                                ),
                                selectInput(
                                    "heatmapPlotDistance",
                                    label = "Distance",
                                    choices = distlist,
                                    selected = "bray"
                                ),
                                selectInput(
                                    "heatmapPlotSample",
                                    label = "Sample Axis",
                                    choices = c("Sample"),
                                    selected = "Sample"
                                ),
                                selectInput(
                                    "heatmapPlotTaxa",
                                    label = "Taxa Axis",
                                    choices = c("OTU"),
                                    selected = "OTU"
                                )
                            )
                        )
                ),

                tabItem(tabName = "scatterTab",
                        fluidRow(
                            box(
                                width = 12,
                                helpText("Coming soon!")
                            )
                        )
                ),
                tabItem(tabName = "replicateTab",
                        fluidRow(
                          box(
                            width = 12,
                            helpText("You can EITHER (1) import a CSV file OR (2) group replicates manually")
                          ),
                          box(
                            width = 12,
                            tabsetPanel(id = 'replicate_tabs',
                              tabPanel("Import group assignment file",
                                       fileInput('replicate_input_file', 'Choose replicate CSV file',
                                                 accept=c('text/csv',
                                                          'text/comma-separated-values,text/plain',
                                                          '.csv')),
                                       actionButton('import_replicates_button', 'Import group assignments')
                              ),
                              tabPanel("Manually place replicate samples into groups",
                                       DT::dataTableOutput("replicateSelectionDT"),
                                       helpText("Define your groups of replicates"),
                                       DT::dataTableOutput('currentReplicateTableDT'),
                                       textInput('replicateName', 'Group Name', width = 300),
                                       textOutput('replicateNameError'),
                                       actionButton('storeReplicate', 'Store replicate-group assignment')
                              ),
                              tabPanel("Replicate assignments",
                                       DT::dataTableOutput('replicateAssignmentTableDT'),
                                       actionButton('finalizeReplicates', 'Finalize Replicate-Group Assignments'),
                                       downloadLink('exportReplicates', 'Export Replicate-Group assignments')
                                       ),
                              tabPanel("Samples with replicate information",
                                       withSpinner(DT::dataTableOutput('phyWithReplicatesSamplesDT'))
                                       )
                            )
                          )
                        )
                ),
                tabItem(tabName = "permanovaTab",
                        fluidRow(
                          box(
                            width = 12,
                            helpText("Permanova analysis of selected replicates; currently using UNFILTERED DATA transformed to 10000 reads"),
                            selectInput(
                              "permanovaDistanceMetric",
                              label = "Distance",
                              choices = distlist,
                              selected = "bray"
                              ),
                            actionButton("runPermanova", "Run Permanova analysis")
                          ),
                          box(
                            width = 12,
                            helpText("Permanova summary"),
                            verbatimTextOutput("permanovaTextResults")
                          ),
                          box(
                            width = 12,
                            helpText("Permutation homogeneity summary"),
                            verbatimTextOutput("betaPermResults")
                          ),
                          box(
                              width = 12,
                              height = "800px",
                              helpText("Ordination colored by replicate group"),
                              selectInput(
                                  "replicateOrdinationMethod",
                                  label = "Method",
                                  choices = ordlist,
                                  selected = "PCoA"
                              ),
                              plotlyOutput("replicateOrdinationPlot", width = "auto", height = "625px")
                          )
                        )
                ),
                tabItem(tabName = "pairwisePermanovaTab",
                        fluidRow(
                            box(width = 12,
                                helpText = "Pairwise Permanova analysis",
                                selectInput('pairwise_permanova_1', 'Group 1', c("none"), selected = "none"),
                                selectInput('pairwise_permanova_2', 'Group 2', c("none"), selected = "none"),
                                selectInput(
                                    "pairwisePermanovaDistanceMetric",
                                    label = "Distance",
                                    choices = distlist,
                                    selected = "bray"
                                ),
                                actionButton("runPairwisePermanova", "Run Pairwise Permanova analysis")
                            ),
                            box(
                                width = 12,
                                helpText("Permanova summary"),
                                verbatimTextOutput("pairwisePermanovaTextResults")
                            ),
                            box(
                                width = 12,
                                helpText("Permutation homogeneity summary"),
                                verbatimTextOutput("pairwiseBetaPermResults")
                            ),
                            box(
                                width = 12,
                                height = "800px",
                                helpText("Ordination colored by replicate group"),
                                selectInput(
                                    "pairwiseReplicateOrdinationMethod",
                                    label = "Method",
                                    choices = ordlist,
                                    selected = "PCoA"
                                ),
                                plotlyOutput("pairwiseReplicateOrdinationPlot", width = "auto", height = "625px")
                            )
                        )
                    ),
                tabItem(tabName = "differentialTab",
                        fluidRow(
                          box(width = 12,
                              tabsetPanel(id='differential_tabs',
                                          tabPanel("Differential Abundance",
                                                   helpText("Differential Abundance; currently using UNFILTERED DATA"),
                                                   selectInput('diff_numerator', 'Group to compare (numerator)', c("none"), selected = "none"),
                                                   selectInput('diff_denominator', 'Group to compare against (denominator)', c("none"), selected = "none"),
                                                   numericInput("differential_abundance_p", "Adjusted P value", 0.05, min=0, max=1, step=0.01, width = 150),
                                                   actionButton("runDifferentialAbundance", "Run Differential Abundance"),
                                                   withSpinner(DT::dataTableOutput("differential_abundance_table"))
                                          ),
                                          tabPanel("Volcano Plot",
                                                   plotlyOutput("differential_abundance_volcano"),
                                                   #verbatimTextOutput("volcano_plot_click"),
                                                   DT::dataTableOutput("volcano_plot_point_table")
                                                   
                                                   ),
                                          tabPanel("Taxonomy plot",
                                                   fluidRow(
                                                     box(
                                                       width = 12,
                                                       helpText("Significantly different ASVs organized by genus"),
                                                       plotlyOutput("genus_taxonomy_plot"),
                                                       DT::dataTableOutput("genus_plot_point_table")
                                                     )
                                                   )),
                                          tabPanel("Counts Plot",
                                                   selectInput('diff_abundance_counts_normalized', 'Normalize counts?', 
                                                               c('TRUE', 'FALSE'), selected = TRUE),
                                                   selectInput('diff_abundance_counts_outlier', 'Outlier replacement?', 
                                                               c('TRUE', 'FALSE'), selected = TRUE),
                                                   selectInput('diff_abundance_counts_transform', 'Log scale for Y axis?', 
                                                               c('TRUE', 'FALSE'), selected = TRUE),
                                                   
                                                   plotOutput("counts_plot"),
                                                   plotOutput("raw_counts_plot"))
                                          
                                        )
                              )
                          
                        )
                ),

                tabItem(tabName = "miningTab",
                        fluidRow(
                            box(
                                width = 12,
                                helpText("Coming soon!")
                            )
                        )
                )
            )
        )
    )
}
