#' @import ggplot2
#' @import shiny
#' @importFrom shinyjs useShinyjs hidden toggle
#' @export
MS2IDgui <- function(annot){
    if(!missing(annot)){
        if(class(annot) == "Annot")
            .GlobalEnv$.jmb.rawAnnotSh <- annot
        else
            stop("'annot' argument is expected to be an Annot object")
    }
    #on.exit(rm(.jmb.rawAnnotSh, envir=.GlobalEnv))
    precDigits <- 4
    ui <- navbarPage(
        "MS2ID GUI",
        id = "tabs",
        theme = bslib::bs_theme(version = 4, bootswatch = "flatly"),
        #tabPanel("Home"),
        #insert here ui side of tabPanelIdentification.R, (in invisible2git)
        tabPanel(
            "Results",
            shinyjs::useShinyjs(),
            fluidPage(
                title = 'MS2ID GUI',
                fluidRow(
                    column(
                        4,
                        tabsetPanel(
                            type = "tabs",
                            tabPanel(
                                "",
                                class = "p-3 border border-top-0 rounded-bottom",
                                icon = icon("home", lib = "font-awesome"),
                                h4("Intro"), hr(),
                                htmlOutput("info"), br(),
                                h4("Settings"), hr(),
                                fileInput("lotFile", "Upload a results file",
                                          multiple = FALSE, accept = c(".rds")
                                          ),
                                ),
                            tabPanel(
                                "QRY/REF", br(),
                                h4(HTML(paste0(
                                    .customColor("QRY", "qry"), "/",
                                    .customColor("REF", "ref"), " MSn spectra")
                                    )),
                                shinyjs::hidden(
                                    div(
                                        id = "hiddenQRYREF",
                                        #content to show/hide
                                        HTML("<b><font color=\"#DF0054\">Pick a
                                             row to show the plot.</font></b>"))
                                ),
                                plotly::plotlyOutput('plotId'), br(),
                                htmlOutput("infoTabId")
                                ),
                            tabPanel(
                                "Cons", br(),
                                h4('QRY consensus formation'),
                                shinyjs::hidden(
                                    div(
                                        id = "hiddenCons",
                                        #content to show/hide
                                        HTML("<b><font color=\"#DF0054\">Pick a row to show the plot.</font></b>")),
                                    div(
                                        id = "hiddenIsCons",
                                        #content to show/hide
                                        HTML("<b><font color=\"#DF0054\">The selected spectrum is NOT a consensus spectrum.</font></b>"))
                                ),
                                plotly::plotlyOutput('plotCons'),
                                htmlOutput("infoTabCons")
                                ),
                            tabPanel(
                                "", br(),
                                h4("Annotation arguments"),

                                htmlOutput("infoTab"),
                                icon = icon("info", lib = "font-awesome")
                                )
                            ),
                        tags$hr(),
                        ),
                    column(
                        8,
                        htmlOutput("subheaderTable"),
                        tags$hr(),
                        fluidRow(
                            column(4, selectInput('ffile', 'mzML file', NULL)),
                            column(2, selectInput('UNKprec', 'precursor MZ ',
                                                  NULL, width=120)),
                            column(6, br(), br(),
                                   checkboxInput('reduntID',
                                                 'Redundant annotations ',
                                                 value = FALSE))
                        ),
                        DT::DTOutput('tbl'),
                        ),
                    ),
                )
            ),
        tabPanel("About")
        )

    server <- function(input, output, session) {
        observe({
            shinyjs::toggle(id = c("hiddenQRYREF"),
                            condition = !length(input$tbl_rows_selected))
            shinyjs::toggle(id = c("hiddenCons"),
                            condition = !length(input$tbl_rows_selected))
            if(!length(input$tbl_rows_selected)){
                shinyjs::hide(id = c("hiddenIsCons"))
                }
        })

        options(shiny.maxRequestSize=30*1024^2)
        #parent environment values
        colVisibles <- NA
        #insert here server side of tabPanelIdentification.R, (in invisible2git)

        output$info <- renderText({.getFormText()})

        rawdata <- reactive({
            if(!is.null(input$lotFile)){
                lot <- readRDS(input$lotFile$datapath)
            }else if(exists(".jmb.rawAnnotSh", envir=.GlobalEnv))
                lot <- .GlobalEnv$.jmb.rawAnnotSh
            else
                req(FALSE)
            dt <- list(
                mtdt = .export2df(lot, summarizeHits = !input$reduntID),
                refSpctr = refSpectra(lot),
                qrySpctr = qrySpectra(lot),
                infoAnnot = infoAnnotation(lot)
                )
            dt$mtdt$massNum <- paste0(
                .customColor(dt$mtdt$QRYmassNum, "qry"), '/', dt$mtdt$cmnMasses,
                '/', .customColor(dt$mtdt$REFmassNum, "ref")
                )
            dt$mtdt$polarity <- paste0(
                .customColor(dt$mtdt$QRYpolarity, "qry"), '/',
                .customColor(dt$mtdt$REFpolarity, "ref")
                )
            dt$mtdt$collisionEnergy <- paste0(
                .customColor(dt$mtdt$QRYcollisionEnergy, "qry"), '/',
                .customColor(dt$mtdt$REFcollisionEnergy_txt, "ref")
            )
            if(all(is.na(dt$mtdt$QRYacquisitionNum))){
                dt$mtdt$QRYacquisitionNum <- dt$mtdt$QRYacquisitionNum_CONS
                }
            #round value to match with input select
            dt$mtdt$cosine <- round(dt$mtdt$cosine, 3)
            dt$mtdt <- .mergeCONS(dt$mtdt, fontsize = 1)
            #ORDER & SUBSET VISIBLE columns
            visiblVar <- c(
                "idQRYspect", "idREFspect","idREFcomp",
                INCRMETRIC, DECRMETRIC, "massNum", "propAdduct",
                "REFname", "REFformula", "collisionEnergy", "ppmPrecMass",
                "REFexactmass", "REFadduct", "REFpredicted", "REFinstrument",
                "REFID_db.comp", "REFID_db.spectra"
            )
            dt$mtdt <- dt$mtdt %>%
                select(-QRYmassNum, -cmnMasses, -REFmassNum, -QRYpolarity,
                       -REFpolarity, -QRYcollisionEnergy,
                       -REFcollisionEnergy) %>%
                #remove SOME empty columns
                select(propAdduct |
                           where(~ !(all(is.na(.)) | all(. == "")))) %>%
                select(sort(current_vars())) %>%
                relocate(visiblVar[visiblVar %in% names(dt$mtdt)]) %>%
                relocate(QRYacquisitionNum, .before = QRYrtime)
            #set columns visibility default
            colVisibles <<- names(dt$mtdt) %in% visiblVar
            return(dt)
        })

        metadataShowed <- reactive({
            req(input$UNKprec)
            rdMtdt <- rawdata()$mtdt
            rdMtdt$QRYprecursorMz <- round(rdMtdt$QRYprecursorMz, precDigits)
            rdMtdt <- rdMtdt[
                isolate(rdMtdt$QRYdataOrigin == input$ffile) &
                    rdMtdt$QRYprecursorMz == input$UNKprec, ]
            #aesthetics modific.
            if("REFexactmass" %in% names(rdMtdt))
                rdMtdt$REFexactmass <- round(rdMtdt$REFexactmass, 3)
            return(rdMtdt)
        })

        getSelId <- reactive({
            rowSelected <- input$tbl_rows_selected
            if(!length(rowSelected)){
                req(FALSE)
            }
            mtdtShw <- isolate(metadataShowed()) %>%
                slice(rowSelected) %>%
                select(idREFspect, idQRYspect, idREFcomp, QRYacquisitionNum,
                       QRYrol)
            return(mtdtShw)
            })

        output$tbl <- DT::renderDT(
            DT::datatable(
                metadataShowed(),
                caption = paste(
                    'file:', isolate(input$ffile), ', QRY precursor m/z:',
                    isolate(as.numeric(input$UNKprec))
                    ),
                colnames = c('QRYprecMz' = 'QRYprecursorMz',
                             'REFprecMz' = 'REFprecursorMz',
                             'QRYacqNum' = 'QRYacquisitionNum'),
                extensions = 'Buttons',
                options = list(
                    pageLength = 8,
                    info = FALSE,
                    lengthMenu = list(c(8,15,-1), c("8","15","All")),
                    columnDefs = list(list(visible=FALSE,
                                           targets=(which(!colVisibles)-1) )
                                      ),
                    stateLoadParams = DT::JS(
                        "function (settings, data) {return false;}"),
                    stateSave = TRUE,
                    scrollX = TRUE,
                    dom = 'Bfrt<"bottom"l>ip',
                    buttons = I('colvis')
                    ),
                selection = list(mode = 'single'),
                rownames= FALSE,
                escape = FALSE
            ) %>%
                DT::formatStyle(
                    'idQRYspect',
                    target = 'row',
                    backgroundColor = DT::styleEqual(
                        unique(metadataShowed()$idQRYspect),
                        head(
                            rep(c("#edf2fb","#ccdbfd"),
                                length(unique(metadataShowed()$idQRYspect))
                                ),
                            length(unique(metadataShowed()$idQRYspect))
                            )
                        )
                    )
            )

        output$subheaderTable <- renderText(
            paste(h3(paste("Annotations in", input$lotFile$name)))
        )

        output$infoTabId <- renderText({
            mtdtShw <- getSelId()
            rd <- isolate(rawdata()$mtdt) %>%
                filter(idREFspect == mtdtShw$idREFspect,
                       idQRYspect == mtdtShw$idQRYspect,
                       idREFcomp == mtdtShw$idREFcomp) %>%
                select(cosine, massNum, propAdduct, QRYprecursorMz, REFname,
                       REFinchikey, REFexactmass, REFadduct, REFprecursorMz,
                       REFpredicted, REFID_db.comp)
            .getFormText("compare", as.list(rd))
            #changed metadataShowed$REFinchikey[rowSelected] for metadata$REFinchikey
        })

        output$infoTab <- renderText({
            rdInfoAnnot <- isolate(rawdata()$infoAnnot)
            rdInfoAnnot$QRYdir <- rdInfoAnnot$... <- NULL
            rdInfoAnnot[rdInfoAnnot==""] <- "default"
            maxCharact <- 40
            toTrim <- rdInfoAnnot$QRYdata
            if(nchar(toTrim) > maxCharact + 3){
                rdInfoAnnot$QRYdata <- paste0(
                    "...",
                    substr(toTrim,
                           nchar(toTrim) - maxCharact + 1, nchar(toTrim)))
            }
            rdInfoAnnot$MS2ID <- basename(rdInfoAnnot$MS2ID)
            for(idArg in seq_along(rdInfoAnnot)){
                rdInfoAnnot$summ[[idArg]] <-
                    paste0("<b>", names(rdInfoAnnot)[idArg],
                           "</b>: ", rdInfoAnnot[[idArg]])
            }
            paste(rdInfoAnnot$summ, collapse = ".<br> ")
        })

        output$infoTabCons <- renderText({
            mtdtShow <- getSelId()
            .getFormText("consens",
                         list(idQRY= mtdtShow$idQRYspect,
                              acqNum = mtdtShow$QRYacquisitionNum))
        })

        output$plotId <- plotly::renderPlotly({
            mtdtShow <- getSelId()
            rd <- rawdata()
            df1 <- .getSpectra2plot(rd$refSpctr, mtdtShow$idREFspect)
            df2 <- .getSpectra2plot(rd$qrySpctr, mtdtShow$idQRYspect)
            #only mz with hits
            mMz <- .matchMz(df2$x, df1$x, isolate(rawdata()$infoAnnot$massErr))
            hits1 <- !is.na(mMz)
            hits2 <- rep(FALSE, nrow(df2))
            hits2[mMz[hits1]] <- TRUE
            suppressWarnings(
                p <- ggplot() +
                    .get_gl(df1[!hits1,]) +
                    coord_cartesian(ylim = c(-100, 100)) +
                    .get_gl(df2[!hits2,], up = FALSE) +
                    .get_gl(df1[hits1,], up = TRUE, hit=TRUE) +
                    .get_gl(df2[hits2,], up = FALSE, hit=TRUE) +
                    labs(x = "m/z", y = "% intensity") +
                    theme_bw() + #black & white theme
                    .draw_precursor(df2, rd$mtdt, mtdtShow, nature = "qry") +
                    .draw_precursor(df1, rd$mtdt, mtdtShow, nature = "ref")
                )
            plotly::ggplotly(p, tooltip = c("text"))
        })

        output$plotCons <- plotly::renderPlotly({
            mtdtShow <- getSelId()
            shinyjs::toggle(id = c("hiddenIsCons"),
                            condition = (mtdtShow$QRYrol != 4))
            req(mtdtShow$QRYrol == 4)
            rd <- rawdata()
            df2 <- .getSpectra2plot(rd$qrySpctr, mtdtShow$idQRYspect)
            srcId <- unlist(rd$qrySpctr[rd$qrySpctr$id == mtdtShow$idQRYspect]$sourceSpect)
            vsd <- rd$qrySpctr[rd$qrySpctr$id %in% srcId]
            if("basePeakIntensity" %in% Spectra::spectraVariables(vsd)){
                vsdOrdr <- order(vsd$basePeakIntensity, decreasing = F)
            }else{
                vsdOrdr <- seq_along(nrow(vsd))
            }

            my_colors <- RColorBrewer::brewer.pal(length(vsd), "Set3")
            p <- ggplot() +
                geom_linerange(data = df2,
                               aes(x = x, ymax = y, ymin = 0,
                                   text=paste('</br>m/z: ',x,'</br>i: ',
                                              round(y,2))), colour="#00786C") +
                labs(x = "m/z", y = "% intensity") +
                theme_bw()#black & white theme

            for(src in vsdOrdr){
                df <- .getSpectra2plot(vsd[src])
                p <- p +
                    geom_linerange(data = df,
                                   aes(x = x, ymax = -y, ymin = 0,
                                       alpha=0.7,
                                       text=paste('</br>m/z: ', x,
                                                  '</br>i: ', round(y,2))
                                       ),
                                   colour = my_colors[src])
            }
            plotly::ggplotly(p, tooltip = c("text"))
        })

        observeEvent(input$lotFile, ignoreNULL = F, {
            freezeReactiveValue(input, "ffile")
            updateSelectInput(session, 'ffile',
                              choices = unique(rawdata()$mtdt$QRYdataOrigin)
            )
            if(exists(".jmb.rawAnnotSh", envir=.GlobalEnv)){
                req(input$lotFile)
                rm(.jmb.rawAnnotSh, envir=.GlobalEnv)
                }
        })

        observeEvent(input$ffile, {
            freezeReactiveValue(input, "UNKprec")
            fileSel <- rawdata()$mtdt$QRYdataOrigin == input$ffile
            precSel <- unique(rawdata()$mtdt$QRYprecursorMz[fileSel])
            updateSelectInput(session, inputId = 'UNKprec',
                              choices = round(sort(precSel), precDigits)
                )
            })

        observeEvent(input$UNKprec, {
            #save columns visibility before the inminnet table refresh
            if(!is.null(isolate(input$tbl_state$columns))) {
                colVisibles <<- isolate(
                    vapply(isolate(input$tbl_state$columns),
                           function(x) x$visible, FUN.VALUE = T)
                )}
        })
        }
    shinyApp(ui, server,  options = list(launch.browser = TRUE))
}
