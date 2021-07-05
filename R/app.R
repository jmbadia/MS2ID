#' @import ggplot2
#' @import shiny
#' @export
#prova amb el boto a executar un annotation amb els arxius query i ms2id.db
MS2ID_gui <- function(){
    ui <- navbarPage(
        "MS2ID GUI",
        id = "tabs",
        theme = bslib::bs_theme(version = 4, bootswatch = "flatly"),
        tabPanel("Home"),
        #insert here ui side of tabPanelIdentification.R, (in invisible2git)
        tabPanel(
            "Results",
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
                                selectInput('ffile', 'mzML file', NULL),
                                selectInput('UNKprec', 'QRY precursor MZ ',
                                            NULL),
                                checkboxInput(
                                    'reduntID',
                                    'Show redundant identifications ',
                                    value = FALSE),
                                ),
                            tabPanel(
                                "QRY/REF", br(),
                                h4("QRY/REF MS2 spectra"),
                                plotly::plotlyOutput('plotId'), br(),
                                htmlOutput("infoTabId")
                                ),
                            tabPanel(
                                "Cons", br(),
                                h4('QRY consensus formation'),
                                plotly::plotlyOutput('plotCons'),
                                htmlOutput("infoTabCons")
                                ),
                            tabPanel(
                                "", htmlOutput("infoTab"),
                                icon = icon("info", lib = "font-awesome")
                                )
                            ),
                        tags$hr(),
                        ),
                    column(
                        8,
                        htmlOutput("subheaderTable"),
                        tags$hr(),
                        DT::DTOutput('tbl'),
                        ),
                    ),
                )
            )
        )

    server <- function(input, output, session) {
        options(shiny.maxRequestSize=30*1024^2)
        #parent environment values
        colVisibles <- NA
        #insert here server side of tabPanelIdentification.R, (in invisible2git)

        output$info <- renderText({.getFormText()})

        rawdata <- reactive({
            req(input$lotFile)
            lot <- readRDS(input$lotFile$datapath)
            dt <- list(
                mtdt = .export2df(lot), refSpctr = refSpectra(lot),
                qrySpctr = qrySpectra(lot),
                infoAnnot = infoAnnotation(lot)
                )
            #html colors for QRY and REF
            fcQ <- '<font color=\"#00786C\">'
            fcR <- '<font color=\"#d64c1d\">'
            ef <- '</font>'
            dt$mtdt$massNum <- paste0(
                fcQ, dt$mtdt$QRYmassNum, ef, '/', dt$mtdt$cmnMasses, '/', fcR,
                dt$mtdt$REFmassNum, ef
                )
            dt$mtdt$polarity <- paste0(
                fcQ, dt$mtdt$QRYpolarity, ef, '/', fcR,
                dt$mtdt$REFpolarity, ef
                )
            dt$mtdt$CE <- paste0(
                fcQ, dt$mtdt$QRYCE, ef, '/', fcR, dt$mtdt$REFCE, ef)
            if(all(is.na(dt$mtdt$QRYacquisitionNum))){
                dt$mtdt$QRYacquisitionNum <- dt$mtdt$QRYacquisitionNum_CONS
                }
            #round value to match with input select
            dt$mtdt$cosine <- round(dt$mtdt$cosine, 3)

            dt$mtdt <- dt$mtdt %>%
                select(-QRYacquisitionNum_CONS, -QRYmassNum, -cmnMasses,
                       -REFmassNum) %>%
                relocate(CE, .after=QRYrtime) %>%
                relocate(massNum, .after=REFformula) %>%
                relocate(QRYacquisitionNum, .before=QRYrtime) %>%
                relocate(idQRYspect, .before=QRYacquisitionNum) %>%
                #remove SOME empty columns
                select(propAdduct | where(~ !(all(is.na(.)) | all(. == ""))))

            #ORDER & SUBSET VISIBLE columns
            visiblVar <- c(
                "idQRYspect", "QRYrtime", "REFMmi", "propAdduct", "REFadduct",
                "REFprecursorMz", INCRMETRIC, DECRMETRIC, "REFname",
                "REFformula", "massNum", "CE", "QRYprecInt", "REFnature",
                "REFinstrument", "idREFspect","idREFcomp", "REFID_db.comp",
                "REFID_db.spect")

            #set columns visibility default
            colVisibles <<- names(dt$mtdt) %in% visiblVar
            return(dt)
        })

        metadataShowed <- reactive({
            req(input$UNKprec)
            #rdmetadata <- rawdata()$metadata
            rdmetadata <- rawdata()$mtdt
            #round value to match with input select
            rdmetadata$QRYprecursorMz <- round(rdmetadata$QRYprecursorMz, 4)
            rdmetadata <- rdmetadata[
                isolate(rdmetadata$QRYdataOrigin == input$ffile) &
                    rdmetadata$QRYprecursorMz==input$UNKprec,]
            #select QRYacq according if they are consensus or not

            #if NOT checked keep only the best score of every scan-REFmetabolite
            if(!input$reduntID){
                rdmetadata <- rdmetadata %>%
                    group_by(idQRYspect, idREFcomp) %>%
                    top_n(1, cosine) %>%
                    distinct(idQRYspect, idREFcomp, .keep_all = T) %>%
                    ungroup() %>% arrange(QRYprecursorMz, idQRYspect) %>%
                    as.data.frame()
            }
            #round leftovers
            rdmetadata$REFMmi <- round(rdmetadata$REFMmi, 4)
            rdmetadata$QRYrtime <- round(rdmetadata$QRYrtime, 2)
            rdmetadata$QRYacquisitionNum <- paste0(
                "<font size=1>", rdmetadata$QRYacquisitionNum,"</font>")
            return(rdmetadata)
        })

        getSelId <- reactive({
            rowSelected <- input$tbl_rows_selected
            if(!length(rowSelected)){
                req(FALSE)
            }
            mtdtShw <- isolate(metadataShowed()) %>%
                slice(rowSelected) %>%
                select(idREFspect, idQRYspect, idREFcomp, QRYacquisitionNum)
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
            paste(h3(paste("Identifications in", input$lotFile$name)))
        )

        output$infoTabId <- renderText({
            mtdtShw <- getSelId()
            rd <- isolate(rawdata()$mtdt) %>%
                filter(idREFspect == mtdtShw$idREFspect,
                       idQRYspect == mtdtShw$idQRYspect,
                       idREFcomp == mtdtShw$idREFcomp) %>%
                select(cosine, massNum, propAdduct, QRYprecursorMz, REFname,
                       REFinchikey, REFMmi, REFadduct, REFprecursorMz,
                       REFnature, REFID_db.comp)
            .getFormText("compare", as.list(rd))
            #changed metadataShowed$REFinchikey[rowSelected] for metadata$REFinchikey
        })

        output$infoTab <- renderText({
            rdInfoAnnot <- rawdata()$infoAnnot
            rdInfoAnnot$QRYdir <- rdInfoAnnot$... <- NULL
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
            .getFormText("consens", as.list(mtdtShow))
        })

        output$plotId <- plotly::renderPlotly({
            mtdtShow <- getSelId()
            rd <- rawdata()
            df1 <- .getSpectra2plot(rd$refSpctr, mtdtShow$idREFspect)
            df2 <- .getSpectra2plot(rd$qrySpctr, mtdtShow$idQRYspect)
            suppressWarnings(
                p <- ggplot() +
                    .get_gl(df1) +
                    coord_cartesian(ylim = c(-100, 100)) +
                    .get_gl(df2, F) +
                    labs(x = "m/z", y = "% intensity") +
                    theme_bw() + #black & white theme
                    .draw_precursor(df2, rd$mtdt, mtdtShow, nature = "qry") +
                    .draw_precursor(df1, rd$mtdt, mtdtShow, nature = "ref")
                )
            plotly::ggplotly(p, tooltip = c("text"))
        })

        output$plotCons <- plotly::renderPlotly({
            mtdtShow <- getSelId()
            rd <- rawdata()
            df2 <- .getSpectra2plot(rd$qrySpctr, mtdtShow$idQRYspect)

            srcId <- unlist(unkSpectra$sourceSpect)
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


        observeEvent(input$lotFile, {
            updateSelectInput(session, 'ffile',
                              choices = unique(rawdata()$mtdt$QRYdataOrigin)
            )
        })

        observeEvent(input$ffile, {
            freezeReactiveValue(input, "UNKprec")
            fileSel <- rawdata()$mtdt$QRYdataOrigin == input$ffile
            precSel <- unique(rawdata()$mtdt$QRYprecursorMz[fileSel])
            updateSelectInput(session, inputId = 'UNKprec',
                              choices = round(precSel, 4)
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
    shinyApp(ui,server)
}
